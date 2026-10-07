// 2R2C 열회로 모델로 일반 서모스탯(히스테리시스)과 위상 평면 슬라이딩 제어를 비교한다.
// 빌드: g++ -std=c++17 -O2 sim.cpp -o sim
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

// ---- 물리 파라미터 (기숙사 1실 + 대형 라디에이터) ----
const double DT = 10.0;        // 샘플링·적분 간격 [s]
const double P = 1500.0;       // 히터 출력 [W]
const double C_h = 7.2e4;      // 히터 열용량 [J/K]
const double R_h = 0.05;       // 히터→실내 열저항 [K/W]  (tau_h = 60분)
const double C_r = 5.4e5;      // 실내 열용량 [J/K]
const double R_o = 0.02;       // 실내→외기 열저항 [K/W]  (tau = 3시간)
const double T_TARGET = 22.0;  // 목표 온도 [C]
const double HOURS = 24.0;

// ---- 제어 파라미터 ----
const double BAND = 0.5;          // 히스테리시스 불감대 [K]
const double LOCKOUT = 300.0;     // 최소 유지 시간 [s]
const double ALPHA = 0.3;         // EMA 계수
const int SLOPE_WIN = 36;         // 기울기 창: 36샘플 = 6분
const double C_SLIDE = 1.0 / 20;  // 슬라이딩 계수 c [1/min] (예측 지평 20분)
const double DELTA_S = 0.02;      // 슬라이딩 히스테리시스 [K/min]

double outdoor(double t) {  // 환절기 외기: 평균 8C, 일교차 +-4C, 15시 최고
    return 8.0 + 4.0 * std::cos(2 * M_PI * (t / 3600.0 - 15.0) / 24.0);
}

struct Sample { double t, T_out, T_in, T_meas, e, edot; int u; };

enum Mode { HYSTERESIS, SLIDING };

std::vector<Sample> run(Mode mode) {
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.05);

    double T_in = 15.0, T_h = 15.0, T_f = T_in;
    int u = 0;
    double last_switch = -LOCKOUT;
    std::vector<double> hist;  // 필터링된 온도 기록 (기울기 계산용)
    std::vector<Sample> out;

    for (double t = 0; t <= HOURS * 3600; t += DT) {
        // 센서: 잡음 + DHT22 분해능 0.1C 양자화, 그 뒤 EMA
        double T_meas = std::round((T_in + noise(rng)) * 10) / 10;
        T_f = ALPHA * T_meas + (1 - ALPHA) * T_f;
        hist.push_back(T_f);

        // 최근 SLOPE_WIN 샘플의 최소제곱 기울기 [K/min]
        double slope = 0;
        int n = std::min<int>(hist.size(), SLOPE_WIN);
        if (n >= 2) {
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int i = 0; i < n; ++i) {
                double x = i * DT / 60.0, y = hist[hist.size() - n + i];
                sx += x; sy += y; sxx += x * x; sxy += x * y;
            }
            slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        }

        double e = T_TARGET - T_f, edot = -slope;
        int want = u;
        if (mode == HYSTERESIS) {
            if (T_f < T_TARGET - BAND) want = 1;
            if (T_f > T_TARGET + BAND) want = 0;
        } else {
            double S = C_SLIDE * e + edot;  // S > 0  <=>  1/c 분 뒤 예측 오차 > 0
            if (S > DELTA_S) want = 1;
            if (S < -DELTA_S) want = 0;
        }
        if (want != u && t - last_switch >= LOCKOUT) { u = want; last_switch = t; }

        double T_out = outdoor(t);
        out.push_back({t, T_out, T_in, T_meas, e, edot, u});

        // 2R2C 모델 오일러 적분
        double q_hr = (T_h - T_in) / R_h;
        T_h += DT * (P * u - q_hr) / C_h;
        T_in += DT * (q_hr - (T_in - T_out) / R_o) / C_r;
    }
    return out;
}

void report(const char* name, const std::vector<Sample>& s) {
    double energy = 0, max_T = -1e9, discomfort = 0;
    int switches = 0;
    bool reached = false;
    for (size_t i = 0; i < s.size(); ++i) {
        energy += P * s[i].u * DT;
        if (i > 0 && s[i].u != s[i - 1].u) ++switches;
        if (s[i].T_in >= T_TARGET) reached = true;
        if (!reached) continue;  // 초기 예열 구간은 쾌적도·오버슈팅에서 제외
        max_T = std::max(max_T, s[i].T_in);
        discomfort += std::max(0.0, std::fabs(s[i].T_in - T_TARGET) - BAND) * DT / 60.0;
    }
    std::printf("%-12s %8.2f %10.2f %12.1f %8d\n", name, max_T - T_TARGET,
                energy / 3.6e6, discomfort, switches);
}

int main() {
    auto hyst = run(HYSTERESIS), slide = run(SLIDING);

    std::printf("%-12s %8s %10s %12s %8s\n", "controller", "over[K]", "energy[kWh]",
                "discomf[K*min]", "switches");
    report("hysteresis", hyst);
    report("sliding", slide);

    FILE* f = std::fopen("results.csv", "w");
    std::fprintf(f, "t_min,T_out,hyst_T,hyst_u,slide_T,slide_u,slide_e,slide_edot\n");
    for (size_t i = 0; i < hyst.size(); ++i)
        std::fprintf(f, "%.2f,%.3f,%.3f,%d,%.3f,%d,%.4f,%.5f\n", hyst[i].t / 60,
                     hyst[i].T_out, hyst[i].T_in, hyst[i].u, slide[i].T_in, slide[i].u,
                     slide[i].e, slide[i].edot);
    std::fclose(f);
}
