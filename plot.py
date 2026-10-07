"""sim이 만든 results.csv를 읽어 result.png로 그린다."""
import csv

import matplotlib.pyplot as plt

T_TARGET, BAND, C_SLIDE = 22.0, 0.5, 1.0 / 20  # sim.cpp와 같은 값
SLIDE, HYST, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"

with open("results.csv") as f:
    rows = list(csv.DictReader(f))
col = {k: [float(r[k]) for r in rows] for k in rows[0]}
hours = [t / 60 for t in col["t_min"]]

plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                     "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(11, 7))
gs = fig.add_gridspec(2, 2, width_ratios=[2, 1], height_ratios=[3, 1])

ax = fig.add_subplot(gs[0, 0])
ax.axhspan(T_TARGET - BAND, T_TARGET + BAND, color="#e8e7e2", lw=0, label="target ±0.5 °C")
ax.plot(hours, col["hyst_T"], color=HYST, lw=1.5, label="Hysteresis thermostat")
ax.plot(hours, col["slide_T"], color=SLIDE, lw=1.5, label="Sliding control")
ax.set(ylabel="Indoor temperature (°C)", xlim=(0, 24), ylim=(20.5, 23.6),
       title="Indoor temperature (2R2C model, τ_heater = 60 min)")
ax.legend(frameon=False, loc="lower right")
ax.grid(axis="y", color="#e8e7e2", lw=0.8)

ax = fig.add_subplot(gs[1, 0], sharex=ax)
ax.step(hours, [u + 1.5 for u in col["hyst_u"]], where="post", color=HYST, lw=1.2)
ax.step(hours, col["slide_u"], where="post", color=SLIDE, lw=1.2)
ax.set(yticks=[0.5, 2.0], yticklabels=["Sliding", "Hysteresis"], xlabel="Time (h)",
       title="Heater relay (ON = up)")
ax.tick_params(axis="y", length=0)

ax = fig.add_subplot(gs[:, 1])
e, edot = col["slide_e"], col["slide_edot"]
ax.plot(e, edot, color=SLIDE, lw=0.8)
ax.plot(e[0], edot[0], "o", color=SLIDE, ms=8)
ax.annotate("start", (e[0], edot[0]), textcoords="offset points", xytext=(-8, 8),
            ha="right", color=MUTED)
xs = [-1, 8]
ax.plot(xs, [-C_SLIDE * x for x in xs], "--", color=MUTED, lw=1, label="S = c·e + ė = 0")
ax.axhline(0, color="#d0cfc8", lw=0.8)
ax.axvline(0, color="#d0cfc8", lw=0.8)
ax.set(xlabel="e = T_target − T_in (°C)", ylabel="ė (°C/min)", xlim=xs, ylim=(-0.12, 0.1),
       title="Phase plane (sliding control)")
ax.legend(frameon=False, loc="upper right")

fig.tight_layout()
fig.savefig("result.png", dpi=150)
