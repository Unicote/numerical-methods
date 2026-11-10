import json
import sys

import matplotlib.pyplot as plt

figure = json.load(sys.stdin)
fig, ax = plt.subplots(figsize=(10, 6.25))

for s in figure["series"]:
    color = s["color"] or None
    if s["points"]:
        size = (2 * (s["width"] or 5)) ** 2
        ax.scatter(s["x"], s["y"], s=size, color=color, edgecolors="white", zorder=3, label=s["name"])
    else:
        style = "--" if s["dashed"] else "-"
        ax.plot(s["x"], s["y"], style, color=color, linewidth=s["width"] or 2, label=s["name"])

ax.set_title(figure["title"])
ax.set_xlabel(figure["xlabel"])
ax.set_ylabel(figure["ylabel"])
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
fig.tight_layout()
fig.savefig(figure["path"])
