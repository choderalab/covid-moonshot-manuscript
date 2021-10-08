import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib.colors import LinearSegmentedColormap


# fig, ax = plt.subplots(figsize=(6, 1.5))
# fig.subplots_adjust(bottom=0.5)

fig, ax = plt.subplots(figsize=(1.75, 6))
# fig.subplots_adjust(bottom=0.5)

c = ["magenta", "white"]
cmap = LinearSegmentedColormap.from_list("custom_cmap", c)
norm = mpl.colors.Normalize(vmin=0, vmax=1)

# cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              cax=ax, orientation='horizontal')

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax, orientation='vertical')

cbar.set_ticks([0,1])
# cbar.ax.set_xticklabels(['Low', 'High'], size=18)
cbar.ax.set_yticklabels(['Low', 'High'], size=18)
cbar.set_label(label='Conservation',size=22)

plt.tight_layout()

plt.savefig("conservation_cbar.png", dpi=300)