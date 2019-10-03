import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

patch_0 = mpatches.Patch(color='#ff0000', label='-2.0')
patch_1 = mpatches.Patch(color='#ff9933', label='-1.5')
patch_2 = mpatches.Patch(color='#ffcc00', label='-1')
patch_3 = mpatches.Patch(color='#ffff00', label='-0.5')
patch_4 = mpatches.Patch(color='#ccff33', label='0.5')
patch_5 = mpatches.Patch(color='#99ff33', label='1')
patch_6 = mpatches.Patch(color='#66ff33', label='1.5')
patch_7 = mpatches.Patch(color='#33cc33', label='2')

plt.legend(title="Node Coloring",handles=[patch_0, patch_1, patch_2, patch_3, patch_4, patch_5, patch_6, patch_7])
plt.show()
