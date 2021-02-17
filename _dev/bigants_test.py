import inspect
import os

from bicon import data_preprocessing
from bicon import BiCoN
from bicon import results_analysis

# print(inspect.signature(results_analysis.show_clustermap))

path_expr, path_net = 'files/gse30219_lung.csv', 'files/biogrid.human.entrez.tsv'

GE, G, labels, _ = data_preprocessing(path_expr, path_net)

L_g_min = 10
L_g_max = 15

# --- Run model
model = BiCoN(GE, G, L_g_min, L_g_max)
solution, scores = model.run_search()

# --- Vizualize the results

results = results_analysis(solution, labels)
# Networks
results.show_networks(GE, G)
# Heatmap
results.show_clustermap(GE, G, true_labels=labels, output="tmp/clustermap.png")
# Convergence
results.convergence_plot(scores)


