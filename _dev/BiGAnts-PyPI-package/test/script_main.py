#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd


#matplotlib.use('agg')
from bigants import data_preprocessing
from bigants import BiGAnts

#from ants import BiGAnts


#Set the paths for the data
path_expr,path_net ='../bigants/data/gse30219_lung.csv', '../bigants/data/biogrid.human.entrez.tsv'

#Load and process the data
GE,G,labels, _= data_preprocessing(path_expr, path_net,log2 = False, size = 2000)

# Set the size of subnetworks
L_g_min = 10
L_g_max = 15

model = BiGAnts(GE,G,L_g_min,L_g_max)

'''
Run the optimisation on default parameters or specify additionally:
        K - number of ants (less ants - less space exploration. Usually set between 20 and 50, default - 20)        
        n_proc = number of processes that should be used (default 1)
        a - pheromone significance (default 1 - does not need to be changed)
        b - heuristic information significance (default 1 - does not need to be changed)
        evaporation - the rate at which pheromone evaporates (default 0.5)
        th - similarity threshold (default 0.5 - does not need to be changed)
        eps - conservative convergence criteria: score_max - score_min < eps (default- 0.02)
        times - allows faster convergence criteria: stop if the maximum so far was reached more than x times (default 6)
        clusters - # of clusters, right now does not work for more than 2
        cost_limit - defines the radius of the search for ants (default 5)
        opt - given if the best score is known apriori (in case of simulated data for instance)
        show_pher - set true if plotting of pheromone heatmap at every iteration is desirable (strickly NOT recommended for # of genes > 2000)
        show_plot - set true if convergence plots should be shown
        save - set an output file name  if the convergence plot should be saved in the end
        show_nets - set true if the selected network should be shown at each iteration       
'''

solution,sc= model.run_search(max_iter = 2)

#Coming back to the original IDs
patients1 = [str(labels[x]) for x in solution[1][0]]
patients1 = "|".join(patients1)

patients2 = [str(labels[x]) for x in solution[1][1]]
patients2 = "|".join(patients2)

genes1 = [str(labels[x]) for x in solution[0][0]]
genes1 = "|".join(genes1)


genes2 = [str(labels[x]) for x in solution[0][1]]
genes2 = "|".join(genes2)

#Saving the solution
pd.DataFrame([[genes1,genes2,patients1,patients2]],columns = ["genes1","genes2","patients1","patients2"]).to_csv("results.csv")

