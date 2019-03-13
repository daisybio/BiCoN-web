
# -*- coding: utf-8 -*-

import time
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
flatten = lambda l: [item for sublist in l for item in sublist]
import polls.weighted_aco_lib as lib
import imp
import seaborn as sns; sns.set(color_codes=True)
import mygene
from polls.tasks import ants_2
#print("Enter 1 for lung cancer data analysis :")
#print("Enter 2 for breast cancer data analysis :")
#s = input()
#print("data is loading")
def run_algo(s,L_g_min,L_g_max):
	if s ==1:
	    # UPLOAD LUNG CANCER DATA 
	    path_expr = "lung_cancer_expr.csv"
	    path_ppi = "biogrid.human.entrez.tsv"
	    #path_expr = "/home/quirin/testproject/polls/data/lung_cancer_expr.csv"
	    #path_ppi = "/home/quirin/testproject/polls/data/biogrid.human.entrez.tsv"
	    #path_expr = "data/lung_cancer_expr.csv"
	    #path_ppi = "data/biogrid.human.entrez.tsv"
	    
	    
	    col = "cancer_type"
	    size = 2000
	    log2 = True
	    B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(path_expr, path_ppi, col,log2 = True, gene_list = None, size = size, sample= None)
	
	
	else:
	    #UPLOAD BREAST CANCER 
	    imp.reload(lib)
	    path_expr = "breast_cancer_expr.csv"
	    path_ppi = "biogrid.human.entrez.tsv"
	    #path_expr = "/home/quirin/testproject/polls/data/lung_cancer_expr.csv"
	    #path_ppi = "/home/quirin/testproject/polls/data/biogrid.human.entrez.tsv"
	    #path_expr = "data/breast_cancer_expr.csv"
	    #path_ppi = "data/biogrid.human.entrez.tsv"
	    
	    
	    col = "prognosis"
	    size = 2000
	    log2 = True
	    B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(path_expr, path_ppi, col,log2 = True, gene_list = None, size = size, sample= None)
	
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = 30 # number of ants
	eps = 0.02 # stopping criteria: score_max-score_av<eps
	b = 1 #HI significance
	evaporation  = 0.5
	a = 1 #pheramone significance
	times = 5 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	
	
	
	
	start = time.time()
	solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = False, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	print(lib.jac_matrix(solution[1],[group1,group2]))
	
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	print(means1)
	print(means2)
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	for label in values:
	    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	                            label=label, linewidth=0)
	g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	                            borderaxespad=0.)
	ax = g.ax_heatmap
	ax.set_xlabel("Genes")
	ax.set_ylabel("Patients")
	
	#plotting convergence
	fig, ax = plt.subplots(figsize=(10, 7))
	bplot1 = ax.boxplot(conv/2,
	                         vert=True,  # vertical box alignment
	                         patch_artist=True)  # will be used to label x-ticks
	
	
	plt.xlabel("iterations")
	plt.ylabel("score per subnetwork")
	plt.show()
	plt.savefig("/home/quirin/testproject/polls/static/algo_output.png")
	return(plt.gcf())



def algo_output(s,L_g_min,L_g_max):
	if s ==1:
	    # UPLOAD LUNG CANCER DATA 
	    path_expr = "lung_cancer_expr.csv"
	    path_ppi = "biogrid.human.entrez.tsv"
	    #path_expr = "/home/quirin/testproject/polls/data/lung_cancer_expr.csv"
	    #path_ppi = "/home/quirin/testproject/polls/data/biogrid.human.entrez.tsv"
	    #path_expr = "data/lung_cancer_expr.csv"
	    #path_ppi = "data/biogrid.human.entrez.tsv"
	    
	    
	    col = "cancer_type"
	    size = 2000
	    log2 = True
	    B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(path_expr, path_ppi, col,log2 = True, gene_list = None, size = size, sample= None)
	
	
	else:
	    #UPLOAD BREAST CANCER 
	    imp.reload(lib)
	    path_expr = "breast_cancer_expr.csv"
	    path_ppi = "biogrid.human.entrez.tsv"
	    #path_expr = "/home/quirin/testproject/polls/data/lung_cancer_expr.csv"
	    #path_ppi = "/home/quirin/testproject/polls/data/biogrid.human.entrez.tsv"
	    #path_expr = "data/breast_cancer_expr.csv"
	    #path_ppi = "data/biogrid.human.entrez.tsv"
	    
	    
	    col = "prognosis"
	    size = 2000
	    log2 = True
	    B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(path_expr, path_ppi, col,log2 = True, gene_list = None, size = size, sample= None)
	
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = 30 # number of ants
	eps = 0.02 # stopping criteria: score_max-score_av<eps
	b = 1 #HI significance
	evaporation  = 0.5
	a = 1 #pheramone significance
	times = 2 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	
	
	
	
	start = time.time()
	solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = False, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	print(lib.jac_matrix(solution[1],[group1,group2]))
	
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	ret2 = means1 + means2
	ret3 = new_genes1 + new_genes2
	adjlist = []
	for line99 in nx.generate_edgelist(G_small,data=False):	
		lineSplit = line99.split()
		adjlist.append([lineSplit[0],lineSplit[1]])
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	#g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	#for label in values:
	#    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	#                            label=label, linewidth=0)
	#g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	#                            borderaxespad=0.)
	#ax = g.ax_heatmap
	#ax.set_xlabel("Genes")
	#ax.set_ylabel("Patients")
	
	#plotting convergence
	#fig, ax = plt.subplots(figsize=(10, 7))
	#bplot1 = ax.boxplot(conv/2,
	 #                        vert=True,  # vertical box alignment
	 #                        patch_artist=True)  # will be used to label x-ticks
	#
	
	#plt.xlabel("iterations")
	#plt.ylabel("score per subnetwork")
	#plt.show()
	#plt.savefig("/home/quirin/testproject/polls/static/algo_output.png")
	return(GE_small.T,row_colors,col_colors,G_small, ret2, ret3, adjlist,new_genes1)



def algo_output_ownfile(s,L_g_min,L_g_max,fh,prot_fh,nbr_iter):
	col = "cancer_type"
	size = 2000
	log2 = True
	#B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing_ownfile(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = 30 # number of ants
	eps = 0.02 # stopping criteria: score_max-score_av<eps
	b = 1 #HI significance
	evaporation  = 0.5
	a = 1 #pheramone significance
	times =int(nbr_iter) #max amount of iterations
	#times =45 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	
	
	
	
	start = time.time()
	solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = False, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	print(lib.jac_matrix(solution[1],[group1,group2]))
	
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	ret2 = means1 + means2
	ret3 = new_genes1 + new_genes2
	adjlist = []
	for line99 in nx.generate_edgelist(G_small,data=False):	
		lineSplit = line99.split()
		adjlist.append([lineSplit[0],lineSplit[1]])
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	plt.savefig("/home/quirin/testproject/polls/static/ntw.png")
	#g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	#for label in values:
	#    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	#                            label=label, linewidth=0)
	#g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	#                            borderaxespad=0.)
	#ax = g.ax_heatmap
	#ax.set_xlabel("Genes")
	#ax.set_ylabel("Patients")
	
	#plotting convergence
	#fig, ax = plt.subplots(figsize=(10, 7))
	plt.clf()
	plt.boxplot(conv/2,
	                        vert=True,  # vertical box alignment
	                        patch_artist=True)  # will be used to label x-ticks
	
	plt.xlabel("iterations")
	plt.ylabel("score per subnetwork")
	#plt.show(bplot1)
	plt.savefig("/home/quirin/testproject/polls/static/conv.png")
	return(GE_small.T,row_colors,col_colors,G_small, ret2, ret3, adjlist,new_genes1)



def algo_output_ownfile_2(s,L_g_min,L_g_max,fh,prot_fh,nbr_iter,nbr_ants,evap):
	col = "cancer_type"
	size = 2000
	log2 = True
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Your files are being processed...")
	#B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing_ownfile(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Starting model run. Progress of the algorithm is shown below...")	
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = int(nbr_ants) # number of ants
	eps = 0.02 # stopping criteria: score_max-score_av<eps
	b = 1 #HI significance
	evaporation  = float(evap)
	a = 1 #pheramone significance
	times =int(nbr_iter) #max amount of iterations
	#times =45 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	print(a)
	print(type(a))
	print(type(b))
	print(type(n))
	print(type(m))
	print(type(H))
	print(H)
	print(type(GE))	
	print(GE)
	print(type(G))
	print(G)
	print(type(clusters))
	print(clusters)
	

	
	
	
	start = time.time()
	#solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	result_99= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = True, save 	= None, show_nets = False)
	print(result_99.get())
	solution,t_best,sc,conv= result_99.get()
	#solution,t_best,sc,conv= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	print(lib.jac_matrix(solution[1],[group1,group2]))
	
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	ret2 = means1 + means2
	ret3 = new_genes1 + new_genes2
	adjlist = []
	for line99 in nx.generate_edgelist(G_small,data=False):	
		lineSplit = line99.split()
		adjlist.append([lineSplit[0],lineSplit[1]])
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	plt.savefig("/home/quirin/testproject/polls/static/ntw.png")
	#g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	#for label in values:
	#    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	#                            label=label, linewidth=0)
	#g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	#                            borderaxespad=0.)
	#ax = g.ax_heatmap
	#ax.set_xlabel("Genes")
	#ax.set_ylabel("Patients")
	
	#plotting convergence
	#fig, ax = plt.subplots(figsize=(10, 7))
	plt.clf()
	plt.boxplot(conv/2,
	                        vert=True,  # vertical box alignment
	                        patch_artist=True)  # will be used to label x-ticks
	
	plt.xlabel("iterations")
	plt.ylabel("score per subnetwork")
	#plt.show(bplot1)
	plt.savefig("/home/quirin/testproject/polls/static/conv.png")
	return(GE_small.T,row_colors,col_colors,G_small, ret2, ret3, adjlist,new_genes1)



def algo_output_ownfile_3(s,L_g_min,L_g_max,expr_str,ppi_str,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig):
	col = "cancer_type"
	size = 2000
	log2 = True
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Your files are being processed...")
	#B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing_strings(expr_str, ppi_str, col,log2 = True, gene_list = None, size = size, sample= None)
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Starting model run. Progress of the algorithm is shown below...")	
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = int(nbr_ants) # number of ants
	eps = float(epsilon) # stopping criteria: score_max-score_av<eps
	b = float(hi_sig) #HI significance
	evaporation  = float(evap)
	a = float(pher_sig) #pheramone significance
	times =int(nbr_iter) #max amount of iterations
	#times =45 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	print(a)
	print(type(a))
	print(type(b))
	print(type(n))
	print(type(m))
	print(type(H))
	print(H)
	print(type(GE))	
	print(GE)
	print(type(G))
	print(G)
	print(type(clusters))
	print(clusters)
	

	
	
	
	start = time.time()
	#solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	result_99= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = True, save 	= None, show_nets = False)
	print(result_99.get())
	solution,t_best,sc,conv= result_99.get()
	#solution,t_best,sc,conv= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	print(lib.jac_matrix(solution[1],[group1,group2]))
	
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	ret2 = means1 + means2
	ret3 = new_genes1 + new_genes2
	adjlist = []
	for line99 in nx.generate_edgelist(G_small,data=False):	
		lineSplit = line99.split()
		adjlist.append([lineSplit[0],lineSplit[1]])
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	plt.savefig("/home/quirin/testproject/polls/static/ntw.png")
	#g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	#for label in values:
	#    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	#                            label=label, linewidth=0)
	#g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	#                            borderaxespad=0.)
	#ax = g.ax_heatmap
	#ax.set_xlabel("Genes")
	#ax.set_ylabel("Patients")
	
	#plotting convergence
	#fig, ax = plt.subplots(figsize=(10, 7))
	plt.clf()
	plt.boxplot(conv/2,
	                        vert=True,  # vertical box alignment
	                        patch_artist=True)  # will be used to label x-ticks
	
	plt.xlabel("iterations")
	plt.ylabel("score per subnetwork")
	#plt.show(bplot1)
	plt.savefig("/home/quirin/testproject/polls/static/conv.png")
	return(GE_small.T,row_colors,col_colors,G_small, ret2, ret3, adjlist,new_genes1)
