#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx
import csv

def data_preprocessing(path_expr, path_net,log2, size = 2000):
    """
    Raw data processing for further analysis
    
    Attributes:
    -----------
    non-default:
    path_expr - path for gene expression
    ATTENTION: expression data format: genes as rows (first column - gene ids), patients as columns
    path_ppi - path for ppi
    log2 - log2 transform (if needed)
    size -   specify size of the gene set  for standard deviation preselection. USually optimal values are between 2000 and 5000
    """
    
    expr = open_file(path_expr)
    
    expr = expr.set_index(expr.columns[0])
    patients_new = list(set(expr.columns))

        
    net = open_file(path_net, header = None)
    nodes_ppi = list(set(net[0]).union(set(net[1])))
    genes_ge = list(expr.index)
    new_genes_ge = set([str(x) for x in genes_ge])
    new_genes_ppi = set([str(x) for x in nodes_ppi])
    intersec_genes = list(set.intersection(new_genes_ge, new_genes_ppi))
    expr.index = [str(x) for x in genes_ge]
    expr = expr.loc[intersec_genes]
    if log2:
        expr = np.log2(expr)        
    
    if size!= None: #std selection
        std_genes = expr.std(axis = 1)
        std_genes, intersec_genes = zip(*sorted(zip(std_genes, intersec_genes)))
        genes_for_expr = list(intersec_genes)[len(std_genes)-size:]
    else:
        genes_for_expr = intersec_genes
        
    expr = expr.loc[genes_for_expr]

    expr = pd.DataFrame(stats.zscore(expr) ,columns = expr.columns, index = expr.index)
    
    labels = dict()
    rev_labels = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in genes_for_expr:
       labels[node] = g
       rev_labels[g] = node
       node = node+1
    for p in patients_new:
       labels[node] = p
       rev_labels[p] = node
       node = node+1
    n,m = expr.shape    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(genes_for_expr) and node2 in set(genes_for_expr):    
            G.add_edge(rev_labels[node1],rev_labels[node2])
    expr.index = np.arange(n)
    expr.columns =  np.arange(n,n+m)
    return expr,G,labels, rev_labels

# allows to determine the delimeter automatically given the path or directly the object
def open_file(file_name, **kwards):
    if isinstance(file_name, str):
        with open(file_name, 'r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
    else:
        file_name.seek(0)
        dialect = csv.Sniffer().sniff(file_name.read(1024))
    file_name.seek(0)

    file = pd.read_csv(file_name,sep = dialect.delimiter, **kwards)
    return file
    