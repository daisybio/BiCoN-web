### this file contains methods for all computationally intensive tasks needed for algorithm runs:
### - algo_output_task - this is needed for running the algorithm based on PPI and expression data. it outputs
###			arrays etc with the algorithm results.
### - script_output_task - this is needed for processing the outputs of algo_output_task to formats used
###  			for data vizualisation. it writes heatmap, ppi graph, survival plot and metadata to 
###			files and outputs links to those files (and an array with metadata). it takes a parameter
###			"session_id" that is included in the path to result files (e.g. "ppi_[SESSION_ID].json").
###			this parameter can be set to "none" if you do not want to use sessions. then it uses
###			static paths (e.g. "ppi.json") instead.
### - import_ndex - this tasks imports PPI files from NDEx based on the UUID and parses them to the correct
###			input format for the algorithm tasks.
### - check_input_files - this task checks given expression and PPI files if they contain data and returns an error
###			string if they do not.
### - preprocess_file - this task preprocesses an input expression data file and tries to find a column with
###			pre-defined clusters. if found, it is renamed to "disease_type".
### - preprocess_file_2 - the same as preprecess_file, but it returns a number of pre-defined clusters additionally.
### - preprocess_ppi_file - preprocesses the PPI file. it finds every row that contains two tab-separated integers (protein IDs)
###			and appends them to the output file.
### - preprocess_clinical_file - converts clinical file to TSV format.
### - list_metadata_from_file - reads metadata from a file and returns 3 arrays with variable names and their frequency
###			in cluster 1 and 2.
### - run_enrichment - runs an enrichment analysis usen given terms on a list of genes.
### - read_enrichment - reads results of enrichment analysis for one patient cluster and outputs dictionary with results
### - read_enrichment_2 - the same as read_enrichment, but reads terms that appear only in cluster 1 or only in cluster 2


import string
from io import StringIO
from multiprocessing import Pool
import time
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite

flatten = lambda l: [item for sublist in l for item in sublist]
import seaborn as sns;

sns.set(color_codes=True)
# from multiprocessing import Pool
from numpy import linalg as LA
# from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import clustering.weighted_aco_lib as lib
from shutil import copyfile
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline
from datetime import datetime
from networkx.readwrite import json_graph
import json
import imp
import seaborn as sns;

sns.set(color_codes=True)
import mygene
from django.contrib.auth.models import User
from django.utils.crypto import get_random_string
import os
from celery import shared_task
import numpy as np
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test
import math
import gseapy as gp
import seaborn as sns
import pandas as pd
from numpy import array
import matplotlib.patches as mpatches
import networkx as nx
from pybiomart import Dataset
from ndex2.nice_cx_network import NiceCXNetwork
import ndex2.client as nc
import ndex2


# make an empty progress file and include session ID in the path
@shared_task(name="make_empty_figure")
def make_empty_figure(session_id):
    fig = plt.figure(figsize=(10, 8))
    if (session_id == "none"):
        plt.savefig("/code/clustering/static/progress.png")
        plt.close(fig)
    else:
        plt.savefig("/code/clustering/static/userfiles/progress_" + session_id + ".png")
        plt.close(fig)


# empty the log file (session ID in path)
@shared_task(name="empty_log_file")
def empty_log_file(session_id):
    if (session_id == "none"):
        text_file = open("/code/clustering/static/output_console.txt", "w")
        text_file.write("")
        text_file.close()
        if (os.path.isfile("/code/clustering/static/userfiles/output_console.txt")):
            text_file = open("/code/clustering/static/output_console.txt", "w")
            text_file.write("")
            text_file.close()
    else:
        text_file = open("/code/clustering/static/output_console_" + session_id + ".txt", "w")
        text_file.write("")
        text_file.close()
        if (os.path.isfile("/code/clustering/static/userfiles/output_console_" + session_id + ".txt")):
            text_file = open("/code/clustering/static/userfiles/output_console_" + session_id + ".txt", "w")
            text_file.write("")
            text_file.close()


########################################################
#### writing and processing metadata, loading 
#### images etc...   ###################################
#### the *actual* algorithm will be further
#### down the page       ###############################
########################################################

# preprocess file with metadata
@shared_task(name="preprocess_clinical_file")
def preprocess_clinical_file(clinical_str):
    if (len(clinical_str.split("\n")[0].split("\t")) > 2):
        return (clinical_str)
    elif ("," in clinical_str):
        # replace comma by tab if file is CSV and not TSV
        if ("\t" not in clinical_str.split("\n")[0]):
            clinical_str = clinical_str.replace(",", "\t")
    return clinical_str


# preprocess PPI file
@shared_task(name="preprocess_ppi_file")
def preprocess_ppi_file(ppistr):
    ppistr_split = ppistr.split("\n")
    # check if file is csv or tsv and convert it to tsv format
    if ("\t" not in ppistr_split[2]):
        ppistr = ppistr.replace(",", "\t")
        ppistr_split = ppistr.split("\n")
    ppistr_split_new = []
    len_first_line = 0
    len_second_line = 0
    for elem in ppistr_split[0].split("\t"):
        if (elem != ""):
            len_first_line = len_first_line + 1
    for elem in ppistr_split[1].split("\t"):
        if (elem != ""):
            len_second_line = len_second_line + 1
    # delete first row if it contains title columns
    if (len_second_line > len_first_line):
        del ppistr_split[0]
    # take only the right two columns
    for line in ppistr_split:
        if (len(line.split("\t")) > 1):
            if (len(line.split("\t")) > 2):
                line_length = len(line.split("\t"))
                line = "\t".join([line.split("\t")[line_length - 2], line.split("\t")[line_length - 1]])
            # check if line contains two integers with protein IDs
            if (str(line.split("\t")[0]).isdigit() and str(line.split("\t")[1].strip().replace("\n", "")).isdigit()):
                ppistr_split_new.append(line)
    ppistr = "\n".join(ppistr_split_new)
    return (ppistr)


# Method to convert expression data file to TSV format and find and rename (for later recongition) column with disease type information	
@shared_task(name="preprocess_file")
def preprocess_file(expr_str):
    expr_str = expr_str.replace("cancer_type", "disease_type")
    if (len(expr_str.split("\n")[0].split("\t")) > 2):
        expr_str_split = expr_str.split("\n")
        # replace column name for disease type
        if ("disease_type" not in expr_str_split[0]):
            if ("subtype" in expr_str_split[0]):
                expr_str = expr_str.replace("subtype", "disease_type")
        # remove name of first column (left upper corner)
        expr_str_first_colname = expr_str_split[0].split("\t")[0]
        expr_str = expr_str.replace(expr_str_first_colname, "", 1)
        expr_stringio = StringIO(expr_str)
        exprdf = pd.read_csv(expr_stringio, sep='\t')
        # check for column with two unique entries (pre-defined clusters)
        for column_name, column in exprdf.transpose().iterrows():
            if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
                if (len(column.unique()) == 2):
                    expr_str = expr_str.replace(column_name, "disease_type")
    elif ("," in expr_str):
        # replace comma by tab if file is CSV and not TSV
        if ("\t" not in expr_str.split("\n")[0]):
            expr_str_split = expr_str.split("\n")
            # replace "subtype" by "cancer type"
            if ("disease_type" not in expr_str):
                if ("subtype" in expr_str):
                    expr_str = expr_str.replace("subtype", "disease_type")
            expr_str_first_colname = expr_str_split[0].split(",")[0]
            expr_str = expr_str.replace(expr_str_first_colname, "", 1)
            expr_str = expr_str.replace(",", "\t")
            expr_str_split = expr_str.split("\n")
            # remove entries after given length if expression data file is too big
            if (len(expr_str_split) > 300):
                expr_str = "\n".join(expr_str_split[:200])
            else:
                expr_str = "\n".join(expr_str_split)
            expr_stringio = StringIO(expr_str)
            expr_str = expr_str.replace("MCI", "CTL")
            exprdf = pd.read_csv(expr_stringio, sep='\t')
            # find column with two unique entries that represents disease type
            for column_name, column in exprdf.transpose().iterrows():
                if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
                    if (len(column.unique()) == 2):
                        expr_str = expr_str.replace(column_name, "disease_type")
            #### uncomment the following lines for automatically selecting the two biggest clusters of patients if more than 2 clusters were given
            done1 = "false"
            for column_name, column in exprdf.transpose().iterrows():
                if (not column_name.isdigit()):
                    if (len(column.unique()) < 6):
                        nbr_col = len(column.unique())
                        expr_str = expr_str.replace(column_name, "disease_type")
                        expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
                        if (len(column.unique()) > 2 and done1 == "false"):
                            expr_str_split_2 = []
                            expr_str_split_2.append(expr_str_split[0])
                            type1 = column.value_counts().index.tolist()[0]
                            type2 = column.value_counts().index.tolist()[1]
                            for i in range(0, len(list(column)) - 1):
                                if (list(column)[i] == type1 or list(column)[i] == type2):
                                    expr_str_split_2.append(expr_str_split[i + 1])
                            expr_str = "\n".join(expr_str_split_2)
                            done1 = "true"
            ########################

            expr_stringio = StringIO(expr_str)
            exprdf = pd.read_csv(expr_stringio, sep='\t')
            return (expr_str)


# the same as preprocess_file, but returns number of pre-defined clusters
@shared_task(name="preprocess_file_2")
def preprocess_file_2(expr_str):
    expr_str = expr_str.replace("cancer_type", "disease_type")
    nbr_col = 1
    if (len(expr_str.split("\n")[0].split("\t")) > 2):
        expr_str_split = expr_str.split("\n")
        # replace column name for disease type
        if ("disease_type" not in expr_str_split[0]):
            if ("subtype" in expr_str_split[0]):
                expr_str = expr_str.replace("subtype", "disease_type")
        # remove name of first column (left upper corner)
        expr_str_first_colname = expr_str_split[0].split("\t")[0]
        expr_str = expr_str.replace(expr_str_first_colname, "", 1)
        expr_stringio = StringIO(expr_str)
        exprdf = pd.read_csv(expr_stringio, sep='\t')
        done1 = "false"
        # check for column with two unique entries (pre-defined clusters)
        for column_name, column in exprdf.transpose().iterrows():
            if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
                if (len(column.unique()) == 2):
                    expr_str = expr_str.replace(column_name, "disease_type")
                    nbr_col = 2
                    done1 = "true"
        # check for column with less than 6 unique entries (pre-defined clusters)
        for column_name, column in exprdf.transpose().iterrows():
            if (not column_name.isdigit()):
                if (len(column.unique()) < 6):
                    nbr_col = len(column.unique())
                    expr_str = expr_str.replace(column_name, "disease_type")
                    expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
                    if (len(column.unique()) > 2 and done1 == "false"):
                        expr_str_split_2 = []
                        expr_str_split_2.append(expr_str_split[0])
                        for i in range(0, len(list(column)) - 1):
                            expr_str_split_2.append(expr_str_split[i + 1])
                        expr_str = "\n".join(expr_str_split_2)
                        done1 = "true"
        return (expr_str, nbr_col)
    elif ("," in expr_str):
        # replace comma by tab if file is CSV and not TSV
        if ("\t" not in expr_str.split("\n")[0]):
            expr_str_split = expr_str.split("\n")
            # replace "subtype" by "cancer type"
            if ("disease_type" not in expr_str):
                if ("subtype" in expr_str):
                    expr_str = expr_str.replace("subtype", "disease_type")
            expr_str_first_colname = expr_str_split[0].split(",")[0]
            expr_str = expr_str.replace(expr_str_first_colname, "", 1)
            expr_str = expr_str.replace(",", "\t")
            expr_str_split = expr_str.split("\n")
            # remove entries after given length if expression data file is too big
            if (len(expr_str_split) > 300):
                expr_str = "\n".join(expr_str_split[:200])
            else:
                expr_str = "\n".join(expr_str_split)
            expr_stringio = StringIO(expr_str)
            expr_str = expr_str.replace("MCI", "CTL")
            exprdf = pd.read_csv(expr_stringio, sep='\t')
            # find column with two unique entries that represents disease type
            for column_name, column in exprdf.transpose().iterrows():
                if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
                    if (len(column.unique()) == 2):
                        # print(column_name)
                        expr_str = expr_str.replace(column_name, "disease_type")
                        nbr_col = 2
            #### uncomment the following lines for automatically selecting the two biggest clusters of patients if more than 2 clusters were given
            done1 = "false"
            for column_name, column in exprdf.transpose().iterrows():
                if (not column_name.isdigit()):
                    if (len(column.unique()) < 6):
                        nbr_col = len(column.unique())
                        expr_str = expr_str.replace(column_name, "disease_type")
                        expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
                        if (len(column.unique()) > 2 and done1 == "false"):
                            expr_str_split_2 = []
                            expr_str_split_2.append(expr_str_split[0])
                            for i in range(0, len(list(column)) - 1):
                                # if(list(column)[i] == type1 or list(column)[i] == type2):
                                expr_str_split_2.append(expr_str_split[i + 1])
                            expr_str = "\n".join(expr_str_split_2)
                            done1 = "true"
            ########################

            expr_stringio = StringIO(expr_str)
            exprdf = pd.read_csv(expr_stringio, sep='\t')
            # return(expr_str,nbr_col)
            return (expr_str, nbr_col)


@shared_task(name="add_loading_image")
def add_loading_image(session_id):
    if (session_id == "none"):
        if (os.path.isfile("/code/clustering/static/loading.gif")):
            copyfile("/code/clustering/static/loading.gif", "/code/clustering/static/loading_1.gif")
        else:
            print("loading image not found")
    else:
        if (os.path.isfile("/code/clustering/static/loading.gif")):
            copyfile("/code/clustering/static/loading.gif",
                     "/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")
        else:
            print("loading image not found")


@shared_task(name="remove_loading_image")
def remove_loading_image(session_id):
    if (session_id == "none"):
        if (os.path.isfile("/code/clustering/static/loading_1.gif")):
            os.unlink("/code/clustering/static/loading_1.gif")
    else:
        if (os.path.isfile("/code/clustering/static/loading_1_" + session_id + ".gif")):
            os.unlink("/code/clustering/static/loading_1_" + session_id + ".gif")
        if (os.path.isfile("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")):
            os.unlink("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")


##################################################################
################## used for metadata display #####################
##################################################################


# read metadata from file and return array with data
@shared_task(name="list_metadata_from_file")
def list_metadata_from_file(path):
    # used for reading metadata
    fh1 = open(path)
    lines = fh1.read()
    if (lines == "NA"):
        return ({}, {}, {})
    # remove html from metadata file and replace table elements by tab
    # if no data in file, remove empty dictionaries
    if (len(lines.split('\n')) < 3):
        return ({}, {}, {})
    # read content from lines
    line0 = lines.split('\n')[0].split('\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    ret = []
    dict3 = {}
    dict1 = {}
    dict2 = {}
    dict0 = {}
    ctr = 0
    dict3['params'] = line0
    dict3['gr1'] = line1
    dict3['gr2'] = line2
    dict3['all'] = zip(dict3['params'], dict3['gr1'], dict3['gr2'])
    # dict 0 is parameter names, dict1 is values for group 1, dict2 is values for group 2
    for i in range(0, len(line0) - 1):
        dict0[i] = line0[i]
        dict1[dict0[i]] = line1[i]
        dict2[dict0[i]] = line2[i]
        ctr = ctr + 1
    return (dict0, dict1, dict2)


##################################################################################################
######### running the algorithm - part 1 #########################################################
##################################################################################################


## method for more than 2 pre-defined clusters, uses session IDs
@shared_task(name="algo_output_task")
def algo_output_task(s, L_g_min, L_g_max, expr_str, ppi_str, nbr_iter, nbr_ants, evap, epsilon, hi_sig, pher_sig,
                     session_id, size, clusters_param):
    col = "disease_type"
    log2 = True
    expr_stringio = StringIO(expr_str)
    exprdf = pd.read_csv(expr_stringio, sep='\t')
    # check if string contains negative numbers.
    if ("-" in expr_str.split("\n")[2]):
        print("expression data are logarithmized")
        log2_2 = False
    else:
        print("expression data not logarithmized")
        log2_2 = True
    ### this checks whether the expression data contain negative numbers
    for i in range(2, 4):
        if (log2_2 and i > len(exprdf.columns)):
            # check only first 1000 lines of column 2 and 3
            for j in range(1, min(len(exprdf.index) - 1, 1000)):
                if (log2_2 and str(exprdf.columns[i]) != "disease_type"):
                    # make integer from negative number (e.g. -1.0 -> 10), check if it is a number and check if number is negative
                    if (exprdf.iloc[[j], [i]].to_string().__contains__('-') and str(exprdf.iloc[j][i]).replace("-", "",
                                                                                                               1).replace(
                            ".", "", 1).isdigit()):
                        print("expression data are logarithmized")
                        log2_2 = False
    if (session_id == "none"):
        with open(("/code/clustering/static/userfiles/output_console.txt"), "w") as text_file:
            text_file.write("Your files are being processed...")
        with open(("/code/clustering/static/output_console.txt"), "w") as text_file:
            text_file.write("Your files are being processed...")
    else:
        with open(("/code/clustering/static/userfiles/output_console_" + session_id + ".txt"), "w") as text_file:
            text_file.write("Your files are being processed...")
        with open(("/code/clustering/static/output_console.txt"), "w") as text_file:
            text_file.write("Your files are being processed...")
    # B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
    if (clusters_param == 2):
        B, G, H, n, m, GE, A_g, group1, group2, labels_B, rev_labels_B, val1, val2, group1_ids, group2_ids = lib.aco_preprocessing_strings(
            expr_str, ppi_str, col, log2=log2_2, gene_list=None, size=int(size), sample=None)
    else:
        B, G, H, n, m, GE, A_g, labels_B, rev_labels_B = lib.aco_preprocessing_strings_2(expr_str, ppi_str, col,
                                                                                         log2=log2_2, gene_list=None,
                                                                                         size=size, sample=None)
    if (session_id == "none"):
        with open(("/code/clustering/static/userfiles/output_console.txt"), "w") as text_file:
            text_file.write("Starting model run...")
    else:
        with open(("/code/clustering/static/userfiles/output_console_" + session_id + ".txt"), "w") as text_file:
            text_file.write("Starting model run...")
    print("How many genes you want per cluster (minimum):")
    # L_g_min = int(input())
    print("How many genes you want per cluster (maximum):")
    # L_g_max = int(input())
    imp.reload(lib)

    # =============================================================================
    # #GENERAL PARAMETERS:
    # =============================================================================
    clusters = clusters_param  # other options are currently unavailable
    K = int(nbr_ants)  # number of ants
    eps = float(epsilon)  # stopping criteria: score_max-score_av<eps
    b = float(hi_sig)  # HI significance
    evaporation = float(evap)
    a = float(pher_sig)  # pheramone significance
    times = int(nbr_iter)  # max amount of iterations
    # =============================================================================
    # #NETWORK SIZE PARAMETERS:
    # =============================================================================
    cost_limit = 20  # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3

    th = 1  # the coefficient to define the search radipus which is supposed to be bigger than
    # mean(heruistic_information[patient]+th*std(heruistic_information[patient])
    # bigger th - less genes are considered (can lead to empty paths if th is too high)
    # will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000

    if (session_id == "none"):
        with open(("/code/clustering/static/output_console.txt"), "w") as text_file:
            text_file.write("Progress of the algorithm is shown below...")
    else:
        with open(("/code/clustering/static/userfiles/output_console_" + session_id + ".txt"), "w") as text_file:
            text_file.write("Progress of the algorithm is shown below...")
    start = time.time()
    # session id is "none" if it is not given
    solution, t_best, sc, conv = lib.ants_new(a, b, n, m, H, GE, G, 2, cost_limit, K, evaporation, th, L_g_min, L_g_max,
                                              eps, times, session_id, opt=None, pts=False, show_pher=False,
                                              show_plot=True, print_runs=False, save=None, show_nets=False)
    end = time.time()
    n_proc = os.getenv("NBR_PROCESSES", '4')
    lib.ants_manager(a, b, n, m, H, GE, G, 2, cost_limit, K, evaporation, th, L_g_min, L_g_max, eps, times, session_id,
                     n_proc, opt=None, pts=False, show_pher=True, show_plot=True, save=None, show_nets=False)
    print("######################################################################")
    print("RESULTS ANALYSIS")
    print("total time " + str(round((end - start) / 60, 2)) + " minutes")
    print("jaccard indexes:")
    jac_1_ret = ""
    jac_2_ret = ""
    if (clusters_param == 2):
        jacindices = lib.jac_matrix(solution[1], [group1, group2])
        print(jacindices)
        jac_1_ret = jacindices[0]
        jac_2_ret = jacindices[1]
        if lib.jac(group1, solution[1][0]) > lib.jac(group1, solution[1][1]):
            values = [val1, val2]
        else:
            values = [val2, val1]
    # mapping to gene names (for now with API)
    mg = mygene.MyGeneInfo()
    new_genes = solution[0][0] + solution[0][1]
    new_genes_entrez = [labels_B[x] for x in new_genes]
    out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
    mapping = dict()
    for line in out:
        if ("symbol" in line):
            mapping[rev_labels_B[line["query"]]] = line["symbol"]
    ###m plotting networks
    new_genes1 = [mapping[key] for key in mapping if key in solution[0][0]]
    new_genes2 = [mapping[key] for key in mapping if key in solution[0][1]]

    genes1, genes2 = solution[0]
    patients1, patients2 = solution[1]
    patients1_ids = []
    patients2_ids = []
    for elem in patients1:
        if (elem in labels_B):
            patients1_ids.append(labels_B[elem])
    for elem in patients2:
        if (elem in labels_B):
            patients2_ids.append(labels_B[elem])
    means1 = [np.mean(GE[patients1].loc[gene]) - np.mean(GE[patients2].loc[gene]) for gene in genes1]
    means2 = [np.mean(GE[patients1].loc[gene]) - np.mean(GE[patients2].loc[gene]) for gene in genes2]

    G_small = nx.subgraph(G, genes1 + genes2)
    G_small = nx.relabel_nodes(G_small, mapping)

    plt.figure(figsize=(15, 15))
    cmap = plt.cm.RdYlGn
    vmin = -2
    vmax = 2
    pos = nx.spring_layout(G_small)
    ec = nx.draw_networkx_edges(G_small, pos)
    if (clusters_param == 2):
        nc1 = nx.draw_networkx_nodes(G_small, nodelist=new_genes1, pos=pos, node_color=means1, node_size=600, alpha=1.0,
                                     vmin=vmin, vmax=vmax, node_shape="^", cmap=cmap, label=values[0])
        nc2 = nx.draw_networkx_nodes(G_small, nodelist=new_genes2, pos=pos, node_color=means2, node_size=600, alpha=1.0,
                                     vmin=vmin, vmax=vmax, node_shape="o", cmap=cmap, label=values[1])
    else:
        nc1 = nx.draw_networkx_nodes(G_small, nodelist=new_genes1, pos=pos, node_color=means1, node_size=600, alpha=1.0,
                                     vmin=vmin, vmax=vmax, node_shape="^", cmap=cmap)
        nc2 = nx.draw_networkx_nodes(G_small, nodelist=new_genes2, pos=pos, node_color=means2, node_size=600, alpha=1.0,
                                     vmin=vmin, vmax=vmax, node_shape="o", cmap=cmap)
    nx.draw_networkx_labels(G_small, pos, font_size=15, font_weight='bold')
    ret2 = means1 + means2
    ret3 = new_genes1 + new_genes2
    adjlist = []
    for line in nx.generate_edgelist(G_small, data=False):
        lineSplit = line.split()
        adjlist.append([lineSplit[0], lineSplit[1]])

    plt.legend(frameon=True)
    try:
        plt.colorbar(nc1)
    except:
        print("no colorbar found")
    plt.axis('off')
    ### plotting expression data
    plt.rc('font', size=30)  # controls default text sizes
    plt.rc('axes', titlesize=20)  # fontsize of the axes title
    plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=15)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=10)  # fontsize of the tick labels
    plt.rc('legend', fontsize=30)

    grouping_p = []
    grouping_g = []
    p_num = list(GE.columns)
    GE_small = GE.T[genes1 + genes2]
    GE_small.rename(columns=mapping, inplace=True)
    GE_small = GE_small.T
    g_num = list(GE_small.index)
    if (clusters_param == 2):
        for g in g_num:
            if g in new_genes1:
                grouping_g.append(values[0])
            elif g in new_genes2:
                grouping_g.append(values[1])
            else:
                grouping_g.append(3)
        for p in p_num:
            if p in solution[1][0]:
                grouping_p.append(values[0])
            else:
                grouping_p.append(values[1])
        grouping_p = pd.DataFrame(grouping_p, index=p_num)
        grouping_g = pd.DataFrame(grouping_g, index=g_num)
        species = grouping_g[grouping_g[0] != 3][0]
        lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
        col_colors = species.map(lut)
        species = grouping_p[0]
        lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
        row_colors = species.map(lut)
    else:
        col_colors = ""
        row_colors = ""
        for g in g_num:
            if g in new_genes1:
                grouping_g.append("cluster1")
            elif g in new_genes2:
                grouping_g.append("cluster2")
            else:
                grouping_g.append(3)
        grouping_g = pd.DataFrame(grouping_g, index=g_num)
        species = grouping_g[grouping_g[0] != 3][0]
        lut = {"cluster1": '#4FB6D3', "cluster2": '#22863E'}
        col_colors = species.map(lut)
    if (session_id == "none"):
        plt.savefig("/code/clustering/static/userfiles/ntw.png")
        plt.savefig("/code/clustering/static/ntw.png")
    else:
        plt.savefig("/code/clustering/static/userfiles/ntw_" + session_id + ".png")
    plt.clf()
    plt.boxplot(conv / 2, vert=True, patch_artist=True)  # vertical box alignment  # will be used to label x-ticks
    plt.xlabel("iterations")
    plt.ylabel("score per subnetwork")
    if (session_id == "none"):
        plt.savefig("/code/clustering/static/userfiles/conv.png")
        plt.savefig("/code/clustering/static/conv.png")
    else:
        plt.savefig("/code/clustering/static/userfiles/conv_" + session_id + ".png")
    return (GE_small.T, row_colors, col_colors, G_small, ret2, ret3, adjlist, new_genes1, patients1_ids, patients2_ids,
            jac_1_ret, jac_2_ret)


##########################################################
#### running the algorithm - part 2 ######################
##########################################################


### Processing of algorithm output with using session ID
@shared_task(name="script_output_task")
def script_output_task(T, row_colors1, col_colors1, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids,
                       clinicalstr, jac_1_ret, jac_2_ret, survival_col, clinicaldf, session_id):
    # define colors depending on z-score differences of genes in graph
    def color_for_graph(v):
        cmap_custom = {-4: 'rgb(255, 0, 0)', -3: 'rgb(255, 153, 51)', -2: 'rgb(255, 204, 0)', -1: 'rgb(255, 255, 0)',
                       0: 'rgb(204, 255, 51)', 1: 'rgb(153, 255, 51)', 2: 'rgb(102, 255, 51)', 3: 'rgb(51, 204, 51)'}
        v = v * 2
        v_int = int(v)
        if (v < -4):
            v_int = -4
        if (v > 3):
            v_int = 3
        return (cmap_custom[v_int])

    nodecolors = []
    genes = {}
    G_list = list(G2.nodes())
    ctr = 0
    G = nx.Graph()
    # make node objects for genes
    for G_tmp in genes_all:
        genes.update({G_tmp: 0})
        # circle/square nodes based on cluster of genes
        tp = "circle"
        if (G_tmp in genes1):
            tp = "square"
        # create node objects with color property based on z-score difference
        G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]), color=color_for_graph(means[ctr]), type=tp, label=G_tmp)
        nodecolors.append(color_for_graph(means[ctr]))
        ctr = ctr + 1
    ctr = 0
    # make edge objects for PPI
    for edg in adjlist:
        G.add_edge(edg[0], edg[1], id=ctr, color="rgb(0,0,0)")
        ctr = ctr + 1
    pos = nx.spring_layout(G)
    x_pos = {}
    y_pos = {}
    # take y and x positions from the given networkx layout
    for k in pos:
        x_pos[k] = pos[k][0]
        y_pos[k] = pos[k][1]
    # write json object of genes and interactions
    nx.set_node_attributes(G, x_pos, 'x')
    nx.set_node_attributes(G, x_pos, 'x')
    nx.set_node_attributes(G, y_pos, 'y')
    nx.set_node_attributes(G, 10, 'size')
    # replace some json object names by correct form
    jsn = json_graph.node_link_data(G)
    jsn2 = str(json.dumps(jsn))
    jsn33 = jsn2.replace('links', 'edges')
    jsn44 = jsn33.replace('Label', 'label')
    jsn55 = jsn44.replace('bels', 'bel')
    jsn3 = jsn55.replace('\"directed\": false, \"multigraph\": false, \"graph\": {},', '')
    if (session_id == "none"):
        json_path = "/code/clustering/static/userfiles/ppi.json"
    else:
        json_path = "/code/clustering/static/userfiles/ppi_" + session_id + ".json"
    with open(json_path, "w") as text_file:
        text_file.write(jsn3)
    colordict = {0: '#BB0000', 1: '#0000BB'}
    # make heatmap (include pre-defined clusters if they were given)
    if (isinstance(col_colors1, str)):
        g = sns.clustermap(T, figsize=(13, 13))
    else:
        if (isinstance(row_colors1, str)):
            g = sns.clustermap(T, figsize=(13, 13), col_colors=col_colors1)
        else:
            g = sns.clustermap(T, figsize=(13, 13), col_colors=col_colors1, row_colors=row_colors1)
    ax = g.ax_heatmap
    ax.set_xlabel("Genes")
    ax.set_ylabel("Patients")
    # do not include session ID in path if no session ID was given
    if (session_id == "none"):
        path_heatmap = "/code/clustering/static/userfiles/heatmap.png"
    else:
        path_heatmap = "/code/clustering/static/userfiles/heatmap_" + session_id + ".png"
    plt.savefig(path_heatmap)
    plt.clf()
    # Array PatientData is for storing survival information
    patientData = {}
    # write lists of genes in files, needed for enrichment analysis
    if (session_id == "none"):
        path_genelist = "/code/clustering/static/userfiles/genelist.txt"
        path_genelist_1 = "/code/clustering/static/userfiles/genelist_1.txt"
        path_genelist_2 = "/code/clustering/static/userfiles/genelist_2.txt"
    else:

        path_genelist = "/code/clustering/static/userfiles/genelist_" + session_id + ".txt"
        path_genelist_1 = "/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"
        path_genelist_2 = "/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"
    with open(path_genelist, "w") as text_file_4:
        for i in G_list:
            text_file_4.write(str(i) + "\n")
    text_file_4.close()
    with open(path_genelist_1, "w") as text_file_5:
        for i in G_list:
            if (i in genes1):
                text_file_5.write(str(i) + "\n")
    text_file_5.close()
    with open(path_genelist_2, "w") as text_file_6:
        for i in G_list:
            if (i not in genes1):
                text_file_6.write(str(i) + "\n")
    text_file_6.close()
    if (session_id == "none"):
        path_metadata = "/code/clustering/static/userfiles/metadata.txt"
    else:
        path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
    # if no metadata given, write an empty metadata file
    p_val = ""
    if (clinicalstr == "empty"):
        print("no clinical data")
        ret_metadata = []
        text_file_3 = open(path_metadata, "w")
        text_file_3.write("NA")
        text_file_3.close()
        plot_div = ""
        if (session_id == "none"):
            output_plot_path = "/code/clustering/static/userfiles/output_plotly.html"
        else:
            output_plot_path = "/code/clustering/static/userfiles/output_plotly_" + session_id + ".html"
        with open(output_plot_path, "w") as text_file_2:
            text_file_2.write("")
        output_plot_path = "empty"
        # fill empty metadata arrays
        ret_metadata = []
        ret_metadata_1 = {}
        ret_metadata_2 = {}
        ret_metadata_3 = {}
        ret_metadata.append(ret_metadata_1)
        ret_metadata.append(ret_metadata_2)
        ret_metadata.append(ret_metadata_3)
        patientids_metadata = []
    else:
        # read clinical data line by line, read title column in separate array
        clinicalLines = clinicalstr.split("\n")
        title_col = clinicalLines[0].split(",")
        # assign some default value to survival col nbr (for the case that no survival column exists)
        survival_col_nbr = 64
        # replace all type of NA in dataframe by standard pandas-NA
        clinicaldf.replace(['NaN', 'nan', '?', '--'], ['NA', 'NA', 'NA', 'NA'], inplace=True)
        clinicaldf.replace(['NTL'], ['NA'], inplace=True)
        clinicaldf.replace(['na'], ['NA'], inplace=True)
        clinicaldf_col_names = list(clinicaldf.columns)
        patientids_metadata = [str(i) for i in clinicaldf.iloc[:, 0].values.tolist()]
        # get patient ids either from first column or from index, add one empty element at the beginning of the column names
        if ("Unnamed" in "\t".join(list(clinicaldf.columns))):
            clinicaldf_col_names_temp = ['empty']
            clinicaldf_col_names_new = clinicaldf_col_names_temp + clinicaldf_col_names
            clinicaldf.columns = list(clinicaldf_col_names_new[:-1])
        if ("GSM" not in patientids_metadata[0]):
            patientids_metadata = list(clinicaldf.index)
        param_names = []
        param_values = []
        param_cols = []
    # if clinical data were uploaded, more than 1 patient exists and patient IDs from metadata and expression data overlap
    if not (clinicalstr == "empty" or ((len(group1_ids) + len(group2_ids)) < 1) or set(patientids_metadata).isdisjoint(
            group1_ids)):
        patients_0 = []
        patients_1 = []
        group1_has = []
        group2_has = []
        # iterate over columns of metadata, get all unique entries
        for column_name, column in clinicaldf.transpose().iterrows():
            column.fillna("NA", inplace=True)
            coluniq = column.unique()
            # replace all instances of NA by pandas-standard NA
            for elem in coluniq:
                if ": " in str(elem):
                    elem = elem.split(": ")[1]
            for elem in column:
                if ": " in str(elem):
                    elem = elem.split(": ")[1]
            for elem in column:
                if (str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
                    elem = "NA"
                elif (elem == float('nan')):
                    elem = "NA"
            for elem in coluniq:
                if (str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
                    elem = "NA"
                elif (elem == float('nan')):
                    elem = "NA"
            if (len(coluniq) == 2 or len(coluniq) == 3):
                patients_temp_0 = []
                patients_temp_1 = []
                # replace simple binary entries like good and bad prognosis by standard 0 and 1 for calculation
                column = column.replace('Good Prognosis', '1')
                column = column.replace('Bad Prognosis', '0')
                column = column.replace('yes', '1')
                column = column.replace('no', '0')
                column = column.replace('P', '1')
                column = column.replace('N', '0')
                column = column.replace('0A', 'NA')
                column = column.replace('relapse (event=1; no event=0): 0', '0')
                column = column.replace('relapse (event=1; no event=0): 1', '0')
                column = column.replace('relapse (event=1; no event=0): na', 'NA')
                column = column.replace('status: ALIVE', '1')
                column = column.replace('status: DEAD', '0')
                column = column.replace('status: NTL', 'NA')
                column = column.replace('ALIVE', '1')
                column = column.replace('DEAD', '0')
                column = column.replace('NTL', 'NA')
                if ("gender" in column_name):
                    column = column.replace('M', '1')
                    column = column.replace('F', '0')
                    column_name = "Gender: Male"
                coluniq2 = column.unique()
                coluniq3 = [str(w) for w in coluniq2]
                # check if sorted column now contains only 0,1 and NA
                if (sorted(coluniq3) == ['0', '1', 'NA'] or sorted(coluniq3) == ['0', '1']):
                    # get column values as array
                    col_as_list = [str(i) for i in column]
                    # append values to patient list
                    for i in range(0, len(col_as_list) - 1):
                        if (col_as_list[i] == '0'):
                            patients_temp_0.append(patientids_metadata[i])
                        elif (col_as_list[i] == '1'):
                            patients_temp_1.append(patientids_metadata[i])
                    # append patient list for metadata variable to overall patient list
                    patients_0.append(patients_temp_0)
                    patients_1.append(patients_temp_1)
                    # add column name to parameter names
                    param_names.append(column_name)
                    param_cols.append(ctr)
                    all_patients = patients_temp_0 + patients_temp_1
                    current_patients_group_1 = []
                    current_patients_group_2 = []
                    # check which patients in both clusters are represented in current column
                    for i in range(0, len(all_patients) - 1):
                        if (all_patients[i] in group1_ids):
                            current_patients_group_1.append(all_patients[i])
                        elif (all_patients[i] in group2_ids):
                            current_patients_group_2.append(all_patients[i])
                    # append to array that lists the available patients for all variables
                    group1_has.append(current_patients_group_1)
                    group2_has.append(current_patients_group_2)
                elif (":" in coluniq[0]):
                    patients_temp_0 = []
                    patients_temp_1 = []
                    coluniq_split = []
                    for elem in coluniq:
                        if (":" in elem):
                            coluniq_split.append(elem.split(":")[1].replace(" ", ""))
                    if (sorted(coluniq_split) == ['0', '1', 'NA'] or sorted(coluniq_split) == ['0', '1']):
                        col_as_list_tmp = [str(i) for i in column]
                        col_as_list = [i.split(":")[1].replace(" ", "") for i in col_as_list_tmp]
                        # append values to patient list
                        for i in range(0, len(col_as_list) - 1):
                            if (col_as_list[i] == '0'):
                                patients_temp_0.append(patientids_metadata[i])
                            elif (col_as_list[i] == '1'):
                                patients_temp_1.append(patientids_metadata[i])
                        # append patient list for metadata variable to overall patient list
                        patients_0.append(patients_temp_0)
                        patients_1.append(patients_temp_1)
                        # add column name to parameter names
                        param_names.append(coluniq[0].split(":")[0])
                        param_cols.append(ctr)
                        all_patients = patients_temp_0 + patients_temp_1
                        current_patients_group_1 = []
                        current_patients_group_2 = []
                        # check which patients in both clusters are represented in current column
                        for i in range(0, len(all_patients) - 1):
                            if (all_patients[i] in group1_ids):
                                current_patients_group_1.append(all_patients[i])
                            elif (all_patients[i] in group2_ids):
                                current_patients_group_2.append(all_patients[i])
                        # append to array that lists the available patients for all variables
                        group1_has.append(current_patients_group_1)
                        group2_has.append(current_patients_group_2)
            ctr = ctr + 1
        jaccards_1 = []
        jaccards_2 = []
        param_names_final = []
        nbr_patients = float(len(group1_ids) + len(group2_ids))
        # calculate fractions of patients for which metadata variables are 0 or 1
        for i in range(0, len(param_names) - 1):
            if ((float(len(patients_0[i]) + len(patients_1[i])) / nbr_patients) > 0.8):
                if not (lib.jac(group1_has[i], patients_0[i]) == 0.0 and lib.jac(group2_has[i], patients_1[i]) == 0.0):
                    param_names_final.append(param_names[i])
                    jaccards_1.append(lib.jac(group1_has[i], patients_0[i]))
                    jaccards_2.append(lib.jac(group2_has[i], patients_1[i]))
        # get list of patient ids
        if ("GSM" not in list(clinicaldf.iloc[:, 0])[1]):
            patient_id_list = list(clinicaldf.index)
        else:
            patient_id_list = list(clinicaldf.iloc[:, 0])
        # check if there is a column with survival data
        if (survival_col in list(clinicaldf.columns)):
            survival_col_nbr = list(clinicaldf.columns).index(survival_col)
            print("column with survival data found")
            clinicaldf.iloc[:, survival_col_nbr].fillna("NA", inplace=True)
            survivalcol_list = list(clinicaldf.iloc[:, survival_col_nbr])
            clinicaldf.iloc[:, survival_col_nbr].fillna("NA", inplace=True)
            survivalcol_list = list(clinicaldf.iloc[:, survival_col_nbr])
        # give empty survival lists if no data given
        else:
            survivalcol_list = []
            patient_id_list = []
        # replace NA by standard NA for all entries in survival column

        # iterate over patient IDs
        for i in range(0, len(patient_id_list)):
            # check if survival column contains number. divide by 12 if it is given in months
            if (survivalcol_list[i] != "--" and survivalcol_list[i] != "name:ch1" and survivalcol_list[i] != "NA" and
                    survivalcol_list[i].replace('.', '', 1).isdigit()):
                if ("month" in survival_col or "MONTH" in survival_col):
                    survivalcol_list_temp = float(survivalcol_list[i]) / 12.0
                    patientData.update({patient_id_list[i]: survivalcol_list_temp})
                else:
                    patientData.update({patient_id_list[i]: survivalcol_list[i]})
        ret_metadata = []
        survival_1 = []
        survival_2 = []
        ctr_surv_1 = 0.001
        sum_surv_1 = 0
        ctr_surv_2 = 0.001
        sum_surv_2 = 0
        # make arrays with survival time of patients in both groups
        for key in patientData:
            if key in group1_ids:
                survival_1.append(float(patientData[key]))
            elif key in group2_ids:
                survival_2.append(float(patientData[key]))
        # calculate p-value for survival times
        if (survival_col in list(clinicaldf.columns) and len(survival_1) > 0 and len(survival_2) > 0):
            surv_results = logrank_test(survival_1, survival_2)
            p_val = surv_results.p_value
        else:
            p_val = ""
        # count survival times in both arrays
        for elem in survival_1:
            sum_surv_1 = sum_surv_1 + float(elem)
            ctr_surv_1 = ctr_surv_1 + 1
        for elem in survival_2:
            sum_surv_2 = sum_surv_2 + float(elem)
            ctr_surv_2 = ctr_surv_2 + 1
        # replace some abbreviated clinical terms by proper description
        param_names = [elem.replace("bm event:ch1", "Breast Metastasis") for elem in param_names]
        param_names = [elem.replace("lm event:ch1", "Lung Metastasis") for elem in param_names]
        param_names = [elem.replace("met event:ch1", "Metastasis") for elem in param_names]
        param_names = [elem.replace("relapse (event=1; no event=0):ch1", "Relapse") for elem in param_names]
        errstr = ""
        if (len(survival_1) == 0):
            errstr = "Unfortunately, no survival data could be computed."
        text_file_4 = open((path_metadata), "w")
        text_file_4.write("\t".join(param_names) + "\n")
        jaccards_1_str = [str(i)[:4] for i in jaccards_1]
        jaccards_2_str = [str(i)[:4] for i in jaccards_2]
        text_file_4.write("\t".join(jaccards_1_str) + "\n")
        text_file_4.write("\t".join(jaccards_2_str) + "\n")
        text_file_4.write("")
        text_file_4.close()
        # write metadata to dicts
        ret_metadata = []
        ret_metadata_1 = {}
        ret_metadata_2 = {}
        ret_metadata_3 = {}
        # for i in range(0, len(param_names)):
        for i in range(0, len(jaccards_1)):
            ret_metadata_1[i] = param_names[i]
            if (len(str(jaccards_1[i])) > 4):
                ret_metadata_2[i] = str(jaccards_1[i])[:4]
            else:
                ret_metadata_2[i] = jaccards_1[i]
            if (len(str(jaccards_2[i])) > 4):
                ret_metadata_3[i] = str(jaccards_2[i])[:4]
            else:
                ret_metadata_3[i] = jaccards_2[i]
        ret_metadata.append(ret_metadata_1)
        ret_metadata.append(ret_metadata_2)
        ret_metadata.append(ret_metadata_3)
        survival_perc_1 = {0: 1}
        survival_perc_2 = {0: 1}
        # calculate data for kaplan meyer plot
        # check for every number of years between 1 and 10, how many patients are alive
        for i in range(1, 10):
            tmp1 = 1.0
            tmp2 = 1.0
            # iterate over list with survival times
            for k in survival_1:
                if (float(k) < float(i)):
                    tmp1 = tmp1 - (1.0 / len(survival_1))
            for k in survival_2:
                if (float(k) < float(i)):
                    tmp2 = tmp2 - (1.0 / len(survival_2))
            survival_perc_1.update({i: tmp1})
            survival_perc_2.update({i: tmp2})
        # write kaplan meyer plots
        trace1 = go.Scatter(
            x=list(survival_perc_1.keys()),
            y=list(survival_perc_1.values()),
            mode='lines+markers',
            name="'Group 1'",
            hoverinfo='name',
            line=dict(
                shape='hv'))
        trace2 = go.Scatter(
            x=list(survival_perc_2.keys()),
            y=list(survival_perc_2.values()),
            mode='lines+markers',
            name="'Group 2'",
            hoverinfo='name',
            line=dict(
                shape='hv'))
        surv_data_for_graph = [trace1, trace2]
        layout = dict(showlegend=False,
                      xaxis=dict(
                          title='Time in years'),
                      yaxis=dict(
                          title='percentage of patients'))
        fig = dict(data=surv_data_for_graph, layout=layout)
        plot_div = plotly.offline.plot(fig, auto_open=False, output_type='div')
        if (session_id == "none"):
            output_plot_path = "/code/clustering/static/userfiles/output_plotly.html"
        else:
            output_plot_path = "/code/clustering/static/userfiles/output_plotly_" + session_id + ".html"
        if (survival_col not in list(clinicaldf.columns) or (len(survival_1) == 0)):
            with open(output_plot_path, "w") as text_file_2:
                text_file_2.write("")
            output_plot_path = "empty"
        elif (errstr == ""):
            with open(output_plot_path, "w") as text_file_2:
                text_file_2.write(plot_div)
        else:
            with open(output_plot_path, "w") as text_file_2:
                text_file_2.write(errstr)
            output_plot_path = "empty"
    else:
        # return empty arrays for metadata / empty file for survival plot  if no metadata were found
        ret_metadata = []
        text_file_3 = open(path_metadata, "w")
        text_file_3.write("NA")
        text_file_3.close()
        plot_div = ""
        if (session_id == "none"):
            output_plot_path = "/code/clustering/static/userfiles/output_plotly.html"
        else:
            output_plot_path = "/code/clustering/static/userfiles/output_plotly_" + session_id + ".html"
        with open(output_plot_path, "w") as text_file_2:
            text_file_2.write("")
        output_plot_path = "empty"
        # fill empty metadata arrays
        ret_metadata = []
        ret_metadata_1 = {}
        ret_metadata_2 = {}
        ret_metadata_3 = {}
        ret_metadata.append(ret_metadata_1)
        ret_metadata.append(ret_metadata_2)
        ret_metadata.append(ret_metadata_3)
    if (session_id == "none"):
        with open("/code/clustering/static/output_console.txt", "w") as text_file:
            text_file.write("Done!")
    else:
        with open("/code/clustering/static/userfiles/output_console_" + session_id + ".txt", "w") as text_file:
            text_file.write("Done!")
    return (ret_metadata, path_heatmap, path_metadata, output_plot_path, json_path, p_val)


## enrichment stuff ##
# run enrichment analysis
@shared_task(name="run_enrichment")
def run_enrichment(path, pval_enr, out_dir, terms):
    fh1 = open(path)
    gene_list = []
    lines = fh1.readlines()
    # read gene list line for line from file
    for line in lines:
        line.replace("\\n", "")
        gene_list.append(line)
    print("running enrichment analysis")
    enr = gp.enrichr(gene_list=gene_list,
                     description='test_name',
                     gene_sets=terms,
                     outdir=out_dir,
                     cutoff=float(pval_enr)  # test dataset, use lower value of range(0,1)
                     )
    return (enr.results)


# read terms in given cluster
@shared_task(name="read_enrichment")
def read_enrichment(path, pval_enr):
    result_file = open(path)
    ret_dict = []
    ctr = 0
    for line in result_file:
        tmp = {}
        lineSplit = line.split("\t")
        if (ctr > 0):
            # check p-value
            if (float(lineSplit[3]) < float(pval_enr)):
                # append genes, enrichment term, p-value etc to list
                for i in range(0, 5):
                    tmp[i] = lineSplit[i]
                tmp[5] = lineSplit[9]
                ret_dict.append(tmp)
        ctr = ctr + 1
    return (ret_dict)


# read terms only in cluster 1 / cluster 2
@shared_task(name="read_enrichment_2")
def read_enrichment_2(path1, path2, pval_enr):
    result_file_1 = open(path1)
    result_file_2 = open(path2)
    temp_dict = {}
    temp_dict_2 = {}
    ret_dict = []
    ret_dict_2 = []
    ctr = 0
    # file 1 and array 1 is results for genes in cluster 1
    for line in result_file_1:
        tmp = {}
        lineSplit = line.split("\t")
        if (ctr > 0):
            # check p-value
            if (float(lineSplit[3]) < float(pval_enr)):
                # append genes, enrichment term, p-value etc to list
                for i in range(0, 5):
                    tmp[i] = lineSplit[i]
                tmp[5] = lineSplit[9]
                temp_dict.update({lineSplit[1]: tmp})
        ctr = ctr + 1
    # file 2 and array 2 is results for genes in cluster 2
    ctr2 = 0
    for line in result_file_2:
        tmp = {}
        lineSplit = line.split("\t")
        if (ctr2 > 0):
            if (float(lineSplit[3]) < float(pval_enr)):
                for i in range(0, 5):
                    tmp[i] = lineSplit[i]
                tmp[5] = lineSplit[9]
                temp_dict_2.update({lineSplit[1]: tmp})
        ctr2 = ctr2 + 1
    # check terms in list 1 but not in list 2
    for key in temp_dict:
        if (key not in temp_dict_2):
            ret_dict.append(temp_dict[key])
    # check terms in list 2 but not in list 1
    for key in temp_dict_2:
        if (key not in temp_dict):
            ret_dict_2.append(temp_dict_2[key])

    return (ret_dict, ret_dict_2)


############################################################
#### check input files / convert stuff #####################
############################################################

@shared_task(name="check_input_files")
def check_input_files(ppistr, exprstr):
    errstr = ""
    ppi_stringio = StringIO(ppistr)
    ppidf = pd.read_csv(ppi_stringio, sep='\t')
    expr_stringio = StringIO(exprstr)
    exprdf = pd.read_csv(expr_stringio, sep='\t')
    # check if PPI file has at least 2 columns
    if (len(ppidf.columns) < 2):
        errstr = errstr + "Input file must contain two columns with interaction partners.\n"
        # errstr = errstr + "\n \n To avoid this error, go to <a href=\"infopage.html\">the infopage</a> and make sure that your input data has the specified format."
        return (errstr)
    contains_numbers = "false"
    # check if PPI file contains lines with two protein IDs
    for i in range(len(ppidf.index)):
        if (contains_numbers == "false"):
            curr_elem = str(ppidf.iloc[[i], 0].values[0])
            curr_elem_2 = str(ppidf.iloc[[i], 1].values[0])
            if (curr_elem.isdigit() and curr_elem_2.isdigit()):
                contains_numbers = "true"
    if (contains_numbers == "false"):
        errstr = errstr + "\n" + "Input file must contain lines with Entrez IDs of interaction partners.\n"
    for column_name, column in exprdf.iterrows():
        coluniq = column.unique().tolist()
        if ("-" in coluniq):
            errstr = "Expression data contain special characters.\n"
    return (errstr)


@shared_task(name="convert_gene_list")
def convert_gene_list(adjlist, filename):
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    # conv is a list of genes with ENSEMBL Id, gene name and Entrez ID
    conv_genelist = conv['Gene name'].tolist()
    retstr = ""
    # convert list of gene IDs from gene names to NCBI ID
    for elem in adjlist:
        prot_1 = elem[0]
        prot_2 = elem[1]
        # read which element of list the gene is and read NCBI ID
        gene_nbr_1 = conv.index[conv['Gene name'] == prot_1]
        gene_nbr_2 = conv.index[conv['Gene name'] == prot_2]
        if (str(gene_nbr_1).isdigit() and str(gene_nbr_2).isdigit()):
            gene_nbr_1_2 = conv.loc[gene_nbr_1, 'NCBI gene ID'].values[0]
            gene_nbr_2_2 = conv.loc[gene_nbr_2, 'NCBI gene ID'].values[0]
            # write genes into tab separated string
            retstr = retstr + str(gene_nbr_1_2).split(".")[0] + "\t" + str(gene_nbr_2_2).split(".")[0] + "\n"
    with open(filename, "w") as text_file:
        text_file.write(retstr)


##########################################################
#### ndex import #########################################
##########################################################

# statically
# this task takes a NDEx file as a string and converts it to a two-column array with interaction partners
@shared_task(name="read_ndex_file_4")
def read_ndex_file_4(fn):
    lines6 = ""
    # read edges and nodes into arrays
    if ("edges" in fn.split("nodes")[1]):
        lines5 = fn.split("{\"nodes\":[")
        lines3 = lines5[1].split("{\"edges\":[")[0]
        # remove "cyTableColumn" from array containing edges
        if ("cyTableColumn" in lines5[1]):
            lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"cyTableColumn\":[")[0]
            lines4 = lines4[:-4]
        # take protein name from networkAttributes or nodeAttributes if it is defined there.
        elif ("networkAttributes" in lines5[1]):
            lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[0]
            lines4 = lines4[:-4]
            if ("nodeAttributes" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[
                1] and "UniprotName" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1]):
                lines6_temp = \
                lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1].split("{\"nodeAttributes\":[")[
                    1]
                lines6 = lines6_temp.split("{\"edgeAttributes\":[")[0]
        else:
            lines4 = lines5[1].split("{\"edges\":[")[1]
    # check if edge-array comes before node-array in file
    elif ("edges" in fn.split("nodes")[0]):
        lines5 = fn.split("{\"nodes\":[")
        lines3 = lines5[1].split("]},")[0] + "]]]"
        lines4 = lines5[0].split("{\"edges\":[")[1][:-4]
    # lines3 contains the nodes, lines4 the edges, lines6 contains nodeAttributes (information from the ndex file usable for the conversion from node IDs to gene IDs)
    # remove signs to allow automatic json to array conversion
    lines3.replace("@", "")
    lines3.replace("uniprot:", "uniprot")
    lines3.replace("signor:", "signor")
    lines3.replace(" ", "")
    lines3.replace("ncbigene:", "")
    lines3.replace("\\n", "")
    lines33 = lines3[:-3].replace("}]", "")
    node_line = lines33.replace("ncbigene:", "")
    nodelinesplit = node_line.split(", ")
    dictlist = []
    # node dict is later filled with keys (node IDs) and the values are NCBI gene IDs
    node_dict = {}
    if not (node_line.endswith("}")):
        node_line = node_line + "}"
    node_line_2 = "[" + node_line + "]"
    tmp2 = json.loads(node_line_2)
    node_dict_2 = {}
    # iterate over lines in nodeAttributes
    if not (lines6 == ""):
        lines6 = "[" + lines6
        # get array with nodeAttributes for current line
        tmp4 = json.loads(lines6[:-4])
        # if node element has attribute "GeneName_A", then the NCBI ID is given in the nodeAttributes
        for item in tmp4:
            if (item['n'] == "GeneName_A"):
                # use node ID and NCBI ID
                node_dict_2[item['po']] = item['v']
    # print(str(item['po']) + " " + str(item['v']))
    # print(node_dict_2)
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    conv_genelist = conv['Gene name'].tolist()
    for item in tmp2:
        dictlist.append(item)
        # write conversion from node ID to gene ID in dictionary, based on nodeAttributes from the data
        if ('r' in item):
            if (any(c.islower() for c in item['r'])):
                gene_name = item['n']
                if (gene_name in conv_genelist):
                    gene_nbr = conv.index[conv['Gene name'] == gene_name]
                    gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                    node_dict[item['@id']] = gene_nbr1
            # print(item)
            else:
                node_dict[item['@id']] = item['r']
        # print(item)
        else:
            if (item['n'].isdigit()):
                # if gene ID is in node attributes
                # print(item)
                node_dict[item['@id']] = item['n']
            elif (item['n'] in node_dict_2):
                # otherwise use conversion table to convert gene ID to NCBI ID
                gene_name = node_dict_2[item['n']]
                # print(gene_name)
                if (gene_name in conv_genelist):
                    gene_nbr = conv.index[conv['Gene name'] == gene_name]
                    gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                    node_dict[item['@id']] = gene_nbr1
        # print(gene_nbr1)
    # print(node_dict)
    # remove signs from string to allow json conversion
    lines4.replace("@", "")
    lines4.replace("uniprot:", "uniprot")
    lines4.replace("signor:", "signor")
    lines4.replace(" ", "")
    lines4 = lines4.replace("]", "")
    edge_line = lines4.rstrip()
    edge_line_2 = "[" + edge_line + "]"
    edgelinesplit = edge_line.split(", ")
    edgelist = []
    tmp4 = json.loads(edge_line_2)
    # get dictionary with gene names and NCBI IDs (entrezgene_id)
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    ret = []
    # convert node IDs in edges to NCBI IDs
    for item in tmp4:
        # print(item)
        if (item['s'] in node_dict and item['t'] in node_dict):
            source = node_dict[item['s']]
            target = node_dict[item['t']]
            # print(source)
            # print(target)
            if (source != target and not (math.isnan(float(source))) and not (math.isnan(float(target)))):
                baz = [str(int(source)), str(int(target))]
                ret.append("\t".join(baz))
    # print("\n".join(ret))
    return ("\n".join(ret))


# from the web
@shared_task(name="import_ndex")
def import_ndex(name):
    # import NDEx from server based on UUID (network contains lists with nodes and edges)
    nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=name)
    tmp4 = []
    # node_dict is for the conversion from Node ID to NCBI ID
    node_dict = {}
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    # get list of genes for later conversion
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    conv_genelist = conv['Gene name'].tolist()
    # iterate over all nodes in network
    for node_id, node in nice_cx_network.get_nodes():
        current_node = {}
        # node has ID and "name"
        current_node['id'] = node_id
        current_node['n'] = node.get('n')
        # gene names are stored in the GeneName variable in this network
        if (name == "9c38ce6e-c564-11e8-aaa6-0ac135e8bacf"):
            # get GeneName for node
            curr_gene_name = nice_cx_network.get_node_attribute_value(node_id, 'GeneName_A')
            if (curr_gene_name in conv_genelist):
                # get index in gene conversion list, convert Gene Name to NCBI ID
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                # check if NCBI ID was found
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))
        else:
            # if gene name is stored in node name
            curr_gene_name = current_node['n']
            if (curr_gene_name in conv_genelist):
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))
        tmp4.append(current_node)
    edgelist = []
    ret = ""
    # iterate over edges
    for edge_id, edge in nice_cx_network.get_edges():
        source = edge.get('s')
        target = edge.get('t')
        if (source in node_dict and target in node_dict and source != target):
            # convert source and target to NCBI IDs and write into string
            curr_edge_str = str(node_dict[source]) + "\t" + str(node_dict[target]) + "\n"
            edgelist.append([node_dict[source], node_dict[target]])
            ret = ret + curr_edge_str
    # return tab separated string
    return (ret)
