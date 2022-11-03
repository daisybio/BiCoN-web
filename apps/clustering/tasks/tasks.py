import json
from io import StringIO, BytesIO
from os import path

import celery.states
import matplotlib.pyplot as plt
import mygene
import networkx as nx
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import seaborn as sns
from bicon import BiCoN
from bicon import data_preprocessing
from bicon import results_analysis
from celery import shared_task, current_task, states
from celery.exceptions import Ignore
from django.conf import settings
from django.core.files import File
from django.utils import timezone
from lifelines.statistics import logrank_test
from networkx.readwrite import json_graph

from .ndex_processing import import_ndex
from apps.clustering.models import Job

flatten = lambda l: [item for sublist in l for item in sublist]
sns.set(color_codes=True)

"""
Note: This file would benefit from some HEAVY refactoring. I rewrote most parts on the way when they were broken or too
slow. I would suggest, that in stead of getting some features to work (e.g. adding enrichment) rewriting the part
might be faster (during development of the code aswell as execution). 

--- Original comments below: Most of them are obsolete but left here for completness sake ---
this file contains methods for all computationally intensive tasks needed for algorithm runs:
- algo_output_task - this is needed for running the algorithm based on PPI and expression data. it outputs
        arrays etc with the algorithm results.
- script_output_task - this is needed for processing the outputs of algo_output_task to formats used
        for data vizualisation. it writes heatmap, ppi graph, survival plot and metadata to
        files and outputs links to those files (and an array with metadata). it takes a parameter
        "session_id" that is included in the path to result files (e.g. "ppi_[SESSION_ID].json").
        this parameter can be set to "none" if you do not want to use sessions. then it uses
        static paths (e.g. "ppi.json") instead.
- import_ndex - this tasks imports PPI files from NDEx based on the UUID and parses them to the correct
        input format for the algorithm tasks.
- check_input_files - this task checks given expression and PPI files if they contain data and returns an error
        string if they do not.
- preprocess_file - this task preprocesses an input expression data file and tries to find a column with
        pre-defined clusters. if found, it is renamed to "disease_type".
- preprocess_file_2 - the same as preprecess_file, but it returns a number of pre-defined clusters additionally.
- preprocess_ppi_file - preprocesses the PPI file. it finds every row that contains two tab-separated integers (protein IDs)
        and appends them to the output file.
- preprocess_clinical_file - converts clinical file to TSV format.
- list_metadata_from_file - reads metadata from a file and returns 3 arrays with variable names and their frequency
        in cluster 1 and 2.
- run_enrichment - runs an enrichment analysis usen given terms on a list of genes.
- read_enrichment - reads results of enrichment analysis for one patient cluster and outputs dictionary with results
- read_enrichment_2 - the same as read_enrichment, but reads terms that appear only in cluster 1 or only in cluster 2
"""


class ClusteringTaskBase(celery.Task):
    def on_failure(self, exc, task_id, args, kwargs, einfo):
        job_finished(task_id, states.FAILURE)


class NetexConfig():
    NODE_TYPE = 'gene'
    IDENTIFIER = 'hugo'
    NODE_SHAPE = 'circle'
    COLOR_MAP = {-4: 'rgb(255, 0, 0)', -3: 'rgb(255, 153, 51)', -2: 'rgb(255, 204, 0)', -1: 'rgb(255, 255, 0)',
                       0: 'rgb(204, 255, 51)', 1: 'rgb(153, 255, 51)', 2: 'rgb(102, 255, 51)', 3: 'rgb(51, 204, 51)'}
    EDGE_GROUP = 'custom'
    EDGE_GROUP_NAME = 'Custom Group'
    EDGE_COLOR = 'black'
    

def convert_json_to_netex(json_data) -> dict:
    """Netex expects minimal input format like:
    
    for nodeGroups like 
    {"genes": {"type": "gene", "color": "blue", "name": "Genes", "shape": "circle"} }

    for edgeGroups like
    {"custom": {"color": "grey", "name": "Default Edge Group"} }
    
    for nodes like 
    [{'label':'xxx', 'id':'unique_y', 'group': 'genes'}]

    for edges like
    [{'from': '2', 'to': '3', 'group': 'custom'}]

    Args:
        json_data ([str]): 
            { "nodes": 
                [{"Name": "CFI", "d": 1.1, "color": "rgb(51, 204, 51)", "type": "square", "label": "CFI", "x": -0.2, "y": -0.8, "size": 10, "id": "CFI"}, ...]

            "edges": 
                [{"id": 0, "color": "rgb(0,0,0)", "source": "CFI", "target": "C4BPA"}, ...]
            }     

    Returns:
        dict: [netex input data]
    """

    def _crop_bicon_value(v):
        v = v * 2
        v_int = int(v)
        if (v < -4):
            v_int = -4
        if (v > 3):
            v_int = 3
        return v_int/2

    def _update_bicon_node(node):
        update_data = {
            'group': str(_crop_bicon_value(node['d'])),
            'name': node['label'],
        }
        node.update(update_data)
        return node
    
    def _update_bicon_edge(edge):
        update_data = {
            'from': edge['source'],
            'to': edge['target'],
            'group': NetexConfig.EDGE_GROUP
        }
        return update_data

    netex_config = {}
    netex_network = {}
    netex_network['nodes'] = [_update_bicon_node(node) for node in json_data['nodes']]
    netex_network['edges'] = [_update_bicon_edge(edge) for edge in json_data['edges']]
    # generate node group for each occuring color
    print(netex_network['nodes'])
    unique_color_values = {float(node['group']) for node in netex_network['nodes']}
    netex_config['nodeGroups'] = {
        str(color_v): {'type': NetexConfig.NODE_TYPE, 'color': NetexConfig.COLOR_MAP[color_v*2], 'groupName': str(color_v), 'shape': NetexConfig.NODE_SHAPE} for color_v in unique_color_values
        }
    # just one edge group needed
    netex_config['edgeGroups'] = {NetexConfig.EDGE_GROUP: {"color": NetexConfig.EDGE_COLOR, "groupName": NetexConfig.EDGE_GROUP_NAME} }
    # additional configuration
    netex_config['idientifier'] = NetexConfig.IDENTIFIER
    return netex_config, netex_network


def job_finished(task_id, job_status):
    job = Job.objects.get(job_id=task_id)

    job.status = job_status
    job.finished_time = timezone.now()
    job.save()


def parse_expression_data(option, expr_raw_str=None):
    dataset_path = path.join(settings.PROJECT_ROOT, 'apps/clustering/datasets')

    # Parse predefined data
    if option == 'lung-cancer':
        with open(path.join(dataset_path, 'lung_cancer_expr_nonorm.csv')) as expr_file:
            expr_str = expr_file.read()

        # This is for METADATA, recheck at it later
        clinical_df = pd.read_csv(path.join(dataset_path, 'lung_cancer_clinical.csv'))
        with open(path.join(dataset_path, 'lung_cancer_clinical.csv')) as clinical_file:
            clinical_str = clinical_file.read()

        survival_col_name = "disease free survival in months:ch1"
        nbr_groups = 2

    # Parse predefined data
    elif option == 'brest-cancer':
        with open(path.join(dataset_path, 'breast_cancer_expr.csv')) as expr_file:
            expr_str = expr_file.read()

        # This is for METADATA look at it later
        clinical_df = pd.read_csv(path.join(dataset_path, 'breast_cancer_clinical.csv'))
        with open(path.join(dataset_path, 'breast_cancer_clinical.csv')) as clinical_file:
            clinical_str = clinical_file.read()

        survival_col_name = "mfs (yr):ch1"
        nbr_groups = 2

    elif option == 'custom':
        expr_str = expr_raw_str
        clinical_df = None
        survival_col_name = None

    return expr_str, clinical_df, survival_col_name


def parse_ppi_data(option, ppi_raw_str=None):
    if option == "apid":
        ppi_str = import_ndex("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
    elif option == "string":
        ppi_str = import_ndex("275bd84e-3d18-11e8-a935-0ac135e8bacf")
    elif option == "biogrid":
        ppi_str = import_ndex("becec556-86d4-11e7-a10d-0ac135e8bacf")
    elif option == "hprd":
        ppi_str = import_ndex("1093e665-86da-11e7-a10d-0ac135e8bacf")
    elif option == 'custom':
        ppi_str = ppi_raw_str
    return ppi_str


# ##################################################################
# ################## used for metadata display #####################
# ##################################################################
#
#
# # read metadata from file and return array with data
# @shared_task(name="list_metadata_from_file")
# def list_metadata_from_file(path):
#     # used for reading metadata
#     fh1 = open(path)
#     lines = fh1.read()
#     if (lines == "NA"):
#         return ({}, {}, {})
#     # remove html from metadata file and replace table elements by tab
#     # if no data in file, remove empty dictionaries
#     if (len(lines.split('\n')) < 3):
#         return ({}, {}, {})
#     # read content from lines
#     line0 = lines.split('\n')[0].split('\t')
#     line1 = lines.split('\n')[1].split('\t')
#     line2 = lines.split('\n')[2].split('\t')
#     ret = []
#     dict3 = {}
#     dict1 = {}
#     dict2 = {}
#     dict0 = {}
#     ctr = 0
#     dict3['params'] = line0
#     dict3['gr1'] = line1
#     dict3['gr2'] = line2
#     dict3['all'] = zip(dict3['params'], dict3['gr1'], dict3['gr2'])
#     # dict 0 is parameter names, dict1 is values for group 1, dict2 is values for group 2
#     for i in range(0, len(line0) - 1):
#         dict0[i] = line0[i]
#         dict1[dict0[i]] = line1[i]
#         dict2[dict0[i]] = line2[i]
#         ctr = ctr + 1
#     return (dict0, dict1, dict2)


##################################################################################################
# Run algorithm and then make plots
##################################################################################################

@shared_task(name="run_algorithm", base=ClusteringTaskBase)
def run_algorithm(job, expr_data_selection, expr_data_str, ppi_network_selection, ppi_network_str, L_g_min, L_g_max,
                  log2, apply_z_transformation, size=2000, n_proc=1, a=1, b=1, k=20, evaporation=0.5, th=1, eps=0.02,
                  times=6, cost_limit=5, max_iter=200, opt=None, show_plot=False,
                  save=None, verbose=False, use_metadata=False, survival_col_name=None, clinical_df=None):
    # ========== Preprocess data and run algorithm ==========

    task_id = str(job.job_id)

    job.status = 'RUNNING'
    job.save()

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'submitted',
              'progress_percent': '0'
              }
    )

    # --- Step 1: Parse the strings or files and create StringIO (file object) again
    try:
        expr_str, clinical_df_demo_data, survival_col_name_demo_data = parse_expression_data(expr_data_selection,
                                                                                             expr_data_str)
        expression_file = StringIO(expr_str)
    except Exception as ex:
        current_task.update_state(
            state='ERROR',
            meta={'progress_step': 'submitted',
                  'progress_percent': '0',
                  'error_message': '\n'.join(ex.args)
                  })
        job_finished(task_id, states.FAILURE)
        raise Ignore()

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'expression_data',
              'progress_percent': '10'
              }
    )
    try:
        ppi_file = StringIO(parse_ppi_data(ppi_network_selection, ppi_network_str))
    except Exception as ex:
        current_task.update_state(
            state='ERROR',
            meta={'progress_step': 'expression_data',
                  'progress_percent': '10',
                  'error_message': '\n'.join(ex.args)
                  })
        job_finished(task_id, states.FAILURE)
        raise Ignore()

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'ppi_data',
              'progress_percent': '30'
              }
    )

    # If both demo data variables are populated use them
    if survival_col_name_demo_data:
        # Check for DF
        if isinstance(clinical_df_demo_data, pd.DataFrame):
            if not clinical_df_demo_data.empty:
                clinical_df = clinical_df_demo_data
                survival_col_name = survival_col_name_demo_data
                use_metadata = True

    # Can be removed when other part is refractored. Then clinical_df can stay none
    if not use_metadata:
        clinical_df = pd.DataFrame()

    # --- Step 2: Try and preprocess files. Catch assertions from the preprocessing function
    try:
        ge, g, labels, rev_labels = data_preprocessing(expression_file, ppi_file, log2,
                                                       zscores=apply_z_transformation, size=size)
    except AssertionError as er:
        current_task.update_state(
            state='ERROR',
            meta={'progress_step': 'ppi_data',
                  'progress_percent': '30',
                  'error_message': '\n'.join(er.args)
                  })
        job_finished(task_id, states.FAILURE)
        raise Ignore()

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'validate_preprocess_data',
              'progress_percent': '40'
              }
    )

    # --- Step 3: Run the clustering algorithm (BiCoN)
    model = BiCoN(ge, g, L_g_min, L_g_max)

    print(f'Execute model.run')
    # solution, scores = model.run_search(max_iter=1)
    # max_iter = 1  # TODO REMOVE LATER
    try:
        solution, scores = model.run_search(n_proc=n_proc, a=a, b=b, K=k, evaporation=evaporation, th=th, eps=eps,
                                            times=times, cost_limit=cost_limit,
                                            max_iter=max_iter, opt=opt, show_plot=show_plot,
                                            save=save, verbose=verbose, logging=False)
    except AssertionError as er:
        current_task.update_state(
            state='ERROR',
            meta={'progress_step': 'validate_preprocess_data',
                  'progress_percent': '40',
                  'error_message': '\n'.join(er.args)
                  })

        job_finished(task_id, states.FAILURE)
        raise Ignore()

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'run_clustering',
              'progress_percent': '90'
              }
    )

    # mapping to gene names (for now with API)
    mg = mygene.MyGeneInfo()
    new_genes = solution[0][0] + solution[0][1]
    new_genes_entrez = [labels[x] for x in new_genes]
    out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
    mapping = dict()
    for line in out:
        if ("symbol" in line):
            mapping[rev_labels[line["query"]]] = line["symbol"]

    # ========== Plotting networks ==========
    results = results_analysis(solution, labels, convert=True, origID="entrezgene")

    # --- Save results
    with StringIO() as output:
        results.save(output=output)
        job.result_csv.save(f'result_{task_id}.csv', File(output))

    # --- Heatmap
    with BytesIO() as output:
        results.show_clustermap(ge, g, output=output)
        job.heatmap_png.save(f'heatmap_{task_id}.png', File(output))

    # --- Convergence
    with BytesIO() as output:
        results.convergence_plot(scores, output=output)
        job.convergence_png.save(f'conv_{task_id}.png', File(output))

    # ToDo: Refactor the lines below
    new_genes1 = [mapping[key] for key in mapping if key in solution[0][0]]
    new_genes2 = [mapping[key] for key in mapping if key in solution[0][1]]

    genes1_algo_id, genes2_algo_id = solution[0]
    patients1_algo_id, patients2_algo_id = solution[1]
    patients1_ids = []
    patients2_ids = []
    for elem in patients1_algo_id:
        if (elem in labels):
            patients1_ids.append(labels[elem])
    for elem in patients2_algo_id:
        if (elem in labels):
            patients2_ids.append(labels[elem])
    means1 = [np.mean(ge[patients1_algo_id].loc[gene]) - np.mean(ge[patients2_algo_id].loc[gene]) for gene in
              genes1_algo_id]
    means2 = [np.mean(ge[patients1_algo_id].loc[gene]) - np.mean(ge[patients2_algo_id].loc[gene]) for gene in
              genes2_algo_id]

    G_small = nx.subgraph(g, genes1_algo_id + genes2_algo_id)
    G_small = nx.relabel_nodes(G_small, mapping)

    plt.figure(figsize=(15, 15))
    cmap = plt.cm.RdYlGn
    vmin = -2
    vmax = 2
    pos = nx.spring_layout(G_small)

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
    p_num = list(ge.columns)
    expr_small = ge.T[genes1_algo_id + genes2_algo_id]
    expr_small.rename(columns=mapping, inplace=True)
    expr_small = expr_small.T
    g_num = list(expr_small.index)

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

    with BytesIO() as output:
        plt.savefig(output)
        job.ppi_png.save(f'ntw_{job.job_id}.png', File(output))

    plt.clf()

    try:
        script_output_task(expr_small.T, row_colors, col_colors, G_small, ret2, ret3, adjlist, new_genes1, patients1_ids,
                        patients2_ids, clinical_df, survival_col_name, job)
    except Exception as e:
        print(e)

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'visualize_data',
              'progress_percent': '100'
              }
    )

    job.finished_time = timezone.now()
    job.status = celery.states.SUCCESS
    job.save()

    current_task.update_state(
        state='SUCCESS'
    )


##########################################################
#### running the algorithm - part 2 ######################
##########################################################

# @shared_task(name="script_output_task")
# ToDo: Refactor this method
def script_output_task(T, row_colors1, col_colors1, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids,
                       clinical_df, survival_col, job):
    session_id = str(job.job_id)

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
        return cmap_custom[v_int]

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


    # json_path = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/ppi_' + session_id + '.json')
    # with open(json_path, "w") as text_file:
    #     text_file.write(jsn3)
    with StringIO() as output:
        output.write(jsn3)
        job.ppi_json.save(f'ppi_{session_id}.json', File(output))

    # convert data to netex input format
    netex_config, netex_network = convert_json_to_netex(json.loads(jsn3))

    with StringIO() as output:
        output.write(json.dumps({'config': netex_config, 'network': netex_network}))
        job.netex_json.save(f'netex_{session_id}.json', File(output))

    colordict = {0: '#BB0000', 1: '#0000BB'}

    # # make heatmap (include pre-defined clusters if they were given)
    # if (isinstance(col_colors1, str)):
    #     g = sns.clustermap(T, figsize=(13, 13))
    # else:
    #     if (isinstance(row_colors1, str)):
    #         g = sns.clustermap(T, figsize=(13, 13), col_colors=col_colors1)
    #     else:
    #         g = sns.clustermap(T, figsize=(13, 13), col_colors=col_colors1, row_colors=row_colors1)
    # ax = g.ax_heatmap
    # ax.set_xlabel("Genes")
    # ax.set_ylabel("Patients")
    #
    # with BytesIO() as output:
    #     plt.savefig(output)
    #     job.heatmap_png.save(f'heatmap_{session_id}.png', File(output))

    plt.clf()
    # Array PatientData is for storing survival information
    patientData = {}

    # TODO reimplement for enrichment analysis
    # # write lists of genes in files, needed for enrichment analysis
    # if (session_id == "none"):
    #     path_genelist = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist.txt')
    #     path_genelist_1 = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist_1.txt')
    #     path_genelist_2 = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist_2.txt')
    # else:
    #
    #     path_genelist = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist_' + session_id + '.txt')
    #     path_genelist_1 = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist_1_' + session_id + '.txt')
    #     path_genelist_2 = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/genelist_2_' + session_id + '.txt')
    # with open(path_genelist, "w") as text_file_4:
    #     for i in G_list:
    #         text_file_4.write(str(i) + "\n")
    # text_file_4.close()
    # with open(path_genelist_1, "w") as text_file_5:
    #     for i in G_list:
    #         if (i in genes1):
    #             text_file_5.write(str(i) + "\n")
    # text_file_5.close()
    # with open(path_genelist_2, "w") as text_file_6:
    #     for i in G_list:
    #         if (i not in genes1):
    #             text_file_6.write(str(i) + "\n")
    # text_file_6.close()
    # if (session_id == "none"):
    #     path_metadata = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/metadata.txt')
    # else:
    #     path_metadata = path.join(settings.MEDIA_ROOT, 'clustering/userfiles/metadata_' + session_id + '.txt')
    # if no metadata given, write an empty metadata file
    p_val = ""

    # fill empty metadata arrays
    ret_metadata = []
    ret_metadata_1 = {}
    ret_metadata_2 = {}
    ret_metadata_3 = {}
    ret_metadata.append(ret_metadata_1)
    ret_metadata.append(ret_metadata_2)
    ret_metadata.append(ret_metadata_3)
    patientids_metadata = []

    if not clinical_df.empty:
        # assign some default value to survival col nbr (for the case that no survival column exists)
        survival_col_nbr = 64
        # replace all type of NA in dataframe by standard pandas-NA
        clinical_df.replace(['NaN', 'nan', '?', '--'], ['NA', 'NA', 'NA', 'NA'], inplace=True)
        clinical_df.replace(['NTL'], ['NA'], inplace=True)
        clinical_df.replace(['na'], ['NA'], inplace=True)
        clinicaldf_col_names = list(clinical_df.columns)
        patientids_metadata = [str(i) for i in clinical_df.iloc[:, 0].values.tolist()]
        # get patient ids either from first column or from index, add one empty element at the beginning of the column names
        if ("Unnamed" in "\t".join(list(clinical_df.columns))):
            clinicaldf_col_names_temp = ['empty']
            clinicaldf_col_names_new = clinicaldf_col_names_temp + clinicaldf_col_names
            clinical_df.columns = list(clinicaldf_col_names_new[:-1])
        if ("GSM" not in patientids_metadata[0]):
            patientids_metadata = list(clinical_df.index)
        param_names = []
        param_values = []
        param_cols = []

    # if clinical data were uploaded, more than 1 patient exists and patient IDs from metadata and expression data overlap
    if not (((len(group1_ids) + len(group2_ids)) < 1) or set(patientids_metadata).isdisjoint(
            group1_ids)):
        patients_0 = []
        patients_1 = []
        group1_has = []
        group2_has = []
        # iterate over columns of metadata, get all unique entries
        for column_name, column in clinical_df.transpose().iterrows():
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
                if not (jac(group1_has[i], patients_0[i]) == 0.0 and jac(group2_has[i], patients_1[i]) == 0.0):
                    param_names_final.append(param_names[i])
                    jaccards_1.append(jac(group1_has[i], patients_0[i]))
                    jaccards_2.append(jac(group2_has[i], patients_1[i]))
        # get list of patient ids
        if ("GSM" not in list(clinical_df.iloc[:, 0])[1]):
            patient_id_list = list(clinical_df.index)
        else:
            patient_id_list = list(clinical_df.iloc[:, 0])
        # check if there is a column with survival data
        if (survival_col in list(clinical_df.columns)):
            survival_col_nbr = list(clinical_df.columns).index(survival_col)
            print("column with survival data found")
            clinical_df.iloc[:, survival_col_nbr].fillna("NA", inplace=True)
            survivalcol_list = list(clinical_df.iloc[:, survival_col_nbr])
            clinical_df.iloc[:, survival_col_nbr].fillna("NA", inplace=True)
            survivalcol_list = list(clinical_df.iloc[:, survival_col_nbr])
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
        if (survival_col in list(clinical_df.columns) and len(survival_1) > 0 and len(survival_2) > 0):
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
        # text_file_4 = open((path_metadata), "w")      # READD FOR ENRICHMENT
        # text_file_4.write("\t".join(param_names) + "\n")
        # jaccards_1_str = [str(i)[:4] for i in jaccards_1]
        # jaccards_2_str = [str(i)[:4] for i in jaccards_2]
        # text_file_4.write("\t".join(jaccards_1_str) + "\n")
        # text_file_4.write("\t".join(jaccards_2_str) + "\n")
        # text_file_4.write("")
        # text_file_4.close()
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
                shape='hv',
                color='#FF7F0E'))
        trace2 = go.Scatter(
            x=list(survival_perc_2.keys()),
            y=list(survival_perc_2.values()),
            mode='lines+markers',
            name="'Group 2'",
            hoverinfo='name',
            line=dict(
                shape='hv',
                color='#1F77B4'))
        surv_data_for_graph = [trace1, trace2]
        layout = dict(showlegend=False,
                      xaxis=dict(
                          title='Time in years'),
                      yaxis=dict(
                          title='Percentage of patients'),
                      template='plotly_white')
        fig = dict(data=surv_data_for_graph, layout=layout)
        plot_div = plotly.offline.plot(fig, auto_open=False, output_type='div')

        if survival_col in list(clinical_df.columns) and len(survival_1) != 0 and errstr == "":
            with StringIO() as output:
                output.write(plot_div)
                job.survival_plotly.save(f'plotly_{session_id}.html', File(output))

    return ret_metadata, p_val


# ## enrichment stuff ##
# # run enrichment analysis
# @shared_task(name="run_enrichment")
# def run_enrichment(path, pval_enr, out_dir, terms):
#     fh1 = open(path)
#     gene_list = []
#     lines = fh1.readlines()
#     # read gene list line for line from file
#     for line in lines:
#         line.replace("\\n", "")
#         gene_list.append(line)
#     print("running enrichment analysis")
#     enr = gp.enrichr(gene_list=gene_list,
#                      description='test_name',
#                      gene_sets=terms,
#                      outdir=out_dir,
#                      cutoff=float(pval_enr)  # test dataset, use lower value of range(0,1)
#                      )
#     return (enr.results)
#
#
# # read terms in given cluster
# @shared_task(name="read_enrichment")
# def read_enrichment(path, pval_enr):
#     result_file = open(path)
#     ret_dict = []
#     ctr = 0
#     for line in result_file:
#         tmp = {}
#         lineSplit = line.split("\t")
#         if (ctr > 0):
#             # check p-value
#             if (float(lineSplit[3]) < float(pval_enr)):
#                 # append genes, enrichment term, p-value etc to list
#                 for i in range(0, 5):
#                     tmp[i] = lineSplit[i]
#                 tmp[5] = lineSplit[9]
#                 ret_dict.append(tmp)
#         ctr = ctr + 1
#     return (ret_dict)
#
#
# # read terms only in cluster 1 / cluster 2
# @shared_task(name="read_enrichment_2")
# def read_enrichment_2(path1, path2, pval_enr):
#     result_file_1 = open(path1)
#     result_file_2 = open(path2)
#     temp_dict = {}
#     temp_dict_2 = {}
#     ret_dict = []
#     ret_dict_2 = []
#     ctr = 0
#     # file 1 and array 1 is results for genes in cluster 1
#     for line in result_file_1:
#         tmp = {}
#         lineSplit = line.split("\t")
#         if (ctr > 0):
#             # check p-value
#             if (float(lineSplit[3]) < float(pval_enr)):
#                 # append genes, enrichment term, p-value etc to list
#                 for i in range(0, 5):
#                     tmp[i] = lineSplit[i]
#                 tmp[5] = lineSplit[9]
#                 temp_dict.update({lineSplit[1]: tmp})
#         ctr = ctr + 1
#     # file 2 and array 2 is results for genes in cluster 2
#     ctr2 = 0
#     for line in result_file_2:
#         tmp = {}
#         lineSplit = line.split("\t")
#         if (ctr2 > 0):
#             if (float(lineSplit[3]) < float(pval_enr)):
#                 for i in range(0, 5):
#                     tmp[i] = lineSplit[i]
#                 tmp[5] = lineSplit[9]
#                 temp_dict_2.update({lineSplit[1]: tmp})
#         ctr2 = ctr2 + 1
#     # check terms in list 1 but not in list 2
#     for key in temp_dict:
#         if (key not in temp_dict_2):
#             ret_dict.append(temp_dict[key])
#     # check terms in list 2 but not in list 1
#     for key in temp_dict_2:
#         if (key not in temp_dict):
#             ret_dict_2.append(temp_dict_2[key])
#
#     return (ret_dict, ret_dict_2)
#
#
# ############################################################
# #### check input files / convert stuff #####################
# ############################################################
#
# @shared_task(name="check_input_files")
# def check_input_files(ppistr, exprstr):
#     errstr = ""
#     ppi_stringio = StringIO(ppistr)
#     ppidf = pd.read_csv(ppi_stringio, sep='\t')
#     expr_stringio = StringIO(exprstr)
#     exprdf = pd.read_csv(expr_stringio, sep='\t')
#     # check if PPI file has at least 2 columns
#     if (len(ppidf.columns) < 2):
#         errstr = errstr + "Input file must contain two columns with interaction partners.\n"
#         # errstr = errstr + "\n \n To avoid this error, go to <a href=\"infopage.html\">the infopage</a> and make sure that your input data has the specified format."
#         return (errstr)
#     contains_numbers = "false"
#     # check if PPI file contains lines with two protein IDs
#     for i in range(len(ppidf.index)):
#         if (contains_numbers == "false"):
#             curr_elem = str(ppidf.iloc[[i], 0].values[0])
#             curr_elem_2 = str(ppidf.iloc[[i], 1].values[0])
#             if (curr_elem.isdigit() and curr_elem_2.isdigit()):
#                 contains_numbers = "true"
#     if (contains_numbers == "false"):
#         errstr = errstr + "\n" + "Input file must contain lines with Entrez IDs of interaction partners.\n"
#     for column_name, column in exprdf.iterrows():
#         coluniq = column.unique().tolist()
#         if ("-" in coluniq):
#             errstr = "Expression data contain special characters.\n"
#     return (errstr)
#
#
# @shared_task(name="convert_gene_list")
# def convert_gene_list(adjlist, filename):
#     dataset = Dataset(name='hsapiens_gene_ensembl',
#                       host='http://www.ensembl.org')
#     conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
#     # conv is a list of genes with ENSEMBL Id, gene name and Entrez ID
#     conv_genelist = conv['Gene name'].tolist()
#     retstr = ""
#     # convert list of gene IDs from gene names to NCBI ID
#     for elem in adjlist:
#         prot_1 = elem[0]
#         prot_2 = elem[1]
#         # read which element of list the gene is and read NCBI ID
#         gene_nbr_1 = conv.index[conv['Gene name'] == prot_1]
#         gene_nbr_2 = conv.index[conv['Gene name'] == prot_2]
#         if (str(gene_nbr_1).isdigit() and str(gene_nbr_2).isdigit()):
#             gene_nbr_1_2 = conv.loc[gene_nbr_1, 'NCBI gene ID'].values[0]
#             gene_nbr_2_2 = conv.loc[gene_nbr_2, 'NCBI gene ID'].values[0]
#             # write genes into tab separated string
#             retstr = retstr + str(gene_nbr_1_2).split(".")[0] + "\t" + str(gene_nbr_2_2).split(".")[0] + "\n"
#     with open(filename, "w") as text_file:
#         text_file.write(retstr)


def jac(x, y):
    if len(x) > 0 and len(y) > 0:
        return len(set(x).intersection(set(y))) / len((set(x).union(set(y))))
    else:
        return 0
