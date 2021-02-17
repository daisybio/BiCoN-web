import json
import math
import traceback
from io import StringIO, BytesIO
from os import path

import celery.states
import matplotlib.pyplot as plt
import mygene
import ndex2
import networkx as nx
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import seaborn as sns
from celery import shared_task, current_task, states
from celery.exceptions import Ignore
from django.conf import settings
from django.core.files import File
from django.utils import timezone
from lifelines.statistics import logrank_test
from networkx.readwrite import json_graph
from pybiomart import Dataset

from bicon import data_preprocessing
from bicon import BiCoN
from bicon import results_analysis

from apps.clustering.models import Job

flatten = lambda l: [item for sublist in l for item in sublist]
sns.set(color_codes=True)

"""
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


def job_finished(task_id, job_status):
    job = Job.objects.get(job_id=task_id)

    job.status = job_status
    job.finished_time = timezone.now()
    job.save()


# def check_input_files():
#     err_str = check_input_files.delay(ppi_str, expr_str).get()
#     if err_str:
#         request.session['errors'] = err_str
#         # Todo change error page to redirect to analysis? or results?
#         # return render(request, 'clustering/errorpage.html', {'errors': err_str})
#         return HttpResponseBadRequest(err_str)

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


# @shared_task(name="make_empty_figure")
# def make_empty_figure(session_id):
#     """
#     Create static picutures of the progress?!
#     TODO Remove? Use bootstrap and js
#     :param session_id:
#     :return:
#     """
#     fig = plt.figure(figsize=(10, 8))
#     if (session_id == "none"):
#         plt.savefig(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/progress.png'))
#         plt.close(fig)
#     else:
#         plt.savefig(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/progress_" + session_id + ".png'))
#         plt.close(fig)


# # empty the log file (session ID in path)
# @shared_task(name="empty_log_file")
# def empty_log_file(session_id):
#     """
#     Some logfiles, presumably just to check the status, not detailed logs?!
#     TODO: Remove and use celary to give status?
#     :param session_id:
#     :return:
#     """
#     if (session_id == "none"):
#         text_file = open(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt", "w'))
#         text_file.write("")
#         text_file.close()
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt'))):
#             text_file = open(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt", "w'))
#             text_file.write("")
#             text_file.close()
#     else:
#         text_file = open(
#             path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt", "w'))
#         text_file.write("")
#         text_file.close()
#         if (
#                 os.path.isfile(
#                     path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt'))):
#             text_file = open(
#                 path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt", "w'))
#             text_file.write("")
#             text_file.close()


########################################################
#### writing and processing metadata, loading 
#### images etc...   ###################################
#### the *actual* algorithm will be further
#### down the page       ###############################
########################################################

# @shared_task(name="preprocess_clinical_file")
# def preprocess_clinical_file(clinical_str):
#     """
#     Create TSV from CSV (or do nothing)
#     TODO: MAKE EFFICIENT; REMOVE?!
#     :param clinical_str:
#     :return:
#     """
#     if (len(clinical_str.split("\n")[0].split("\t")) > 2):
#         return (clinical_str)
#     elif ("," in clinical_str):
#         # replace comma by tab if file is CSV and not TSV
#         if ("\t" not in clinical_str.split("\n")[0]):
#             clinical_str = clinical_str.replace(",", "\t")
#     return clinical_str


# # preprocess PPI file
# @shared_task(name="preprocess_ppi_file")
# def preprocess_ppi_file(ppi_str):
#     """
#     Converts csv into tsv (see above, can be made more efficient)
#     Removes first row if it contains title columns ?!
#     Take only the right two columns
#     Misc??
#     TODO Refractor with pandas
#     :param ppi_str:
#     :return:
#     """
#     ppistr_split = ppi_str.split("\n")
#     # check if file is csv or tsv and convert it to tsv format
#     if ("\t" not in ppistr_split[2]):
#         ppi_str = ppi_str.replace(",", "\t")
#         ppistr_split = ppi_str.split("\n")
#     ppistr_split_new = []
#     len_first_line = 0
#     len_second_line = 0
#     for elem in ppistr_split[0].split("\t"):
#         if (elem != ""):
#             len_first_line = len_first_line + 1
#     for elem in ppistr_split[1].split("\t"):
#         if (elem != ""):
#             len_second_line = len_second_line + 1
#     # delete first row if it contains title columns
#     if (len_second_line > len_first_line):
#         del ppistr_split[0]
#     # take only the right two columns
#     for line in ppistr_split:
#         if (len(line.split("\t")) > 1):
#             if (len(line.split("\t")) > 2):
#                 line_length = len(line.split("\t"))
#                 line = "\t".join([line.split("\t")[line_length - 2], line.split("\t")[line_length - 1]])
#             # check if line contains two integers with protein IDs
#             if (str(line.split("\t")[0]).isdigit() and str(line.split("\t")[1].strip().replace("\n", "")).isdigit()):
#                 ppistr_split_new.append(line)
#     ppi_str = "\n".join(ppistr_split_new)
#     return (ppi_str)


# # Method to convert expression data file to TSV format and find and rename (for later recongition) column with disease type information
# @shared_task(name="preprocess_file")
# def preprocess_file(expr_str):
#     expr_str = expr_str.replace("cancer_type", "disease_type")
#     if (len(expr_str.split("\n")[0].split("\t")) > 2):
#         expr_str_split = expr_str.split("\n")
#         # replace column name for disease type
#         if ("disease_type" not in expr_str_split[0]):
#             if ("subtype" in expr_str_split[0]):
#                 expr_str = expr_str.replace("subtype", "disease_type")
#         # remove name of first column (left upper corner)
#         expr_str_first_colname = expr_str_split[0].split("\t")[0]
#         expr_str = expr_str.replace(expr_str_first_colname, "", 1)
#         expr_stringio = StringIO(expr_str)
#         exprdf = pd.read_csv(expr_stringio, sep='\t')
#         # check for column with two unique entries (pre-defined clusters)
#         for column_name, column in exprdf.transpose().iterrows():
#             if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
#                 if (len(column.unique()) == 2):
#                     expr_str = expr_str.replace(column_name, "disease_type")
#     elif ("," in expr_str):
#         # replace comma by tab if file is CSV and not TSV
#         if ("\t" not in expr_str.split("\n")[0]):
#             expr_str_split = expr_str.split("\n")
#             # replace "subtype" by "cancer type"
#             if ("disease_type" not in expr_str):
#                 if ("subtype" in expr_str):
#                     expr_str = expr_str.replace("subtype", "disease_type")
#             expr_str_first_colname = expr_str_split[0].split(",")[0]
#             expr_str = expr_str.replace(expr_str_first_colname, "", 1)
#             expr_str = expr_str.replace(",", "\t")
#             expr_str_split = expr_str.split("\n")
#             # remove entries after given length if expression data file is too big
#             if (len(expr_str_split) > 300):
#                 expr_str = "\n".join(expr_str_split[:200])
#             else:
#                 expr_str = "\n".join(expr_str_split)
#             expr_stringio = StringIO(expr_str)
#             expr_str = expr_str.replace("MCI", "CTL")
#             exprdf = pd.read_csv(expr_stringio, sep='\t')
#             # find column with two unique entries that represents disease type
#             for column_name, column in exprdf.transpose().iterrows():
#                 if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
#                     if (len(column.unique()) == 2):
#                         expr_str = expr_str.replace(column_name, "disease_type")
#             #### uncomment the following lines for automatically selecting the two biggest clusters of patients if more than 2 clusters were given
#             done1 = "false"
#             for column_name, column in exprdf.transpose().iterrows():
#                 if (not column_name.isdigit()):
#                     if (len(column.unique()) < 6):
#                         nbr_col = len(column.unique())
#                         expr_str = expr_str.replace(column_name, "disease_type")
#                         expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
#                         if (len(column.unique()) > 2 and done1 == "false"):
#                             expr_str_split_2 = []
#                             expr_str_split_2.append(expr_str_split[0])
#                             type1 = column.value_counts().index.tolist()[0]
#                             type2 = column.value_counts().index.tolist()[1]
#                             for i in range(0, len(list(column)) - 1):
#                                 if (list(column)[i] == type1 or list(column)[i] == type2):
#                                     expr_str_split_2.append(expr_str_split[i + 1])
#                             expr_str = "\n".join(expr_str_split_2)
#                             done1 = "true"
#             ########################
#
#             expr_stringio = StringIO(expr_str)
#             exprdf = pd.read_csv(expr_stringio, sep='\t')
#             return (expr_str)


# # the same as preprocess_file, but returns number of pre-defined clusters
# @shared_task(name="preprocess_file_2")
# def preprocess_file_2(expr_str):
#     expr_str = expr_str.replace("cancer_type", "disease_type")
#     nbr_col = 1
#     if (len(expr_str.split("\n")[0].split("\t")) > 2):
#         expr_str_split = expr_str.split("\n")
#         # replace column name for disease type
#         if ("disease_type" not in expr_str_split[0]):
#             if ("subtype" in expr_str_split[0]):
#                 expr_str = expr_str.replace("subtype", "disease_type")
#         # remove name of first column (left upper corner)
#         expr_str_first_colname = expr_str_split[0].split("\t")[0]
#         expr_str = expr_str.replace(expr_str_first_colname, "", 1)
#         expr_stringio = StringIO(expr_str)
#         exprdf = pd.read_csv(expr_stringio, sep='\t')
#         done1 = "false"
#         # check for column with two unique entries (pre-defined clusters)
#         for column_name, column in exprdf.transpose().iterrows():
#             if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
#                 if (len(column.unique()) == 2):
#                     expr_str = expr_str.replace(column_name, "disease_type")
#                     nbr_col = 2
#                     done1 = "true"
#         # check for column with less than 6 unique entries (pre-defined clusters)
#         for column_name, column in exprdf.transpose().iterrows():
#             if (not column_name.isdigit()):
#                 if (len(column.unique()) < 6):
#                     nbr_col = len(column.unique())
#                     expr_str = expr_str.replace(column_name, "disease_type")
#                     expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
#                     if (len(column.unique()) > 2 and done1 == "false"):
#                         expr_str_split_2 = []
#                         expr_str_split_2.append(expr_str_split[0])
#                         for i in range(0, len(list(column)) - 1):
#                             expr_str_split_2.append(expr_str_split[i + 1])
#                         expr_str = "\n".join(expr_str_split_2)
#                         done1 = "true"
#         return (expr_str, nbr_col)
#     elif ("," in expr_str):
#         # replace comma by tab if file is CSV and not TSV
#         if ("\t" not in expr_str.split("\n")[0]):
#             expr_str_split = expr_str.split("\n")
#             # replace "subtype" by "cancer type"
#             if ("disease_type" not in expr_str):
#                 if ("subtype" in expr_str):
#                     expr_str = expr_str.replace("subtype", "disease_type")
#             expr_str_first_colname = expr_str_split[0].split(",")[0]
#             expr_str = expr_str.replace(expr_str_first_colname, "", 1)
#             expr_str = expr_str.replace(",", "\t")
#             expr_str_split = expr_str.split("\n")
#             # remove entries after given length if expression data file is too big
#             if (len(expr_str_split) > 300):
#                 expr_str = "\n".join(expr_str_split[:200])
#             else:
#                 expr_str = "\n".join(expr_str_split)
#             expr_stringio = StringIO(expr_str)
#             expr_str = expr_str.replace("MCI", "CTL")
#             exprdf = pd.read_csv(expr_stringio, sep='\t')
#             # find column with two unique entries that represents disease type
#             for column_name, column in exprdf.transpose().iterrows():
#                 if ((not column_name.isdigit()) and (not (column_name == "disease_type"))):
#                     if (len(column.unique()) == 2):
#                         # print(column_name)
#                         expr_str = expr_str.replace(column_name, "disease_type")
#                         nbr_col = 2
#             #### uncomment the following lines for automatically selecting the two biggest clusters of patients if more than 2 clusters were given
#             done1 = "false"
#             for column_name, column in exprdf.transpose().iterrows():
#                 if (not column_name.isdigit()):
#                     if (len(column.unique()) < 6):
#                         nbr_col = len(column.unique())
#                         expr_str = expr_str.replace(column_name, "disease_type")
#                         expr_str_split[0] = expr_str_split[0].replace(column_name, "disease_type")
#                         if (len(column.unique()) > 2 and done1 == "false"):
#                             expr_str_split_2 = []
#                             expr_str_split_2.append(expr_str_split[0])
#                             for i in range(0, len(list(column)) - 1):
#                                 # if(list(column)[i] == type1 or list(column)[i] == type2):
#                                 expr_str_split_2.append(expr_str_split[i + 1])
#                             expr_str = "\n".join(expr_str_split_2)
#                             done1 = "true"
#             ########################
#
#             return expr_str, nbr_col
#
#
# @shared_task(name="add_loading_image")
# def add_loading_image(session_id):
#     if (session_id == "none"):
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading.gif'))):
#             copyfile(path.join(settings.MEDIA_ROOT,
#                                'clustering/userfiles/loading.gif", "/code/clustering/static/loading_1.gif'))
#         else:
#             print("loading image not found")
#     else:
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading.gif'))):
#             copyfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading.gif'),
#                      path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1_" + session_id + ".gif'))
#         else:
#             print("loading image not found")
#
#
# @shared_task(name="remove_loading_image")
# def remove_loading_image(session_id):
#     if (session_id == "none"):
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1.gif'))):
#             os.unlink(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1.gif'))
#     else:
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1_" + session_id + ".gif'))):
#             os.unlink(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1_" + session_id + ".gif'))
#         if (os.path.isfile(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1_" + session_id + ".gif'))):
#             os.unlink(path.join(settings.MEDIA_ROOT, 'clustering/userfiles/loading_1_" + session_id + ".gif'))
#
#
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
# def run_algorithm(lg_min, lg_max, expr_str, ppi_str, nbr_iter, nbr_ants, evap,
#                   epsilon, hi_sig, pher_sig, session_id, gene_set_size, nbr_groups, clinical_str, survival_col_name,
#

# (T, row_colors, col_colors, G2, means, genes_all, adj_list, genes1_algo_id, group1_ids, group2_ids, jac_1,
# jac_2) = algo_output_task(1, lg_min, lg_max, expr_str, ppi_str, nbr_iter, nbr_ants, evap,
# epsilon, hi_sig, pher_sig, session_id, gene_set_size, nbr_groups, job)
#
# (ret_metadata, path_metadata, p_val) = script_output_task(T,
# row_colors,
# col_colors,
# G2,
# means,
# genes_all,
# adj_list,
# genes1_algo_id,
# group1_ids,
# group2_ids,
# clinical_str,
# jac_1,
# jac_2,
# survival_col_name,
# clinical_df,
# session_id,
# job)
#
# job.finished_time = timezone.now()
# job.status = celery.states.SUCCESS
# job.save()
@shared_task(name="run_algorithm", base=ClusteringTaskBase)
def run_algorithm(job, expr_data_selection, expr_data_str, ppi_network_selection, ppi_network_str, L_g_min, L_g_max,
                  log2, apply_z_transformation, size=2000, n_proc=1, a=1, b=1, k=20, evaporation=0.5, th=1, eps=0.02,
                  times=6, clusters=2, cost_limit=5, max_iter=200, opt=None, show_pher=False, show_plot=False,
                  save=None, show_nets=False, use_metadata=False, survival_col_name=None, clinical_df=None):
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
                                            times=times, clusters=clusters, cost_limit=cost_limit,
                                            max_iter=max_iter, opt=opt, show_pher=show_pher, show_plot=show_plot,
                                            save=save, show_nets=show_nets)
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
    ec = nx.draw_networkx_edges(G_small, pos)

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

    script_output_task(expr_small.T, row_colors, col_colors, G_small, ret2, ret3, adjlist, new_genes1, patients1_ids,
                       patients2_ids, clinical_df, survival_col_name, job)

    current_task.update_state(
        state='RUNNING',
        meta={'progress_step': 'visualize_data',
              'progress_percent': '100'
              }
    )

    job.finished_time = timezone.now()
    job.status = celery.states.SUCCESS
    job.save()


#

##################################################################################################
######### running the algorithm - part 1 #########################################################
##################################################################################################

# @shared_task(name="algo_output_task")
# TODO Add log2

#
# def algo_output_task(s, L_g_min, L_g_max, expr_str, ppi_str, nbr_iter, nbr_ants, evap, epsilon, hi_sig, pher_sig,
#                      session_id, size, clusters_param, job, log2):
#     """
#     method for more than 2 pre-defined clusters, uses session IDs
#     :param s:
#     :param L_g_min:
#     :param L_g_max:
#     :param expr_str:
#     :param ppi_str:
#     :param nbr_iter:
#     :param nbr_ants:
#     :param evap:
#     :param epsilon:
#     :param hi_sig:
#     :param pher_sig:
#     :param session_id:
#     :param size:
#     :param clusters_param:
#     :return:
#     """
#
#     col = "disease_type"
#     not_log = True
#     expr_stringio = StringIO(expr_str)
#     expr_df = pd.read_csv(expr_stringio, sep='\t')
#
#     # check if string contains negative numbers.
#     # USED FOR: setting not_log (contains negative numbers => numbers are logarithmic
#     if ("-" in expr_str.split("\n")[2]):
#         print("expression data are logarithmized")
#         not_log = False
#     else:
#         print("expression data not logarithmized")
#
#     # this checks whether the expression data contain negative numbers
#     # USED FOR: setting not_log (contains negative numbers => numbers are logarithmic
#     for i in range(2, 4):
#         if (not_log and i > len(expr_df.columns)):
#             # check only first 1000 lines of column 2 and 3
#             for j in range(1, min(len(expr_df.index) - 1, 1000)):
#                 if (not_log and str(expr_df.columns[i]) != "disease_type"):
#                     # make integer from negative number (e.g. -1.0 -> 10), check if it is a number and check if number is negative
#                     if (expr_df.iloc[[j], [i]].to_string().__contains__('-') and str(expr_df.iloc[j][i]).replace("-",
#                                                                                                                  "",
#                                                                                                                  1).replace(
#                         ".", "", 1).isdigit()):
#                         print("expression data are logarithmized")
#                         not_log = False
#
#     # USED FOR: don't know, remove?
#     """
#     if (session_id == "none"):
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt"), "w')) as text_file:
#             text_file.write("Your files are being processed...")
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt"), "w')) as text_file:
#             text_file.write("Your files are being processed...")
#     else:
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt"), "w')) as text_file:
#             text_file.write("Your files are being processed...")
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt"), "w')) as text_file:
#             text_file.write("Your files are being processed...")
#     """
#
#     # Preprocess everything?
#     # RUN THE ACTUAL TASK
#     if (clusters_param == 2):
#         # B, G, H, n, m, GE, A_g, group1, group2, labels_B, rev_labels_B, val1, val2, group1_ids, group2_ids = lib.aco_preprocessing_strings(
#         #     expr_str, ppi_str, col, log2=not_log, gene_list=None, size=int(size), sample=None)
#
#         expression_file = StringIO(expr_str)
#         ppi_file = StringIO(ppi_str)
#
#         expr, G, labels, rev_labels = data_preprocessing(expression_file, ppi_str, log2, int(size))
#     else:
#         B, G, H, n, m, GE, A_g, labels_B, rev_labels_B = lib.aco_preprocessing_strings_2(expr_str, ppi_str, col,
#                                                                                          log2=not_log, gene_list=None,
#                                                                                          size=size, sample=None)
#
#     # USED FOR: don't know? remove?
#     """
#     if (session_id == "none"):
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt"), "w')) as text_file:
#             text_file.write("Starting model run...")
#     else:
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt"), "w')) as text_file:
#             text_file.write("Starting model run...")
#     """
#     print("How many genes you want per cluster (minimum):")
#     # L_g_min = int(input())
#     print("How many genes you want per cluster (maximum):")
#     # L_g_max = int(input())
#     imp.reload(lib)
#
#     # =============================================================================
#     # GENERAL PARAMETERS:
#     # =============================================================================
#     clusters = clusters_param  # other options are currently unavailable
#     K = int(nbr_ants)  # number of ants
#     eps = float(epsilon)  # stopping criteria: score_max-score_av<eps
#     b = float(hi_sig)  # HI significance
#     evaporation = float(evap)
#     a = float(pher_sig)  # pheramone significance
#     times = int(nbr_iter)  # max amount of iterations
#     # =============================================================================
#     # NETWORK SIZE PARAMETERS:
#     # =============================================================================
#     cost_limit = 20  # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
#
#     th = 1  # the coefficient to define the search radipus which is supposed to be bigger than
#     # mean(heruistic_information[patient]+th*std(heruistic_information[patient])
#     # bigger th - less genes are considered (can lead to empty paths if th is too high)
#     # will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
#
#     # USED FOR: don't know? remove?
#     """
#     if (session_id == "none"):
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console.txt"), "w')) as text_file:
#             text_file.write("Progress of the algorithm is shown below...")
#     else:
#         with open((path.join(settings.MEDIA_ROOT, 'clustering/userfiles/output_console_" + session_id + ".txt"), "w')) as text_file:
#             text_file.write("Progress of the algorithm is shown below...")
#     """
#
#     start = time.time()
#
#     # session id is "none" if it is not given
#     solution, t_best, sc, conv = lib.ants_new(a, b, n, m, H, GE, G, 2, cost_limit, K, evaporation, th, L_g_min, L_g_max,
#                                               eps, times, session_id, opt=None, pts=False, show_pher=False,
#                                               show_plot=True, print_runs=False, save=None, show_nets=False)
#     end = time.time()
#     n_proc = os.getenv("NBR_PROCESSES", '4')
#     lib.ants_manager(a, b, n, m, H, GE, G, 2, cost_limit, K, evaporation, th, L_g_min, L_g_max, eps, times, session_id,
#                      n_proc, opt=None, pts=False, show_pher=True, show_plot=True, save=None, show_nets=False)
#
#     print("######################################################################")
#     print("RESULTS ANALYSIS")
#     print("total time " + str(round((end - start) / 60, 2)) + " minutes")
#     print("jaccard indexes:")
#     jac_1_ret = ""
#     jac_2_ret = ""
#     if (clusters_param == 2):
#         jacindices = lib.jac_matrix(solution[1], [group1, group2])
#         print(jacindices)
#         jac_1_ret = jacindices[0]
#         jac_2_ret = jacindices[1]
#         if lib.jac(group1, solution[1][0]) > lib.jac(group1, solution[1][1]):
#             values = [val1, val2]
#         else:
#             values = [val2, val1]
#     # mapping to gene names (for now with API)
#     mg = mygene.MyGeneInfo()
#     new_genes = solution[0][0] + solution[0][1]
#     new_genes_entrez = [labels_B[x] for x in new_genes]
#     out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
#     mapping = dict()
#     for line in out:
#         if ("symbol" in line):
#             mapping[rev_labels_B[line["query"]]] = line["symbol"]
#     ###m plotting networks
#     new_genes1 = [mapping[key] for key in mapping if key in solution[0][0]]
#     new_genes2 = [mapping[key] for key in mapping if key in solution[0][1]]
#
#     genes1, genes2 = solution[0]
#     patients1, patients2 = solution[1]
#     patients1_ids = []
#     patients2_ids = []
#     for elem in patients1:
#         if (elem in labels_B):
#             patients1_ids.append(labels_B[elem])
#     for elem in patients2:
#         if (elem in labels_B):
#             patients2_ids.append(labels_B[elem])
#     means1 = [np.mean(GE[patients1].loc[gene]) - np.mean(GE[patients2].loc[gene]) for gene in genes1]
#     means2 = [np.mean(GE[patients1].loc[gene]) - np.mean(GE[patients2].loc[gene]) for gene in genes2]
#
#     G_small = nx.subgraph(G, genes1 + genes2)
#     G_small = nx.relabel_nodes(G_small, mapping)
#
#     plt.figure(figsize=(15, 15))
#     cmap = plt.cm.RdYlGn
#     vmin = -2
#     vmax = 2
#     pos = nx.spring_layout(G_small)
#     ec = nx.draw_networkx_edges(G_small, pos)
#     if (clusters_param == 2):
#         nc1 = nx.draw_networkx_nodes(G_small, nodelist=new_genes1, pos=pos, node_color=means1, node_size=600, alpha=1.0,
#                                      vmin=vmin, vmax=vmax, node_shape="^", cmap=cmap, label=values[0])
#         nc2 = nx.draw_networkx_nodes(G_small, nodelist=new_genes2, pos=pos, node_color=means2, node_size=600, alpha=1.0,
#                                      vmin=vmin, vmax=vmax, node_shape="o", cmap=cmap, label=values[1])
#     else:
#         nc1 = nx.draw_networkx_nodes(G_small, nodelist=new_genes1, pos=pos, node_color=means1, node_size=600, alpha=1.0,
#                                      vmin=vmin, vmax=vmax, node_shape="^", cmap=cmap)
#         nc2 = nx.draw_networkx_nodes(G_small, nodelist=new_genes2, pos=pos, node_color=means2, node_size=600, alpha=1.0,
#                                      vmin=vmin, vmax=vmax, node_shape="o", cmap=cmap)
#     nx.draw_networkx_labels(G_small, pos, font_size=15, font_weight='bold')
#     ret2 = means1 + means2
#     ret3 = new_genes1 + new_genes2
#     adjlist = []
#     for line in nx.generate_edgelist(G_small, data=False):
#         lineSplit = line.split()
#         adjlist.append([lineSplit[0], lineSplit[1]])
#
#     plt.legend(frameon=True)
#     try:
#         plt.colorbar(nc1)
#     except:
#         print("no colorbar found")
#     plt.axis('off')
#     ### plotting expression data
#     plt.rc('font', size=30)  # controls default text sizes
#     plt.rc('axes', titlesize=20)  # fontsize of the axes title
#     plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
#     plt.rc('xtick', labelsize=15)  # fontsize of the tick labels
#     plt.rc('ytick', labelsize=10)  # fontsize of the tick labels
#     plt.rc('legend', fontsize=30)
#
#     grouping_p = []
#     grouping_g = []
#     p_num = list(GE.columns)
#     GE_small = GE.T[genes1 + genes2]
#     GE_small.rename(columns=mapping, inplace=True)
#     GE_small = GE_small.T
#     g_num = list(GE_small.index)
#     if (clusters_param == 2):
#         for g in g_num:
#             if g in new_genes1:
#                 grouping_g.append(values[0])
#             elif g in new_genes2:
#                 grouping_g.append(values[1])
#             else:
#                 grouping_g.append(3)
#         for p in p_num:
#             if p in solution[1][0]:
#                 grouping_p.append(values[0])
#             else:
#                 grouping_p.append(values[1])
#         grouping_p = pd.DataFrame(grouping_p, index=p_num)
#         grouping_g = pd.DataFrame(grouping_g, index=g_num)
#         species = grouping_g[grouping_g[0] != 3][0]
#         lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
#         col_colors = species.map(lut)
#         species = grouping_p[0]
#         lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
#         row_colors = species.map(lut)
#     else:
#         col_colors = ""
#         row_colors = ""
#         for g in g_num:
#             if g in new_genes1:
#                 grouping_g.append("cluster1")
#             elif g in new_genes2:
#                 grouping_g.append("cluster2")
#             else:
#                 grouping_g.append(3)
#         grouping_g = pd.DataFrame(grouping_g, index=g_num)
#         species = grouping_g[grouping_g[0] != 3][0]
#         lut = {"cluster1": '#4FB6D3', "cluster2": '#22863E'}
#         col_colors = species.map(lut)
#
#     with BytesIO() as output:
#         plt.savefig(output)
#         job.ppi_png.save(f'ntw_{session_id}.png', File(output))
#
#     plt.clf()
#     plt.boxplot(conv / 2, vert=True, patch_artist=True)  # vertical box alignment  # will be used to label x-ticks
#     plt.xlabel("iterations")
#     plt.ylabel("score per subnetwork")
#
#     with BytesIO() as output:
#         plt.savefig(output)
#         job.convergence_png.save(f'conv_{session_id}.png', File(output))
#
#     return (GE_small.T, row_colors, col_colors, G_small, ret2, ret3, adjlist, new_genes1, patients1_ids, patients2_ids,
#             jac_1_ret, jac_2_ret)


##########################################################
#### running the algorithm - part 2 ######################
##########################################################


### Processing of algorithm output with using session ID
# @shared_task(name="script_output_task")
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
                    lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1].split(
                        "{\"nodeAttributes\":[")[
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
        # node has ID and "name"
        current_node = {'id': node_id, 'n': node.get('n')}
        # gene names are stored in the GeneName variable in this network
        if name == "9c38ce6e-c564-11e8-aaa6-0ac135e8bacf":
            # get GeneName for node
            curr_gene_name = nice_cx_network.get_node_attribute_value(node_id, 'GeneName_A')
            if curr_gene_name in conv_genelist:
                # get index in gene conversion list, convert Gene Name to NCBI ID
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene (formerly Entrezgene) ID'].values[0]
                # check if NCBI ID was found
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))

        else:
            # if gene name is stored in node name
            curr_gene_name = current_node['n']
            if curr_gene_name in conv_genelist:
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene (formerly Entrezgene) ID'].values[0]
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))

    edgelist = []
    ret = ""
    # iterate over edges
    for edge_id, edge in nice_cx_network.get_edges():
        source = edge.get('s')
        target = edge.get('t')
        if source in node_dict and target in node_dict and source != target:
            # convert source and target to NCBI IDs and write into string
            curr_edge_str = str(node_dict[source]) + "\t" + str(node_dict[target]) + "\n"
            edgelist.append([node_dict[source], node_dict[target]])
            ret = ret + curr_edge_str
    # return tab separated string
    return ret


def jac(x, y):
    if len(x) > 0 and len(y) > 0:
        return len(set(x).intersection(set(y))) / len((set(x).union(set(y))))
    else:
        return 0
