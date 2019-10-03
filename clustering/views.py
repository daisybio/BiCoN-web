from django.shortcuts import render
from io import StringIO
from django.http import HttpResponse
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
from django.contrib.auth.hashers import check_password
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.views import generic
from django.core.cache import cache
from datetime import datetime
from django.contrib.auth.models import User
from django.shortcuts import render_to_response, render
from django.template import RequestContext
from django.contrib.auth import authenticate, login
from shutil import copyfile
import shutil
#### own imports
import clustering
from clustering.models import GraphForm
from django.contrib.auth import authenticate, login, logout
from clustering.tasks import make_empty_figure, algo_output_task, empty_log_file, add_loading_image, \
    remove_loading_image, import_ndex, read_enrichment, read_ndex_file_4, read_enrichment_2, convert_gene_list, \
    check_input_files, script_output_task, preprocess_file, list_metadata_from_file, preprocess_clinical_file, \
    preprocess_ppi_file, preprocess_file_2, run_enrichment
from django.core.cache import cache
import os.path
from django.core.mail import send_mail

### *ACTUAL* imports (that have dependencies other than django and my own stuff) ####
import pandas as pd
import matplotlib.pyplot as plt


def logout_2(request):
    if (request.method == "POST"):
        logout(request)
        # return redirect('polls/logout.html')
        return redirect('clustering/logout.html')
    else:
        return render(request, 'clustering/logout.html')


def login_2(request):
    if ('username' in request.POST and 'password' in request.POST):
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            # return render(request,'polls/login.html')
            return redirect('clustering/clustering.html')
        else:
            text = "Username or password are incorrect"
            return render(request, 'clustering/login.html', {'text': text})
    # return redirect('polls/clustering.html')
    else:
        return render(request, 'clustering/login.html')


# return redirect('polls/clustering.html')

def signup(request):
    if ('username' in request.POST and 'password' in request.POST):
        username = request.POST['username']
        password = request.POST['password']
        email = ""
        if ('email' in request.POST):
            email = request.POST['email']
        if (User.objects.filter(username=username).exists()):
            text = "Username already exists. Please choose another username!"
            return render(request, 'clustering/signup.html', {'text': text})
        else:
            user = User.objects.create_user(username, email, password)
            user.save()
            current_site = get_current_site(request)
            mail_subject = 'Activate your account.'
            # message = render_to_string('acc_active_email.html', {
            # 'user': user,
            # 'domain': current_site.domain,
            # 'uid':urlsafe_base64_encode(force_bytes(user.pk)),
            # 'token':account_activation_token.make_token(user),
            # })
            userdir = "user_uploaded_files/" + username
            if not (os.path.isdir(userdir)):
                os.mkdir(userdir)
            to_email = email
            # email_message = EmailMessage(
            # mail_subject, message, to=[to_email]
            # )
            # email_message.send()
            # send_mail('Account activation', 'Your account was activated.', 'sender@example.com', ['receiver1@example.com',])
            text = "Account is being created. You will receive a confirmation e-mail soon!"
            return render(request, 'clustering/signup.html', {'text': text, 'new_user': "true"})
    # return redirect('polls/clustering.html')
    else:
        text = "Please input username and password!"
        return render(request, 'clustering/signup.html', {'text': text})


# return redirect('polls/clustering.html')


def delete_user(request):
    if ('username' in request.POST and 'password' in request.POST and request.user.is_authenticated):
        username = request.POST['username']
        password = request.POST['password']
        print(str(request.user))
        print(username)
        print(str(request.user.password))
        print(password)
        if (str(request.user) == username and check_password(password, request.user.password)):
            print("password found")
            u = User.objects.get(username=username)
            u.delete()
            user_dir = "user_uploaded_files/" + username
            if (os.path.isdir(user_dir)):
                shutil.rmtree(user_dir)
            text = "Your account was deleted."
            return render(request, 'clustering/delete_user.html', {'text': text, 'deleted': "true"})
        return render(request, 'clustering/delete_user.html', {'text': ""})
    # return redirect('polls/clustering.html')
    else:
        text = "Please input username and password!"
        return render(request, 'clustering/delete_user.html', {'text': text})


# return redirect('polls/clustering.html')


def errorpage(request):
    errors = ""
    if ('errors' in request.session):
        errors = request.session['errors']
    errors_from_cache = cache.get('errors', 'has expired')
    if not (errors_from_cache == ""):
        errors = errors_from_cache
        cache.set('errors', '')
    return render(request, 'clustering/errorpage.html', {'errors': errors})


#########################################################################
## This method displays everything on one page with no Session IDs     ##
#########################################################################


# def clustering_6_new(request):
def clustering_no_sessions(request):
    ret_metadata1 = {}
    ret_metadata2 = {}
    ret_metadata3 = {}
    metadata_dict = []
    enrichment_dict = []
    pval_enr = 0.5
    list_of_files = ""
    list_of_files_2 = ""
    save_data = request.POST.get("save_data", None)
    gene_set_size = request.POST.get("gene_set_size", 2000)
    nbr_iter = request.POST.get("nbr_iter", 45)
    nbr_ants = request.POST.get("nbr_ants", 30)
    evap = request.POST.get("evap", 0.3)
    epsilon = request.POST.get("stopcr", 0.02)
    hi_sig = request.POST.get("hisig", 1)
    pher_sig = request.POST.get("pher", 1)
    savedata_param = "false"
    if ('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
        if (request.POST['input_own_file'] and request.POST['display_old_results']):
            # configure loading page
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            make_empty_figure.delay("none")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            filename1 = request.POST.get("input_own_file")
            # get name of selected file, and path/name of other stored result files from same run
            path_json = filename1
            path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
            path_metadata = filename1.split("_json.json")[0] + "metadata.txt"
            path_plotly = filename1.split("_json.json")[0] + "plotly.html"
            path_genelist = filename1.split("_json.json")[0] + "_genelist.txt"
            path_genelist_1 = filename1.split("_json.json")[0] + "_genelist_1.txt"
            path_genelist_2 = filename1.split("_json.json")[0] + "_genelist_2.txt"
            # get locations to copy old result files to
            json_path = "userfiles/ppi.json"
            path_heatmap_2 = "userfiles/heatmap.png"
            path_metadata_2 = "userfiles/metadata.txt"
            path_plotly_2 = "userfiles/output_plotly.html"
            # copy files to static directory
            copyfile(path_json, ("clustering/static/" + json_path))
            copyfile(path_heatmap, ("clustering/static/" + path_heatmap_2))
            copyfile(path_genelist, ("clustering/static/userfiles/genelist.txt"))
            copyfile(path_genelist_1, ("clustering/static/userfiles/genelist_1.txt"))
            copyfile(path_genelist_2, ("clustering/static/userfiles/genelist_2.txt"))
            output_plot_path_2 = ""
            ret_metadata_1 = ""
            ret_metadata_2 = ""
            ret_metadata_3 = ""
            # check if plotly file exists and copy
            if (os.path.isfile(path_plotly)):
                copyfile(path_plotly, ("clustering/static/" + path_plotly_2))
                output_plot_path_2 = path_plotly_2
            # print("plot copied to")
            # print(path_plotly)
            # print(output_plot_path_2)
            # read metadata (must copy file to shared volume for processing via celery)
            if (os.path.isfile(path_metadata)):
                # print("found metadata")
                # print(path_metadata)
                copyfile(path_metadata, ("/code/clustering/static/metadata.txt"))
                filename_for_old_metadata = "/code/clustering/static/metadata.txt"
                # print(filename_for_old_metadata)
                metd = list_metadata_from_file.apply_async(args=[filename_for_old_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            # print(ret_metadata1)
            cache.clear()
            # set session ID in cache
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            make_empty_figure.apply_async(args=["none"], countdown=10)
            empty_log_file.apply_async(args=["none"], countdown=10)
            # list old files
            list_of_files = ""
            list_of_files_2 = ""
            if request.user.is_authenticated:
                username = str(request.user)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
            return render(request, 'clustering/clustering_no_sessions.html',
                          {'form': "", 'images': "", 'plot_div': "", 'script': "", 'path_heatmap': path_heatmap_2,
                           'output_plot_path': output_plot_path_2, 'json_path': json_path,
                           'list_of_files': list_of_files, 'ret_dat': "", 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'list_of_files_2': list_of_files_2})

    elif (('myfile' in request.FILES or 'predef_file' in request.POST) and (
            'protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
        # check if input files exist
        input_valid = "false"
        if ('myfile' in request.FILES and 'protfile' in request.FILES):
            if (request.FILES['myfile'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'protfile' in request.FILES):
            if (request.POST['predef_file'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        # if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
        if (input_valid == "true"):
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                if (request.POST['L_g_min'] != "" and request.POST['L_g_max'] != ""):
                    lgmin = int(request.POST['L_g_min'])
                    lgmax = int(request.POST['L_g_max'])
                else:
                    lgmin = 10
                    lgmax = 20
            else:
                lgmin = 10
                lgmax = 20
            # if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
            # if(save_data in ["save_data"]):
            #	if request.user.is_authenticated:
            #		print("saving data is true")
            # lgmin = int(request.POST['L_g_min'])
            # lgmax = int(request.POST['L_g_max'])
            # assign standard result size
            if (request.POST['L_g_min'] == ""):
                lgmin = 10
            if (request.POST['L_g_max'] == ""):
                lgmax = 20
            clinicalstr = ""
            clinicaldf = ""
            # configure loading page
            add_loading_image.delay("none")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("Your request is being processed...")
                text_file.close()
            make_empty_figure.delay("none")
            clinicalstr = "empty"
            clinicaldf = ""
            survival_col_name = ""
            # read expression file
            if ('myfile' in request.FILES):
                exprstr = request.FILES['myfile'].read().decode('utf-8')
                result10 = preprocess_file_2.delay(exprstr)
                (exprstr, nbr_groups) = result10.get()
            # read predefined expression file and clinical data
            elif ('predef_file' in request.POST and 'cancer_type' in request.POST):
                cancer_type = request.POST.get("cancer_type")
                if (cancer_type == "1"):
                    # print("babababababa")
                    fh1 = open("clustering/data/lung_cancer_expr.csv")
                    exprstr = fh1.read()
                    clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
                    fh4 = open("clustering/data/lung_cancer_clinical.csv")
                    clinicalstr = fh4.read()
                    fh4.flush()
                    fh4.close()
                    survival_col_name = "disease free survival in months:ch1"
                    nbr_groups = 2
                else:
                    fh1 = open("clustering/data/breast_cancer_expr.csv")
                    exprstr = fh1.read()
                    clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
                    fh4 = open("clustering/data/breast_cancer_clinical.csv")
                    clinicalstr = fh4.read()
                    fh4.flush()
                    fh4.close()
                    survival_col_name = "mfs (yr):ch1"
                    nbr_groups = 2
            # read PPI file
            if ('protfile' in request.FILES):
                ppistr = request.FILES['protfile'].read().decode('utf-8')
                result3 = preprocess_ppi_file.delay(ppistr)
                ppistr = result3.get()
                result4 = check_input_files.delay(ppistr, exprstr)
                errstr = result4.get()
                if (errstr != ""):
                    request.session['errors'] = errstr
                    return render(request, 'clustering/errorpage.html', {'errors': errstr})
            # read ndex file from web
            elif ('ndex_name_2' in request.POST):
                ndex_file_id = request.POST.get("ndex_name_2")
                if (ndex_file_id == "1"):
                    result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "2"):
                    # result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
                    result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "3"):
                    result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "4"):
                    result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
                    ppistr = result_ndex.get()
            # read metadata if given
            if ('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
                if (request.FILES['patientdata']):
                    clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
                    clinicalstr_first_line = clinicalstr.split("\n")[0]
                    if (len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
                        # print("is tsv")
                        clinicalstr = clinicalstr.replace("\t", ",")
                    clinical_stringio = StringIO(clinicalstr)
                    clinicaldf = pd.read_csv(clinical_stringio)
                    if ('survival_col' in request.POST):
                        if (request.POST['survival_col']):
                            survival_col_name = request.POST['survival_col']
            session_id = ""
            # assign standard value to gene set size
            if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                gene_set_size = 2000
            # run algorithm and read results
            try:
                result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap, epsilon,
                                                 hi_sig, pher_sig, "none", gene_set_size, nbr_groups)
                (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                 jac_2) = result1.get()
                # make plots and process results
                result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                   group1_ids, group2_ids, clinicalstr, jac_1, jac_2, survival_col_name,
                                                   clinicaldf, "none")
                (plot_div, ret_metadata, p_val) = result2.get()
            except:
                return render(request, 'clustering/errorpage.html',
                              {'errors': "An error occurred during the algorithm run.",
                               'hide_standard_message': "true"})
            output_plot_path = "output_plotly.html"
            json_path = "ppi.json"
            path_metadata = "/code/clustering/static/metadata.txt"
            path_heatmap = "heatmap.png"
            # json_path = "ppi_" + session_id + ".json"
            # path_heatmap = "heatmap_" + session_id + ".png"
            if (save_data in ["save_data"]):
                if request.user.is_authenticated:
                    # print("saving data in views.py")
                    username = str(request.user)
                    if not (survival_col_name == ""):
                        if ("month" in survival_col_name):
                            clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN_MONTH", 1)
                        else:
                            clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN", 1)
                    # save input data
                    GraphForm.save_user_data_3(exprstr, ppistr, clinicalstr, username)
                    curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
                    # save output data
                    copyfile(("/code/clustering/static/heatmap.png"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_heatmap.png"))
                    copyfile(("/code/clustering/static/ppi.json"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_json.json"))
                    copyfile("/code/clustering/static/metadata.txt",
                             ("user_uploaded_files/" + username + "/" + curr_time + "metadata.txt"))
                    copyfile(("/code/clustering/static/output_plotly.html"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "plotly.html"))
                    copyfile(("/code/clustering/static/genelist.txt"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_genelist.txt"))
                    copyfile(("/code/clustering/static/genelist_1.txt"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_1.txt"))
                    copyfile(("/code/clustering/static/genelist_2.txt"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_2.txt"))
            # read metadata
            ret_metadata1 = ret_metadata[0]
            ret_metadata2 = ret_metadata[1]
            ret_metadata3 = ret_metadata[2]
            # empty enrichment data from cache
            enrichment_dict = cache.get('enrichment_dict', "")
            if not (enrichment_dict == ""):
                cache.set("enrichment_dict", "")
                cache.set("enrichment_dict_2", "")
                cache.set("enrichment_dict_3", "")
                cache.set("enrichment_dict_4", "")
                cache.set("enrichment_dict_5", "")
            # paths for showing results
            # write list of genes to downloadable file
            convert_gene_list.delay(adjlist, "/code/clustering/static/genelist_temp.txt")
            # save uploaded files if specified
            # render list of previously uploaded files if user is logged in (needed if user submits another request)
            if request.user.is_authenticated:
                username = str(request.user)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
            # remove the loading-gif and progress image, clear cache
            remove_loading_image.delay("none")
            # cache.clear()
            make_empty_figure.apply_async(args=["none"], countdown=10)
            empty_log_file.apply_async(args=["none"], countdown=10)
            if (os.path.isfile("clustering/static/loading_1.gif")):
                os.unlink("clustering/static/loading_1.gif")
            cache.clear()
            # copy static files from shared directory to static-file-dir on web container
            copyfile(("/code/clustering/static/heatmap.png"), ("clustering/static/userfiles/heatmap.png"))
            copyfile(("/code/clustering/static/ppi.json"), ("clustering/static/userfiles/ppi.json"))
            copyfile(("/code/clustering/static/output_plotly.html"), ("clustering/static/userfiles/output_plotly.html"))
            copyfile(("/code/clustering/static/genelist.txt"), ("clustering/static/userfiles/genelist.txt"))
            copyfile(("/code/clustering/static/genelist_1.txt"), ("clustering/static/userfiles/genelist_1.txt"))
            copyfile(("/code/clustering/static/genelist_2.txt"), ("clustering/static/userfiles/genelist_2.txt"))
            # save session ID and metadata in cache
            cache.set('session_id', session_id)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            cache.set('json_path', "ppi.json")
            cache.set('p_val', p_val)
            cache.set('analysis_running', 'analysis_running')
            if (clinicalstr == "empty"):
                output_plot_path = "empty"
            return render(request, 'clustering/clustering_no_sessions.html',
                          {'list_of_files': list_of_files, 'ret_dat': ret_metadata, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'list_of_files_2': list_of_files_2, 'pval': p_val})
    if ('redo_analysis' in request.POST and request.user.is_authenticated):
        if (request.POST['redo_analysis']):
            with open("clustering/static/output_console.txt", "w") as text_file:
                text_file.write("Your request is being processed...")
            add_loading_image.delay("none")
            filename1 = request.POST.get("input_own_file_redo")
            fh1 = open(filename1)
            exprstr = fh1.read()
            filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
            fh2 = open(filename2)
            filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
            has_clin_data = "false"
            if (os.path.isfile(filename3)):
                fh3 = open(filename3)
                has_clin_data = "true"
                clinicalstr = fh3.read()
            ppistr = fh2.read()
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                make_empty_figure.delay("none")
                lgmin = int(request.POST['L_g_min'])
                lgmax = int(request.POST['L_g_max'])
                ## assign standard result size
                # if(request.POST['L_g_min'] == ""):
                #	lgmin = 10
                # if(request.POST['L_g_max'] == ""):
                #	lgmax = 20
                if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                    gene_set_size = 2000
                if not (has_clin_data == "true"):
                    clinicalstr = "empty"
                    ret_metadata = ""
                try:
                    result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap,
                                                     epsilon, hi_sig, pher_sig, "none", gene_set_size)
                    (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                     jac_2) = result1.get()
                    result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                       group1_ids, group2_ids, clinicalstr, jac_1, jac_2,
                                                       survival_col_name, clinicaldf, "none")
                    # (div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
                    (plot_div, ret_metadata, p_val) = result2.get()
                except:
                    return render(request, 'clustering/errorpage.html',
                                  {'errors': "An error occurred during the algorithm run.",
                                   'hide_standard_message': "true"})

                # plot2 = "test.png"
                cache.clear()
                make_empty_figure.apply_async(args=["none"], countdown=10)
                empty_log_file.apply_async(args=["none"], countdown=10)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
                remove_loading_image.delay("none")
                # copy static files from shared directory to static-file-dir on web container
                copyfile(("/code/clustering/static/heatmap.png"), ("clustering/static/userfiles/heatmap.png"))
                copyfile(("/code/clustering/static/ppi.json"), ("clustering/static/userfiles/ppi.json"))
                copyfile(("/code/clustering/static/output_plotly.html"),
                         ("clustering/static/userfiles/output_plotly.html"))
                copyfile(("/code/clustering/static/genelist.txt"), ("clustering/static/userfiles/genelist.txt"))
                copyfile(("/code/clustering/static/genelist_1.txt"), ("clustering/static/userfiles/genelist_1.txt"))
                copyfile(("/code/clustering/static/genelist_2.txt"), ("clustering/static/userfiles/genelist_2.txt"))
                return render(request, 'clustering/clustering_no_sessions.html',
                              {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                               'ret_dat': ret_metadata})
    elif ('enrichment_type' in request.POST):
        enr_type = request.POST.get("enrichment_type")
        group_for_enr = "both"
        if ('pval_enr' in request.POST):
            pval_enr = request.POST.get('pval_enr')
            print(pval_enr)
        if ('group_for_enr' in request.POST):
            group_for_enr = request.POST.get('group_for_enr')
        # print(pval_enr)
        enrichment_dict = {}
        enrichment_dict_2 = {}
        enrichment_dict_3 = {}
        enrichment_dict_4 = {}
        enrichment_dict_5 = {}
        if not (os.path.isfile("genelist.txt") and os.path.isfile("genelist_1.txt") and os.path.isfile(
                "genelist_2.txt")):
            return render(request, 'clustering/clustering_no_sessions.html', {'errstr': ""})
        try:
            if (enr_type == "kegg_enrichment"):
                result1 = run_enrichment.delay("genelist.txt", pval_enr, "/code/clustering/data/test/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                result2 = run_enrichment.delay("genelist_1.txt", pval_enr, "/code/clustering/data/test2/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                result3 = run_enrichment.delay("genelist_2.txt", pval_enr, "/code/clustering/data/test3/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_enrichment"):
                result1 = run_enrichment.delay("genelist.txt", pval_enr, "/code/clustering/data/test/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay("genelist_1.txt", pval_enr, "/code/clustering/data/test2/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay("genelist_2.txt", pval_enr, "/code/clustering/data/test3/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                # result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_molecular"):
                result1 = run_enrichment.delay("genelist.txt", pval_enr, "/code/clustering/data/test/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay("genelist_1.txt", pval_enr, "/code/clustering/data/test2/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay("genelist_2.txt", pval_enr, "/code/clustering/data/test3/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "reactome_enrichment"):
                # result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
                result2 = run_enrichment.delay("genelist_2.txt", pval_enr,
                                               "/code/clustering/data/test2/enrichr_reactome",
                                               ['Reactome_2013', 'Reactome_2016'])
                # result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
                # enr_results = result1.get()
                enr_results_2 = result2.get()
                # enr_results_3 = result3.get()
                print("enr")
                # result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",
                    pval_enr)
                # result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
                # enrichment_dict = result4.get()
                enrichment_dict = {}
                enrichment_dict_2 = result5.get()
                # enrichment_dict_3 = result6.get()
                enrichment_dict_3 = {}
            return render(request, 'clustering/clustering_no_sessions.html',
                          {'list_of_files': list_of_files, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'enrichment_dict': enrichment_dict, 'enrichment_dict_2': enrichment_dict_2,
                           'enrichment_dict_3': enrichment_dict_3, 'enrichment_dict_4': enrichment_dict_4,
                           'enrichment_dict_5': enrichment_dict_5, 'enrichment_open': "true"})
        except:
            return render(request, 'clustering/errorpage.html',
                          {'errors': "An error occurred during the algorithm run.", 'hide_standard_message': "true"})
    else:
        # remove_loading_image.delay()
        ret_metadata = ""
        if not (request.user.is_authenticated):
            cache.clear()
            metd = list_metadata_from_file.apply_async(args=["/code/clustering/static/metadata.txt"], countdown=0)
            (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            metadata_dict = [ret_metadata1, ret_metadata2, ret_metadata3]
        else:
            username = str(request.user)
            list_of_files = GraphForm.list_user_data_2(username)
            list_of_files_2 = GraphForm.list_user_data(username)
            cache.clear()
            metd = list_metadata_from_file.apply_async(args=["/code/clustering/static/metadata.txt"], countdown=0)
            (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            metadata_dict = [ret_metadata1, ret_metadata2, ret_metadata3]
        return render(request, 'clustering/clustering_no_sessions.html',
                      {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2, 'ret_metadata': ret_metadata,
                       'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                       'metadata_dict': metadata_dict, 'enrichment_dict': enrichment_dict})


#########################################################################
#### version of the page with separate input form and result display ####
#########################################################################

# def clustering_6_4_part_2(request):
def clustering_step_1(request):
    # check if user has clicked 'return' button on result page, then request.POST['newAnalysis'] is "true"
    if ('newAnalysis' in request.POST):
        print(request.POST.get("newAnalysis"))
        print(request.POST['newAnalysis'])
        request.POST._mutable = True
        done_from_cache = cache.get("done", "")
        print(done_from_cache)
        if (request.POST['newAnalysis'] != "false"):
            if ('done' in request.session):
                if (request.session['done'] == "true"):
                    # set done parameter to false if user has clicked return on result page
                    request.session['done'] = "False"
            if (done_from_cache == "done" or done_from_cache == 'done'):
                # print("done from cache")
                # set done parameter to false if user has clicked return on result page
                cache.set('done', "False")
            # remove parameter from request.POST to allow later switching to result page
            request.POST['newAnalysis'] = "false"
    if not (os.path.isdir("/code/clustering/static/userfiles")):
        os.mkdir("/code/clustering/static/userfiles")
    done_from_cache = cache.get("done", "")
    print(done_from_cache)
    # if(done_from_cache == "done"):
    #	print("done in cache")
    #	return(clustering_6_part_3_2(request))
    ret_metadata1 = {}
    ret_metadata2 = {}
    ret_metadata3 = {}
    metadata_dict = []
    enrichment_dict = []
    pval_enr = 0.5
    list_of_files = ""
    list_of_files_2 = ""
    # assign standard parameters
    save_data = request.POST.get("save_data", None)
    gene_set_size = request.POST.get("gene_set_size", 2000)
    nbr_iter = request.POST.get("nbr_iter", 45)
    nbr_ants = request.POST.get("nbr_ants", 30)
    evap = request.POST.get("evap", 0.3)
    epsilon = request.POST.get("stopcr", 0.02)
    hi_sig = request.POST.get("hisig", 1)
    pher_sig = request.POST.get("pher", 1)
    # check if the user has uploaded or chosen an expression and a PPI file
    # if(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
    # check if these files are not empty and exist
    #	if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
    if ('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
        if (request.POST['input_own_file'] and request.POST['display_old_results']):
            # configure loading page
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            make_empty_figure.delay("none")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            filename1 = request.POST.get("input_own_file")
            # get name of selected file, and path/name of other stored result files from same run
            path_json = filename1
            path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
            path_metadata = filename1.split("_json.json")[0] + "metadata.txt"
            path_plotly = filename1.split("_json.json")[0] + "plotly.html"
            path_genelist = filename1.split("_json.json")[0] + "_genelist.txt"
            path_genelist_1 = filename1.split("_json.json")[0] + "_genelist_1.txt"
            path_genelist_2 = filename1.split("_json.json")[0] + "_genelist_2.txt"
            path_ntw_2 = filename1.split("_json.json")[0] + "_ntw.png"
            # get locations to copy old result files to
            session_id = ""
            session_id_from_cache = cache.get("session_id", "none")
            if (session_id_from_cache == "none" or session_id_from_cache == ""):
                # start session for storing result data
                session_id = request.session._get_or_create_session_key()
            else:
                session_id = session_id_from_cache

            json_path = "userfiles/ppi_" + session_id + ".json"
            path_heatmap_2 = "userfiles/heatmap_" + session_id + ".png"
            path_metadata_2 = "userfiles/metadata_" + session_id + ".txt"
            path_plotly_2 = "userfiles/output_plotly_" + session_id + ".html"
            # copy files to static directory
            copyfile(path_json, ("clustering/static/" + json_path))
            copyfile(path_heatmap, ("clustering/static/" + path_heatmap_2))
            copyfile(path_genelist, ("clustering/static/userfiles/genelist_" + session_id + ".txt"))
            copyfile(path_genelist_1, ("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
            copyfile(path_genelist_2, ("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
            if (os.path.isfile(path_ntw_2)):
                copyfile(path_ntw_2, ("clustering/static/userfiles/ntw_" + session_id + ".png"))
            output_plot_path_2 = ""
            ret_metadata_1 = ""
            ret_metadata_2 = ""
            ret_metadata_3 = ""
            has_survival_plot = "false"
            # check if plotly file exists and copy
            if (os.path.isfile(path_plotly)):
                copyfile(path_plotly, ("clustering/static/" + path_plotly_2))
                output_plot_path_2 = path_plotly_2
                has_survival_plot = "true"
            # read metadata (must copy file to shared volume for processing via celery)
            if (os.path.isfile(path_metadata)):
                copyfile(path_metadata, ("/code/clustering/static/userfiles/metadata_" + session_id + ".txt"))
                filename_for_old_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
                metd = list_metadata_from_file.apply_async(args=[filename_for_old_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            cache.clear()
            # set session ID in cache
            cache.set('session_id', session_id)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('has_survival_plot', has_survival_plot)
            cache.set('ret_metadata3', ret_metadata3)
            cache.set('done', "done")
            make_empty_figure.apply_async(args=["none"], countdown=10)
            empty_log_file.apply_async(args=["none"], countdown=10)
            # list old files
            list_of_files = ""
            list_of_files_2 = ""
            network_path = "userfiles/ntw_" + session_id + ".png"
            if request.user.is_authenticated:
                username = str(request.user)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
            return render(request, 'clustering/clustering_step_2.html',
                          {'form': "", 'images': "", 'plot_div': "", 'script': "", 'path_heatmap': path_heatmap_2,
                           'output_plot_path': output_plot_path_2, 'json_path': json_path,
                           'list_of_files': list_of_files, 'ret_dat': "", 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'list_of_files_2': list_of_files_2, 'has_survival_plot': has_survival_plot,
                           'network_path': network_path})

    elif (('myfile' in request.FILES or 'predef_file' in request.POST) and (
            'protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
        # check if input files exist
        input_valid = "false"
        if ('myfile' in request.FILES and 'protfile' in request.FILES):
            if (request.FILES['myfile'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'protfile' in request.FILES):
            if (request.POST['predef_file'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        # if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
        if (input_valid == "true"):
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                # if(save_data in ["save_data"]):
                #	if request.user.is_authenticated:
                #		print("saving data is true")
                lgmin = int(request.POST['L_g_min'])
                lgmax = int(request.POST['L_g_max'])
                # start session for storing result data
                session_id = ""
                session_id_from_cache = cache.get("session_id", "none")
                print("session id from cache:" + str(session_id_from_cache))
                if (session_id_from_cache == "none" or session_id_from_cache == ""):
                    # start session for storing result data
                    session_id = request.session._get_or_create_session_key()
                else:
                    session_id = session_id_from_cache

                # assign standard result size
                if (request.POST['L_g_min'] == ""):
                    lgmin = 10
                if (request.POST['L_g_max'] == ""):
                    lgmax = 20
                clinicalstr = ""
                clinicaldf = ""
                # configure loading page
                add_loading_image.delay("none")
                add_loading_image.delay(session_id)
                with open("/code/clustering/static/output_console.txt", "w") as text_file:
                    text_file.write("Your request is being processed...")
                    text_file.close()
                make_empty_figure.delay("none")
                clinicalstr = "empty"
                clinicaldf = ""
                survival_col_name = ""
                # read expression file
                if ('myfile' in request.FILES):
                    exprstr = request.FILES['myfile'].read().decode('utf-8')
                    result10 = preprocess_file_2.delay(exprstr)
                    (exprstr, nbr_groups) = result10.get()
                # result10 = preprocess_file.delay(exprstr)
                # exprstr = result10.get()
                # read predefined expression file and clinical data
                elif ('predef_file' in request.POST and 'cancer_type' in request.POST):
                    cancer_type = request.POST.get("cancer_type")
                    if (cancer_type == "1"):
                        fh1 = open("clustering/data/lung_cancer_expr.csv")
                        exprstr = fh1.read()
                        clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
                        fh4 = open("clustering/data/lung_cancer_clinical.csv")
                        clinicalstr = fh4.read()
                        fh4.flush()
                        fh4.close()
                        survival_col_name = "disease free survival in months:ch1"
                        nbr_groups = 2
                    else:
                        fh1 = open("clustering/data/breast_cancer_expr.csv")
                        exprstr = fh1.read()
                        clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
                        fh4 = open("clustering/data/breast_cancer_clinical.csv")
                        clinicalstr = fh4.read()
                        fh4.flush()
                        fh4.close()
                        survival_col_name = "mfs (yr):ch1"
                        nbr_groups = 2
                # read PPI file
                if ('protfile' in request.FILES):
                    ppistr = request.FILES['protfile'].read().decode('utf-8')
                    result3 = preprocess_ppi_file.delay(ppistr)
                    ppistr = result3.get()
                    result4 = check_input_files.delay(ppistr, exprstr)
                    errstr = result4.get()
                    if (errstr != ""):
                        request.session['errors'] = errstr
                        return render(request, 'clustering/errorpage.html', {'errors': errstr})
                # read ndex file from web
                elif ('ndex_name_2' in request.POST):
                    ndex_file_id = request.POST.get("ndex_name_2")
                    if (ndex_file_id == "1"):
                        result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
                        ppistr = result_ndex.get()
                    elif (ndex_file_id == "2"):
                        result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
                        ppistr = result_ndex.get()
                    elif (ndex_file_id == "3"):
                        result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
                        ppistr = result_ndex.get()
                    elif (ndex_file_id == "4"):
                        result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
                        ppistr = result_ndex.get()
                # read metadata if given
                if ('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
                    if (request.FILES['patientdata']):
                        clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
                        clinicalstr_first_line = clinicalstr.split("\n")[1]
                        if (len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
                            clinicalstr = clinicalstr.replace("\t", ",")
                        clinical_stringio = StringIO(clinicalstr)
                        clinicaldf = pd.read_csv(clinical_stringio)
                        if ('survival_col' in request.POST):
                            if (request.POST['survival_col']):
                                survival_col_name = request.POST['survival_col']

                # assign standard value to gene set size
                if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                    gene_set_size = 2000
                # run algorithm and read results
                try:
                    result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap,
                                                     epsilon, hi_sig, pher_sig, session_id, gene_set_size, nbr_groups)
                    (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                     jac_2) = result1.get()
                    # make plots and process results
                    result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                       group1_ids, group2_ids, clinicalstr, jac_1, jac_2,
                                                       survival_col_name, clinicaldf, session_id)
                    (ret_metadata, path_heatmap, path_metadata, output_plot_path, json_path, p_val) = result2.get()
                except:
                    return render(request, 'clustering/errorpage.html',
                                  {'errors': "An error occurred during the algorithm run.",
                                   'hide_standard_message': "true"})
                has_survival_plot = "true"
                if (output_plot_path == "empty"):
                    has_survival_plot = "false"
                output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
                json_path = "userfiles/ppi_" + session_id + ".json"
                path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
                path_heatmap = "userfiles/heatmap_" + session_id + ".png"
                if (save_data in ["save_data"]):
                    if request.user.is_authenticated:
                        print("saving data")
                        username = str(request.user)
                        if not (survival_col_name == ""):
                            if ("month" in survival_col_name):
                                clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN_MONTH", 1)
                            else:
                                clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN", 1)
                        # save input data
                        GraphForm.save_user_data_3(exprstr, ppistr, clinicalstr, username)
                        curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
                        # save output data
                        copyfile(("/code/clustering/static/" + path_heatmap),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_heatmap.png"))
                        copyfile(("/code/clustering/static/" + json_path),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_json.json"))
                        copyfile(("/code/clustering/static/userfiles/ntw_" + session_id + ".png"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_ntw.png"))
                        copyfile(path_metadata, ("user_uploaded_files/" + username + "/" + curr_time + "metadata.txt"))
                        if (has_survival_plot == "true"):
                            copyfile(("/code/clustering/static/" + output_plot_path),
                                     ("user_uploaded_files/" + username + "/" + curr_time + "plotly.html"))
                        copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist.txt"))
                        copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_1.txt"))
                        copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_2.txt"))
                print(ret_metadata1)
                # read metadata
                ret_metadata1 = ret_metadata[0]
                ret_metadata2 = ret_metadata[1]
                ret_metadata3 = ret_metadata[2]
                # empty enrichment data from cache
                enrichment_dict = cache.get('enrichment_dict', "")
                if not (enrichment_dict == ""):
                    cache.set("enrichment_dict", "")
                    cache.set("enrichment_dict_2", "")
                    cache.set("enrichment_dict_3", "")
                    cache.set("enrichment_dict_4", "")
                    cache.set("enrichment_dict_5", "")
                # paths for showing results
                # write list of genes to downloadable file
                convert_gene_list.delay(adjlist, "/code/clustering/static/genelist_temp.txt")
                # save uploaded files if specified
                # render list of previously uploaded files if user is logged in (needed if user submits another request)
                # remove the loading-gif and progress image, clear cache
                remove_loading_image.delay("none")
                cache.clear()
                make_empty_figure.apply_async(args=["none"], countdown=10)
                empty_log_file.apply_async(args=["none"], countdown=10)
                if (os.path.isfile("clustering/static/loading_1.gif")):
                    os.unlink("clustering/static/loading_1.gif")
                if request.user.is_authenticated:
                    request.session['done'] = "true"
                    cache.set("done", "done")
                    cache.set('done', "done")
                    username = str(request.user)
                    list_of_files = GraphForm.list_user_data_2(username)
                    list_of_files_2 = GraphForm.list_user_data(username)
                cache.set("has_survival_plot", has_survival_plot)
                make_empty_figure.apply_async(args=[session_id], countdown=10)
                empty_log_file.apply_async(args=[session_id], countdown=10)
                # copy static files from shared directory to static-file-dir on web container
                copyfile(("/code/clustering/static/" + path_heatmap), ("clustering/static/" + path_heatmap))
                copyfile(("/code/clustering/static/" + json_path), ("clustering/static/" + json_path))
                copyfile(("/code/clustering/static/" + output_plot_path), ("clustering/static/" + output_plot_path))
                copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                         ("clustering/static/userfiles/genelist_" + session_id + ".txt"))
                copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                         ("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
                copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                         ("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
                path_heatmap = "userfiles/heatmap_" + session_id + ".png"
                network_file = "userfiles/ntw_" + session_id + ".png"
                conv_file = "userfiles/conv_" + session_id + ".png"
                # save session ID and metadata in cache
                cache.set('session_id', session_id)
                cache.set('ret_metadata1', ret_metadata1)
                cache.set('ret_metadata2', ret_metadata2)
                cache.set('ret_metadata3', ret_metadata3)
                cache.set('json_path', json_path)
                cache.set('p_val', p_val)
                # set "done" parameter
                cache.set('done', "done")
                cache.set('analysis_running', 'analysis_running')
                cache.set('has_survival_plot', has_survival_plot)
                if (clinicalstr == "empty"):
                    output_plot_path = "empty"
                return render(request, 'clustering/clustering_step_2.html',
                              {'form': "", 'images': "", 'plot_div': "", 'script': "", 'plot2': "",
                               'path_heatmap': path_heatmap, 'json_path': json_path,
                               'output_plot_path': output_plot_path, 'list_of_files': list_of_files,
                               'ret_dat': ret_metadata, 'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                               'ret_metadata3': ret_metadata3, 'list_of_files_2': list_of_files_2, 'pval': p_val,
                               'has_survival_plot': has_survival_plot, 'network_file': network_file,
                               'conv_file': conv_file}, status=301)
    elif ('redo_analysis' in request.POST and request.user.is_authenticated):
        if (request.POST['redo_analysis']):
            # configure loading page
            analysis_running = cache.get('analysis_running', 'none')
            # set analysis running parameter to allow display of loading images
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("Your request is being processed...")
            add_loading_image.delay("none")
            session_id = ""
            session_id_from_cache = cache.get("session_id", "none")
            if (session_id_from_cache == "none" or session_id_from_cache == ""):
                # start session for storing result data
                session_id = request.session._get_or_create_session_key()
            else:
                session_id = session_id_from_cache

            add_loading_image.delay(session_id)

            # get expression file
            filename1 = request.POST.get("input_own_file_redo")
            fh1 = open(filename1)
            exprstr = fh1.read()
            # get other filenames and files
            filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
            fh2 = open(filename2)
            ppistr = fh2.read()
            filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
            # see if clinical data were given, and load them
            has_clin_data = "false"
            clinicaldf = ""
            survival_col_name = ""
            if ('survival_col' in request.POST):
                if (request.POST['survival_col']):
                    survival_col_name = request.POST['survival_col']
            if (os.path.isfile(filename3)):
                fh3 = open(filename3)
                has_clin_data = "true"
                clinicalstr = fh3.read()
                if ("SURVIVAL_COLUMN_MONTH" in clinicalstr):
                    survival_col_name = "SURVIVAL_COLUMN_MONTH"
                elif ("SURVIVAL_COLUMN" in clinicalstr):
                    survival_col_name = "SURVIVAL_COLUMN"
                clinical_stringio = StringIO(clinicalstr)
                clinicaldf = pd.read_csv(clinical_stringio)
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                make_empty_figure.delay("none")
                make_empty_figure.delay(session_id)
                lgmin = int(request.POST['L_g_min'])
                lgmax = int(request.POST['L_g_max'])
                # assign standard result size
                if (request.POST['L_g_min'] == ""):
                    lgmin = 10
                if (request.POST['L_g_max'] == ""):
                    lgmax = 20
                if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                    gene_set_size = 2000
                # check if clinical data exist
                if not (has_clin_data == "true"):
                    survival_col_name = ""
                    clinicalstr = "empty"
                    ret_metadata = ""
                result10 = preprocess_file_2.delay(exprstr)
                (exprstr_2, nbr_groups) = result10.get()
                # run algorithm
                try:
                    result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap,
                                                     epsilon, hi_sig, pher_sig, session_id, gene_set_size, nbr_groups)
                    (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                     jac_2) = result1.get()
                    result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                       group1_ids, group2_ids, clinicalstr, jac_1, jac_2,
                                                       survival_col_name, clinicaldf, session_id)
                    (ret_metadata, path_heatmap, path_metadata, output_plot_path, json_path, p_val) = result2.get()
                except:
                    return render(request, 'clustering/errorpage.html',
                                  {'errors': "An error occurred during the algorithm run.",
                                   'hide_standard_message': "true"})

                has_survival_plot = "true"
                if (output_plot_path == "empty"):
                    has_survival_plot = "false"
                ret_metadata1 = ret_metadata[0]
                ret_metadata2 = ret_metadata[1]
                ret_metadata3 = ret_metadata[2]
                enrichment_dict = cache.get('enrichment_dict', "")
                if not (enrichment_dict == ""):
                    cache.set("enrichment_dict", "")
                    cache.set("enrichment_dict_2", "")
                    cache.set("enrichment_dict_3", "")
                    cache.set("enrichment_dict_4", "")
                    cache.set("enrichment_dict_5", "")
                output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
                json_path = "userfiles/ppi_" + session_id + ".json"
                path_heatmap = "userfiles/heatmap_" + session_id + ".png"
                path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
                if (has_clin_data == "true"):
                    if (os.path.isfile(path_metadata)):
                        metd = list_metadata_from_file.apply_async(args=[path_metadata], countdown=0)
                        (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
                cache.clear()
                make_empty_figure.apply_async(args=["none"], countdown=10)
                make_empty_figure.apply_async(args=[session_id], countdown=10)
                empty_log_file.apply_async(args=["none"], countdown=10)
                empty_log_file.apply_async(args=[session_id], countdown=10)
                list_of_files = ""
                list_of_files_2 = ""
                if (request.user.is_authenticated):
                    username = str(request.user)
                    list_of_files = GraphForm.list_user_data_2(username)
                    list_of_files_2 = GraphForm.list_user_data(username)
                cache.set('session_id', session_id)
                cache.set('ret_metadata1', ret_metadata1)
                cache.set('ret_metadata2', ret_metadata2)
                cache.set('ret_metadata3', ret_metadata3)
                cache.set('p_val', p_val)
                cache.set('done', "done")
                cache.set("has_survival_plot", has_survival_plot)
                network_file = "userfiles/ntw_" + session_id + ".png"
                conv_file = "userfiles/conv_" + session_id + ".png"
                remove_loading_image.delay("none")
                remove_loading_image.delay(session_id)
                return render(request, 'clustering/clustering_step_2.html',
                              {'path_heatmap': path_heatmap, 'output_plot_path': output_plot_path,
                               'json_path': json_path, 'list_of_files': list_of_files, 'ret_dat': "",
                               'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                               'ret_metadata3': ret_metadata3, 'list_of_files_2': list_of_files_2,
                               'has_survival_plot': has_survival_plot, 'network_file': network_file,
                               'conv_file': conv_file})
    elif ('enrichment_type' in request.POST):
        enr_type = request.POST.get("enrichment_type")
        group_for_enr = "both"
        # get p-value from POST data
        if ('pval_enr' in request.POST):
            pval_enr = request.POST.get('pval_enr')
            print(pval_enr)
        if ('group_for_enr' in request.POST):
            group_for_enr = request.POST.get('group_for_enr')
        enrichment_dict = {}
        enrichment_dict_2 = {}
        enrichment_dict_3 = {}
        enrichment_dict_4 = {}
        enrichment_dict_5 = {}
        session_id_from_cache = cache.get('session_id', "none")
        if (session_id_from_cache == "none"):
            session_id = request.session._get_or_create_session_key()
        else:
            session_id = session_id_from_cache
        has_survival_plot = ""
        surv_from_cache = cache.get('has_survival_plot', 'none')
        print(surv_from_cache)
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        try:
            # one if loop for each enrichment type due to complicated naming of result files
            if (enr_type == "kegg_enrichment"):
                # run enrichment and write to directories
                result1 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                result2 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test2/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                result3 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test3/enrichr_kegg",
                                               ['KEGG_2016', 'KEGG_2013'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                # get enrichment results as dict
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt", pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_enrichment"):
                result1 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test2/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",
                                               pval_enr, "clustering/data/test3/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_molecular"):
                result1 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",
                                               pval_enr, "/code/clustering/data/test/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",
                                               pval_enr, "/code/clustering/data/test2/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",
                                               pval_enr, "/code/clustering/data/test3/enrichr_go",
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay(
                    "/code/clustering/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result6 = read_enrichment.delay(
                    "/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                result7 = read_enrichment_2.delay(
                    "/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    "/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",
                    pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "reactome_enrichment"):
                result2 = run_enrichment.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",
                                               pval_enr, "/code/clustering/data/test2/enrichr_reactome",
                                               ['Reactome_2013', 'Reactome_2016'])
                enr_results_2 = result2.get()
                print("enr")
                result5 = read_enrichment.delay(
                    "/code/clustering/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",
                    pval_enr)
                enrichment_dict = {}
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = {}
            if ('enr' not in request.POST):
                mutable = request.POST._mutable
                request.POST._mutable = True
                request.POST['enr'] = "true"
                request.POST._mutable = mutable
            if ('enr' in request.POST):
                print("enr in request")
            # get current session id and result files for session
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            json_path = "userfiles/ppi_" + session_id + ".json"
            output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
            conv_file = "userfiles/conv_" + session_id + ".png"
            ret_metadata1 = cache.get('ret_metadata1', "")
            ret_metadata2 = cache.get('ret_metadata2', "")
            ret_metadata3 = cache.get('ret_metadata3', "")
            if (ret_metadata1 == "" and os.path.isfile(path_metadata)):
                metd = list_metadata_from_file.apply_async(args=[path_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            if ('ppi_path' in request.POST and 'heatmap_path' in request.POST and 'plot_path' in request.POST):
                print(request.POST.get('ppi_path'))
                path_heatmap = request.POST.get('heatmap_path')
                json_path = request.POST.get('ppi_path')
                output_plot_path = request.POST.get('plot_path')
            return render(request, 'clustering/clustering_step_2.html',
                          {'list_of_files': list_of_files, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3, 'path_heatmap': path_heatmap,
                           'json_path': json_path, 'output_plot_path': output_plot_path,
                           'enrichment_dict': enrichment_dict, 'enrichment_dict_2': enrichment_dict_2,
                           'enrichment_dict_3': enrichment_dict_3, 'enrichment_dict_4': enrichment_dict_4,
                           'enrichment_dict_5': enrichment_dict_5, 'enrichment_open': "true",
                           'has_survival_plot': has_survival_plot, 'conv_file': conv_file})
        except:
            return render(request, 'clustering/errorpage.html',
                          {'errors': "An error occurred during the algorithm run.", 'hide_standard_message': "true"})
    else:
        ret_metadata = ""
        metdata_dict = ""
        done_from_cache = cache.get("done", "")
        analysis_running = cache.get('analysis_running', 'none')
        session_id_from_cache = cache.get('session_id', 'has_expired')
        surv_from_cache = cache.get('has_survival_plot', 'none')
        has_survival_plot = "false"
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        if (session_id_from_cache == 'has_expired' or session_id_from_cache == ""):
            session_id = request.session._get_or_create_session_key()
            cache.set('session_id', session_id)
            print("session ID: " + str(session_id))
        else:
            session_id = session_id_from_cache
        # if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
        if (analysis_running == 'none'):
            print("analysis not running")
            if (os.path.isfile("clustering/static/loading_1.gif")):
                os.unlink("clustering/static/loading_1.gif")
            if (os.path.isfile("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")):
                os.unlink("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")
            if (os.path.isfile("clustering/static/progress.png")):
                os.unlink("clustering/static/progress.png")
            if (os.path.isfile("/code/clustering/static/progress.png")):
                os.unlink("/code/clustering/static/progress.png")
            if (os.path.isfile("/code/clustering/static/userfiles/progress_" + session_id + ".png")):
                os.unlink("/code/clustering/static/userfiles/progress_" + session_id + ".png")
            remove_loading_image.delay("none")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            with open("clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
        # set session ID (always), read metadata (if loading page new because user has hit reload)
        if not (request.user.is_authenticated):
            cache.clear()
            cache.set('analysis_running', analysis_running)
            if (session_id_from_cache != 'has_expired'):
                cache.set('session_id', session_id_from_cache)
            cache.set('done', done_from_cache)
            if not (session_id == ""):
                if (os.path.isfile("/code/clustering/static/userfiles/metadata_" + session_id + ".txt")):
                    metd = list_metadata_from_file.apply_async(
                        args=["/code/clustering/static/userfiles/metadata_" + session_id + ".txt"], countdown=0)
                    (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
                    metadata_dict = [ret_metadata1, ret_metadata2, ret_metadata3]
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
        # get session ID and user files if user is logged in
        else:
            username = str(request.user)
            list_of_files = GraphForm.list_user_data_2(username)
            list_of_files_2 = GraphForm.list_user_data(username)
            if not (session_id_from_cache == 'has expired'):
                session_id = session_id_from_cache
            else:
                session_id = request.session._get_or_create_session_key()
            if not (session_id == ""):
                if (os.path.isfile("/code/clustering/static/userfiles/metadata_" + session_id + ".txt")):
                    metd = list_metadata_from_file.apply_async(
                        args=["/code/clustering/static/userfiles/metadata_" + session_id + ".txt"], countdown=0)
                    (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
                    metadata_dict = [ret_metadata1, ret_metadata2, ret_metadata3]
            cache.clear()
            cache.set('analysis_running', analysis_running)
            cache.set('session_id', session_id)
            cache.set('done', done_from_cache)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
        # print(ret_metadata1)
        # if('done' in request.session):
        #	if(request.session['done'] == "true"):
        #		session_id = request.session._get_or_create_session_key()
        #		path_heatmap = "userfiles/heatmap_" + session_id + ".png"
        #		json_path = "userfiles/ppi_" + session_id + ".json"
        #		path_metadata = "userfiles/metadata_" + session_id + ".txt"
        #		output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
        #		conv_file = "userfiles/conv_" + session_id + ".png"
        #		return render(request,'clustering/clustering_step_2.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'has_survival_plot':has_survival_plot,'conv_file':conv_file})
        # this redirects the user to the result page if the done parameter is true (e.g. reloading page after an analysis)
        if (done_from_cache == "done"):
            session_id_from_cache = cache.get("session_id", 'has expired')
            if (session_id_from_cache == 'has expired'):
                session_id = request.session._get_or_create_session_key()
            else:
                session_id = session_id_from_cache
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            json_path = "userfiles/ppi_" + session_id + ".json"
            path_metadata = "userfiles/metadata_" + session_id + ".txt"
            output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
            conv_file = "userfiles/conv_" + session_id + ".png"
            return render(request, 'clustering/clustering_step_2.html',
                          {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                           'ret_metadata': ret_metadata, 'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                           'ret_metadata3': ret_metadata3, 'path_heatmap': path_heatmap, 'json_path': json_path,
                           'output_plot_path': output_plot_path, 'metadata_dict': metadata_dict,
                           'enrichment_dict': enrichment_dict, 'has_survival_plot': has_survival_plot,
                           'conv_file': conv_file}, status=301)
        if (session_id_from_cache == 'has_expired'):
            session_id = request.session._get_or_create_session_key()
        else:
            session_id = session_id_from_cache
        cache.set('session_id', session_id)
        output_console_file = "userfiles/output_console_" + session_id + ".txt"
        progress_file = "userfiles/progress_" + session_id + ".png"
        loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
        print("session id: " + str(session_id))
        return render(request, 'clustering/clustering_step_1.html',
                      {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2, 'ret_metadata': ret_metadata,
                       'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                       'metadata_dict': metadata_dict, 'enrichment_dict': enrichment_dict,
                       'progress_file': progress_file, 'output_console_file': output_console_file,
                       'loading_gif_1': loading_gif_1})


#########################################################################
#### version of the page that displays input form and result together ###
#########################################################################


def clustering(request):
    # the parameter analysis_running is true when an analysis has been run while the current cache exists. If it is false and an empty request is submitted (which is when an user first accesses the
    # page), and for some reason the output-console file for the progress page is filled with text, it gets emptied and the loading-gif removed.
    analysis_running = cache.get('analysis_running', 'none')
    # print(analysis_running)
    session_id_from_cache = cache.get('session_id', 'has expired')
    # print("session ID: " +str(session_id_from_cache))
    nbr_processes = os.getenv("NBR_PROCESSES", '4')
    print(nbr_processes)
    ret_metadata1 = {}
    ret_metadata2 = {}
    ret_metadata3 = {}
    metadata_dict = []
    enrichment_dict = []
    # assign standard parameters to input variables if none are given
    pval_enr = 0.5
    list_of_files = ""
    list_of_files_2 = ""
    save_data = request.POST.get("save_data", None)
    gene_set_size = request.POST.get("gene_set_size", 2000)
    nbr_iter = request.POST.get("nbr_iter", 45)
    nbr_ants = request.POST.get("nbr_ants", 30)
    evap = request.POST.get("evap", 0.3)
    epsilon = request.POST.get("stopcr", 0.02)
    hi_sig = request.POST.get("hisig", 1)
    pher_sig = request.POST.get("pher", 1)
    if not (os.path.isdir("/code/clustering/static/userfiles")):
        os.mkdir("/code/clustering/static/userfiles")
    # if the user wants to display old results
    if ('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
        if (request.POST['input_own_file'] and request.POST['display_old_results']):
            # configure loading page
            make_empty_figure.delay("none")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            filename1 = request.POST.get("input_own_file")
            # get name of selected file, and path/name of other stored result files from same run
            path_json = filename1
            path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
            path_metadata = filename1.split("_json.json")[0] + "metadata.txt"
            path_plotly = filename1.split("_json.json")[0] + "plotly.html"
            path_genelist = filename1.split("_json.json")[0] + "_genelist.txt"
            path_genelist_1 = filename1.split("_json.json")[0] + "_genelist_1.txt"
            path_genelist_2 = filename1.split("_json.json")[0] + "_genelist_2.txt"
            path_ntw_2 = filename1.split("_json.json")[0] + "_ntw.png"
            # get locations to copy old result files to
            session_id = request.session._get_or_create_session_key()
            json_path = "userfiles/ppi_" + session_id + ".json"
            path_heatmap_2 = "userfiles/heatmap_" + session_id + ".png"
            path_metadata_2 = "userfiles/metadata_" + session_id + ".txt"
            path_plotly_2 = "userfiles/output_plotly_" + session_id + ".html"
            # copy files to static directory
            copyfile(path_json, ("/code/clustering/static/" + json_path))
            copyfile(path_heatmap, ("/code/clustering/static/" + path_heatmap_2))
            if (os.path.isfile(path_ntw_2)):
                copyfile(path_ntw_2, ("/code/clustering/static/userfiles/ntw_" + session_id + ".png"))
            copyfile(path_genelist, ("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"))
            copyfile(path_genelist_1, ("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
            copyfile(path_genelist_2, ("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
            output_plot_path_2 = ""
            has_survival_plot = "false"
            # check if plotly file exists and copy
            if (os.path.isfile(path_plotly)):
                copyfile(path_plotly, ("/code/clustering/static/" + path_plotly_2))
                output_plot_path_2 = path_plotly_2
                has_survival_plot = "true"
            ret_metadata_1 = ""
            ret_metadata_2 = ""
            ret_metadata_3 = ""
            # read metadata (must copy file to shared volume for processing via celery)
            if (os.path.isfile(path_metadata)):
                copyfile(path_metadata, ("/code/clustering/static/userfiles/metadata_" + session_id + ".txt"))
                filename_for_old_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
                metd = list_metadata_from_file.apply_async(args=[filename_for_old_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            cache.clear()
            cache.set('has_survival_plot', has_survival_plot)
            # set session ID in cache
            cache.set('session_id', session_id)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            make_empty_figure.apply_async(args=["none"], countdown=10)
            empty_log_file.apply_async(args=["none"], countdown=10)
            # list old files
            list_of_files = ""
            list_of_files_2 = ""
            if request.user.is_authenticated:
                username = str(request.user)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
            output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
            progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
            loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
            network_path = "userfiles/ntw_" + session_id + ".png"
            return render(request, 'clustering/clustering.html',
                          {'form': "", 'images': "", 'plot_div': "", 'script': "", 'path_heatmap': path_heatmap_2,
                           'output_plot_path': output_plot_path_2, 'json_path': json_path,
                           'list_of_files': list_of_files, 'ret_dat': "", 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'list_of_files_2': list_of_files_2, 'has_survival_plot': has_survival_plot,
                           'output_console_file': output_console_file, 'progress_file': progress_file,
                           'loading_gif_1': loading_gif_1, 'network_path': network_path})

    elif (('myfile' in request.FILES or 'predef_file' in request.POST) and (
            'protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
        # check if input files exist
        input_valid = "false"
        if ('myfile' in request.FILES and 'protfile' in request.FILES):
            if (request.FILES['myfile'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'protfile' in request.FILES):
            if (request.POST['predef_file'] and request.FILES['protfile']):
                input_valid = "true"
        elif ('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
            if (request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
                input_valid = "true"
        if (input_valid == "true"):
            # set analysis running parameter to allow display of loading image
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            # get number of genes per cluster/assign standard parameter
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                if (request.POST['L_g_min'] != "" and request.POST['L_g_max'] != ""):
                    lgmin = int(request.POST['L_g_min'])
                    lgmax = int(request.POST['L_g_max'])
                else:
                    lgmin = 10
                    lgmax = 20
            else:
                lgmin = 10
                lgmax = 20
            ## assign standard result size
            session_id = ""
            session_id_from_cache = cache.get("session_id", "none")
            if (session_id_from_cache == "none" or session_id_from_cache == ""):
                # start session for storing result data
                session_id = request.session._get_or_create_session_key()
            else:
                session_id = session_id_from_cache

            # configure loading page
            add_loading_image.delay("none")
            add_loading_image.delay(session_id)
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("Your request is being processed...")
                text_file.close()
            make_empty_figure.delay("none")
            make_empty_figure.delay(session_id)
            clinicalstr = "empty"
            clinicaldf = ""
            survival_col_name = ""
            # read expression file
            if ('myfile' in request.FILES):
                exprstr = request.FILES['myfile'].read().decode('utf-8')
                result10 = preprocess_file_2.delay(exprstr)
                (exprstr, nbr_groups) = result10.get()
            # read predefined expression file and clinical data
            elif ('predef_file' in request.POST and 'cancer_type' in request.POST):
                cancer_type = request.POST.get("cancer_type")
                if (cancer_type == "1"):
                    fh1 = open("clustering/data/lung_cancer_expr.csv")
                    exprstr = fh1.read()
                    clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
                    fh4 = open("clustering/data/lung_cancer_clinical.csv")
                    clinicalstr = fh4.read()
                    fh4.flush()
                    fh4.close()
                    survival_col_name = "disease free survival in months:ch1"
                    nbr_groups = 2
                else:
                    fh1 = open("clustering/data/breast_cancer_expr.csv")
                    exprstr = fh1.read()
                    clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
                    fh4 = open("clustering/data/breast_cancer_clinical.csv")
                    clinicalstr = fh4.read()
                    fh4.flush()
                    fh4.close()
                    survival_col_name = "mfs (yr):ch1"
                    nbr_groups = 2
            # read PPI file
            if ('protfile' in request.FILES):
                ppistr = request.FILES['protfile'].read().decode('utf-8')
                result3 = preprocess_ppi_file.delay(ppistr)
                ppistr = result3.get()
                result4 = check_input_files.delay(ppistr, exprstr)
                errstr = result4.get()
                if (errstr != ""):
                    request.session['errors'] = errstr
                    return render(request, 'clustering/errorpage.html', {'errors': errstr})
            # read ndex file from web
            elif ('ndex_name_2' in request.POST):
                ndex_file_id = request.POST.get("ndex_name_2")
                if (ndex_file_id == "1"):
                    result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "2"):
                    result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "3"):
                    result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
                    ppistr = result_ndex.get()
                elif (ndex_file_id == "4"):
                    result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
                    ppistr = result_ndex.get()
            # read metadata if given
            if ('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
                if (request.FILES['patientdata']):
                    clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
                    clinicalstr_first_line = clinicalstr.split("\n")[1]
                    if (len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
                        print("converting metadata to csv format")
                        clinicalstr = clinicalstr.replace("\t", ",")
                    clinical_stringio = StringIO(clinicalstr)
                    clinicaldf = pd.read_csv(clinical_stringio)
                    if ('survival_col' in request.POST):
                        if (request.POST['survival_col']):
                            survival_col_name = request.POST['survival_col']
            # assign standard value to gene set size
            if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                gene_set_size = 2000
            # run algorithm and read results
            try:
                result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap, epsilon,
                                                 hi_sig, pher_sig, session_id, gene_set_size, nbr_groups)
                (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                 jac_2) = result1.get()
                # make plots and process results
                result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                   group1_ids, group2_ids, clinicalstr, jac_1, jac_2, survival_col_name,
                                                   clinicaldf, session_id)
                (ret_metadata, path_heatmap, path_metadata, output_plot_path, json_path, p_val) = result2.get()
            except:
                return render(request, 'clustering/errorpage.html',
                              {'errors': "An error occurred during the algorithm run.",
                               'hide_standard_message': "true"})
            # set variable for display of survival plot to true if survival plot exists
            has_survival_plot = "true"
            if (output_plot_path == "empty"):
                cache.set("has_survival_plot", "false")
                has_survival_plot = "false"
            else:
                cache.set("has_survival_plot", "true")
            # set paths for result files
            output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
            json_path = "userfiles/ppi_" + session_id + ".json"
            path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            # save uploaded files if specified
            if (save_data in ["save_data"]):
                if request.user.is_authenticated:
                    username = str(request.user)
                    # replace name of survival column by "survival column"
                    if not (survival_col_name == ""):
                        if ("month" in survival_col_name):
                            # write "month" to indicate that survival time is given in months
                            clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN_MONTH", 1)
                        else:
                            clinicalstr = clinicalstr.replace(survival_col_name, "SURVIVAL_COLUMN", 1)
                    # save input data
                    GraphForm.save_user_data_3(exprstr, ppistr, clinicalstr, username)
                    curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
                    # save output data
                    copyfile(("/code/clustering/static/" + path_heatmap),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_heatmap.png"))
                    copyfile(("/code/clustering/static/" + json_path),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_json.json"))
                    copyfile(("/code/clustering/static/userfiles/ntw_" + session_id + ".png"),
                             ("user_uploaded_files/" + username + "/" + curr_time + "_ntw.png"))
                    copyfile(path_metadata, ("user_uploaded_files/" + username + "/" + curr_time + "metadata.txt"))
                    if (has_survival_plot == "true"):
                        copyfile(("/code/clustering/static/" + output_plot_path),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "plotly.html"))
                    if (os.path.isfile("/code/clustering/static/userfiles/genelist_" + session_id + ".txt")):
                        copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist.txt"))
                        copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_1.txt"))
                        copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                 ("user_uploaded_files/" + username + "/" + curr_time + "_genelist_2.txt"))
            # read metadata
            ret_metadata1 = ret_metadata[0]
            ret_metadata2 = ret_metadata[1]
            ret_metadata3 = ret_metadata[2]
            # empty enrichment data from cache
            enrichment_dict = cache.get('enrichment_dict', "")
            if not (enrichment_dict == ""):
                cache.set("enrichment_dict", "")
                cache.set("enrichment_dict_2", "")
                cache.set("enrichment_dict_3", "")
                cache.set("enrichment_dict_4", "")
                cache.set("enrichment_dict_5", "")
            # write list of genes to downloadable file
            convert_gene_list.delay(adjlist, "/code/clustering/static/userfiles/genelist_temp.txt")
            # render list of previously uploaded files if user is logged in (needed if user submits another request)
            if request.user.is_authenticated:
                username = str(request.user)
                list_of_files = GraphForm.list_user_data_2(username)
                list_of_files_2 = GraphForm.list_user_data(username)
            # remove the loading-gif and progress image, clear cache
            remove_loading_image.delay("none")
            remove_loading_image.delay(session_id)
            cache.clear()
            make_empty_figure.apply_async(args=["none"], countdown=10)
            make_empty_figure.apply_async(args=[session_id], countdown=10)
            empty_log_file.apply_async(args=["none"], countdown=10)
            empty_log_file.apply_async(args=[session_id], countdown=10)
            # remove loading gif
            if (os.path.isfile("clustering/static/loading_1.gif")):
                os.unlink("clustering/static/loading_1.gif")
            # copy static files from shared directory to static-file-dir on web container
            copyfile(("/code/clustering/static/" + path_heatmap), ("clustering/static/" + path_heatmap))
            copyfile(("/code/clustering/static/" + json_path), ("clustering/static/" + json_path))
            copyfile(("/code/clustering/static/" + output_plot_path), ("clustering/static/" + output_plot_path))
            copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                     ("clustering/static/userfiles/genelist_" + session_id + ".txt"))
            copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                     ("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
            copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                     ("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            # save session ID and metadata in cache
            cache.set('session_id', session_id)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            cache.set('json_path', json_path)
            cache.set('p_val', p_val)
            cache.set('analysis_running', 'analysis_running')
            cache.set('has_survival_plot', has_survival_plot)
            output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
            progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
            loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
            network_file = "userfiles/ntw_" + session_id + ".png"
            conv_file = "userfiles/conv_" + session_id + ".png"
            if (clinicalstr == "empty"):
                output_plot_path = "empty"
            return render(request, 'clustering/clustering.html',
                          {'path_heatmap': path_heatmap, 'output_plot_path': output_plot_path, 'json_path': json_path,
                           'list_of_files': list_of_files, 'ret_dat': ret_metadata, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'list_of_files_2': list_of_files_2, 'p_val': p_val, 'has_survival_plot': has_survival_plot,
                           'output_console_file': output_console_file, 'progress_file': progress_file,
                           'loading_gif_1': loading_gif_1, 'network_file': network_file, 'conv_file': conv_file})
    if ('redo_analysis' in request.POST and request.user.is_authenticated):
        if (request.POST['redo_analysis']):
            # configure loading page
            if (session_id_from_cache == 'has expired' or session_id_from_cache == ""):
                session_id = request.session._get_or_create_session_key()
            else:
                session_id = session_id_from_cache
            analysis_running = cache.get('analysis_running', 'none')
            if (analysis_running == 'none'):
                cache.set('analysis_running', 'analysis_running')
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("Your request is being processed...")
            add_loading_image.delay("none")
            add_loading_image.delay(session_id)
            # get expression file
            filename1 = request.POST.get("input_own_file_redo")
            fh1 = open(filename1)
            exprstr = fh1.read()
            # get other filenames and files
            filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
            fh2 = open(filename2)
            filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
            # see if clinical data were given, and load them
            has_clin_data = "false"
            clinicaldf = ""
            survival_col_name = ""
            if ('survival_col' in request.POST):
                if (request.POST['survival_col']):
                    survival_col_name = request.POST['survival_col']
            if (os.path.isfile(filename3)):
                fh3 = open(filename3)
                has_clin_data = "true"
                clinicalstr = fh3.read()
                # get survival column name (standard name assigned)
                if ("SURVIVAL_COLUMN_MONTH" in clinicalstr):
                    survival_col_name = "SURVIVAL_COLUMN_MONTH"
                elif ("SURVIVAL_COLUMN" in clinicalstr):
                    survival_col_name = "SURVIVAL_COLUMN"
                clinical_stringio = StringIO(clinicalstr)
                clinicaldf = pd.read_csv(clinical_stringio)
            ppistr = fh2.read()
            if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
                make_empty_figure.delay("none")
                make_empty_figure.delay(session_id)
                lgmin = int(request.POST['L_g_min'])
                lgmax = int(request.POST['L_g_max'])
                if (request.POST['L_g_min'] == ""):
                    lgmin = 10
                if (request.POST['L_g_max'] == ""):
                    lgmax = 20
                ## assign standard result size
                if (gene_set_size == "" or not str(gene_set_size).isdigit()):
                    gene_set_size = 2000
                session_id_from_cache = cache.get('session_id', 'has expired')
                # run algorithm
                if not (has_clin_data == "true"):
                    survival_col_name = ""
                    clinicalstr = "empty"
                    ret_metadata = ""
                result10 = preprocess_file_2.delay(exprstr)
                (exprstr_2, nbr_groups) = result10.get()
                # check if clinical data exist
                try:
                    # run algorithm
                    result1 = algo_output_task.delay(1, lgmin, lgmax, exprstr, ppistr, nbr_iter, nbr_ants, evap,
                                                     epsilon, hi_sig, pher_sig, session_id, gene_set_size, nbr_groups)
                    (T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1, group1_ids, group2_ids, jac_1,
                     jac_2) = result1.get()
                    # make plots and process results
                    result2 = script_output_task.delay(T, row_colors, col_colors, G2, means, genes_all, adjlist, genes1,
                                                       group1_ids, group2_ids, clinicalstr, jac_1, jac_2,
                                                       survival_col_name, clinicaldf, session_id)
                    (ret_metadata, path_heatmap, path_metadata, output_plot_path, json_path, p_val) = result2.get()
                except:
                    return render(request, 'clustering/errorpage.html',
                                  {'errors': "An error occurred during the algorithm run.",
                                   'hide_standard_message': "true"})

                # get metadata
                ret_metadata1 = ret_metadata[0]
                ret_metadata2 = ret_metadata[1]
                ret_metadata3 = ret_metadata[2]
                # set enrichment
                enrichment_dict = cache.get('enrichment_dict', "")
                if not (enrichment_dict == ""):
                    cache.set("enrichment_dict", "")
                    cache.set("enrichment_dict_2", "")
                    cache.set("enrichment_dict_3", "")
                    cache.set("enrichment_dict_4", "")
                    cache.set("enrichment_dict_5", "")
                output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
                json_path = "userfiles/ppi_" + session_id + ".json"
                path_heatmap = "userfiles/heatmap_" + session_id + ".png"
                path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
                if (has_clin_data == "true"):
                    if (os.path.isfile(path_metadata)):
                        metd = list_metadata_from_file.apply_async(args=[path_metadata], countdown=0)
                        (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
                cache.clear()
                make_empty_figure.apply_async(args=["none"], countdown=10)
                make_empty_figure.apply_async(args=[session_id], countdown=10)
                empty_log_file.apply_async(args=["none"], countdown=10)
                empty_log_file.apply_async(args=[session_id], countdown=10)
                list_of_files = ""
                list_of_files_2 = ""
                if (request.user.is_authenticated):
                    username = str(request.user)
                    list_of_files = GraphForm.list_user_data_2(username)
                    list_of_files_2 = GraphForm.list_user_data(username)
                cache.set('session_id', session_id)
                cache.set('ret_metadata1', ret_metadata1)
                cache.set('ret_metadata2', ret_metadata2)
                cache.set('ret_metadata3', ret_metadata3)
                cache.set('p_val', p_val)
                has_survival_plot = ""
                if (output_plot_path == "empty"):
                    cache.set("has_survival_plot", "false")
                    has_survival_plot = "false"
                else:
                    cache.set("has_survival_plot", "true")
                    has_survival_plot = "true"
                remove_loading_image.delay("none")
                remove_loading_image.delay(session_id)
                output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
                progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
                loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
                network_file = "userfiles/ntw_" + session_id + ".png"
                conv_file = "userfiles/conv_" + session_id + ".png"
                return render(request, 'clustering/clustering.html',
                              {'path_heatmap': path_heatmap, 'output_plot_path': output_plot_path,
                               'json_path': json_path, 'list_of_files': list_of_files, 'ret_dat': "",
                               'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                               'ret_metadata3': ret_metadata3, 'list_of_files_2': list_of_files_2,
                               'has_survival_plot': has_survival_plot, 'output_console_file': output_console_file,
                               'progress_file': progress_file, 'loading_gif_1': loading_gif_1,
                               'network_file': network_file, 'conv_file': conv_file})
    elif ('enrichment_type' in request.POST):
        enr_type = request.POST.get("enrichment_type")
        group_for_enr = "both"
        # get p-value cutoff
        if ('pval_enr' in request.POST):
            pval_enr = request.POST.get('pval_enr')
        if ('group_for_enr' in request.POST):
            group_for_enr = request.POST.get('group_for_enr')
        analysis_running = cache.get('analysis_running', 'none')
        # set analysis running parameter to allow display of "loading"-gif + text
        if (analysis_running == 'none'):
            cache.set('analysis_running', 'analysis_running')
        has_survival_plot = ""
        surv_from_cache = cache.get('has_survival_plot', 'none')
        print(surv_from_cache)
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        enrichment_dict = {}
        enrichment_dict_2 = {}
        enrichment_dict_3 = {}
        enrichment_dict_4 = {}
        enrichment_dict_5 = {}
        list_of_files = ""
        list_of_files_2 = ""
        if request.user.is_authenticated:
            username = str(request.user)
            list_of_files = GraphForm.list_user_data_2(username)
            list_of_files_2 = GraphForm.list_user_data(username)
        # get session ID
        session_id_from_cache = cache.get('session_id', 'has expired')
        if not (session_id_from_cache == 'has expired'):
            session_id = session_id_from_cache
        else:
            session_id = request.session._get_or_create_session_key()
        conv_file = "userfiles/conv_" + session_id + ".png"
        try:
            genelist = "/code/clustering/static/userfiles/genelist_" + session_id + ".txt"
            genelist1 = "/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"
            genelist2 = "/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"
            if not (os.path.isfile(genelist) and os.path.isfile(genelist1) and os.path.isfile(genelist2)):
                print("gene list not found")
                return render(request, 'clustering/clustering.html', {'errstr': ""})
            if (enr_type == "kegg_enrichment"):

                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                # both groups
                enrichment_dict = result4.get()
                # group 1, group 2
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                # results "only in group 1/2"
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_enrichment"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_molecular"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "reactome_enrichment"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_reactome/" + session_id + "/Reactome_2016.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
            # give links to result files from last analysis
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            json_path = "userfiles/ppi_" + session_id + ".json"
            ret_metadata1 = cache.get('ret_metadata1', "")
            ret_metadata2 = cache.get('ret_metadata2', "")
            ret_metadata3 = cache.get('ret_metadata3', "")
            path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
            if (ret_metadata1 == "" and os.path.isfile(path_metadata)):
                metd = list_metadata_from_file.apply_async(args=[path_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
            # write enrichment results to cache
            if not (session_id_from_cache == 'has expired'):
                session_id = session_id_from_cache
                cache.set("enrichment_dict", enrichment_dict)
                cache.set("enrichment_dict_2", enrichment_dict_2)
                cache.set("enrichment_dict_3", enrichment_dict_3)
                cache.set("enrichment_dict_4", enrichment_dict_4)
                cache.set("enrichment_dict_5", enrichment_dict_5)
            output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
            progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
            loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
            network_path = "userfiles/ntw_" + session_id_from_cache + ".png"
            return render(request, 'clustering/clustering.html',
                          {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                           'path_heatmap': path_heatmap, 'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                           'ret_metadata3': ret_metadata3,
                           'json_path': (json_path + "?foo=bar"), 'output_plot_path': output_plot_path,
                           'network_path': network_path, 'metadata_dict': metadata_dict,
                           'enrichment_dict': enrichment_dict, 'enrichment_dict_2': enrichment_dict_2,
                           'enrichment_dict_3': enrichment_dict_3, 'enrichment_dict_4': enrichment_dict_4,
                           'enrichment_dict_5': enrichment_dict_5, 'enrichment_open': "true",
                           'has_survival_plot': has_survival_plot, 'output_console_file': output_console_file,
                           'progress_file': progress_file, 'loading_gif_1': loading_gif_1, 'conv_file': conv_file})
        except:
            return render(request, 'clustering/errorpage.html',
                          {'errors': "An error occurred during the algorithm run.", 'hide_standard_message': "true"})
        if ('enr' not in request.POST):
            mutable = request.POST._mutable
            request.POST._mutable = True
            request.POST['enr'] = "true"
            request.POST._mutable = mutable
        has_survival_plot = ""
        surv_from_cache = cache.get('has_survival_plot', 'none')
        print(surv_from_cache)
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        # get session id from cache
        session_id_from_cache = cache.get('session_id', 'has expired')
        if not (session_id_from_cache == 'has expired' or session_id_from_cache == ""):
            path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
            json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
            path_metadata = "userfiles/metadata_" + session_id_from_cache + ".txt"
            output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
            ret_metadata1 = cache.get('ret_metadata1', 'none')
            ret_metadata2 = cache.get('ret_metadata2', 'none')
            ret_metadata3 = cache.get('ret_metadata3', 'none')
        else:
            session_id = request.session._get_or_create_session_key()
            # create new session if session id does not exist in cache
            path_heatmap = "userfiles/heatmap_" + session_id + ".png"
            json_path = "userfiles/ppi_" + session_id + ".json"
            output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
        network_file = ""
        if (os.path.isfile("/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png")):
            network_file = "/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png"
        output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
        progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
        loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
        network_path = "userfiles/ntw_" + session_id_from_cache + ".png"
        return render(request, 'clustering/clustering.html',
                      {'list_of_files': list_of_files, 'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2,
                       'ret_metadata3': ret_metadata3, 'network_path': network_path, 'path_heatmap': path_heatmap,
                       'json_path': json_path, 'output_plot_path': output_plot_path, 'enrichment_dict': enrichment_dict,
                       'enrichment_dict_2': enrichment_dict_2, 'enrichment_dict_3': enrichment_dict_3,
                       'enrichment_dict_4': enrichment_dict_4, 'enrichment_dict_5': enrichment_dict_5,
                       'enrichment_open': "true", 'has_survival_plot': has_survival_plot,
                       'output_console_file': output_console_file, 'progress_file': progress_file,
                       'loading_gif_1': loading_gif_1, 'network_file': network_file, 'conv_file': conv_file})
    else:
        # check if analysis is running (and user has hit reload)
        analysis_running = cache.get('analysis_running', 'none')
        session_id_from_cache = cache.get('session_id', 'has_expired')
        # start new session and set session ID in cache if none exists
        if (session_id_from_cache == 'has_expired' or session_id_from_cache == ""):
            session_id = request.session._get_or_create_session_key()
            cache.set('session_id', session_id)
            print("session ID: " + str(session_id))
        else:
            session_id = session_id_from_cache
        # if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
        if (analysis_running == 'none'):
            print("no analysis running")
            if (os.path.isfile("clustering/static/loading_1.gif")):
                os.unlink("clustering/static/loading_1.gif")
            if (os.path.isfile("clustering/static/loading_1_" + session_id + ".gif")):
                os.unlink("clustering/static/loading_1_" + session_id + ".gif")
            if (os.path.isfile("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")):
                os.unlink("/code/clustering/static/userfiles/loading_1_" + session_id + ".gif")
            if (os.path.isfile("clustering/static/progress.png")):
                os.unlink("clustering/static/progress.png")
            if (os.path.isfile("/code/clustering/static/progress.png")):
                os.unlink("/code/clustering/static/progress.png")
            if (os.path.isfile("/code/clustering/static/userfiles/progress_" + session_id + ".png")):
                os.unlink("/code/clustering/static/userfiles/progress_" + session_id + ".png")
            remove_loading_image.delay("none")
            remove_loading_image.delay(session_id)
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            with open("clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
        ret_metadata = ""
        # check if session already exists for current user (e.g. when user has hit the reload button)
        # get session ID and create list of previously saved files if user is authenticated
        if (request.user.is_authenticated):
            username = str(request.user)
            list_of_files = GraphForm.list_user_data_2(username)
            list_of_files_2 = GraphForm.list_user_data(username)
        ret_metadata1 = ""
        ret_metadata2 = ""
        ret_metadata3 = ""
        has_survival_plot = ""
        # check if last analysis (currently displayed) has survival plot
        surv_from_cache = cache.get('has_survival_plot', 'none')
        print("survival from cache" + str(surv_from_cache))
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        # display results from most recent analysis
        if not (session_id_from_cache == 'has_expired' or session_id_from_cache == ""):
            # take result files from storage
            path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
            json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
            if (os.path.isfile(json_path)):
                network_path = "userfiles/ntw_" + session_id_from_cache + ".png"
            else:
                network_path = "none"
            # has_survival_plot = "false"
            output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
            conv_file = "userfiles/conv_" + session_id + ".png"
            # take metadata from cache
            ret_metadata1 = cache.get('ret_metadata1', "")
            ret_metadata2 = cache.get('ret_metadata2', "")
            ret_metadata3 = cache.get('ret_metadata3', "")
            p_val = cache.get('p_val', "")
            # get enrichment results (in cache if user has hit reload button after running analysis)
            enrichment_dict = cache.get('enrichment_dict', "")
            enrichment_dict_2 = {}
            enrichment_dict_3 = {}
            enrichment_dict_4 = {}
            enrichment_dict_5 = {}
            if not (enrichment_dict == ""):
                enrichment_dict_2 = cache.get('enrichment_dict_2', "")
                enrichment_dict_3 = cache.get('enrichment_dict_3', "")
                enrichment_dict_4 = cache.get('enrichment_dict_4', "")
                enrichment_dict_5 = cache.get('enrichment_dict_5', "")
            cache.clear()
            # if(surv_from_cache == "false"):
            #	has_survival_plot = "false"
            #	cache.set("has_survival_plot","false")
            # set cache variables after clearing cache
            cache.set('session_id', session_id_from_cache)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            cache.set('p_val', p_val)
            cache.set('analysis_running', analysis_running)
            cache.set('has_survival_plot', has_survival_plot)
            surv_from_cache = cache.get('has_survival_plot', 'none')
            print("survival from cache" + str(surv_from_cache))
            if (surv_from_cache == "true"):
                print("survival from cache is true")
            output_console_file = "userfiles/output_console_" + session_id_from_cache + ".txt"
            progress_file = "userfiles/progress_" + session_id_from_cache + ".png"
            loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
            return render(request, 'clustering/clustering.html',
                          {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                           'ret_metadata': ret_metadata, 'path_heatmap': path_heatmap, 'json_path': json_path,
                           'output_plot_path': output_plot_path, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'metadata_dict': metadata_dict, 'enrichment_dict': enrichment_dict,
                           'enrichment_dict_2': enrichment_dict_2, 'enrichment_dict_3': enrichment_dict_3,
                           'enrichment_dict_4': enrichment_dict_4, 'enrichment_dict_5': enrichment_dict_5,
                           'p_val': p_val, 'has_survival_plot': has_survival_plot,
                           'output_console_file': output_console_file, 'progress_file': progress_file,
                           'loading_gif_1': loading_gif_1, 'network_path': network_path, 'conv_file': conv_file})
        cache.clear()
        # set cache variables after clearing cache
        cache.set('session_id', session_id)
        cache.set('analysis_running', analysis_running)
        cache.set('has_survival_plot', has_survival_plot)
        surv_from_cache = cache.get('has_survival_plot', 'none')
        print("survival from cache" + str(surv_from_cache))
        if (surv_from_cache == "true"):
            print("survival from cache is true")
        conv_file = "userfiles/conv_" + session_id + ".png"
        output_console_file = "userfiles/output_console_" + session_id + ".txt"
        progress_file = "userfiles/progress_" + session_id + ".png"
        loading_gif_1 = "userfiles/loading_1_" + session_id + ".gif"
        network_path = "none"
        if (os.path.isfile("/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png")):
            network_path = "/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png"
        return render(request, 'clustering/clustering.html',
                      {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2, 'ret_metadata': ret_metadata,
                       'ret_metadata1': ret_metadata1, 'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                       'metadata_dict': metadata_dict, 'enrichment_dict': "", 'has_survival_plot': has_survival_plot,
                       'output_console_file': output_console_file, 'progress_file': progress_file,
                       'loading_gif_1': loading_gif_1, 'network_path': network_path, 'conv_file': conv_file})


def errorpage(request):
    errors = ""
    if ('errors' in request.session):
        errors = request.session['errors']
    return render(request, 'clustering/errorpage.html', {'errors': errors})


# method for generating the result page belonging to clustering_6_part_4_2
# def clustering_6_part_3_2(request):
def clustering_step_2(request):
    ret_metadata1 = {}
    ret_metadata2 = {}
    ret_metadata3 = {}
    metadata_dict = []
    enrichment_dict = []
    pval_enr = 0.5
    list_of_files = ""
    list_of_files_2 = ""
    save_data = request.POST.get("save_data", None)
    gene_set_size = request.POST.get("gene_set_size", 2000)
    nbr_iter = request.POST.get("nbr_iter", 45)
    nbr_ants = request.POST.get("nbr_ants", 30)
    evap = request.POST.get("evap", 0.3)
    epsilon = request.POST.get("stopcr", 0.02)
    hi_sig = request.POST.get("hisig", 1)
    pher_sig = request.POST.get("pher", 1)
    if request.user.is_authenticated:
        username = str(request.user)
        list_of_files = GraphForm.list_user_data_2(username)
    if ('enrichment_type' in request.POST):
        session_id_from_cache = cache.get('session_id', 'has expired')
        if not (session_id_from_cache == 'has expired'):
            session_id = session_id_from_cache
        else:
            session_id = request.session._get_or_create_session_key()
        try:
            genelist = "/code/clustering/static/userfiles/genelist_" + session_id + ".txt"
            genelist1 = "/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"
            genelist2 = "/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"
            if not (os.path.isfile(genelist) and os.path.isfile(genelist1) and os.path.isfile(genelist2)):
                return render(request, 'clustering/clustering_step_2.html', {'errstr': ""})
            if (enr_type == "kegg_enrichment"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_kegg/" + session_id),
                                               ['KEGG_2016', 'KEGG_2013'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                # both groups
                enrichment_dict = result4.get()
                # group 1, group 2
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                # results "only in group 1/2"
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_enrichment"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "go_molecular"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_go/" + session_id),
                                               ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018',
                                                'GO_Molecular_Function_2018'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result7 = read_enrichment_2.delay((
                                                              "/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                  (
                                                              "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                  pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
                (enrichment_dict_4, enrichment_dict_5) = result7.get()
            elif (enr_type == "reactome_enrichment"):
                result1 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                result2 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test2/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                result3 = run_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),
                                               pval_enr, ("/code/clustering/data/test3/enrichr_reactome/" + session_id),
                                               ['Reactome_2013', 'Reactome_2016'])
                enr_results = result1.get()
                enr_results_2 = result2.get()
                enr_results_3 = result3.get()
                print("enr")
                result4 = read_enrichment.delay((
                                                            "/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result5 = read_enrichment.delay((
                                                            "/code/clustering/data/test2/enrichr_reactome/" + session_id + "/Reactome_2016.test_name.enrichr.reports.txt"),
                                                pval_enr)
                result6 = read_enrichment.delay((
                                                            "/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),
                                                pval_enr)
                enrichment_dict = result4.get()
                enrichment_dict_2 = result5.get()
                enrichment_dict_3 = result6.get()
            # give links to result files from last analysis
            path_heatmap = "heatmap_" + session_id + ".png"
            json_path = "ppi_" + session_id + ".json"
            path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
            conv_file = "userfiles/conv_" + session_id + ".png"
            if (os.path.isfile(path_metadata)):
                metd = list_metadata_from_file.apply_async(args=[path_metadata], countdown=0)
                (ret_metadata1, ret_metadata2, ret_metadata3) = metd.get()
            output_plot_path = "output_plotly_" + session_id + ".html"
            # write enrichment results to cache
            if not (session_id_from_cache == 'has expired'):
                session_id = session_id_from_cache
                cache.set("enrichment_dict", enrichment_dict)
                cache.set("enrichment_dict_2", enrichment_dict_2)
                cache.set("enrichment_dict_3", enrichment_dict_3)
                cache.set("enrichment_dict_4", enrichment_dict_4)
                cache.set("enrichment_dict_5", enrichment_dict_5)
            return render(request, 'polls/clustering_step_2.html',
                          {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                           'path_heatmap': ("userfiles/" + path_heatmap), 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'json_path': ("userfiles/" + json_path + "?foo=bar"),
                           'output_plot_path': ("userfiles/" + output_plot_path), 'metadata_dict': metadata_dict,
                           'enrichment_dict': enrichment_dict, 'enrichment_dict_2': enrichment_dict_2,
                           'enrichment_dict_3': enrichment_dict_3, 'enrichment_dict_4': enrichment_dict_4,
                           'enrichment_dict_5': enrichment_dict_5, 'enrichment_open': "true",
                           'has_survival_plot': "true", 'conv_file': conv_file})
        except:
            return render(request, 'clustering/errorpage.html',
                          {'errors': "An error occurred during the algorithm run.", 'hide_standard_message': "true"})
        return render(request, 'polls/clustering_step_2.html',
                      {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2, 'has_survival_plot': "true"})
    else:
        analysis_running = cache.get('analysis_running', 'none')
        # if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
        if (analysis_running == 'none'):
            if (os.path.isfile("clustering/static/loading_1.gif")):
                os.unlink("clustering/static/loading_1.gif")
            if (os.path.isfile("clustering/static/progress.png")):
                os.unlink("clustering/static/progress.png")
            remove_loading_image.delay("none")
            # print("removed loading image")
            with open("/code/clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
            with open("clustering/static/output_console.txt", "w") as text_file:
                text_file.write("")
        ret_metadata = ""
        # check if session already exists for current user (e.g. when user has hit the reload button)
        session_id_from_cache = cache.get('session_id', 'has_expired')
        # get session ID and create list of previously saved files if user is authenticated
        if (request.user.is_authenticated):
            list_of_files = GraphForm.list_user_data_2(username)
            list_of_files_2 = GraphForm.list_user_data(username)
        ret_metadata1 = ""
        ret_metadata2 = ""
        ret_metadata3 = ""
        # display results from most recent analysis
        if not (session_id_from_cache == "has_expired"):
            # take result files from storage
            path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
            json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
            output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
            # take metadata from cache
            ret_metadata1 = cache.get('ret_metadata1', "")
            ret_metadata2 = cache.get('ret_metadata2', "")
            ret_metadata3 = cache.get('ret_metadata3', "")
            p_val = cache.get('p_val', "")
            # get enrichment results (in cache if user has hit reload button after running analysis)
            enrichment_dict = cache.get('enrichment_dict', "")
            enrichment_dict_2 = {}
            enrichment_dict_3 = {}
            enrichment_dict_4 = {}
            enrichment_dict_5 = {}
            if not (enrichment_dict == ""):
                enrichment_dict_2 = cache.get('enrichment_dict_2', "")
                enrichment_dict_3 = cache.get('enrichment_dict_3', "")
                enrichment_dict_4 = cache.get('enrichment_dict_4', "")
                enrichment_dict_5 = cache.get('enrichment_dict_5', "")
            surv_from_cache = cache.get('has_survival_plot', 'none')
            cache.clear()
            if (surv_from_cache == "false"):
                has_survival_plot = "false"
            elif (surv_from_cache == "true"):
                has_survival_plot = "true"
            cache.set('session_id', session_id_from_cache)
            cache.set('ret_metadata1', ret_metadata1)
            cache.set('ret_metadata2', ret_metadata2)
            cache.set('ret_metadata3', ret_metadata3)
            cache.set('p_val', p_val)
            network_path = "none"
            conv_file = "userfiles/conv_" + session_id + ".png"
            if (os.path.isfile("/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png")):
                network_path = "/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png"
            return render(request, 'clustering/clustering_step_2.html',
                          {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2,
                           'ret_metadata': ret_metadata, 'path_heatmap': path_heatmap, 'json_path': json_path,
                           'output_plot_path': output_plot_path, 'ret_metadata1': ret_metadata1,
                           'ret_metadata2': ret_metadata2, 'ret_metadata3': ret_metadata3,
                           'metadata_dict': metadata_dict, 'enrichment_dict': enrichment_dict,
                           'enrichment_dict_2': enrichment_dict_2, 'enrichment_dict_3': enrichment_dict_3,
                           'enrichment_dict_4': enrichment_dict_4, 'enrichment_dict_5': enrichment_dict_5,
                           'p_val': p_val, 'has_survival_plot': "true", 'network_path': network_path,
                           'conv_file': conv_file})
        session_id = request.session._get_or_create_session_key()
        path99 = "heatmap_" + session_id + ".png"
        json_path = "ppi_" + session_id + ".json"
        path_metadata = "metadata_" + session_id + ".txt"
        output_plot_path = "output_plotly_" + session_id + ".html"
        surv_from_cache = cache.get('has_survival_plot', 'none')
        cache.clear()
        cache.set('analysis_running', analysis_running)
        conv_file = "userfiles/conv_" + session_id + ".png"
        network_path = "none"
        print("survival from cache" + str(surv_from_cache))
        if (surv_from_cache == "false"):
            has_survival_plot = "false"
        elif (surv_from_cache == "true"):
            has_survival_plot = "true"
        if (os.path.isfile("/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png")):
            network_path = "/code/clustering/static/userfiles/ntw_" + session_id_from_cache + ".png"
        return render(request, 'clustering/clustering_step_2.html',
                      {'list_of_files': list_of_files, 'list_of_files_2': list_of_files_2, 'path_heatmap': path99,
                       'output_plot_path': output_plot_path, 'metadata_dict': metadata_dict,
                       'enrichment_dict': enrichment_dict, 'has_survival_plot': "true", 'network_path': network_path,
                       'has_survival_plot': has_survival_plot, 'conv_file': conv_file})


def infopage(request):
    return render(request, 'clustering/infopage.html')


def sources(request):
    return render(request, 'clustering/sources.html')
