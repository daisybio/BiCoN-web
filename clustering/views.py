from django.shortcuts import render
from io import StringIO
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
from django.contrib.auth.hashers import check_password
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.views import generic
from django.core.cache import cache
from datetime import datetime
#from .models import Choice, Question
#from sklearn.ensemble import RandomForestClassifier
from django.contrib.auth.models import User
from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.urls import reverse
from django.contrib.auth import authenticate, login
from shutil import copyfile
import shutil
#### own imports
import clustering
from clustering.models import Upload,UploadForm,GraphForm
from django.contrib.auth import authenticate, login, logout
from clustering.tasks import make_empty_figure,algo_output_task,empty_log_file,add_loading_image,remove_loading_image,show_old_data,import_ndex,script_output_task_9,read_kegg_enrichment,read_ndex_file_4,run_enrichment_2,run_go_enrichment_2,run_reac_enrichment,read_kegg_enrichment_2,convert_gene_list,check_input_files,script_output_task_10,preprocess_file,write_pval,algo_output_task_new,list_metadata_5,preprocess_clinical_file,preprocess_ppi_file,preprocess_file_2,algo_output_task_2,algo_output_task_3
from django.core.cache import cache
import os.path
from django.core.mail import send_mail

### *ACTUAL* imports (that have dependencies other than django and my own stuff) ####
import networkx as nx
#from biomart import BiomartServer
#from bokeh.palettes import Spectral4
from pybiomart import Dataset
import seaborn as sns
import pandas as pd
from numpy import array
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from plotly.offline import plot_mpl
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline


def logout_2(request):
	if(request.method == "POST"):
		logout(request)
		#return redirect('polls/logout.html')
		return redirect('clustering/logout.html')
	else:
		return render(request,'clustering/logout.html')

def login_2(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		user = authenticate(request, username=username, password=password)
		if user is not None:
			login(request, user)
			#return render(request,'polls/login.html')
			return redirect('clustering/clustering_6_part_4.html')
		else:
			text = "Username or password are incorrect"
			return render(request,'clustering/login.html',{'text':text})
			#return redirect('polls/clustering.html')
	else:
		return render(request,'clustering/login.html')
		#return redirect('polls/clustering.html')

def signup(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		email = ""
		if('email' in request.POST):
			email = request.POST['email']
		if(User.objects.filter(username=username).exists()):
			text = "Username already exists. Please choose another username!"
			return render(request,'clustering/signup.html',{'text':text})
		else:
			user = User.objects.create_user(username, email, password)
			user.save()
			current_site = get_current_site(request)
			mail_subject = 'Activate your account.'
			#message = render_to_string('acc_active_email.html', {
			#'user': user,
			#'domain': current_site.domain,
			#'uid':urlsafe_base64_encode(force_bytes(user.pk)),
			#'token':account_activation_token.make_token(user),
			#})
			userdir = "user_uploaded_files/" + username
			if not(os.path.isdir(userdir)):
				os.mkdir(userdir)
			to_email = email
			#email_message = EmailMessage(
			#mail_subject, message, to=[to_email]
			#)
			#email_message.send()
			#send_mail('Account activation', 'Your account was activated.', 'sender@example.com', ['receiver1@example.com',])
			text = "Account is being created. You will receive a confirmation e-mail soon!"
			return render(request,'clustering/signup.html',{'text':text,'new_user':"true"})
			#return redirect('polls/clustering.html')
	else:
		text = "Please input username and password!"
		return render(request,'clustering/signup.html',{'text':text})
		#return redirect('polls/clustering.html')


def delete_user(request):
	if('username' in request.POST and 'password' in request.POST and request.user.is_authenticated):
		username = request.POST['username']
		password = request.POST['password']
		print(str(request.user))
		print(username)
		print(str(request.user.password))
		print(password)
		if(str(request.user) == username and check_password(password,request.user.password)):
			print("password found")
			u = User.objects.get(username=username)
			u.delete()
			user_dir = "user_uploaded_files/" + username
			if (os.path.isdir(user_dir)):
				shutil.rmtree(user_dir)
			text = "Your account was deleted."
			return render(request,'clustering/delete_user.html',{'text':text,'deleted':"true"})
		return render(request,'clustering/delete_user.html',{'text':""})
			#return redirect('polls/clustering.html')
	else:
		text = "Please input username and password!"
		return render(request,'clustering/delete_user.html',{'text':text})
		#return redirect('polls/clustering.html')



def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		errors_from_cache = cache.get('errors', 'has expired')
		if not(errors_from_cache == ""):
			errors = errors_from_cache
			cache.set('errors','')
		return render(request,'clustering/errorpage.html',{'errors':errors})


#########################################################################
## This method displays everything on one page with no Session IDs     ##
#########################################################################



def clustering_6_new(request):
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	gene_set_size = request.POST.get("gene_set_size",2000)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	if('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
		if(request.POST['input_own_file'] and request.POST['display_old_results']):
			# configure loading page
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			make_empty_figure.delay()
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
			#json_path = "ppi_" + session_id + ".json"
			#path_heatmap_2 = "heatmap_" + session_id + ".png"
			path_metadata_2 = "userfiles/metadata.txt"
			path_plotly_2 = "userfiles/output_plotly.html"
			# copy files to static directory
			copyfile(path_json,("clustering/static/" + json_path))	
			copyfile(path_heatmap,("clustering/static/" + path_heatmap_2))
			copyfile(path_genelist,("clustering/static/userfiles/genelist.txt"))
			copyfile(path_genelist_1,("clustering/static/userfiles/genelist_1.txt"))
			copyfile(path_genelist_2,("clustering/static/userfiles/genelist_2.txt"))
			output_plot_path_2 = ""
			ret_metadata_1 = ""
			ret_metadata_2 = ""
			ret_metadata_3 = ""
			# check if plotly file exists and copy
			if(os.path.isfile(path_plotly)):
				copyfile(path_plotly,("clustering/static/" + path_plotly_2))
				output_plot_path_2 = path_plotly_2
				print("plot copied to")
				print(path_plotly)
				print(output_plot_path_2)
			# read metadata (must copy file to shared volume for processing via celery)
			#if(os.path.isfile(path_metadata+"_2")):
			if(os.path.isfile(path_metadata)):
				print("found metadata")
				print(path_metadata)
				copyfile(path_metadata,("/code/clustering/static/metadata.txt"))
				filename_for_old_metadata = "/code/clustering/static/metadata.txt"
				print(filename_for_old_metadata)
				#metd = list_metadata_4.apply_async(args=[filename_for_old_metadata],countdown=0)
				metd = list_metadata_5.apply_async(args=[filename_for_old_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print(ret_metadata1)
			cache.clear()
			# set session ID in cache
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			make_empty_figure.apply_async(countdown=10)
			empty_log_file.apply_async(countdown=10)
			# list old files
			list_of_files = ""
			list_of_files_2 = ""
			if request.user.is_authenticated:
				username = str(request.user)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)              				
			return render(request, 'clustering/clustering_6.html', {'form':"",'images':"",'plot_div':"",'script':"",'path_heatmap':path_heatmap_2,'output_plot_path':output_plot_path_2,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	
	elif(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		# check if input files exist
		input_valid = "false"
		if('myfile' in request.FILES and 'protfile' in request.FILES):
			if(request.FILES['myfile'] and request.FILES['protfile']):
				input_valid = "true"
		elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'protfile' in request.FILES):
			if(request.POST['predef_file'] and request.FILES['protfile']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		#if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
		if(input_valid == "true"):
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
				if(request.POST['L_g_min'] != "" and request.POST['L_g_max'] != ""):
					lgmin = int(request.POST['L_g_min'])
					lgmax = int(request.POST['L_g_max'])
				else:
					lgmin = 10
					lgmax = 20
			else:
				lgmin = 10
				lgmax = 20
			#if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(save_data in ["save_data"]):
				if request.user.is_authenticated:
					print("saving data is true")
			#lgmin = int(request.POST['L_g_min'])
			#lgmax = int(request.POST['L_g_max'])
			# assign standard result size
			if(request.POST['L_g_min'] == ""):
				lgmin = 10
			if(request.POST['L_g_max'] == ""):
				lgmax = 20
			clinicalstr = ""
			clinicaldf = ""
			# configure loading page
			add_loading_image.delay()
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("Your request is being processed...")
				text_file.close()
			make_empty_figure.delay()
			clinicalstr = "empty"
			clinicaldf = ""
			survival_col_name = ""
			# read expression file
			if('myfile' in request.FILES):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result10 = preprocess_file_2.delay(exprstr)
				(exprstr,nbr_col) = result10.get()
				#result10 = preprocess_file.delay(exprstr)
				#exprstr = result10.get()
			# read predefined expression file and clinical data
			elif('predef_file' in request.POST and 'cancer_type' in request.POST):
				cancer_type = request.POST.get("cancer_type")
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("clustering/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
					fh4 = open("clustering/data/lung_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("clustering/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
					fh4 = open("clustering/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
			# read PPI file
			if('protfile' in request.FILES):
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				result3 = preprocess_ppi_file.delay(ppistr)
				ppistr = result3.get()
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'clustering/errorpage.html',{'errors':errstr})
			# read ndex file from web
			elif('ndex_name_2' in request.POST):
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					#result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
			# read metadata if given
			if('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
				if(request.FILES['patientdata']):
					clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
					clinicalstr_first_line = clinicalstr.split("\n")[0]
					if(len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
						print("is tsv")
						clinicalstr = clinicalstr.replace("\t",",")
					clinical_stringio = StringIO(clinicalstr)
					clinicaldf = pd.read_csv(clinical_stringio)
					if('survival_col' in request.POST):
						if(request.POST['survival_col']):
							survival_col_name = request.POST['survival_col']
			session_id = ""
			# assign standard value to gene set size
			if(gene_set_size == "" or not str(gene_set_size).isdigit()):
				gene_set_size = 2000
				# run algorithm and read results
			result1 = algo_output_task_3.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,gene_set_size,nbr_col)
			#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id)
			(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()
			#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size)
			#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id)
			#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()			
			# make plots and process results	
			result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
			(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
			output_plot_path = "output_plotly.html"
			json_path = "ppi.json"
			path_metadata = "/code/clustering/static/metadata.txt"
			path_heatmap = "heatmap.png"
			#json_path = "ppi_" + session_id + ".json"
			#path_heatmap = "heatmap_" + session_id + ".png"
			if(save_data in ["save_data"]):
				if request.user.is_authenticated:
					print("saving data in views.py")
					username = str(request.user)
					if not (survival_col_name == ""):
						if("month" in survival_col_name):
							clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN_MONTH",1)
						else:
							clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN",1)
					# save input data
					GraphForm.save_user_data_3(exprstr,ppistr,clinicalstr,username)
					curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]	
					# save output data
					copyfile(("/code/clustering/static/heatmap.png"),("user_uploaded_files/"+ username + "/" + curr_time + "_heatmap.png"))	
					copyfile(("/code/clustering/static/ppi.json"),("user_uploaded_files/"+ username + "/" + curr_time + "_json.json"))	
					copyfile( "/code/clustering/static/metadata.txt",("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt"))
					#if(os.path.isfile("metadata.txt" + "_2")):
					#	copyfile(("metadata.txt"+ "_2"),("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt_2"))
					copyfile(("/code/clustering/static/output_plotly.html"),("user_uploaded_files/"+ username + "/" + curr_time + "plotly.html"))
					copyfile(("/code/clustering/static/genelist.txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist.txt"))
					copyfile(("/code/clustering/static/genelist_1.txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_1.txt"))
					copyfile(("/code/clustering/static/genelist_2.txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_2.txt"))
			# read metadata
			ret_metadata1=ret_metadata[0]
			ret_metadata2=ret_metadata[1]
			ret_metadata3=ret_metadata[2]
			# empty enrichment data from cache
			enrichment_dict = cache.get('enrichment_dict',"")
			if not(enrichment_dict == ""):
				cache.set("enrichment_dict","")
				cache.set("enrichment_dict_2","")
				cache.set("enrichment_dict_3","")
				cache.set("enrichment_dict_4","")
				cache.set("enrichment_dict_5","")
			# paths for showing results
			# write list of genes to downloadable file
			convert_gene_list.delay(adjlist,"/code/clustering/static/genelist_temp.txt")
			# save uploaded files if specified
			# render list of previously uploaded files if user is logged in (needed if user submits another request)
			if request.user.is_authenticated:
				username = str(request.user)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)              
			# remove the loading-gif and progress image, clear cache             
			remove_loading_image.delay()
			#cache.clear()				
			make_empty_figure.apply_async(countdown=10)
			empty_log_file.apply_async(countdown=10)
			#if request.user.is_authenticated:
			#	request.session['done'] = "true"
			#remove_loading_image.delay()
			if(os.path.isfile("clustering/static/loading_1.gif")):
				os.unlink("clustering/static/loading_1.gif")
			cache.clear()				
			#make_empty_figure.apply_async(countdown=10)
			#empty_log_file.apply_async(countdown=10)
			# copy static files from shared directory to static-file-dir on web container
			copyfile(("/code/clustering/static/heatmap.png"),("clustering/static/userfiles/heatmap.png"))	
			copyfile(("/code/clustering/static/ppi.json"),("clustering/static/userfiles/ppi.json"))
			copyfile(("/code/clustering/static/output_plotly.html"),("clustering/static/userfiles/output_plotly.html"))
			copyfile(("/code/clustering/static/genelist.txt"),("clustering/static/userfiles/genelist.txt"))
			copyfile(("/code/clustering/static/genelist_1.txt"),("clustering/static/userfiles/genelist_1.txt"))
			copyfile(("/code/clustering/static/genelist_2.txt"),("clustering/static/userfiles/genelist_2.txt"))
			# save session ID and metadata in cache
			cache.set('session_id', session_id)	
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			cache.set('json_path', "ppi.json")	
			cache.set('p_val', p_val)
			cache.set('analysis_running','analysis_running')	
			if(clinicalstr == "empty"):
				output_plot_path = "empty"		
			return render(request, 'clustering/clustering_6.html', {'form':"",'images':"",'plot_div':"",'script':"",'plot2':"", 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val})
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("clustering/static/output_console.txt", "w") as text_file:
				text_file.write("Your request is being processed...")
			add_loading_image.delay()
			filename1 = request.POST.get("input_own_file_redo")
			fh1 = open(filename1)
			exprstr = fh1.read()
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
			has_clin_data = "false"
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
			#path_json = filename1
			#path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				## assign standard result size
				#if(request.POST['L_g_min'] == ""):
				#	lgmin = 10
				#if(request.POST['L_g_max'] == ""):
				#	lgmax = 20
				if(gene_set_size == "" or not str(gene_set_size).isdigit()):
					gene_set_size = 2000
				result1 = algo_output_task_3.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,gene_set_size)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				if not(has_clin_data == "true"):
					clinicalstr = "empty"
					ret_metadata = ""
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				#plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				return render(request, 'clustering/clustering_6.html', {'form':"",'images':"",'div':"",'script':"", 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if not(os.path.isfile("genelist.txt") and os.path.isfile("genelist_1.txt") and os.path.isfile("genelist_2.txt")):
				return render(request,'clustering/clustering_6.html',{'errstr':""})
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist.txt",pval_enr,"/code/clustering/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/clustering/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/clustering/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/clustering/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/clustering/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/clustering/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/clustering/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/clustering/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/clustering/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"/code/clustering/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		#print(request.POST.get("newAnalysis"))
		#print(request.POST['newAnalysis'])
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		return render(request,'clustering/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()	
			metd = list_metadata_5.apply_async(args=["/code/clustering/static/metadata.txt"],countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			metd = list_metadata_5.apply_async(args=["/code/clustering/static/metadata.txt"],countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
		return render(request,'clustering/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})




#########################################################################
#### version of the page with separate input form and result display ####
#########################################################################

def clustering_6_4_part_2(request):
	# check if user has clicked 'return' button on result page, then request.POST['newAnalysis'] is "true"
	if('newAnalysis' in request.POST):
		print(request.POST.get("newAnalysis"))
		print(request.POST['newAnalysis'])
		request.POST._mutable = True
		done_from_cache = cache.get("done","")
		print(done_from_cache)
		if(request.POST['newAnalysis'] != "false"):
			if('done' in request.session):
				if(request.session['done'] == "true"):
					#set done parameter to false if user has clicked return on result page
					request.session['done'] = "False"
			if(done_from_cache == "done" or done_from_cache == 'done'):
				print("done from cache")
				#set done parameter to false if user has clicked return on result page
				cache.set('done',"False")
			# remove parameter from request.POST to allow later switching to result page
			request.POST['newAnalysis'] = "false"
	#if(request.session['done'] == "true"):
	#if('done' in request.session):
	#	return HttpResponseRedirect('polls/clustering_6_part_3.html')
	#	print("done")
	if not(os.path.isdir("/code/clustering/static/userfiles")):
		os.mkdir("/code/clustering/static/userfiles")
	done_from_cache = cache.get("done","")
	print(done_from_cache)
	if(done_from_cache == "done"):
		print("done in cache")
		#return(clustering_6_part_3_2(request))
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
	gene_set_size = request.POST.get("gene_set_size",2000)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	#for key, value in request.POST.items():
	#	print(key)
	#	print(value)
	#if(('myfile' in request.POST) or ('predef_file' in request.POST)):
	#	print("in request")
	#if(('protfile' in request.POST) or (('parse_ndex_file' in request.POST) and ('ndex_name_2' in request.POST))):
	#	print("in request 2")
	# check if the user has uploaded or chosen an expression and a PPI file
	#if(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		# check if these files are not empty and exist
	#	if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
	if('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
		if(request.POST['input_own_file'] and request.POST['display_old_results']):
			# configure loading page
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			make_empty_figure.delay()
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
			session_id = request.session._get_or_create_session_key() 
			json_path = "userfiles/ppi_" + session_id + ".json"
			path_heatmap_2 = "userfiles/heatmap_" + session_id + ".png"
			#json_path = "ppi_" + session_id + ".json"
			#path_heatmap_2 = "heatmap_" + session_id + ".png"
			path_metadata_2 = "userfiles/metadata_" + session_id + ".txt"
			path_plotly_2 = "userfiles/output_plotly_" + session_id + ".html"
			# copy files to static directory
			copyfile(path_json,("clustering/static/" + json_path))	
			copyfile(path_heatmap,("clustering/static/" + path_heatmap_2))
			copyfile(path_genelist,("clustering/static/userfiles/genelist_" + session_id + ".txt"))
			copyfile(path_genelist_1,("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
			copyfile(path_genelist_2,("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
			output_plot_path_2 = ""
			ret_metadata_1 = ""
			ret_metadata_2 = ""
			ret_metadata_3 = ""
			# check if plotly file exists and copy
			if(os.path.isfile(path_plotly)):
				copyfile(path_plotly,("clustering/static/" + path_plotly_2))
				output_plot_path_2 = path_plotly_2
				#print("plot copied to")
				#print(path_plotly)
				#print(output_plot_path_2)
			# read metadata (must copy file to shared volume for processing via celery)
			#if(os.path.isfile(path_metadata+"_2")):
			if(os.path.isfile(path_metadata)):
				#print("found metadata")
				#print(path_metadata)
				copyfile(path_metadata,("/code/clustering/static/userfiles/metadata_" + session_id + ".txt"))
				filename_for_old_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				#copyfile(path_metadata,("/code/clustering/static/metadata_" + session_id + ".txt"))
				#filename_for_old_metadata = "/code/clustering/static/metadata_" + session_id + ".txt"
				#print(filename_for_old_metadata)
				#metd = list_metadata_4.apply_async(args=[filename_for_old_metadata],countdown=0)
				metd = list_metadata_5.apply_async(args=[filename_for_old_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#print(ret_metadata1)
			cache.clear()
			# set session ID in cache
			cache.set('session_id',session_id)
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			cache.set('done',"done")
			make_empty_figure.apply_async(countdown=10)
			empty_log_file.apply_async(countdown=10)
			# list old files
			list_of_files = ""
			list_of_files_2 = ""
			if request.user.is_authenticated:
				username = str(request.user)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)              				
			return render(request, 'clustering/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':"",'script':"",'path_heatmap':path_heatmap_2,'output_plot_path':output_plot_path_2,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'has_survival_plot':"true"})
	
	elif(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		# check if input files exist
		input_valid = "false"
		if('myfile' in request.FILES and 'protfile' in request.FILES):
			if(request.FILES['myfile'] and request.FILES['protfile']):
				input_valid = "true"
		elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'protfile' in request.FILES):
			if(request.POST['predef_file'] and request.FILES['protfile']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		#if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
		if(input_valid == "true"):
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				if(save_data in ["save_data"]):
					if request.user.is_authenticated:
						print("saving data is true")
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				## assign standard result size
				#if(request.POST['L_g_min'] == ""):
				#	lgmin = 10
				#if(request.POST['L_g_max'] == ""):
				#	lgmax = 20
				clinicalstr = ""
				clinicaldf = ""
				# configure loading page
				add_loading_image.delay()
				with open("/code/clustering/static/output_console.txt", "w") as text_file:
					text_file.write("Your request is being processed...")
					text_file.close()
				make_empty_figure.delay()
				clinicalstr = "empty"
				clinicaldf = ""
				survival_col_name = ""
				# read expression file
				if('myfile' in request.FILES):
					exprstr = request.FILES['myfile'].read().decode('utf-8')
					result10 = preprocess_file_2.delay(exprstr)
					(exprstr,nbr_col) = result10.get()
					#result10 = preprocess_file.delay(exprstr)
					#exprstr = result10.get()
				# read predefined expression file and clinical data
				elif('predef_file' in request.POST and 'cancer_type' in request.POST):
					cancer_type = request.POST.get("cancer_type")
					if(cancer_type == "1"):
						#print("babababababa")
						fh1 = open("clustering/data/lung_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
						fh4 = open("clustering/data/lung_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "disease free survival in months:ch1"
					else:
						fh1 = open("clustering/data/breast_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
						fh4 = open("clustering/data/breast_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "mfs (yr):ch1"
				# read PPI file
				if('protfile' in request.FILES):
					ppistr = request.FILES['protfile'].read().decode('utf-8')
					result3 = preprocess_ppi_file.delay(ppistr)
					ppistr = result3.get()
					result4 = check_input_files.delay(ppistr,exprstr)
					errstr = result4.get()
					if(errstr != ""):
						request.session['errors'] = errstr
						return render(request,'clustering/errorpage.html',{'errors':errstr})
				# read ndex file from web
				elif('ndex_name_2' in request.POST):
					ndex_file_id = request.POST.get("ndex_name_2")
					if(ndex_file_id == "1"):
						result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "2"):
						#result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "3"):
						result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "4"):
						result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
				# read metadata if given
				if('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
					if(request.FILES['patientdata']):
						clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
						clinicalstr_first_line = clinicalstr.split("\n")[1]
						if(len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
							#print("converting file to csv")
							clinicalstr = clinicalstr.replace("\t",",")
						clinical_stringio = StringIO(clinicalstr)
						clinicaldf = pd.read_csv(clinical_stringio)
						if('survival_col' in request.POST):
							if(request.POST['survival_col']):
								survival_col_name = request.POST['survival_col']
				#session_id = ""
				# start session for storing result data			
				#session_id = request.session._get_or_create_session_key()
				session_id = ""
				session_id_from_cache = cache.get("session_id","none")
				if(session_id_from_cache == "none" or session_id_from_cache == ""):
					# start session for storing result data			
					session_id = request.session._get_or_create_session_key()
				else:
					session_id = session_id_from_cache
				
				# assign standard value to gene set size
				if(gene_set_size == "" or not str(gene_set_size).isdigit()):
					gene_set_size = 2000
				# run algorithm and read results
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size)
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()	
				result1 = algo_output_task_2.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size,nbr_col)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()			
				# make plots and process results	
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(ret_metadata,path_heatmap,path_metadata,output_plot_path,json_path,p_val) = result2.get()
				has_survival_plot = "true"
				if(output_plot_path == "empty"):
					has_survival_plot = "false"
				output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
				json_path = "userfiles/ppi_" + session_id + ".json"
				path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				#json_path = "ppi_" + session_id + ".json"
				#path_heatmap = "heatmap_" + session_id + ".png"
				if(save_data in ["save_data"]):
					if request.user.is_authenticated:
						print("saving data")
						username = str(request.user)
						if not (survival_col_name == ""):
							if("month" in survival_col_name):
								clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN_MONTH",1)
							else:
								clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN",1)
						# save input data
						GraphForm.save_user_data_3(exprstr,ppistr,clinicalstr,username)
						curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]	
						# save output data
						copyfile(("/code/clustering/static/" + path_heatmap),("user_uploaded_files/"+ username + "/" + curr_time + "_heatmap.png"))	
						copyfile(("/code/clustering/static/" + json_path),("user_uploaded_files/"+ username + "/" + curr_time + "_json.json"))	
						copyfile( path_metadata,("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt"))
						if(os.path.isfile(path_metadata + "_2")):
							copyfile((path_metadata+ "_2"),("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt_2"))
						copyfile(("/code/clustering/static/" + output_plot_path),("user_uploaded_files/"+ username + "/" + curr_time + "plotly.html"))
						copyfile(("/code/clustering/static/genelist_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist.txt"))
						copyfile(("/code/clustering/static/genelist_1_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_1.txt"))
						copyfile(("/code/clustering/static/genelist_2_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_2.txt"))
				print(ret_metadata1)
				# read metadata
				ret_metadata1=ret_metadata[0]
				ret_metadata2=ret_metadata[1]
				ret_metadata3=ret_metadata[2]
				# empty enrichment data from cache
				enrichment_dict = cache.get('enrichment_dict',"")
				if not(enrichment_dict == ""):
					cache.set("enrichment_dict","")
					cache.set("enrichment_dict_2","")
					cache.set("enrichment_dict_3","")
					cache.set("enrichment_dict_4","")
					cache.set("enrichment_dict_5","")
				# paths for showing results
				# write list of genes to downloadable file
				convert_gene_list.delay(adjlist,"/code/clustering/static/genelist_temp.txt")
				# save uploaded files if specified
				# render list of previously uploaded files if user is logged in (needed if user submits another request)
				# remove the loading-gif and progress image, clear cache             
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#if request.user.is_authenticated:
				#	request.session['done'] = "true"
				remove_loading_image.delay()
				if(os.path.isfile("clustering/static/loading_1.gif")):
					os.unlink("clustering/static/loading_1.gif")
				cache.clear()	
				if request.user.is_authenticated:
					request.session['done'] = "true"
					cache.set("done","done")
					cache.set('done',"done")
					username = str(request.user)
					list_of_files = GraphForm.list_user_data_2(username)	
					list_of_files_2 = GraphForm.list_user_data(username)              
				cache.set("has_survival_plot",has_survival_plot)	
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				# copy static files from shared directory to static-file-dir on web container
				copyfile(("/code/clustering/static/" + path_heatmap),("clustering/static/" + path_heatmap))	
				copyfile(("/code/clustering/static/" + json_path),("clustering/static/" + json_path))
				copyfile(("/code/clustering/static/" + output_plot_path),("clustering/static/" + output_plot_path))
				copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),("clustering/static/userfiles/genelist_" + session_id + ".txt"))
				copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
				copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				# save session ID and metadata in cache
				cache.set('session_id', session_id)	
				cache.set('ret_metadata1', ret_metadata1)	
				cache.set('ret_metadata2', ret_metadata2)	
				cache.set('ret_metadata3', ret_metadata3)	
				cache.set('json_path', json_path)	
				cache.set('p_val', p_val)
				# set "done" parameter
				cache.set('done',"done")
				cache.set('analysis_running','analysis_running')
				if(clinicalstr == "empty"):
					output_plot_path = "empty"
				#return render(request, 'clustering/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':"",'script':"",'plot2':"",'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val})
				return render(request, 'clustering/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':"",'script':"",'plot2':"",'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val,'has_survival_plot':has_survival_plot},status=301)
	elif('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			# configure loading page
			analysis_running = cache.get('analysis_running', 'none')
			# set analysis running parameter to allow display of loading images
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("Your request is being processed...")
			add_loading_image.delay()
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
			if('survival_col' in request.POST):
				if(request.POST['survival_col']):
					survival_col_name = request.POST['survival_col']
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
				if("SURVIVAL_COLUMN_MONTH" in clinicalstr):
					survival_col_name = "SURVIVAL_COLUMN_MONTH"
				elif("SURVIVAL_COLUMN" in clinicalstr):
					survival_col_name = "SURVIVAL_COLUMN"
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				## assign standard result size
				#if(request.POST['L_g_min'] == ""):
				#	lgmin = 10
				#if(request.POST['L_g_max'] == ""):
				#	lgmax = 20
				if(gene_set_size == "" or not str(gene_set_size).isdigit()):
					gene_set_size = 2000
				if (session_id_from_cache == 'has expired'):
					session_id = request.session._get_or_create_session_key()
				else:
					session_id = session_id_from_cache
				
				# run algorithm
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,gene_set_size)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()
				result1 = algo_output_task_2.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size,nbr_col)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()			
				# check if clinical data exist				
				if not(has_clin_data == "true"):
					survival_col_name = ""
					clinicalstr = "empty"
					ret_metadata = ""
				session_id_from_cache = cache.get('session_id', 'has expired')
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(ret_metadata,path_heatmap,path_metadata,output_plot_path,json_path,p_val) = result2.get()
				has_survival_plot = "true"
				if(output_plot_path == "empty"):
					has_survival_plot = "false"
				ret_metadata1 = ""
				ret_metadata2 = ""
				ret_metadata3 = ""
				ret_metadata1 = ret_metadata[0]
				ret_metadata2 = ret_metadata[1]
				ret_metadata3 = ret_metadata[2]
				enrichment_dict = cache.get('enrichment_dict',"")
				if not(enrichment_dict == ""):
					cache.set("enrichment_dict","")
					cache.set("enrichment_dict_2","")
					cache.set("enrichment_dict_3","")
					cache.set("enrichment_dict_4","")
					cache.set("enrichment_dict_5","")
				#metd = list_metadata_3.apply_async(countdown=0)
				#(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
				#json_path = "test15_" + session_id + ".json"
				#path_heatmap = "test_" + session_id + ".png"
				json_path = "userfiles/ppi_" + session_id + ".json"
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				#json_path = "ppi_" + session_id + ".json"
				#path_heatmap = "heatmap_" + session_id + ".png"
				#metd = list_metadata_3.apply_async(countdown=0)
				#path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				if(has_clin_data == "true"):
					#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
					if(os.path.isfile(path_metadata)):
						metd = list_metadata_5.apply_async(args=[path_metadata],countdown=0)
						(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				#(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				cache.clear()	
				#cache.set('session_id',session_id)		
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
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
				cache.set('done',"done")
				cache.set("has_survival_plot",has_survival_plot)
				remove_loading_image.delay()	
				return render(request, 'clustering/clustering_6_part_3.html', {'path_heatmap':path_heatmap,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'has_survival_plot':has_survival_plot})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		# get p-value from POST data
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		print("enrichment type")
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		session_id_from_cache = cache.get('session_id',"none")
		if(session_id_from_cache == "none"):
			session_id = request.session._get_or_create_session_key() 
		else:
			session_id = session_id_from_cache
		# one if loop for each enrichment type due to complicated naming of result files
		if(enr_type == "kegg_enrichment"):
			# run enrichment and write to directories
			result1 = run_enrichment_2.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",pval_enr,"clustering/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",pval_enr,"clustering/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",pval_enr,"clustering/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			# get enrichment results as dict
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",pval_enr,"clustering/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",pval_enr,"clustering/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",pval_enr,"clustering/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_" + session_id + ".txt",pval_enr,"/code/clustering/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt",pval_enr,"/code/clustering/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",pval_enr,"/code/clustering/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/clustering/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/clustering/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt",pval_enr,"/code/clustering/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/clustering/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		#print(request.POST.get("newAnalysis"))
		#print(request.POST['newAnalysis'])
		if('enr' not in request.POST):
			mutable = request.POST._mutable
			request.POST._mutable = True
			request.POST['enr'] = "true"
			request.POST._mutable = mutable
		if('enr' in request.POST):
			print("enr in request")
		# get current session id and result files for session
		#session_id = request.session._get_or_create_session_key()
		path_heatmap = "userfiles/heatmap_" + session_id + ".png"
		json_path = "userfiles/ppi_" + session_id + ".json"
		output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
		if('ppi_path' in request.POST and 'heatmap_path' in request.POST and 'plot_path' in request.POST):
			print(request.POST.get('ppi_path'))
			path_heatmap = request.POST.get('heatmap_path')
			json_path = request.POST.get('ppi_path')
			output_plot_path = request.POST.get('plot_path')
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true",'has_survival_plot':"true"})

	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		metdata_dict = ""
		done_from_cache = cache.get("done","")
		analysis_running = cache.get('analysis_running', 'none')
		session_id_from_cache = cache.get('session_id', 'has_expired')
		# if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
		if (analysis_running == 'none'):
			print("analysis not running")
			if(os.path.isfile("clustering/static/loading_1.gif")):
				os.unlink("clustering/static/loading_1.gif")
			if(os.path.isfile("clustering/static/progress.png")):
				os.unlink("clustering/static/progress.png")
			if(os.path.isfile("/code/clustering/static/progress.png")):
				os.unlink("/code/clustering/static/progress.png")
			remove_loading_image.delay()	
			#print("removed loading image")
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			with open("clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			#cache.set('analysis_running','analysis_running')
		
		if not (request.user.is_authenticated):
			cache.clear()
			cache.set('analysis_running','analysis_running')
			if(session_id_from_cache != 'has_expired'):
				cache.set('session_id',session_id_from_cache)
			cache.set('done',done_from_cache)
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			if not(session_id_from_cache == 'has expired'):
				session_id = session_id_from_cache
			else:
				session_id = request.session._get_or_create_session_key()
			if not(session_id == ""):
				if(os.path.isfile("/code/clustering/static/userfiles/metadata_" + session_id + ".txt")):
					metd = list_metadata_5.apply_async(args=["/code/clustering/static/userfiles/metadata_" + session_id + ".txt"],countdown=0)
					(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
					metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			cache.clear()
			cache.set('analysis_running','analysis_running')
			cache.set('session_id',session_id)
			cache.set('done',done_from_cache)
			#metd = list_metadata.apply_async(countdown=0)
			#metd = list_metadata_3.apply_async(countdown=0)
			print(ret_metadata1) 
			print("iubaerb")
			#result2 = read_kegg_enrichment.delay("/code/clustering/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			##result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			#enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	
		if('done' in request.session):
			print("clustering 6 part 3")
			print(request.session['done'])
			if(request.session['done'] == "true"): 
				session_id = request.session._get_or_create_session_key()
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				json_path = "userfiles/ppi_" + session_id + ".json"
				path_metadata = "userfiles/metadata_" + session_id + ".txt"
				output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
				return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'has_survival_plot':"true"})		
		# this redirects the user to the result page if the done parameter is true (e.g. reloading page after an analysis)
		if(done_from_cache=="done"):
			print("clustering 6 part 3")
			print(done_from_cache)
			session_id_from_cache = cache.get("session_id","none")
			print(session_id_from_cache)
			if(session_id_from_cache == "none"):
				session_id = request.session._get_or_create_session_key()
			else:
				session_id = session_id_from_cache
			path_heatmap = "userfiles/heatmap_" + session_id + ".png"
			json_path = "userfiles/ppi_" + session_id + ".json"
			path_metadata = "userfiles/metadata_" + session_id + ".txt"
			output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
			return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'has_survival_plot':"true"},status=301)				
			#return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})				
		return render(request,'clustering/clustering_6_part_1.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'has_survival_plot':"true"})




#########################################################################
#### version of the page that displays input form and result together ###
#########################################################################


def clustering_6_4(request):
	# the parameter analysis_running is true when an analysis has been run while the current cache exists. If it is false and an empty request is submitted (which is when an user first accesses the
	# page), and for some reason the output-console file for the progress page is filled with text, it gets emptied and the loading-gif removed.
	analysis_running = cache.get('analysis_running', 'none')
	print(analysis_running)
	session_id_from_cache = cache.get('session_id', 'has expired')
	print(session_id_from_cache)
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
	#print("save data parameter" + str(save_data))
	gene_set_size = request.POST.get("gene_set_size",2000)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	if not(os.path.isdir("/code/clustering/static/userfiles")):
		os.mkdir("/code/clustering/static/userfiles")
	# if the user wants to display old results
	if('input_own_file' in request.POST and 'display_old_results' in request.POST and request.user.is_authenticated):
		if(request.POST['input_own_file'] and request.POST['display_old_results']):
			# configure loading page
			#analysis_running = cache.get('analysis_running', 'none')
			#if (analysis_running == 'none'):
			#	cache.set('analysis_running','analysis_running')
			make_empty_figure.delay()
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
			session_id = request.session._get_or_create_session_key() 
			json_path = "userfiles/ppi_" + session_id + ".json"
			path_heatmap_2 = "userfiles/heatmap_" + session_id + ".png"
			#json_path = "ppi_" + session_id + ".json"
			#path_heatmap_2 = "heatmap_" + session_id + ".png"
			path_metadata_2 = "userfiles/metadata_" + session_id + ".txt"
			path_plotly_2 = "userfiles/output_plotly_" + session_id + ".html"
			# copy files to static directory
			copyfile(path_json,("clustering/static/" + json_path))	
			copyfile(path_heatmap,("clustering/static/" + path_heatmap_2))
			copyfile(path_genelist,("clustering/static/userfiles/genelist_" + session_id + ".txt"))
			copyfile(path_genelist_1,("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
			copyfile(path_genelist_2,("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
			output_plot_path_2 = ""
			ret_metadata_1 = ""
			ret_metadata_2 = ""
			ret_metadata_3 = ""
			# check if plotly file exists and copy
			if(os.path.isfile(path_plotly)):
				copyfile(path_plotly,("clustering/static/" + path_plotly_2))
				output_plot_path_2 = path_plotly_2
				#print("plot copied to")
				#print(path_plotly)
				#print(output_plot_path_2)
			# read metadata (must copy file to shared volume for processing via celery)
			#if(os.path.isfile(path_metadata+"_2")):
			if(os.path.isfile(path_metadata)):
				#print("found metadata")
				#print(path_metadata)
				#copyfile(path_metadata,("/code/clustering/static/metadata_" + session_id + ".txt"))
				#filename_for_old_metadata = "/code/clustering/static/metadata_" + session_id + ".txt"
				copyfile(path_metadata,("/code/clustering/static/userfiles/metadata_" + session_id + ".txt"))
				filename_for_old_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				#print(filename_for_old_metadata)
				#metd = list_metadata_4.apply_async(args=[filename_for_old_metadata],countdown=0)
				metd = list_metadata_5.apply_async(args=[filename_for_old_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#print(ret_metadata1)
			cache.clear()
			# set session ID in cache
			cache.set('session_id',session_id)
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			make_empty_figure.apply_async(countdown=10)
			empty_log_file.apply_async(countdown=10)
			# list old files
			list_of_files = ""
			list_of_files_2 = ""
			if request.user.is_authenticated:
				username = str(request.user)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)              				
			return render(request, 'clustering/clustering_6_part_4.html', {'form':"",'images':"",'plot_div':"",'script':"",'path_heatmap':path_heatmap_2,'output_plot_path':output_plot_path_2,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'has_survival_plot':"true"})
	
	elif(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		# check if input files exist
		input_valid = "false"
		if('myfile' in request.FILES and 'protfile' in request.FILES):
			if(request.FILES['myfile'] and request.FILES['protfile']):
				input_valid = "true"
		elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'protfile' in request.FILES):
			if(request.POST['predef_file'] and request.FILES['protfile']):
				input_valid = "true"
		elif('predef_file' in request.POST and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
			if(request.POST['predef_file'] and request.POST['parse_ndex_file'] and request.POST['ndex_name_2']):
				input_valid = "true"
		#if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
		if(input_valid == "true"):
			# set analysis running parameter to allow display of loading image
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			#if ('L_g_min' in request.POST and 'L_g_max' in request.POST):
			#	if(request.POST['L_g_min'] != "" and request.POST['L_g_max'] != ""):
			#		lgmin = int(request.POST['L_g_min'])
			#		lgmax = int(request.POST['L_g_max'])
			#	else:
			#		lgmin = 10
			#		lgmax = 20
			#else:
			#	lgmin = 10
			#	lgmax = 20
			#if(1 == 1):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				if(save_data in ["save_data"]):
					if request.user.is_authenticated:
						print("saving data is true")
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				## assign standard result size
				#if(request.POST['L_g_min'] == ""):
				#	lgmin = 10
				#if(request.POST['L_g_max'] == ""):
				#	lgmax = 20
				clinicalstr = ""
				clinicaldf = ""
				# configure loading page
				add_loading_image.delay()
				with open("/code/clustering/static/output_console.txt", "w") as text_file:
					text_file.write("Your request is being processed...")
					text_file.close()
				make_empty_figure.delay()
				clinicalstr = "empty"
				clinicaldf = ""
				survival_col_name = ""
				# read expression file
				if('myfile' in request.FILES):
					exprstr = request.FILES['myfile'].read().decode('utf-8')
					#result10 = preprocess_file.delay(exprstr)
					#exprstr = result10.get()
					result10 = preprocess_file_2.delay(exprstr)
					(exprstr,nbr_col) = result10.get()
				# read predefined expression file and clinical data
				elif('predef_file' in request.POST and 'cancer_type' in request.POST):
					cancer_type = request.POST.get("cancer_type")
					if(cancer_type == "1"):
						print("babababababa")
						fh1 = open("clustering/data/lung_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("clustering/data/lung_cancer_clinical.csv")
						fh4 = open("clustering/data/lung_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "disease free survival in months:ch1"
						nbr_col = 2
					else:
						fh1 = open("clustering/data/breast_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("clustering/data/breast_cancer_clinical.csv")
						fh4 = open("clustering/data/breast_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "mfs (yr):ch1"
						nbr_col = 2
				# read PPI file
				if('protfile' in request.FILES):
					ppistr = request.FILES['protfile'].read().decode('utf-8')
					result3 = preprocess_ppi_file.delay(ppistr)
					ppistr = result3.get()
					result4 = check_input_files.delay(ppistr,exprstr)
					errstr = result4.get()
					if(errstr != ""):
						request.session['errors'] = errstr
						return render(request,'clustering/errorpage.html',{'errors':errstr})
				# read ndex file from web
				elif('ndex_name_2' in request.POST):
					ndex_file_id = request.POST.get("ndex_name_2")
					if(ndex_file_id == "1"):
						result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "2"):
						#result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "3"):
						result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "4"):
						result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
				# read metadata if given
				if('analyze_metadata' in request.POST and 'patientdata' in request.FILES):
					if(request.FILES['patientdata']):
						clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
						clinicalstr_first_line = clinicalstr.split("\n")[1]
						if(len(clinicalstr_first_line.split("\t")) > len(clinicalstr_first_line.split(","))):
							print("converting metadata to csv format")
							clinicalstr = clinicalstr.replace("\t",",")
						clinical_stringio = StringIO(clinicalstr)
						clinicaldf = pd.read_csv(clinical_stringio)
						if('survival_col' in request.POST):
							if(request.POST['survival_col']):
								survival_col_name = request.POST['survival_col']
				session_id = ""
				session_id_from_cache = cache.get("session_id","none")
				if(session_id_from_cache == "none" or session_id_from_cache == ""):
					# start session for storing result data			
					session_id = request.session._get_or_create_session_key()
				else:
					session_id = session_id_from_cache
				# assign standard value to gene set size
				if(gene_set_size == "" or not str(gene_set_size).isdigit()):
					gene_set_size = 2000
				# run algorithm and read results
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size)
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id)
				result1 = algo_output_task_2.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size,nbr_col)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()			
				# make plots and process results	
				#print(group1_ids)
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(ret_metadata,path_heatmap,path_metadata,output_plot_path,json_path,p_val) = result2.get()
				has_survival_plot = "true"
				#print(output_plot_path)
				if(output_plot_path == "empty"):
					cache.set("has_survival_plot","false")
					has_survival_plot = "false"
				else:
					cache.set("has_survival_plot","true")
				#print(has_survival_plot)
				#output_plot_path = "output_plotly_" + session_id + ".html"
				#json_path = "ppi_" + session_id + ".json"
				output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
				json_path = "userfiles/ppi_" + session_id + ".json"
				path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				#path_metadata = "/code/clustering/static/metadata_" + session_id + ".txt"
				#path_heatmap = "heatmap_" + session_id + ".png"
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				#json_path = "ppi_" + session_id + ".json"
				#path_heatmap = "heatmap_" + session_id + ".png"
				if(save_data in ["save_data"]):
					if request.user.is_authenticated:
						#print("saving data in views.py")
						username = str(request.user)
						# replace name of survival column by "survival column"
						if not (survival_col_name == ""):
							if("month" in survival_col_name):
								# write "month" to indicate that survival time is given in months
								clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN_MONTH",1)
							else:
								clinicalstr = clinicalstr.replace(survival_col_name,"SURVIVAL_COLUMN",1)
						# save input data
						GraphForm.save_user_data_3(exprstr,ppistr,clinicalstr,username)
						curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]	
						# save output data
						copyfile(("/code/clustering/static/" + path_heatmap),("user_uploaded_files/"+ username + "/" + curr_time + "_heatmap.png"))	
						copyfile(("/code/clustering/static/" + json_path),("user_uploaded_files/"+ username + "/" + curr_time + "_json.json"))	
						copyfile( path_metadata,("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt"))
						#if(os.path.isfile(path_metadata + "_2")):
						#	copyfile((path_metadata+ "_2"),("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt_2"))
						copyfile(("/code/clustering/static/" + output_plot_path),("user_uploaded_files/"+ username + "/" + curr_time + "plotly.html"))
						if(os.path.isfile("/code/clustering/static/userfiles/genelist_" + session_id + ".txt")):
							copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist.txt"))
							copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_1.txt"))
							copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),("user_uploaded_files/"+ username + "/" + curr_time + "_genelist_2.txt"))
				
				# read metadata
				ret_metadata1=ret_metadata[0]
				ret_metadata2=ret_metadata[1]
				ret_metadata3=ret_metadata[2]
				# empty enrichment data from cache
				enrichment_dict = cache.get('enrichment_dict',"")
				if not(enrichment_dict == ""):
					cache.set("enrichment_dict","")
					cache.set("enrichment_dict_2","")
					cache.set("enrichment_dict_3","")
					cache.set("enrichment_dict_4","")
					cache.set("enrichment_dict_5","")
				# paths for showing results
				# write list of genes to downloadable file
				convert_gene_list.delay(adjlist,"/code/clustering/static/userfiles/genelist_temp.txt")
				convert_gene_list.delay(adjlist,"/code/clustering/static/userfiles/genelist_temp.txt")
				# save uploaded files if specified
				# render list of previously uploaded files if user is logged in (needed if user submits another request)
				if request.user.is_authenticated:
					username = str(request.user)
					list_of_files = GraphForm.list_user_data_2(username)	
					list_of_files_2 = GraphForm.list_user_data(username)              
				# remove the loading-gif and progress image, clear cache             
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#if request.user.is_authenticated:
				#	request.session['done'] = "true"
				#remove_loading_image.delay()
				# remove loading gif
				if(os.path.isfile("clustering/static/loading_1.gif")):
					os.unlink("clustering/static/loading_1.gif")
				#cache.clear()				
				#make_empty_figure.apply_async(countdown=10)
				#empty_log_file.apply_async(countdown=10)
				# copy static files from shared directory to static-file-dir on web container
				copyfile(("/code/clustering/static/" + path_heatmap),("clustering/static/" + path_heatmap))	
				copyfile(("/code/clustering/static/" + json_path),("clustering/static/" + json_path))
				copyfile(("/code/clustering/static/" + output_plot_path),("clustering/static/" + output_plot_path))
				#copyfile(("/code/clustering/static/genelist_" + session_id + ".txt"),("clustering/static/userfiles/genelist_" + session_id + ".txt"))
				#copyfile(("/code/clustering/static/genelist_1_" + session_id + ".txt"),("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
				#copyfile(("/code/clustering/static/genelist_2_" + session_id + ".txt"),("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
				copyfile(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),("clustering/static/userfiles/genelist_" + session_id + ".txt"))
				copyfile(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),("clustering/static/userfiles/genelist_1_" + session_id + ".txt"))
				copyfile(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),("clustering/static/userfiles/genelist_2_" + session_id + ".txt"))
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				# save session ID and metadata in cache
				cache.set('session_id', session_id)	
				cache.set('ret_metadata1', ret_metadata1)	
				cache.set('ret_metadata2', ret_metadata2)	
				cache.set('ret_metadata3', ret_metadata3)	
				cache.set('json_path', json_path)	
				cache.set('p_val', p_val)
				cache.set('analysis_running','analysis_running')
				if(clinicalstr == "empty"):
					output_plot_path = "empty"		
				return render(request, 'clustering/clustering_6_part_4.html', {'path_heatmap':path_heatmap,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'p_val':p_val,'has_survival_plot':has_survival_plot})
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			# configure loading page
			analysis_running = cache.get('analysis_running', 'none')
			if (analysis_running == 'none'):
				cache.set('analysis_running','analysis_running')
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("Your request is being processed...")
			add_loading_image.delay()
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
			if('survival_col' in request.POST):
				if(request.POST['survival_col']):
					survival_col_name = request.POST['survival_col']
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
				# get survival column name (standard name assigned)
				if("SURVIVAL_COLUMN_MONTH" in clinicalstr):
					survival_col_name = "SURVIVAL_COLUMN_MONTH"
				elif("SURVIVAL_COLUMN" in clinicalstr):
					survival_col_name = "SURVIVAL_COLUMN"
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				## assign standard result size
				#if(request.POST['L_g_min'] == ""):
				#	lgmin = 10
				#if(request.POST['L_g_max'] == ""):
				#	lgmax = 20
				if(gene_set_size == "" or not str(gene_set_size).isdigit()):
					gene_set_size = 2000
				session_id_from_cache = cache.get('session_id', 'has expired')
				if (session_id_from_cache == 'has expired' or session_id_from_cache == ""):
					session_id = request.session._get_or_create_session_key()
				else:
					session_id = session_id_from_cache
				# run algorithm
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size)
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,gene_set_size)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()
				result1 = algo_output_task_2.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id,gene_set_size,nbr_col)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()			
				# check if clinical data exist				
				if not(has_clin_data == "true"):
					survival_col_name = ""
					clinicalstr = "empty"
					ret_metadata = ""
				
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(ret_metadata,path_heatmap,path_metadata,output_plot_path,json_path,p_val) = result2.get()
				has_survival_plot = "true"
				if(output_plot_path == "empty"):
					cache.set("has_survival_plot","false")
					has_survival_plot = "false"
				else:
					cache.set("has_survival_plot","true")
				ret_metadata1 = ""
				ret_metadata2 = ""
				ret_metadata3 = ""
				ret_metadata1 = ret_metadata[0]
				ret_metadata2 = ret_metadata[1]
				ret_metadata3 = ret_metadata[2]
				enrichment_dict = cache.get('enrichment_dict',"")
				if not(enrichment_dict == ""):
					cache.set("enrichment_dict","")
					cache.set("enrichment_dict_2","")
					cache.set("enrichment_dict_3","")
					cache.set("enrichment_dict_4","")
					cache.set("enrichment_dict_5","")
				#metd = list_metadata_3.apply_async(countdown=0)
				#(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
				#json_path = "test15_" + session_id + ".json"
				#path_heatmap = "test_" + session_id + ".png"
				json_path = "userfiles/ppi_" + session_id + ".json"
				path_heatmap = "userfiles/heatmap_" + session_id + ".png"
				#json_path = "ppi_" + session_id + ".json"
				#path_heatmap = "heatmap_" + session_id + ".png"
				#metd = list_metadata_3.apply_async(countdown=0)
				#path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
				if(has_clin_data == "true"):
					#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
					if(os.path.isfile(path_metadata)):
						metd = list_metadata_5.apply_async(args=[path_metadata],countdown=0)
						(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				#(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				cache.clear()	
				#cache.set('session_id',session_id)		
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
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
				remove_loading_image.delay()	
				return render(request, 'clustering/clustering_6_part_4.html', {'path_heatmap':path_heatmap,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'has_survival_plot':has_survival_plot})
				#return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})
	
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"	
		# get p-value cutoff
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		analysis_running = cache.get('analysis_running', 'none')
		# set analysis running parameter to allow display of "loading"-gif + text
		if (analysis_running == 'none'):
			cache.set('analysis_running','analysis_running')
		surv_from_cache = cache.get('has_survival_plot','none')
		print(surv_from_cache)
		if(surv_from_cache == "false"):
			has_survival_plot = "false"

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
		session_id_from_cache = cache.get('session_id', 'has expired')
		if not(session_id_from_cache == 'has expired'):
			session_id = session_id_from_cache
		else:
			session_id = request.session._get_or_create_session_key()
		if not(session_id == ""):
			#genelist = "/code/clustering/static/genelist_" + session_id + ".txt"
			#genelist1 = "/code/clustering/static/genelist_1_" + session_id + ".txt"
			#genelist2 = "/code/clustering/static/genelist_2_" + session_id + ".txt"
			genelist = "/code/clustering/static/userfiles/genelist_" + session_id + ".txt"
			genelist1 = "/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"
			genelist2 = "/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"
			if not(os.path.isfile(genelist) and os.path.isfile(genelist1) and os.path.isfile(genelist2)):
				return render(request,'clustering/clustering_6_part_4.html',{'errstr':""})
			#kegg_dir = "/code/clustering/data/test/enrichr_kegg/" + session_id
			#kegg_output_dir = kegg_dir + "/KEGG_2013.test_name.enrichr.reports.txt"
			#print(genelist)
			if(enr_type == "kegg_enrichment"):
				result1 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_kegg/" + session_id))
				result2 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_kegg/" + session_id))
				result3 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_kegg/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				# both groups
				enrichment_dict = result4.get()
				# group 1, group 2
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				# results "only in group 1/2"
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_enrichment"):	
				result1 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_go/" + session_id))
				#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()	
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_molecular"):
				result1 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_go/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "reactome_enrichment"):
				result1 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_reactome/" + session_id))
				result2 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_reactome/" + session_id))
				result3 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_reactome/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_reactome/" + session_id + "/Reactome_2016.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
			# give links to result files from last analysis
			#path_heatmap = "test_" + session_id + ".png"
			#json_path = "test15_" + session_id + ".json"
			#path_heatmap = "heatmap_" + session_id + ".png"
			#json_path = "ppi_" + session_id + ".json"
			path_heatmap = "userfiles/heatmap_" + session_id + ".png"
			json_path = "userfiles/ppi_" + session_id + ".json"
			#path_metadata = "/code/clustering/static/metadata_" + session_id + ".txt"
			path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
			#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			if(os.path.isfile(path_metadata)):
				metd = list_metadata_5.apply_async(args=[path_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 			
			#output_plot_path = "output_plotly_" + session_id + ".html"
			output_plot_path = "userfiles/output_plotly_" + session_id + ".html"
			# write enrichment results to cache
			if not(session_id_from_cache == 'has expired'):
				session_id = session_id_from_cache
				cache.set("enrichment_dict",enrichment_dict)
				cache.set("enrichment_dict_2",enrichment_dict_2)
				cache.set("enrichment_dict_3",enrichment_dict_3)
				cache.set("enrichment_dict_4",enrichment_dict_4)
				cache.set("enrichment_dict_5",enrichment_dict_5)
			return render(request,'clustering/clustering_6_part_4.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'path_heatmap':path_heatmap,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,
'json_path':(json_path + "?foo=bar"),'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true",'has_survival_plot':has_survival_plot})
		if('enr' not in request.POST):
			mutable = request.POST._mutable
			request.POST._mutable = True
			request.POST['enr'] = "true"
			request.POST._mutable = mutable
		has_survival_plot = ""
		surv_from_cache = cache.get('has_survival_plot','none')
		print(surv_from_cache)
		if(surv_from_cache == "false"):
			has_survival_plot = "false"
		# get session id from cache
		session_id_from_cache = cache.get('session_id', 'has expired')
		if not(session_id_from_cache == 'has expired' or session_id_from_cache == ""):
			#path_heatmap = "heatmap_" + session_id_from_cache + ".png"
			#json_path = "ppi_" + session_id_from_cache + ".json"
			#path_metadata = "metadata_" + session_id_from_cache + ".txt"
			#output_plot_path = "output_plotly_" + session_id_from_cache + ".html"
			path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
			json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
			path_metadata = "userfiles/metadata_" + session_id_from_cache + ".txt"
			output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
			ret_metadata1 = cache.get('ret_metadata1','none')
			ret_metadata2 = cache.get('ret_metadata2','none')
			ret_metadata3 = cache.get('ret_metadata3','none')
		else:	
			session_id = request.session._get_or_create_session_key()
			# create new session if session id does not exist in cache
			#path_heatmap = "heatmap_" + session_id + ".png"
			#json_path = "ppi_" + session_id + ".json"
			#output_plot_path = "output_plotly_" + session_id + ".html"
			path_heatmap = "userfiles/heatmap_" + session_id + ".png"
			json_path = "userfiles/ppi_" + session_id + ".json"
			output_plot_path = "userfiles/output_plotly_" + session_id + ".html"	
		return render(request,'clustering/clustering_6_part_4.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true",'has_survival_plot':has_survival_plot})
	else:		
		analysis_running = cache.get('analysis_running', 'none')
		# if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
		if (analysis_running == 'none'):
			if(os.path.isfile("clustering/static/loading_1.gif")):
				os.unlink("clustering/static/loading_1.gif")
			if(os.path.isfile("clustering/static/progress.png")):
				os.unlink("clustering/static/progress.png")
			if(os.path.isfile("/code/clustering/static/progress.png")):
				os.unlink("/code/clustering/static/progress.png")
			remove_loading_image.delay()	
			#print("removed loading image")
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			with open("clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			#cache.set('analysis_running','analysis_running')
		ret_metadata = ""
		# check if session already exists for current user (e.g. when user has hit the reload button)
		session_id_from_cache = cache.get('session_id', 'has_expired')
		# get session ID and create list of previously saved files if user is authenticated
		if (request.user.is_authenticated):	
			#session_id = request.session._get_or_create_session_key()
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			#cache.clear()
			#ret_metadata1 = ""
			#ret_metadata2 = ""
			#ret_metadata3 = ""
			#if(os.path.isfile("/code/clustering/static/metadata_" + session_id + ".txt")):
			#	metd = list_metadata_5.apply_async(args=["/code/clustering/static/metadata_" + session_id + ".txt"],countdown=0)
			#	(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			#metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
		ret_metadata1 = ""
		ret_metadata2 = ""
		ret_metadata3 = ""
		has_survival_plot = ""
		surv_from_cache = cache.get('has_survival_plot','none')
		if(surv_from_cache == "false"):
			has_survival_plot = "false"
		
		# display results from most recent analysis
		if not (session_id_from_cache == 'has_expired' or session_id_from_cache == ""):
			#cache.set('session_id',session_id_from_cache)
			# take result files from storage
			path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
			json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
			#path_heatmap = "userfiles/test_" + session_id_from_cache + ".png"
			#json_path = "userfiles/test15_" + session_id_from_cache + ".json"
			output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
			#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			#(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			# take metadata from cache
			ret_metadata1 = cache.get('ret_metadata1',"")
			ret_metadata2 = cache.get('ret_metadata2',"")
			ret_metadata3 = cache.get('ret_metadata3',"")
			p_val = cache.get('p_val',"")
			# get enrichment results (in cache if user has hit reload button after running analysis)
			enrichment_dict = cache.get('enrichment_dict',"")
			enrichment_dict_2 = {}
			enrichment_dict_3 = {}
			enrichment_dict_4 = {}
			enrichment_dict_5 = {}			
			if not(enrichment_dict == ""):
				enrichment_dict_2 = cache.get('enrichment_dict_2',"")
				enrichment_dict_3 = cache.get('enrichment_dict_3',"")
				enrichment_dict_4 = cache.get('enrichment_dict_4',"")
				enrichment_dict_5 = cache.get('enrichment_dict_5',"")
			cache.clear()	
			if(surv_from_cache == "false"):
				has_survival_plot = "false"
				cache.set("has_survival_plot","false")
			cache.set('session_id', session_id_from_cache)	
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			cache.set('p_val', p_val)	
			
			return render(request,'clustering/clustering_6_part_4.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'p_val':p_val,'has_survival_plot':has_survival_plot})
		cache.clear()
		cache.set('analysis_running', analysis_running)
		return render(request,'clustering/clustering_6_part_4.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':"",'has_survival_plot':has_survival_plot})


def clustering_6_part_2_2(request):
		#print("badfbasdbasdbasdb")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#print(exprstr_par)
		print(request.session['ppistr'])
		if('done' in request.session):
			print(request.session['done'])
		#if(request.session.get('done') == "True"):
		#	print("done")
		return HttpResponseRedirect('clustering/clustering_6_part_3.html')
		#return render(request,'polls/clustering_6_part_3.html')
		#return(clustering_6_part_3_2(request))




def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		#print("badfbasdbasdbasdb")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#print(exprstr_par)
		#print(request.session['ppistr'])
		#if('done' in request.session):
		#	print(request.session['done'])
		#if(request.session.get('done') == "True"):
		#	print("done")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		return render(request,'clustering/errorpage.html',{'errors':errors})
		#return(clustering_6_part_3_2(request))


# method for generating the result page belonging to clustering_6_part_4_2
def clustering_6_part_3_2(request):
	#if('session_id' in request.POST):
	#	print(request.POST['session_id'])
	print("basdbasdbasdbasdbs")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	gene_set_size = request.POST.get("gene_set_size",2000)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		session_id_from_cache = cache.get('session_id', 'has expired')
		if not(session_id_from_cache == 'has expired'):
			session_id = session_id_from_cache
		else:
			session_id = request.session._get_or_create_session_key()
		if not(session_id == ""):
			genelist = "/code/clustering/static/userfiles/genelist_" + session_id + ".txt"
			genelist1 = "/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"
			genelist2 = "/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"
			if not(os.path.isfile(genelist) and os.path.isfile(genelist1) and os.path.isfile(genelist2)):
				return render(request,'clustering/clustering_6_part_3.html',{'errstr':""})
			#kegg_dir = "/code/clustering/data/test/enrichr_kegg/" + session_id
			#kegg_output_dir = kegg_dir + "/KEGG_2013.test_name.enrichr.reports.txt"
			#print(genelist)
			if(enr_type == "kegg_enrichment"):
				result1 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_kegg/" + session_id))
				result2 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_kegg/" + session_id))
				result3 = run_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_kegg/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				# both groups
				enrichment_dict = result4.get()
				# group 1, group 2
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				# results "only in group 1/2"
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_enrichment"):	
				result1 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_go/" + session_id))
				#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()	
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_molecular"):
				result1 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_go/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/clustering/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "reactome_enrichment"):
				result1 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test/enrichr_reactome/" + session_id))
				result2 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_1_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test2/enrichr_reactome/" + session_id))
				result3 = run_reac_enrichment.delay(("/code/clustering/static/userfiles/genelist_2_" + session_id + ".txt"),pval_enr,("/code/clustering/data/test3/enrichr_reactome/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/clustering/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/clustering/data/test2/enrichr_reactome/" + session_id + "/Reactome_2016.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/clustering/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
			# give links to result files from last analysis
			#path_heatmap = "test_" + session_id + ".png"
			#json_path = "test15_" + session_id + ".json"
			path_heatmap = "heatmap_" + session_id + ".png"
			json_path = "ppi_" + session_id + ".json"
			path_metadata = "/code/clustering/static/userfiles/metadata_" + session_id + ".txt"
			#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			if(os.path.isfile(path_metadata)):
				metd = list_metadata_5.apply_async(args=[path_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 			
			output_plot_path = "output_plotly_" + session_id + ".html"
			# write enrichment results to cache
			if not(session_id_from_cache == 'has expired'):
				session_id = session_id_from_cache
				cache.set("enrichment_dict",enrichment_dict)
				cache.set("enrichment_dict_2",enrichment_dict_2)
				cache.set("enrichment_dict_3",enrichment_dict_3)
				cache.set("enrichment_dict_4",enrichment_dict_4)
				cache.set("enrichment_dict_5",enrichment_dict_5)
			return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'path_heatmap':("userfiles/"+path_heatmap),'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,
'json_path':("userfiles/"+json_path + "?foo=bar"),'output_plot_path':("userfiles/"+output_plot_path),'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true",'has_survival_plot':"true"})
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'has_survival_plot':"true"})
	else:		
		analysis_running = cache.get('analysis_running', 'none')
		# if no analysis is running, remove loading image and text. this is to make sure after an incomplete analysis no "leftover" text with the status of last run is displayed
		if (analysis_running == 'none'):
			if(os.path.isfile("clustering/static/loading_1.gif")):
				os.unlink("clustering/static/loading_1.gif")
			if(os.path.isfile("clustering/static/progress.png")):
				os.unlink("clustering/static/progress.png")
			remove_loading_image.delay()	
			#print("removed loading image")
			with open("/code/clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			with open("clustering/static/output_console.txt", "w") as text_file:
				text_file.write("")
			#cache.set('analysis_running','analysis_running')
		ret_metadata = ""
		# check if session already exists for current user (e.g. when user has hit the reload button)
		session_id_from_cache = cache.get('session_id', 'has_expired')
		# get session ID and create list of previously saved files if user is authenticated
		if (request.user.is_authenticated):	
			#session_id = request.session._get_or_create_session_key()
			#username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			#cache.clear()
			#ret_metadata1 = ""
			#ret_metadata2 = ""
			#ret_metadata3 = ""
			#if(os.path.isfile("/code/clustering/static/metadata_" + session_id + ".txt")):
			#	metd = list_metadata_5.apply_async(args=["/code/clustering/static/metadata_" + session_id + ".txt"],countdown=0)
			#	(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			#metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
		ret_metadata1 = ""
		ret_metadata2 = ""
		ret_metadata3 = ""
		# display results from most recent analysis
		if not (session_id_from_cache == "has_expired"):
			#cache.set('session_id',session_id_from_cache)
			# take result files from storage
			path_heatmap = "userfiles/heatmap_" + session_id_from_cache + ".png"
			json_path = "userfiles/ppi_" + session_id_from_cache + ".json"
			#path_heatmap = "userfiles/test_" + session_id_from_cache + ".png"
			#json_path = "userfiles/test15_" + session_id_from_cache + ".json"
			output_plot_path = "userfiles/output_plotly_" + session_id_from_cache + ".html"
			#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			#(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			# take metadata from cache
			ret_metadata1 = cache.get('ret_metadata1',"")
			ret_metadata2 = cache.get('ret_metadata2',"")
			ret_metadata3 = cache.get('ret_metadata3',"")
			p_val = cache.get('p_val',"")
			# get enrichment results (in cache if user has hit reload button after running analysis)
			enrichment_dict = cache.get('enrichment_dict',"")
			enrichment_dict_2 = {}
			enrichment_dict_3 = {}
			enrichment_dict_4 = {}
			enrichment_dict_5 = {}			
			if not(enrichment_dict == ""):
				enrichment_dict_2 = cache.get('enrichment_dict_2',"")
				enrichment_dict_3 = cache.get('enrichment_dict_3',"")
				enrichment_dict_4 = cache.get('enrichment_dict_4',"")
				enrichment_dict_5 = cache.get('enrichment_dict_5',"")
			cache.clear()	
			cache.set('session_id', session_id_from_cache)	
			cache.set('ret_metadata1', ret_metadata1)	
			cache.set('ret_metadata2', ret_metadata2)	
			cache.set('ret_metadata3', ret_metadata3)	
			cache.set('p_val', p_val)	
			return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'path_heatmap':path_heatmap,'json_path':json_path,'output_plot_path':output_plot_path,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'p_val':p_val,'has_survival_plot':"true"})
		session_id = request.session._get_or_create_session_key()
		path99 = "heatmap_" + session_id + ".png"
		json_path = "ppi_" + session_id + ".json"
		path_metadata = "metadata_" + session_id + ".txt"
		output_plot_path = "output_plotly_" + session_id + ".html"
		cache.clear()
		cache.set('analysis_running', analysis_running)
		return render(request,'clustering/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'path_heatmap':path99,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'has_survival_plot':"true"})





def infopage(request):			
	return render(request,'clustering/infopage.html')

def sources(request):			
	return render(request,'clustering/sources.html')



