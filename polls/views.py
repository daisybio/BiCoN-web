from django.shortcuts import render
from io import StringIO
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
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
#### own imports
import polls
from polls.models import Upload,UploadForm,GraphForm
from django.contrib.auth import authenticate, login, logout
from polls.script3 import run_algo,algo_output,algo_output_ownfile,algo_output_ownfile_2,algo_output_ownfile_3
#from polls.tasks import make_empty_figure,algo_output_task,script_output_task,empty_log_file,write_to_file_1,add_loading_image,remove_loading_image,script_output_task_2,show_old_data,script_output_task_3,write_metadata_to_file,list_metadata,metadata_to_string,script_output_task_4,read_ndex_file,read_ndex_file_2,read_ndex_file_3,list_metadata_2,list_metadata,script_output_task_7,script_output_task_8,script_output_task_9,list_metadata_3,run_enrichment,read_kegg_enrichment,run_go_enrichment,read_ndex_file_4,run_enrichment_2,run_go_enrichment_2,run_reac_enrichment,import_ndex,read_kegg_enrichment_2,convert_gene_list,check_input_files,script_output_task_10,list_metadata_4
from polls.tasks import make_empty_figure,algo_output_task,empty_log_file,write_to_file_1,add_loading_image,remove_loading_image,script_output_task_2,show_old_data,write_metadata_to_file,metadata_to_string,script_output_task_9,list_metadata_3,read_kegg_enrichment,read_ndex_file_4,run_enrichment_2,run_go_enrichment_2,run_reac_enrichment,read_kegg_enrichment_2,convert_gene_list,check_input_files,script_output_task_10,list_metadata_4,preprocess_file,write_pval,algo_output_task_new

import os.path

### *ACTUAL* imports (that have dependencies other than django and my own stuff) ####
import networkx as nx
from biomart import BiomartServer
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
		return redirect('polls/logout.html')
	else:
		return render(request,'polls/logout.html')

def login_2(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		user = authenticate(request, username=username, password=password)
		if user is not None:
			login(request, user)
			#return render(request,'polls/login.html')
			return redirect('polls/clustering_6.html')
		else:
			return render(request,'polls/login.html')
			#return redirect('polls/clustering.html')
	else:
		return render(request,'polls/login.html')
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
			return render(request,'polls/signup.html',{'text':text})
		else:
			user = User.objects.create_user(username, email, password)
			user.save()
			current_site = get_current_site(request)
			mail_subject = 'Activate your account.'
			message = render_to_string('acc_active_email.html', {
			'user': user,
			'domain': current_site.domain,
			'uid':urlsafe_base64_encode(force_bytes(user.pk)),
			'token':account_activation_token.make_token(user),
			})
			to_email = email
			email_message = EmailMessage(
			mail_subject, message, to=[to_email]
			)
			email_message.send()
			foobar = "user_uploaded_files/" + username
			if not(os.path.isdir(foobar)):
				os.mkdir(foobar)
			text = "Account is being created. You will receive a confirmation e-mail soon!"
			return render(request,'polls/signup.html',{'text':text,'new_user':"true"})
			#return redirect('polls/clustering.html')
	else:
		text = "Please input username and password!"
		return render(request,'polls/signup.html',{'text':text})
		#return redirect('polls/clustering.html')




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
		return render(request,'polls/errorpage.html',{'errors':errors})
		#return(clustering_6_part_3_2(request))


#########################################################################
#################### this is not used yet ###############################
#########################################################################



def clustering_6_new(request):
	#print("in clustering")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
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
	if(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				clinicalstr = ""
				clinicaldf = ""
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				clinicalstr = "empty"
				clinicaldf = ""
				survival_col_name = ""
				if('myfile' in request.FILES):
					exprstr = request.FILES['myfile'].read().decode('utf-8')
					result10 = preprocess_file.delay(exprstr)
					exprstr = result10.get()
				elif('predef_file' in request.POST and 'cancer_type' in request.POST):
					cancer_type = request.POST.get("cancer_type")
					if(cancer_type == "1"):
						print("babababababa")
						fh1 = open("polls/data/lung_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
						fh4 = open("polls/data/lung_cancer_clinical.csv")
						#fh4 = open("polls/data/breast_cancer_clinical.csv")
						is_lungc = "true"
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "disease free survival in months:ch1"
					else:
						fh1 = open("polls/data/breast_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
						fh4 = open("polls/data/breast_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "mfs (yr):ch1"
				if('protfile' in request.FILES):
					ppistr = request.FILES['protfile'].read().decode('utf-8')
					result4 = check_input_files.delay(ppistr,exprstr)
					errstr = result4.get()
					if(errstr != ""):
						request.session['errors'] = errstr
						return render(request,'polls/errorpage.html',{'errors':errstr})
				elif('ndex_name_2' in request.POST):
					ndex_file_id = request.POST.get("ndex_name_2")
					if(ndex_file_id == "1"):
						result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "2"):
						result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "3"):
						result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "4"):
						result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
				if('analyze_metadata' in request.POST):
					if(request.FILES['patientdata']):
						clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
						clinical_stringio = StringIO(clinicalstr)
						clinicaldf = pd.read_csv(clinical_stringio)
						if('survival_col' in request.POST):
							#print("barabsfrbasdb")
							if(request.POST['survival_col']):
								survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#print(adjlist)
				convert_gene_list.delay(adjlist,"/code/polls/static/genelist_temp.txt")
				#print("iasdfasdfsf")
				#print(ret_metadata1)
				write_pval.apply_async([p_val,"/code/polls/static/pvalue.txt"],countdown=0)	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			#print("save data")
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val})
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
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

				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				if not(has_clin_data == "true"):
					clinicalstr = "empty"
					ret_metadata = ""
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})
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
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
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
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	 				
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})




#########################################################################
#### version of the page with separate input form and result display ####
#########################################################################

def clustering_6_4_part_2(request):
	# check if user has clicked 'return' button on result page, then request.POST['newAnalysis'] is "true"
	if('newAnalysis' in request.POST):
		print(request.POST.get("newAnalysis"))
		print(request.POST['newAnalysis'])
		request.POST._mutable = True
		if(request.POST['newAnalysis'] != "false"):
			if('done' in request.session):
				if(request.session['done'] == "true"):
					#set done parameter to false if user has clicked return on result page
					request.session['done'] = "False"
			# remove parameter from request.POST to allow later switching to result page
			request.POST['newAnalysis'] = "false"
	#if(request.session['done'] == "true"):
	if('done' in request.session):
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		print("done")
	if('session_id' in request.POST):
		print(request.POST['session_id'])
	#print("in clustering")
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
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
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
	if(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		# check if these files are not empty and exist
		if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				clinicalstr = ""
				clinicaldf = ""
				add_loading_image.delay()
				with open("/code/polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				clinicalstr = "empty"
				clinicaldf = ""
				survival_col_name = ""
				# read expression file
				if('myfile' in request.FILES):
					exprstr = request.FILES['myfile'].read().decode('utf-8')
					result10 = preprocess_file.delay(exprstr)
					exprstr = result10.get()
				# read predefined expression file and clinical data
				elif('predef_file' in request.POST and 'cancer_type' in request.POST):
					cancer_type = request.POST.get("cancer_type")
					if(cancer_type == "1"):
						print("babababababa")
						fh1 = open("polls/data/lung_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
						fh4 = open("polls/data/lung_cancer_clinical.csv")
						#fh4 = open("polls/data/breast_cancer_clinical.csv")
						is_lungc = "true"
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "disease free survival in months:ch1"
					else:
						fh1 = open("polls/data/breast_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
						fh4 = open("polls/data/breast_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "mfs (yr):ch1"
				# read PPI file
				if('protfile' in request.FILES):
					ppistr = request.FILES['protfile'].read().decode('utf-8')
					result4 = check_input_files.delay(ppistr,exprstr)
					errstr = result4.get()
					if(errstr != ""):
						request.session['errors'] = errstr
						return render(request,'polls/errorpage.html',{'errors':errstr})
				# read ndex file from web
				elif('ndex_name_2' in request.POST):
					ndex_file_id = request.POST.get("ndex_name_2")
					if(ndex_file_id == "1"):
						result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "2"):
						result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "3"):
						result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "4"):
						result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
				# read metadata if given
				if('analyze_metadata' in request.POST):
					if(request.FILES['patientdata']):
						clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
						clinical_stringio = StringIO(clinicalstr)
						clinicaldf = pd.read_csv(clinical_stringio)
						if('survival_col' in request.POST):
							#print("barabsfrbasdb")
							if(request.POST['survival_col']):
								survival_col_name = request.POST['survival_col']
				# run algorithm and read results
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()	
				# start session for storing result data			
				session_id = request.session._get_or_create_session_key()
				# make plots and process results
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(div,script,plot1,plot_div,ret_metadata,path99,path_metadata,output_plot_path,json_path) = result2.get()
				# read metadata
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				# paths for showing results
				output_plot_path = "output_plotly_" + session_id + ".html"
				json_path = "test15_" + session_id + ".json"
				#metd = list_metadata_3.apply_async(countdown=0)
				path_metadata = "/code/polls/static/metadata_" + session_id + ".txt"
				p_val = ""
				#print(adjlist)
				convert_gene_list.delay(adjlist,"/code/polls/static/genelist_temp.txt")
				#print("iasdfasdfsf")
				#print(ret_metadata1)
				write_pval.apply_async([p_val,"/code/polls/static/pvalue.txt"],countdown=0)	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				# save data
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			#print("save data")
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("/code/polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				with open("/code/polls/static/metadata_test.txt","w") as text_file_4:
					text_file_4.write("bla")
				metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#return HttpResponseRedirect(reverse('polls/clustering_6_part_2_2/'),args=(exprstr,ppistr))
				#return clustering_6_part_2_2(request,exprstr,ppistr)
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				if request.user.is_authenticated:
					request.session['done'] = "true"
				# store session id in POST data
				else:
					request.POST._mutable = True
					request.POST['session_id'] = session_id
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#session_id = request.session._get_or_create_session_key()
				#print(session_id)
				#request.session['done'] = "true"
				# copy plots from shared volume to static directory
				path99 = "test_" + session_id + ".png"				
				copyfile(("/code/polls/static/" + path99),("polls/static/" + path99))	
				copyfile(("/code/polls/static/" + json_path),("polls/static/" + json_path))
				copyfile(("/code/polls/static/" + output_plot_path),("polls/static/" + output_plot_path))
				if request.user.is_authenticated:
					request.session['done'] = "true"
				else:
					request.POST._mutable = True
					request.POST['session_id'] = session_id
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val})
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
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
				
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				if not(has_clin_data == "true"):
					clinicalstr = "empty"
					ret_metadata = ""
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)  	 	   
				remove_loading_image.delay()	
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})

	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
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
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
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
		session_id = request.session._get_or_create_session_key()
		path99 = "test_" + session_id + ".png"
		json_path = "test15_" + session_id + ".json"
		output_plot_path = "output_plotly_" + session_id + ".html"
		if('ppi_path' in request.POST and 'heatmap_path' in request.POST and 'plot_path' in request.POST):
			print(request.POST.get('ppi_path'))
			path99 = request.POST.get('heatmap_path')
			json_path = request.POST.get('ppi_path')
			output_plot_path = request.POST.get('plot_path')
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	
		if('done' in request.session):
			print("clustering 6 part 3")
			print(request.session['done'])
			if(request.session['done'] == "true"): 
				session_id = request.session._get_or_create_session_key()
				path99 = "test_" + session_id + ".png"
				json_path = "test15_" + session_id + ".json"
				path_metadata = "metadata_" + session_id + ".txt"
				output_plot_path = "output_plotly_" + session_id + ".html"
				return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})		
		return render(request,'polls/clustering_6_part_1.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})




#########################################################################
#### version of the page that displays input form and result together ###
#########################################################################


def clustering_6_4(request):
	if('newAnalysis' in request.POST):
		print(request.POST.get("newAnalysis"))
		print(request.POST['newAnalysis'])
		request.POST._mutable = True
		if(request.POST['newAnalysis'] != "false"):
			if('done' in request.session):
				if(request.session['done'] == "true"):
					request.session['done'] = "False"
			request.POST['newAnalysis'] = "false"
	#if(request.session['done'] == "true"):
	if('done' in request.session):
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		print("done")
	if('session_id' in request.POST):
		print(request.POST['session_id'])
	#print("in clustering")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
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
	if('input_own_file' in request.POST and 'display_old_results' in request.POST and not ('myfile' in request.FILES) and request.user.is_authenticated):
		if(request.POST['input_own_file'] and request.POST['display_old_results']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			#fh1 = open(filename1)
			#exprstr = fh1.read()
			#filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			#fh2 = open(filename2)
			path_json = filename1
			path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			path_metadata = filename1.split("_json.json")[0] + "metadata.txt"
			path_plotly = filename1.split("_json.json")[0] + "plotly.html"
			#ppistr = fh2.read()
			#list_of_files = GraphForm.list_user_data_2(username)	
			#list_of_files_2 = GraphForm.list_user_data(username)     
			session_id = request.session._get_or_create_session_key() 
			json_path = "test15_" + session_id + ".json"
			#metd = list_metadata_3.apply_async(countdown=0)
			#path_metadata = "/code/polls/static/metadata_" + session_id + ".txt"
			path99 = "test_" + session_id + ".png"
			path_metadata_2 = "metadata_" + session_id + ".txt"
			path_plotly_2 = "output_plotly_" + session_id + ".html"
			copyfile(path_json,("polls/static/" + json_path))	
			copyfile(path_heatmap,("polls/static/" + path99))
			#copyfile( path_metadata,("polls/static/"+path_metadata_2))
			#copyfile(path_plotly ,("polls/static/"+path_plotly_2))
			cache.clear()			
			make_empty_figure.apply_async(countdown=10)
			empty_log_file.apply_async(countdown=10)
			return render(request, 'polls/clustering_6_part_4.html', {'form':"",'images':"",'plot_div':"",'script':"",'plot2':"",'path4':path99,'output_plot_path':"",'json_path':json_path, 'list_of_files':"",'ret_dat':"",'ret_metadata1':"",'ret_metadata2':"",'ret_metadata3':"",'list_of_files_2':""})
	
	elif(('myfile' in request.FILES or 'predef_file' in request.POST) and ('protfile' in request.FILES or ('parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST))):
		if((request.FILES['myfile'] or request.POST['predef_file']) and (request.FILES['protfile'] or (request.POST['parse_ndex_file'] and request.POST['ndex_name_2']))):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				clinicalstr = ""
				clinicaldf = ""
				add_loading_image.delay()
				with open("/code/polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				clinicalstr = "empty"
				clinicaldf = ""
				survival_col_name = ""
				if('myfile' in request.FILES):
					exprstr = request.FILES['myfile'].read().decode('utf-8')
					result10 = preprocess_file.delay(exprstr)
					exprstr = result10.get()
				elif('predef_file' in request.POST and 'cancer_type' in request.POST):
					cancer_type = request.POST.get("cancer_type")
					if(cancer_type == "1"):
						print("babababababa")
						fh1 = open("polls/data/lung_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
						fh4 = open("polls/data/lung_cancer_clinical.csv")
						#fh4 = open("polls/data/breast_cancer_clinical.csv")
						is_lungc = "true"
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "disease free survival in months:ch1"
					else:
						fh1 = open("polls/data/breast_cancer_expr.csv")
						exprstr = fh1.read()
						clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
						fh4 = open("polls/data/breast_cancer_clinical.csv")
						clinicalstr = fh4.read()
						fh4.flush()
						fh4.close()
						survival_col_name = "mfs (yr):ch1"
				if('protfile' in request.FILES):
					ppistr = request.FILES['protfile'].read().decode('utf-8')
					result4 = check_input_files.delay(ppistr,exprstr)
					errstr = result4.get()
					if(errstr != ""):
						request.session['errors'] = errstr
						return render(request,'polls/errorpage.html',{'errors':errstr})
				elif('ndex_name_2' in request.POST):
					ndex_file_id = request.POST.get("ndex_name_2")
					if(ndex_file_id == "1"):
						result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "2"):
						result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "3"):
						result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
					elif(ndex_file_id == "4"):
						result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
						ppistr = result_ndex.get()
				if('analyze_metadata' in request.POST):
					if(request.FILES['patientdata']):
						clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
						clinical_stringio = StringIO(clinicalstr)
						clinicaldf = pd.read_csv(clinical_stringio)
						if('survival_col' in request.POST):
							#print("barabsfrbasdb")
							if(request.POST['survival_col']):
								survival_col_name = request.POST['survival_col']
				session_id = ""
				#result1 = algo_output_task_new.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig,session_id)
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				session_id = request.session._get_or_create_session_key()
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(div,script,plot1,plot_div,ret_metadata,path99,path_metadata,output_plot_path,json_path) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				output_plot_path = "output_plotly_" + session_id + ".html"
				json_path = "test15_" + session_id + ".json"
				#metd = list_metadata_3.apply_async(countdown=0)
				path_metadata = "/code/polls/static/metadata_" + session_id + ".txt"
				
				#print(adjlist)
				convert_gene_list.delay(adjlist,"/code/polls/static/genelist_temp.txt")
				#print("iasdfasdfsf")
				#print(ret_metadata1)
				#write_pval.apply_async([p_val,"/code/polls/static/pvalue.txt"],countdown=0)	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			#print("save data")
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("/code/polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				with open("/code/polls/static/metadata_test.txt","w") as text_file_4:
					text_file_4.write("bla")
				metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#return HttpResponseRedirect(reverse('polls/clustering_6_part_2_2/'),args=(exprstr,ppistr))
				#return clustering_6_part_2_2(request,exprstr,ppistr)
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"				
				path99 = "test_" + session_id + ".png"
				if save_data in ["save_data"]:
					if request.user.is_authenticated:
						savedata_param = "true"
						username = str(request.user)
						#GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
						curr_time = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
						copyfile(("/code/polls/static/" + path99),("user_uploaded_files/"+ username + "/" + curr_time + "_heatmap.png"))	
						copyfile(("/code/polls/static/" + json_path),("user_uploaded_files/"+ username + "/" + curr_time + "_json.json"))	
						copyfile( path_metadata,("user_uploaded_files/"+ username + "/" + curr_time + "metadata.txt"))
						copyfile(("/code/polls/static/" + output_plot_path),("user_uploaded_files/"+ username + "/" + curr_time + "plotly.html"))
				
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				if request.user.is_authenticated:
					request.session['done'] = "true"
				else:
					request.POST._mutable = True
					request.POST['session_id'] = session_id
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				copyfile(("/code/polls/static/" + path99),("polls/static/" + path99))	
				copyfile(("/code/polls/static/" + json_path),("polls/static/" + json_path))
				copyfile(("/code/polls/static/" + output_plot_path),("polls/static/" + output_plot_path))
				#session_id = request.session._get_or_create_session_key()
				#print(session_id)
				#request.session['done'] = "true"
				path99 = "test_" + session_id + ".png"
				#pval = ""
				return render(request, 'polls/clustering_6_part_4.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2,'path4':path99,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
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
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				if not(has_clin_data == "true"):
					clinicalstr = "empty"
					ret_metadata = ""
				session_id = request.session._get_or_create_session_key()
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(div,script,plot1,plot_div,ret_metadata,path99,path_metadata,output_plot_path,json_path) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				output_plot_path = "output_plotly_" + session_id + ".html"
				json_path = "test15_" + session_id + ".json"
				#metd = list_metadata_3.apply_async(countdown=0)
				path_metadata = "/code/polls/static/metadata_" + session_id + ".txt"
				ret_metadata1 = ""
				ret_metadata2 = ""
				ret_metadata3 = ""
				if(has_clin_data == "true"):
					metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
					(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				#result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				#(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)  	 	   
				remove_loading_image.delay()	
				return render(request, 'polls/clustering_6_part_4.html', {'form':"",'images':"",'plot_div':"",'script':"",'plot2':plot2,'path4':path99,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':"",'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
				#return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})
	
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
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
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
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		return render(request,'polls/clustering_6_part_4.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
			if not ('session_id' in request.POST):
				session_id = request.session._get_or_create_session_key()
				request.POST._mutable = True
				request.POST['session_id'] = session_id
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	
		#if('done' in request.session):
		#	print("clustering 6 part 3")
		#	print(request.session['done'])
		#	if(request.session['done'] == "true"): 
		#		session_id = request.session._get_or_create_session_key()
		#		path99 = "test_" + session_id + ".png"
		#		json_path = "test15_" + session_id + ".json"
		#		path_metadata = "metadata_" + session_id + ".txt"
		#		output_plot_path = "output_plotly_" + session_id + ".html"
		#		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})		
		return render(request,'polls/clustering_6_part_4.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})



##This is the version of the current webpage without using the session ID.##
def clustering_6(request):
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	#check which option (redo, new analysis etc) and which files (own file, predefined file) the user inputs


	# this is called if the user wants to redo an analysis.
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
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
				
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				if not(has_clin_data == "true"):
					clinicalstr = "empty"
					ret_metadata = ""
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})       		

	##################################################
	###########  LOOK HERE!                 ##########
	##################################################

	#### this method is called if the user uploads everything himself. very nice for testing!

	elif('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		#check if the needed parameters and files are not empty
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				# prepare the progress bar etc
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				
				#read the files
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result10 = preprocess_file.delay(exprstr)
				exprstr = result10.get()				
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				#dataframe of clinical data for method
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				#check input files format, if not correct display error page
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				survival_col_name = ""
				# check for specified 'survival' column name
				if('survival_col' in request.POST):
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']

				#######################
				#running the algorithm!
				#######################

				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				
				write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
				#get metadata from file
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				convert_gene_list.delay(adjlist,"polls/static/genelist_temp.txt")
				print(ret_metadata1) 	
				write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)					
				# save data if user wants to
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				# list the filenames of previous saved requests				
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				# remove the loading page stuff
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				#write survival plot to a file
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				plot2 = "test.png"
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2,'pval':p_val})
	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST and 'analyze_metadata' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result10 = preprocess_file.delay(exprstr)
				exprstr = result10.get()
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	
					
					write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	
	elif('myfile' in request.FILES and 'ndexfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['ndexfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				result3 = read_ndex_file_2.delay(request.FILES['ndexfile'].read().decode('utf-8'))
				ppistr_1 = result3.get()
				ppistr_2 = request.FILES['protfile'].read().decode('utf-8')
				ppistr = "\n".join([ppistr_2,ppistr_1])
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					survival_col_name = ""
					clinicalstr = "empty"
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					survival_col_name = ""
					clinicalstr = "empty"
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	
					
					write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('predef_file' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['predef_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	
					
					write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST and 'analyze_metadata' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	
					
					write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_file_input' in request.FILES):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['ndex_file_input']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_str = request.FILES['ndex_file_input'].read().decode('utf-8')				
				result_ndex = read_ndex_file_4.delay(ndex_file_str)
				ppistr = result_ndex.get()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					survival_col_name = ""
					clinicalstr = "empty"
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST and 'analyze_metadata' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#clinicalstr = "empty"
				#if(request.user.is_authenticated and savedata_param == "true"):
				#	result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
				#	(div,script,plot1) = result2.get()			
				#else:				
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	

				write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	
	
	
	
	elif('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			#fh1 = open(filename1)
			#exprstr = fh1.read()
			#filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			#fh2 = open(filename2)
			path_json = filename1
			path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			#ppistr = fh2.read()
			list_of_files = GraphForm.list_user_data_2(username)	
			list_of_files_2 = GraphForm.list_user_data(username)      
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				show_old_data(path_json,path_heatmap)             	
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})       		
	elif('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				#with open("polls/static/output_console.txt", "w") as text_file:
				with open("/code/polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				clinicalstr = "empty"
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()	
					
					write_pval.apply_async([p_val,"/code/polls/static/pvalue.txt"],countdown=0)
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			foobar = "/code/user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                		
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})
	elif('predef_file' in request.POST and 'protfile' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file'] and request.POST['protfile']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh2 = open("polls/data/biogrid.human.entrez.tsv")
				#ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				is_lungc = "false"					
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				print(clinicaldf)
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				
				write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
				print(ret_metadata)
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				with open("polls/static/metadata_99.html", "w") as text_file_4:
   					text_file_4.write(str(ret_metadata[0]))
   					text_file_4.write(str(ret_metadata[1]))
   					text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				result_99 = string_99.get()
				return render(request, 'polls/clustering_6.html', {'plot2':plot2, 'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	
	elif('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				is_lungc = "false"					
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				print(clinicaldf)
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata,p_val) = result2.get()
				
				write_pval.apply_async([p_val,"polls/static/pvalue.txt"],countdown=0)
				print(ret_metadata)
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				with open("polls/static/metadata_99.html", "w") as text_file_4:
   					text_file_4.write(str(ret_metadata[0]))
   					text_file_4.write(str(ret_metadata[1]))
   					text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				result_99 = string_99.get()
				return render(request, 'polls/clustering_6.html', {'plot2':plot2, 'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
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
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
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
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
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
		#return clustering_6(request)
		#return HttpResponseRedirect('polls/clustering_6.html')
		path99 = ""
		json_path = ""
		output_plot_path = ""
		if('ppi_path' in request.POST and 'heatmap_path' in request.POST and 'plot_path' in request.POST):
			print(request.POST.get('ppi_path'))
			path99 = request.POST.get('heatmap_path')
			json_path = request.POST.get('ppi_path')
			output_plot_path = request.POST.get('plot_path')
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	 				
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})


def clustering_6_part_2_2(request):
		#print("badfbasdbasdbasdb")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#print(exprstr_par)
		print(request.session['ppistr'])
		if('done' in request.session):
			print(request.session['done'])
		#if(request.session.get('done') == "True"):
		#	print("done")
		return HttpResponseRedirect('polls/clustering_6_part_3.html')
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
		return render(request,'polls/errorpage.html',{'errors':errors})
		#return(clustering_6_part_3_2(request))

def clustering_6_part_3_2(request):
	if('session_id' in request.POST):
		print(request.POST['session_id'])
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
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES and (1 == 0)):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				#return HttpResponseRedirect(reverse('clustering_6_step_2'))
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				session_id = request.session._get_or_create_session_key()
				path99 = "polls/static/test_" + session_id + ".png"
				path_metadata = "polls/static/metadata_" + session_id + ".txt"
				output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'output_plot_path':output_plot_path,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
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
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"/code/polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"/code/polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("/code/polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("/code/polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","/code/polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"/code/polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("/code/polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		session_id = request.session._get_or_create_session_key()
		path99 = "test_" + session_id + ".png"
		json_path = "test15_" + session_id + ".json"
		output_plot_path = "output_plotly_" + session_id + ".html"
		if('ppi_path' in request.POST and 'heatmap_path' in request.POST and 'plot_path' in request.POST):
			print(request.POST.get('ppi_path'))
			path99 = request.POST.get('heatmap_path')
			json_path = request.POST.get('ppi_path')
			output_plot_path = request.POST.get('plot_path')
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'json_path':json_path,'output_plot_path':output_plot_path,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})
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
		if(request.user.is_authenticated):
			session_id = request.session._get_or_create_session_key()
			genelist = "genelist_" + session_id + ".txt"
			genelist1 = "genelist_1_" + session_id + ".txt"
			genelist2 = "genelist_2_" + session_id + ".txt"
			kegg_dir = "/code/polls/data/test/enrichr_kegg/" + session_id
			kegg_dir_2 = "/code/polls/data/test2/enrichr_kegg/" + session_id
			kegg_dir_3 = "/code/polls/data/test3/enrichr_kegg/" + session_id
			kegg_output_dir = kegg_dir + "/KEGG_2013.test_name.enrichr.reports.txt"
			kegg_output_dir_2 = kegg_dir_2 + "/KEGG_2013.test_name.enrichr.reports.txt"
			kegg_output_dir_3 = kegg_dir_3 + "/KEGG_2013.test_name.enrichr.reports.txt"
			if(enr_type == "kegg_enrichment"):
				result1 = run_enrichment_2.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test/enrichr_kegg/" + session_id))
				result2 = run_enrichment_2.delay(("genelist_1_" + session_id + ".txt"),pval_enr,("/code/polls/data/test2/enrichr_kegg/" + session_id))
				result3 = run_enrichment_2.delay(("genelist_2_" + session_id + ".txt"),pval_enr,("/code/polls/data/test3/enrichr_kegg/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/polls/data/test/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/polls/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/polls/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/polls/data/test2/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),("polls/data/test3/enrichr_kegg/" + session_id + "/KEGG_2013.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_enrichment"):	
				result1 = run_go_enrichment_2.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("genelist_1_" + session_id + ".txt"),pval_enr,("/code/polls/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("genelist_2_" + session_id + ".txt"),pval_enr,("/code/polls/data/test3/enrichr_go/" + session_id))
				#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/polls/data/test/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/polls/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/polls/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/codepolls/data/test2/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),("/code/polls/data/test3/enrichr_go/" + session_id + "/GO_Biological_Process_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()	
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "go_molecular"):
				result1 = run_go_enrichment_2.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test/enrichr_go/" + session_id))
				result2 = run_go_enrichment_2.delay(("genelist_1_" + session_id + ".txt"),pval_enr,("/code/polls/data/test2/enrichr_go/" + session_id))
				result3 = run_go_enrichment_2.delay(("genelist_2_" + session_id + ".txt"),pval_enr,("/code/polls/data/test3/enrichr_go/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/polls/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/polls/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/polls/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result7 = read_kegg_enrichment_2.delay(("/code/polls/data/test2/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),("/code/polls/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				(enrichment_dict_4,enrichment_dict_5) = result7.get()
			elif(enr_type == "reactome_enrichment"):
				result1 = run_reac_enrichment.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test/enrichr_reactome/" + session_id))
				result2 = run_reac_enrichment.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test2/enrichr_reactome/" + session_id))
				result3 = run_reac_enrichment.delay(("genelist_" + session_id + ".txt"),pval_enr,("/code/polls/data/test3/enrichr_reactome/" + session_id))
				enr_results = result1.get()
				enr_results_2 = result2.get()
				enr_results_3 = result3.get()
				print("enr")
				result4 = read_kegg_enrichment.delay(("/code/polls/data/test/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				result5 = read_kegg_enrichment.delay(("/code/polls/data/test2/enrichr_reactome/" + session_id + "/Reactome_2016.test_name.enrichr.reports.txt"),pval_enr)
				result6 = read_kegg_enrichment.delay(("/code/polls/data/test3/enrichr_go/" + session_id + "/GO_Molecular_Function_2018.test_name.enrichr.reports.txt"),pval_enr)
				enrichment_dict = result4.get()
				#enrichment_dict = {}
				enrichment_dict_2 = result5.get()
				enrichment_dict_3 = result6.get()
				#enrichment_dict_3 = {}		
				#return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})
			#path99 = "/home/quirin/testproject/polls/static/test_" + session_id + ".png"
			path99 = "test_" + session_id + ".png"
			json_path = "test15_" + session_id + ".json"
			#path_metadata = "polls/static/metadata_" + session_id + ".txt"
			path_metadata = "/code/polls/static/metadata_" + session_id + ".txt"
			metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 			
			#output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
			output_plot_path = "output_plotly_" + session_id + ".html"
			return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'path4':path99,
'json_path':json_path,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5,'enrichment_open':"true"})
		
	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			#metd = list_metadata_3.apply_async(countdown=0)
			session_id = request.session._get_or_create_session_key()
			#path_metadata = "polls/static/metadata_" + session_id + ".txt"
			#metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)			
			#(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			#print(ret_metadata1) 
			print("iubaerb")
			#metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
			print("clustering_6_part_3_2")
		session_id = request.session._get_or_create_session_key()
		path99 = "test_" + session_id + ".png"
		json_path = "test15_" + session_id + ".json"
		path_metadata = "metadata_" + session_id + ".txt"
		output_plot_path = "output_plotly_" + session_id + ".html"
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'path4':path99,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})

			#return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2) 				
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})




def infopage(request):			
	return render(request,'polls/infopage.html')




