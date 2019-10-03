import datetime
from django.db import models
from django.utils import timezone
# Create your models here.
import os, sys
from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic
##from scipy.stats.stats import pearsonr
##import plotly
##import plotly.plotly as py
##import plotly.graph_objs as go
##import plotly.io as pio
##import plotly.offline
#from .models import Choice, Question
from datetime import datetime
##from networkx.readwrite import json_graph
##import json
from django.shortcuts import render_to_response,render
from django.template import RequestContext
##from django.http import HttpResponseRedirect
##from django.urls import reverse
import clustering
#from polls.models import Document
##from clustering.forms import DocumentForm
#from polls.models import Upload,UploadForm
##import numpy as np
##import matplotlib.pyplot as plt
#import mpld3

##import seaborn as sns
##import pandas as pd
##from numpy import array

import matplotlib.patches as mpatches


#import networkx as nx
#from bokeh.io import show, output_notebook, output_file, save
#from bokeh.plotting import figure
#from bokeh.models import Circle, HoverTool, TapTool, BoxSelectTool
#from bokeh.models.graphs import from_networkx
#from bokeh.transform import linear_cmap
#from bokeh.models import ColumnDataSource, LabelSet
#from bokeh.models.graphs import NodesAndLinkedEdges, EdgesAndLinkedNodes
#from biomart import BiomartServer
#from bokeh.embed import components
#from bokeh.palettes import Spectral4
#from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxSelectTool

##from pybiomart import Dataset

from django.forms import ModelForm

#class Upload(models.Model):
#        pic = models.FileField(upload_to="images/")
#        upload_date=models.DateTimeField(auto_now_add =True)

# FileUpload form class.
#class UploadForm(ModelForm):
#        class Meta:
#                model = Upload
#                fields = ('pic',)

class GraphForm(models.Model):
	#def makehref(term):
	#	ret = term + ""
	#	ret = "<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term="+ret+"\">"+ret+"</a>"
	#	return(ret)
	#def is_logged_in(username,password):
	#	return(1)
	
	def save_user_data(fn,prot_fn,username):
		user_dir = "user_uploaded_files/" + username
		if not(os.path.isdir(user_dir)):
			os.mkdir(user_dir)
		fn.seek(0)
		prot_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = user_dir + "/" + filename_1 + "_expr.txt"
		filename_3 = user_dir + "/" + filename_1 + "_prot.txt"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
	def save_user_data_2(fn,prot_fn,clinical_fn,username):
		user_dir = "user_uploaded_files/" + username
		if not(os.path.isdir(user_dir)):
			os.mkdir(user_dir)
		fn.seek(0)
		prot_fn.seek(0)
		clinical_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		str3 = clinical_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = user_dir + "/" + filename_1 + "_expr.txt"
		filename_3 = user_dir + "/" + filename_1 + "_prot.txt"
		filename_4 = user_dir + "/" + filename_1 + "_clin.txt"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
		outfile3 = open(filename_4, "w")
		outfile3.write(str3)
		outfile3.close()
	def save_user_data_3(exprstr,ppistr,clinicalstr,username):
		user_dir = "user_uploaded_files/" + username
		# check if directory with stored files for user exists
		if not(os.path.isdir(user_dir)):
			os.mkdir(user_dir)
		#if not(os.path.isdir(foobar + "/bla")):
		#	os.mkdir(foobar + "/bla")	
		if not(os.path.isdir("/code/user_uploaded_files/" + username)):
			os.mkdir("/code/user_uploaded_files/" + username)	
		if(clinicalstr == ""):
			clinicalstr = "empty"
		print("saving files")
		#fn.seek(0)
		#prot_fn.seek(0)
		#clinical_fn.seek(0)
		#str1 = fn.read().decode('utf-8')
		#str2 = prot_fn.read().decode('utf-8')
		#str3 = clinical_fn.read().decode('utf-8')
		#print(str1)
		# make filenames with current date
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = user_dir + "/" + filename_1 + "_expr.txt"
		filename_3 = user_dir + "/" + filename_1 + "_prot.txt"
		filename_4 = user_dir + "/" + filename_1 + "_clin.txt"
		# write strings to files
		outfile1 = open(filename_2, "w")
		outfile1.write(exprstr)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(ppistr)
		outfile2.close()
		outfile3 = open(filename_4, "w")
		outfile3.write(clinicalstr)
		outfile3.close()


	def save_results(username):
		user_dir = "user_uploaded_files/" + username
		if not(os.path.isdir(user_dir)):
			os.mkdir(user_dir)
		fn.seek(0)
		prot_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = user_dir + "/" + filename_1 + "_json.json"
		filename_3 = user_dir + "/" + filename_1 + "_heatmap.png"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
	
	#def list_user_data(username):
	#	foobar = "user_uploaded_files/" + username
	#	fileslist = os.listdir(foobar)
	#	print(fileslist)
	#	bar = []
	#	for f in fileslist:
	#		if "_expr.txt" in f:
	#			print(f)
	#			ret1 = "user_uploaded_files/" + username + "/" + f
	#			fn_temp = f.split("_expr.txt")[0] + "_prot.txt"
	#			ret2 = "user_uploaded_files/" + username + "/" + fn_temp
	#			bar.append({'f1':ret1,'f2':ret2})
	#			#time_temp = f.split("_")
	#			#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
		return bar
	#def list_user_data(username):
	#	foobar = "user_uploaded_files/" + username
	#	fileslist = os.listdir(foobar)
	#	print(fileslist)
	#	bar = []
	#	for f in fileslist:
	#		if "_expr.txt" in f:
	#			print(f)
	#			ret1 = "user_uploaded_files/" + username + "/" + f
	#			fn_temp = f.split("_expr.txt")[0] + "_prot.txt"
	#			ret2 = "user_uploaded_files/" + username + "/" + fn_temp
	#			bar.append({'f1':ret1,'f2':ret2})
	#			#time_temp = f.split("_")
	#			#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
	#	return bar


	def list_user_data(username):
		# function to list all files with stored input files for a given user.
		print("listing user data")
		user_dir = "user_uploaded_files/" + username
		fileslist = os.listdir(user_dir)
		print(fileslist)
		bar = []
		for f in fileslist:
			if "_expr.txt" in f:
				#print(f)
				# get date from filename
				f_split=f.split("_")
				f_new = f_split[0] + "/" + f_split[1] + "/" + f_split[2] + ", " + f_split[3] + ":" + f_split[4]
				print(f_new)
				ret1 = "user_uploaded_files/" + username + "/" + f
				# get name of PPI file
				fn_temp = f.split("_expr.txt")[0] + "_prot.txt"
				ret2 = "user_uploaded_files/" + username + "/" + fn_temp
				#f1: the date in nice format, f2: the filename (used as input variable in form)
				bar.append({'f1':f_new,'f2':ret1})
				#time_temp = f.split("_")
				#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
		return bar
	

	def list_user_data_2(username):
		# function to list all files with stored result files (heatmap, PPI) for a given user.
		print("listing user data 2")
		user_dir = "user_uploaded_files/" + username
		fileslist = os.listdir(user_dir)
		bar = []
		for f in fileslist:
			if "_json.json" in f:
				# get date from filename
				f_split=f.split("_")				
				f_new = f_split[0] + "/" + f_split[1] + "/" + f_split[2] + ", " + f_split[3] + ":" + f_split[4]
				print(f_new)
				ret1 = "user_uploaded_files/" + username + "/" + f
				# get name of heatmap file
				fn_temp = f.split("_json.json")[0] + "_heatmap.png"
				ret2 = "user_uploaded_files/" + username + "/" + fn_temp
				#f1: the date in nice format, f2: the filename (used as input variable in form)
				bar.append({'f1':f_new,'f2':ret1})
				#time_temp = f.split("_")
				#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
		return bar
	
	
