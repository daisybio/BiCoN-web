##### LOTS OF OLD CODE!
import datetime
from django.db import models
from django.utils import timezone
# Create your models here.

from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic
from scipy.stats.stats import pearsonr
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline
#from .models import Choice, Question


from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.urls import reverse
import polls
#from polls.models import Document
from polls.forms import DocumentForm
#from polls.models import Upload,UploadForm
import numpy as np
import matplotlib.pyplot as plt
import mpld3

import seaborn as sns
import pandas as pd
from numpy import array

import matplotlib.patches as mpatches


import networkx as nx
from bokeh.io import show, output_notebook, output_file, save
from bokeh.plotting import figure
from bokeh.models import Circle, HoverTool, TapTool, BoxSelectTool
from bokeh.models.graphs import from_networkx
from bokeh.transform import linear_cmap
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.models.graphs import NodesAndLinkedEdges, EdgesAndLinkedNodes
from biomart import BiomartServer
from bokeh.embed import components
from bokeh.palettes import Spectral4
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxSelectTool

from pybiomart import Dataset

class Question(models.Model):
    question_text = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')
    def __str__(self):
        return self.question_text
    def was_published_recently(self):
        return self.pub_date >= timezone.now() - datetime.timedelta(days=1)



class Choice(models.Model):
    question = models.ForeignKey(Question, on_delete=models.CASCADE)
    choice_text = models.CharField(max_length=200)
    votes = models.IntegerField(default=0)
    def __str__(self):
        return self.choice_text


class Document(models.Model):
#	docfile = models.FileField(upload_to='documents/%Y/%m/%d')
	docfile = models.FileField(upload_to='/home/quirin/bla')
	class Meta:
		app_label = 'Document'

class DocumentForm(models.Model):
	docfile = models.FileField(
	#label='Select a file',
	help_text='max. 42 megabytes'
	)
	class Meta:
		app_label = 'DocumentForm'

from django.forms import ModelForm

class Upload(models.Model):
        pic = models.FileField(upload_to="images/")
        upload_date=models.DateTimeField(auto_now_add =True)

# FileUpload form class.
class UploadForm(ModelForm):
        class Meta:
                model = Upload
                fields = ('pic',)

class GraphForm(models.Model):
	def handle_upload_3(fn,prot_fn):
		patients = []	
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		for line in fn.readlines():
			lineSplit = line.decode('utf8').split()
			#patients.append(['1','2'])
			#patients1.append(['1','2'])
			#patients.append(['3','4'])
			#patients1.append(['3','4'])
			#group_labels1.append([1,1])
			#group_labels2.append([2,2])
			#ctr = ctr + 1
			geneName = lineSplit[0]
			nbr = lineSplit[2]
			if(float(nbr) > 100.0):
				if(geneName in data1):
					data1[geneName].append(float(lineSplit[1]))
					data1[geneName].append(float(lineSplit[1]))
					data2[geneName].append(float(lineSplit[2]))
					data2[geneName].append(float(lineSplit[2]))
				else:
					data1[geneName] = []
					data2[geneName] = []
					data1[geneName].append(float(lineSplit[1]))
					data1[geneName].append(float(lineSplit[1]))
					data2[geneName].append(float(lineSplit[2]))
					data2[geneName].append(float(lineSplit[2]))
					logstring = logstring + str(data1[geneName][1])
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])
		genes = ()
	
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		for i in data1:
			#print(data1[1000])
			if(i in data2):
				if(1 == 1):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					group1.append(data1[i])
					group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					ensg_id = i.split('.')[0]
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,len(data1[i])):
						diff_curr = diff_curr + data1[i][j] - data2[i][j]
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
	                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		dataArray1 = np.array(group1)
		dataArray2 = np.array(group2)
		endData = np.concatenate((dataArray1, dataArray2), axis=1)
		logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		logstring = logstring + str(df)
		sns.clustermap(df,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		##print(logstring)
		script, div = components(plot)
		return(script,div,plt.gcf())
	def handle_upload_3_2(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		for line in fn.readlines():
			lineSplit = line.decode('utf8').split()
			#patients.append(['1','2'])
			#patients1.append(['1','2'])
			#patients.append(['3','4'])
			#patients1.append(['3','4'])
			#group_labels1.append([1,1])
			#group_labels2.append([2,2])
			#ctr = ctr + 1
			geneName = lineSplit[0]
			nbr = lineSplit[2]
			ctr_bar = ctr_bar + 1
			if(ctr_bar > 2):
				#if(float(nbr) > 100.0):
				if(1==1):
					if(geneName in data1):
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
					else:
						data1[geneName] = []
						data2[geneName] = []
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
						logstring = logstring + str(data1[geneName][1])
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])
		genes = ()
	
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		for i in data1:
			#print(data1[1000])
			if(i in data2):
				if(1 == 1):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					group1.append(data1[i])
					group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					ensg_id = i.split('.')[0]
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,len(data1[i])):
						diff_curr = diff_curr + data1[i][j] - data2[i][j]
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		dataArray1 = np.array(group1)
		dataArray2 = np.array(group2)
		endData = np.concatenate((dataArray1, dataArray2), axis=1)
		logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		logstring = logstring + str(df)
		df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		print(foobar)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df2.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>500.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)
		return(script,div,plt.gcf())
	def handle_upload_3_3(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		for line in fn.readlines():
			lineSplit = line.decode('utf8').split()
			#patients.append(['1','2'])
			#patients1.append(['1','2'])
			#patients.append(['3','4'])
			#patients1.append(['3','4'])
			#group_labels1.append([1,1])
			#group_labels2.append([2,2])
			#ctr = ctr + 1
			geneName = lineSplit[0]
			nbr = lineSplit[2]
			ctr_bar = ctr_bar + 1
			if(ctr_bar > 2):
				#if(float(nbr) > 100.0):
				if(1==1):
					if(geneName in data1):
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
					else:
						data1[geneName] = []
						data2[geneName] = []
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
						logstring = logstring + str(data1[geneName][1])
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					#print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					#print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_3 and prot1 in genes_3):
				G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		print(df2)
		#print(foobar)
		#print(colors)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,2]>500.0]
		#print(df7)
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		#print(logstring)
		script, div = components(plot)
		return(script,div,plt.gcf())
	def handle_upload_3_4(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		genes_4 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		conv_dict_2 = {}
		#print(ensg_id_list)
		df2shape = df2.shape[1]-1
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			#if(float(df2.loc[i].tolist()[2]) > 100.0):
			ensg_id = i.split('.')[0]
			if(ensg_id in conv_dict):
				prot_id = conv_dict[ensg_id]['Gene_name']
				conv_dict_2[prot_id]=i
				genes_4.update({prot_id:0})
		G = nx.Graph()
		nodes = []
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		#dataset.list_attributes()
		genes = {}
		genes_3 = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_4):
				nodes.append(prot1)
				#logstring = logstring + str(genes_3[prot1])
				diff = 0
				diff_curr = 0
				gene_tmp = conv_dict_2[prot1]
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot1, Name=prot1, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:diff_curr_2})
			if(prot2 not in nodes and prot2 in genes_4):
				nodes.append(prot2)
				#G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				#genes.update({prot2:genes_3[prot2]})
				gene_tmp = conv_dict_2[prot2]
				diff = 0
				diff_curr = 0
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot2, Name=prot2, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot2:diff_curr_2})
			#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			#if(prot2 in genes_3 and prot1 in genes_3):
			G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		#logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		z = tuple([(bar-0.1) for bar in x])
		#print(x)
		#print(z)
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='z', y='y', text='Name', source=source,
                #  background_fill_color='white',level='glyph',render_mode='canvas')
		labels2 = LabelSet(x='x', y='y', text='Name',x_offset=-30, source=source,
                  background_fill_color='white',level='glyph',background_fill_alpha=0.0,render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		#logstring = logstring + "\n\nPPI Graph created..."
		#output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		print(df2)
		#print(foobar)
		#print(colors)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		#print(df7)
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		#print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		plt.close()
		return(script,div,plot1)
	def handle_upload_3_survival(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		for line in fn.readlines():
			lineSplit = line.decode('utf8').split()
			#patients.append(['1','2'])
			#patients1.append(['1','2'])
			#patients.append(['3','4'])
			#patients1.append(['3','4'])
			#group_labels1.append([1,1])
			#group_labels2.append([2,2])
			#ctr = ctr + 1
			geneName = lineSplit[0]
			nbr = lineSplit[2]
			ctr_bar = ctr_bar + 1
			if(ctr_bar > 2):
				#if(float(nbr) > 100.0):
				if(1==1):
					if(geneName in data1):
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
					else:
						data1[geneName] = []
						data2[geneName] = []
						data1[geneName].append(float(lineSplit[1]))
						data1[geneName].append(float(lineSplit[1]))
						data2[geneName].append(float(lineSplit[2]))
						data2[geneName].append(float(lineSplit[2]))
						logstring = logstring + str(data1[geneName][1])
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])
		genes = ()
	
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		for i in data1:
			#print(data1[1000])
			if(i in data2):
				if(1 == 1):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					group1.append(data1[i])
					group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					ensg_id = i.split('.')[0]
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,len(data1[i])):
						diff_curr = diff_curr + data1[i][j] - data2[i][j]
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		dataArray1 = np.array(group1)
		dataArray2 = np.array(group2)
		endData = np.concatenate((dataArray1, dataArray2), axis=1)
		logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		logstring = logstring + str(df)
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>500.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		return(script,div,plot_1,plot_div)
	def handle_upload_3_survival_2(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0		
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])		
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		#diff_frame = pd.DataFrame(list(genes.items()), columns=['Gene', 'Difference'])
		genelist_diff = []
		#diff_frame_2 = diff_frame.sort_values('Difference')
		#for i in range(0,len(diff_frame_2)):
		#	curr_tmp = diff_frame_2.iloc[i,0]
		#	print(curr_tmp)
		#	genelist_diff.append({'gene_name':curr_tmp,'difference':diff_frame_2.loc[i,'Difference']})
		for k,v in list(genes.items()):
			#print(k)
			#print(v)
			genelist_diff.append({'gene_name':k,'difference':v})
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		#print(colors)			
		#df22 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		#colors = [colordict[i] for i in foobar]
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		print(genelist_diff)
		return(script,div,plot_1,plot_div,genelist_diff)
	def handle_upload_3_survival_3(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0		
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])		
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		survival_yrs =  df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		genes_corr = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						tmp_3 = df2.loc[i].tolist()
						tmp_4 = [float(x) for x in tmp_3]
						tmp_corr = pearsonr(tmp_4,survival_yrs)
						print(tmp_corr)
						genes_corr.update({prot_id:tmp_corr})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		#diff_frame = pd.DataFrame(list(genes.items()), columns=['Gene', 'Difference'])
		genelist_diff = []
		#diff_frame_2 = diff_frame.sort_values('Difference')
		#for i in range(0,len(diff_frame_2)):
		#	curr_tmp = diff_frame_2.iloc[i,0]
		#	print(curr_tmp)
		#	genelist_diff.append({'gene_name':curr_tmp,'difference':diff_frame_2.loc[i,'Difference']})
		genelist_ret =[]
		for k,v in list(genes.items()):
			#print(k)
			#print(v)
			genelist_diff.append({'gene_name':k,'difference':v})
			genelist_ret.append({'gene_name':k,'difference':v,'correlation':genes_corr[k]})
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		#print(colors)			
		#df22 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		#colors = [colordict[i] for i in foobar]
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		print(genelist_diff)
		plt.close()
		return(script,div,plot_1,plot_div,genelist_ret)
	def handle_upload_3_survival_4(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0		
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])		
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		survival_yrs =  df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		genes_corr = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						tmp_3 = df2.loc[i].tolist()
						tmp_4 = [float(x) for x in tmp_3]
						tmp_corr = pearsonr(tmp_4,survival_yrs)
						print(tmp_corr)
						genes_corr.update({prot_id:tmp_corr})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		#diff_frame = pd.DataFrame(list(genes.items()), columns=['Gene', 'Difference'])
		genelist_diff = []
		#diff_frame_2 = diff_frame.sort_values('Difference')
		#for i in range(0,len(diff_frame_2)):
		#	curr_tmp = diff_frame_2.iloc[i,0]
		#	print(curr_tmp)
		#	genelist_diff.append({'gene_name':curr_tmp,'difference':diff_frame_2.loc[i,'Difference']})
		genelist_ret =[]
		for k,v in list(genes.items()):
			#print(k)
			#print(v)
			genelist_diff.append({'gene_name':k,'difference':v})
			genelist_ret.append({'gene_name':k,'difference':v,'correlation':genes_corr[k]})
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		#print(colors)			
		#df22 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		#colors = [colordict[i] for i in foobar]
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		print(genelist_diff)
		return(script,div,plot_1,plot_div,genelist_ret)
	def handle_upload_5(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		#group_labels1 = []
		#group_labels2 = []
		group1_data = []
		group2_data = []
		
		#patients.append([1,2,3,4])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients2.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#for line in fn.readlines():
		#	lineSplit = line.split()
			#patients.append(['1','2'])
			#patients1.append(['1','2'])
			#patients.append(['3','4'])
			#patients1.append(['3','4'])
			#group_labels1.append([1,1])
			#group_labels2.append([2,2])
			#ctr = ctr + 1
		#	geneName = lineSplit[0]
		#	nbr = lineSplit[2]
		#	ctr_bar = ctr_bar + 1
		#	if(ctr_bar > 2):
		#		#if(float(nbr) > 100.0):
		#		if(1==1):
		#			if(geneName in data1):
		#				data1[geneName].append(float(lineSplit[1]))
		#				data1[geneName].append(float(lineSplit[1]))
		#				data2[geneName].append(float(lineSplit[2]))
		#				data2[geneName].append(float(lineSplit[2]))
		#			else:
		#				data1[geneName] = []
		#				data2[geneName] = []
		#				data1[geneName].append(float(lineSplit[1]))
		#				data1[geneName].append(float(lineSplit[1]))
		#				data2[geneName].append(float(lineSplit[2]))
		#				data2[geneName].append(float(lineSplit[2]))
		#				logstring = logstring + str(data1[geneName][1])
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		#print(ensg_id_list)
		conv.to_csv("Output3.txt", sep='\t')
		with open("Output.txt", "w") as text_file:
			#text_file.write(conv)
			text_file.write("\n".join(ensg_id_list))
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					#print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					#print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_3 and prot1 in genes_3):
				G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		z = tuple([(bar-0.1) for bar in x])
		#print(x)
		#print(z)
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='z', y='y', text='Name', source=source,
                #  background_fill_color='white',level='glyph',render_mode='canvas')
		labels2 = LabelSet(x='x', y='y', text='Name',x_offset=-30, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		#print(foobar)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)	
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		return(script,div,plot_1,plot_div)
	def handle_upload_6(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group1_data = []
		group2_data = []
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		genes = ()
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])	
		print(df2)
		#df16 = df2.iloc[5:,:]
		#df17 = df16.loc[df2[1]>100.0]
		#print(df17)
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		#print(ensg_id_list)
		df2shape = df2.shape[1]-1
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
			#if(1==1):
				ensg_id = i.split('.')[0]
				if(ensg_id in conv_dict):
				#if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					#print(ensg_id)
					#prot_nbr = conv.index[conv['Gene_stable_ID'] == ensg_id].tolist()
					#print(prot_nbr)
					#prot_id = conv.loc[prot_nbr,'Gene_name'].values[0]
					prot_id = conv_dict[ensg_id]['Gene_name']
					#prot_id = conv.loc[ensg_id,'Gene_name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					#if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					#logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		#dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				#logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_3 and prot1 in genes_3):
				G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		z = tuple([(bar-0.1) for bar in x])
		#print(x)
		#print(z)
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='z', y='y', text='Name', source=source,
                #  background_fill_color='white',level='glyph',render_mode='canvas')
		labels2 = LabelSet(x='x', y='y', text='Name',x_offset=-30, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		#print(foobar)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)	
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div')
		return(script,div,plot_1,plot_div)	
	def handle_upload_7(fn,prot_fn):
		patients = []
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		group1 = []
		group2 = []
		group1_data = []
		group2_data = []
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		genes = ()
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		genes_4 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])	
		#print(df2)
		#df16 = df2.iloc[5:,:]
		#df17 = df16.loc[df2[1]>100.0]
		#print(df17)
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		conv_dict_2 = {}
		conv_dict_3 = {}
		#print(ensg_id_list)
		df2shape = df2.shape[1]-1
		#ctr99 = 4
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			#if(float(df2.loc[i].tolist()[2]) > 100.0):
			ensg_id = i.split('.')[0]
			if(ensg_id in conv_dict):
				prot_id = conv_dict[ensg_id]['Gene_name']
				conv_dict_2[prot_id]=i	
				#conv_dict_3[i]=prot_id
				#df2.iloc[ctr99,0] = prot_id
				genes_4.update({prot_id:0})
		G = nx.Graph()
		nodes = []
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		#dataset.list_attributes()
		genes = {}
		genes_3 = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_4):
				nodes.append(prot1)
				#logstring = logstring + str(genes_3[prot1])
				diff = 0
				diff_curr = 0
				gene_tmp = conv_dict_2[prot1]
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot1, Name=prot1, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:diff_curr_2})
			if(prot2 not in nodes and prot2 in genes_4):
				nodes.append(prot2)
				#G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				#genes.update({prot2:genes_3[prot2]})
				gene_tmp = conv_dict_2[prot2]
				diff_curr = 0
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot2, Name=prot2, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot2:diff_curr_2})
				#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_4 and prot1 in genes_4):
				G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		#logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		z = tuple([(bar-0.1) for bar in x])
		#print(x)
		#print(z)
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='z', y='y', text='Name', source=source,
                #  background_fill_color='white',level='glyph',render_mode='canvas')
		labels2 = LabelSet(x='x', y='y', text='Name',x_offset=-30, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		#logstring = logstring + "\n\nPPI Graph created..."
		#   output_file("polls/interactive_graphs.html")
		#   save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		#print(foobar)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,3]>100.0]
		#df7.iloc[:,0].replace(conv_dict_3, inplace=True)
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),robust=True,col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)	
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		plt.close()
		return(script,div,plot_1,plot_div)	


#### old code from views.py


def handle_upload_3(fn,prot_fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	for line in fn.readlines():
		lineSplit = line.decode('utf8').split()
		#patients.append(['1','2'])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients1.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		#ctr = ctr + 1
		geneName = lineSplit[0]
		nbr = lineSplit[2]
		if(float(nbr) > 100.0):
			if(geneName in data1):
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
			else:
				data1[geneName] = []
				data2[geneName] = []
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
				logstring = logstring + str(data1[geneName][1])
	diff = {}
	ctr2 = 0
	#group1.append(group_labels1)
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()

	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	for i in data1:
		#print(data1[1000])
		if(i in data2):
			if(1 == 1):
			#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
				#print(data1[1000])
				ctr2 = ctr2 + 1
				#print(i)
				#if(data1[i][0] > 10000.0):
				#	print(data1[i])
				geneNames.append(ctr2)
				group1.append(data1[i])
				group2.append(data2[i])
				#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
				ensg_id = i.split('.')[0]
				prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
				prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
				#print(prot_id)
				diff_curr = 0
				for j in range(0,len(data1[i])):
					diff_curr = diff_curr + data1[i][j] - data2[i][j]
				if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
	
	G = nx.Graph()
	nodes = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
	#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
	dataset.list_attributes()
	genes = {}
	logstring = logstring + "\n\nReading PPI file...\n"
	for line in prot_fn.readlines():
		lineSplit = line.decode('utf8').split()
		prot1 = lineSplit[0]
		prot2 = lineSplit[1]
		logstring = logstring + prot1 + prot2
		if(prot1 not in nodes and prot1 in genes_3):
			nodes.append(prot1)
			logstring = logstring + str(genes_3[prot1])
			G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
			genes.update({prot1:genes_3[prot1]})
		if(prot2 not in nodes and prot2 in genes_3):
			nodes.append(prot2)
			G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
			genes.update({prot2:genes_3[prot2]})
		logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
		G.add_edge(prot1,prot2)
	output_notebook()
	plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
               TapTool(), 
               BoxSelectTool())
	# create bokeh graph
	graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
	
	# add name to node data
	#graph.node_renderer.data_source.data['d'] = list(G.nodes())
	
	
	graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
	
	# add name to node data
	graph.node_renderer.data_source.data['Name'] = list(G.nodes())
	x,y = zip(*graph.layout_provider.graph_layout.values())
	#x = [i for i in G.nodes]
	#y = [genes[i] for i in G.nodes]
	node_labels = nx.get_node_attributes(G, 'Name')
	logstring = logstring + str(node_labels)
	#source = ColumnDataSource({'x': x, 'y': y,
	 #                          'Name': [node_labels[i] for i in range(len(x))]})
	
	source = ColumnDataSource(data=dict(x=x, y=y, 
                           Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',level='glyph',render_mode='canvas')
	#print(labels2)
	#logstring = logstring + str(graph.layout_provider.graph_layout.values())
	# add club to node data
	#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
	
	# set node size
	#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
	graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
	
	plot.renderers.append(graph)
	plot.renderers.append(labels2)
	#graph.inspection_policy = EdgesAndLinkedNodes()
	plot.add_layout(labels2)
	logstring = logstring + "\n\nPPI Graph created..."
	output_file("polls/interactive_graphs.html")
	save(plot)
	dataArray1 = np.array(group1)
	dataArray2 = np.array(group2)
	endData = np.concatenate((dataArray1, dataArray2), axis=1)
	logstring = logstring + str(endData)
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	logstring = logstring + str(df)
	sns.clustermap(df,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
	#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
	#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
	#plt.yticks(range(data_present_length), geneNames, size='small')
	#plt.xticks(range(len(patients)), patients, size='small')
	plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
	#plt.switch_backend('SVG') 
	plt.savefig("/home/quirin/mysite/static/test.png")
	print(logstring)
	script, div = components(plot)
	return(script,div,plt.gcf())


def handle_upload_3_2(fn,prot_fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	ctr_bar = 0
	for line in fn.readlines():
		lineSplit = line.decode('utf8').split()
		#patients.append(['1','2'])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients1.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		#ctr = ctr + 1
		geneName = lineSplit[0]
		nbr = lineSplit[2]
		ctr_bar = ctr_bar + 1
		if(ctr_bar > 2):
			#if(float(nbr) > 100.0):
			if(1==1):
				if(geneName in data1):
					data1[geneName].append(float(lineSplit[1]))
					data1[geneName].append(float(lineSplit[1]))
					data2[geneName].append(float(lineSplit[2]))
					data2[geneName].append(float(lineSplit[2]))
				else:
					data1[geneName] = []
					data2[geneName] = []
					data1[geneName].append(float(lineSplit[1]))
					data1[geneName].append(float(lineSplit[1]))
					data2[geneName].append(float(lineSplit[2]))
					data2[geneName].append(float(lineSplit[2]))
					logstring = logstring + str(data1[geneName][1])
	diff = {}
	ctr2 = 0
	#group1.append(group_labels1)
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()

	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	for i in data1:
		#print(data1[1000])
		if(i in data2):
			if(1 == 1):
			#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
				#print(data1[1000])
				ctr2 = ctr2 + 1
				#print(i)
				#if(data1[i][0] > 10000.0):
				#	print(data1[i])
				geneNames.append(ctr2)
				group1.append(data1[i])
				group2.append(data2[i])
				#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
				ensg_id = i.split('.')[0]
				prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
				prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
				#print(prot_id)
				diff_curr = 0
				for j in range(0,len(data1[i])):
					diff_curr = diff_curr + data1[i][j] - data2[i][j]
				if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
	
	G = nx.Graph()
	nodes = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
	#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
	dataset.list_attributes()
	genes = {}
	logstring = logstring + "\n\nReading PPI file...\n"
	for line in prot_fn.readlines():
		lineSplit = line.decode('utf8').split()
		prot1 = lineSplit[0]
		prot2 = lineSplit[1]
		logstring = logstring + prot1 + prot2
		if(prot1 not in nodes and prot1 in genes_3):
			nodes.append(prot1)
			logstring = logstring + str(genes_3[prot1])
			G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
			genes.update({prot1:genes_3[prot1]})
		if(prot2 not in nodes and prot2 in genes_3):
			nodes.append(prot2)
			G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
			genes.update({prot2:genes_3[prot2]})
		logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
		G.add_edge(prot1,prot2)
	output_notebook()
	plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
               TapTool(), 
               BoxSelectTool())
	# create bokeh graph
	graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
	
	# add name to node data
	#graph.node_renderer.data_source.data['d'] = list(G.nodes())
	
	
	graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
	
	# add name to node data
	graph.node_renderer.data_source.data['Name'] = list(G.nodes())
	x,y = zip(*graph.layout_provider.graph_layout.values())
	#x = [i for i in G.nodes]
	#y = [genes[i] for i in G.nodes]
	node_labels = nx.get_node_attributes(G, 'Name')
	logstring = logstring + str(node_labels)
	#source = ColumnDataSource({'x': x, 'y': y,
	 #                          'Name': [node_labels[i] for i in range(len(x))]})
	
	source = ColumnDataSource(data=dict(x=x, y=y, 
                           Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',level='glyph',render_mode='canvas')
	#print(labels2)
	#logstring = logstring + str(graph.layout_provider.graph_layout.values())
	# add club to node data
	#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
	
	# set node size
	#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
	graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
	
	plot.renderers.append(graph)
	plot.renderers.append(labels2)
	#graph.inspection_policy = EdgesAndLinkedNodes()
	plot.add_layout(labels2)
	logstring = logstring + "\n\nPPI Graph created..."
	output_file("polls/interactive_graphs.html")
	save(plot)
	dataArray1 = np.array(group1)
	dataArray2 = np.array(group2)
	endData = np.concatenate((dataArray1, dataArray2), axis=1)
	logstring = logstring + str(endData)
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	logstring = logstring + str(df)
	df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
	foobar = df2.loc['SURVIVAL',:]
	colors = [colordict[i] for i in foobar]
	#print(colors)
	df6 = df2.iloc[5:500,1:].apply(pd.to_numeric)
	sns.clustermap(df6,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
	#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
	#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
	#plt.yticks(range(data_present_length), geneNames, size='small')
	#plt.xticks(range(len(patients)), patients, size='small')
	plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
	#plt.switch_backend('SVG') 
	plt.savefig("/home/quirin/mysite/static/test.png")
	print(logstring)
	script, div = components(plot)
	return(script,div,plt.gcf())

def handle_upload_4(fn,prot_fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	for line in fn.readlines():
		lineSplit = line.split()
		#patients.append(['1','2'])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients1.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		#ctr = ctr + 1
		geneName = lineSplit[0]
		nbr = lineSplit[2]
		#if(float(nbr) > 100.0):
		if(1==1):
			if(geneName in data1):
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
			else:
				data1[geneName] = []
				data2[geneName] = []
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
				logstring = logstring + str(data1[geneName][1])
	diff = {}
	ctr2 = 0
	#group1.append(group_labels1)
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()

	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	for i in data1:
		#print(data1[1000])
		if(i in data2):
			if(1 == 1):
			#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
				#print(data1[1000])
				ctr2 = ctr2 + 1
				#print(i)
				#if(data1[i][0] > 10000.0):
				#	print(data1[i])
				geneNames.append(ctr2)
				group1.append(data1[i])
				group2.append(data2[i])
				#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
				ensg_id = i.split('.')[0]
				prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
				prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
				#print(prot_id)
				diff_curr = 0
				for j in range(0,len(data1[i])):
					diff_curr = diff_curr + data1[i][j] - data2[i][j]
				if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
	
	G = nx.Graph()
	nodes = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
	#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
	dataset.list_attributes()
	genes = {}
	logstring = logstring + "\n\nReading PPI file...\n"
	for line in prot_fn.readlines():
		lineSplit = line.split()
		prot1 = lineSplit[0]
		prot2 = lineSplit[1]
		logstring = logstring + prot1 + prot2
		if(prot1 not in nodes and prot1 in genes_3):
			nodes.append(prot1)
			logstring = logstring + str(genes_3[prot1])
			G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
			genes.update({prot1:genes_3[prot1]})
		if(prot2 not in nodes and prot2 in genes_3):
			nodes.append(prot2)
			G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
			genes.update({prot2:genes_3[prot2]})
		logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
		G.add_edge(prot1,prot2)
	output_notebook()
	plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
               TapTool(), 
               BoxSelectTool())
	# create bokeh graph
	graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
	
	# add name to node data
	#graph.node_renderer.data_source.data['d'] = list(G.nodes())
	
	
	graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
	
	# add name to node data
	graph.node_renderer.data_source.data['Name'] = list(G.nodes())
	x,y = zip(*graph.layout_provider.graph_layout.values())
	#x = [i for i in G.nodes]
	#y = [genes[i] for i in G.nodes]
	node_labels = nx.get_node_attributes(G, 'Name')
	logstring = logstring + str(node_labels)
	#source = ColumnDataSource({'x': x, 'y': y,
	 #                          'Name': [node_labels[i] for i in range(len(x))]})
	
	source = ColumnDataSource(data=dict(x=x, y=y, 
                           Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',level='glyph',render_mode='canvas')
	#print(labels2)
	#logstring = logstring + str(graph.layout_provider.graph_layout.values())
	# add club to node data
	#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
	
	# set node size
	#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
	graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
	#ADDED: SELECTION NOW WORKS AS NTENEDED, COPY TO handle..._3 WHEN NEEDED
	graph.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
	graph.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
	graph.inspection_policy = EdgesAndLinkedNodes()
	graph.selection_policy = NodesAndLinkedEdges()
	plot.renderers.append(graph)
	plot.renderers.append(labels2)
	plot.add_layout(labels2)
	logstring = logstring + "\n\nPPI Graph created..."
	output_file("polls/interactive_graphs.html")
	save(plot)
	dataArray1 = np.array(group1)
	dataArray2 = np.array(group2)
	endData = np.concatenate((dataArray1, dataArray2), axis=1)
	logstring = logstring + str(endData)
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	colordict={0:'#BB0000',1:'#0000BB'}
	df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	logstring = logstring + str(df)
	sns.clustermap(df,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
	#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
	#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
	#plt.yticks(range(data_present_length), geneNames, size='small')
	#plt.xticks(range(len(patients)), patients, size='small')
	plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
	#plt.switch_backend('SVG') 
	plt.savefig("/home/quirin/mysite/static/test.png")
	print(logstring)
	script, div = components(plot)
	return(script,div,plt.gcf())





def handle_upload_4(fn,prot_fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	for line in fn.readlines():
		lineSplit = line.split()
		#patients.append(['1','2'])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients1.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		#ctr = ctr + 1
		geneName = lineSplit[0]
		nbr = lineSplit[2]
		data_temp = [float(i) for i in lineSplit[1:]]
		if(float(nbr) > 100.0):
			if(geneName in data1):
				data1[geneName].append(data_temp)
				data1[geneName].append(float(lineSplit[1]))
			else:
				data1[geneName] = []
				data1[geneName].append(data_temp)
	diff = {}
	ctr2 = 0


	
	#group1.append(group_labels1)
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	for i in data1:
		#print(data1[1000])
		if(i in data2):
			if(1 == 1):
			#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
				#print(data1[1000])
				ctr2 = ctr2 + 1
				#print(i)
				#if(data1[i][0] > 10000.0):
				#	print(data1[i])
				geneNames.append(ctr2)
				group1.append(data1[i])
				group2.append(data2[i])
				#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
				ensg_id = i.split('.')[0]
				prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
				prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
				#print(prot_id)
				diff_curr = 0
				for j in range(0,len(data1[i])):
					diff_curr = diff_curr + data1[i][j] - data2[i][j]
				if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
	
	G = nx.Graph()
	nodes = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
	#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
	dataset.list_attributes()
	genes = {}
	logstring = logstring + "\n\nReading PPI file...\n"
	for line in prot_fn.readlines():
		lineSplit = line.split()
		prot1 = lineSplit[0]
		prot2 = lineSplit[1]
		logstring = logstring + prot1 + prot2
		if(prot1 not in nodes and prot1 in genes_3):
			nodes.append(prot1)
			logstring = logstring + str(genes_3[prot1])
			G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
			genes.update({prot1:genes_3[prot1]})
		if(prot2 not in nodes and prot2 in genes_3):
			nodes.append(prot2)
			G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
			genes.update({prot2:genes_3[prot2]})
		logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
		G.add_edge(prot1,prot2)
	output_notebook()
	plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
               TapTool(), 
               BoxSelectTool())
	# create bokeh graph
	graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
	
	# add name to node data
	#graph.node_renderer.data_source.data['d'] = list(G.nodes())
	
	
	graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
	
	# add name to node data
	graph.node_renderer.data_source.data['Name'] = list(G.nodes())
	x,y = zip(*graph.layout_provider.graph_layout.values())
	#x = [i for i in G.nodes]
	#y = [genes[i] for i in G.nodes]
	node_labels = nx.get_node_attributes(G, 'Name')
	logstring = logstring + str(node_labels)
	#source = ColumnDataSource({'x': x, 'y': y,
	 #                          'Name': [node_labels[i] for i in range(len(x))]})
	
	source = ColumnDataSource(data=dict(x=x, y=y, 
                           Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',level='glyph',render_mode='canvas')
	#print(labels2)
	#logstring = logstring + str(graph.layout_provider.graph_layout.values())
	# add club to node data
	#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
	
	# set node size
	#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
	graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
	#ADDED: SELECTION NOW WORKS AS NTENEDED, COPY TO handle..._3 WHEN NEEDED
	graph.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
	graph.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
	graph.inspection_policy = EdgesAndLinkedNodes()
	graph.selection_policy = NodesAndLinkedEdges()
	plot.renderers.append(graph)
	plot.renderers.append(labels2)
	plot.add_layout(labels2)
	logstring = logstring + "\n\nPPI Graph created..."
	output_file("polls/interactive_graphs.html")
	save(plot)
	dataArray1 = np.array(group1)
	dataArray2 = np.array(group2)
	endData = np.concatenate((dataArray1, dataArray2), axis=1)
	logstring = logstring + str(endData)
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	logstring = logstring + str(df)
	sns.clustermap(df,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
	#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
	#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
	#plt.yticks(range(data_present_length), geneNames, size='small')
	#plt.xticks(range(len(patients)), patients, size='small')
	plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
	#plt.switch_backend('SVG') 
	plt.savefig("/home/quirin/mysite/static/test.png")
	print(logstring)
	script, div = components(plot)
	return(script,div,plt.gcf())

def handle_upload_5(fn,prot_fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	#patient_ids = ['3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917','3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917']
	#patientfilename = 'nationwidechildrens.org_clinical_patient_brca.txt'
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	for line in fn.readlines():
		lineSplit = line.split()
		#patients.append(['1','2'])
		#patients1.append(['1','2'])
		#patients.append(['3','4'])
		#patients1.append(['3','4'])
		#group_labels1.append([1,1])
		#group_labels2.append([2,2])
		#ctr = ctr + 1
		geneName = lineSplit[0]
		nbr = lineSplit[2]
		if(float(nbr) > 100.0):
			if(geneName in data1):
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
			else:
				data1[geneName] = []
				data2[geneName] = []
				data1[geneName].append(float(lineSplit[1]))
				data1[geneName].append(float(lineSplit[1]))
				data2[geneName].append(float(lineSplit[2]))
				data2[geneName].append(float(lineSplit[2]))
				logstring = logstring + str(data1[geneName][1])
	diff = {}
	ctr2 = 0
	#group1.append(group_labels1)
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()

	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	for i in data1:
		#print(data1[1000])
		if(i in data2):
			if(1 == 1):
			#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
				#print(data1[1000])
				ctr2 = ctr2 + 1
				#print(i)
				#if(data1[i][0] > 10000.0):
				#	print(data1[i])
				geneNames.append(ctr2)
				group1.append(data1[i])
				group2.append(data2[i])
				#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
				ensg_id = i.split('.')[0]
				prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
				prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
				#print(prot_id)
				diff_curr = 0
				for j in range(0,len(data1[i])):
					diff_curr = diff_curr + data1[i][j] - data2[i][j]
				if(1==1):
					#group1.append(data1[i])
					#group2.append(data2[i])
					logstring = logstring + str(prot_id)
					diff_curr_2 = diff_curr / 2000.0
					diff.update({i:diff_curr_2})
					genes_3.update({prot_id:diff_curr_2})
					#print(prot_id)
	
	G = nx.Graph()
	nodes = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
	#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
	dataset.list_attributes()
	genes = {}
	logstring = logstring + "\n\nReading PPI file...\n"
	for line in prot_fn.readlines():
		lineSplit = line.split()
		prot1 = lineSplit[0]
		prot2 = lineSplit[1]
		logstring = logstring + prot1 + prot2
		if(prot1 not in nodes and prot1 in genes_3):
			nodes.append(prot1)
			logstring = logstring + str(genes_3[prot1])
			G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
			genes.update({prot1:genes_3[prot1]})
		if(prot2 not in nodes and prot2 in genes_3):
			nodes.append(prot2)
			G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
			genes.update({prot2:genes_3[prot2]})
		logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
		G.add_edge(prot1,prot2)
	output_notebook()
	plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
               TapTool(), 
               BoxSelectTool())
	# create bokeh graph
	graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
	
	# add name to node data
	#graph.node_renderer.data_source.data['d'] = list(G.nodes())
	
	
	graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
	
	# add name to node data
	graph.node_renderer.data_source.data['Name'] = list(G.nodes())
	x,y = zip(*graph.layout_provider.graph_layout.values())
	#x = [i for i in G.nodes]
	#y = [genes[i] for i in G.nodes]
	node_labels = nx.get_node_attributes(G, 'Name')
	logstring = logstring + str(node_labels)
	#source = ColumnDataSource({'x': x, 'y': y,
	 #                          'Name': [node_labels[i] for i in range(len(x))]})
	
	source = ColumnDataSource(data=dict(x=x, y=y, 
                           Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                  background_fill_color='white',level='glyph',render_mode='canvas')
	#print(labels2)
	#logstring = logstring + str(graph.layout_provider.graph_layout.values())
	# add club to node data
	#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
	
	# set node size
	#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
	graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
	#ADDED: SELECTION NOW WORKS AS NTENEDED, COPY TO handle..._3 WHEN NEEDED
	graph.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
	graph.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
	graph.inspection_policy = EdgesAndLinkedNodes()
	graph.selection_policy = NodesAndLinkedEdges()
	plot.renderers.append(graph)
	plot.renderers.append(labels2)
	plot.add_layout(labels2)
	logstring = logstring + "\n\nPPI Graph created..."
	output_file("polls/interactive_graphs.html")
	save(plot)
	dataArray1 = np.array(group1)
	dataArray2 = np.array(group2)
	endData = np.concatenate((dataArray1, dataArray2), axis=1)
	logstring = logstring + str(endData)
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	df9 = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	colordict={0:'#BB0000',1:'#0000BB'}
	logstring = logstring + str(df9)
	df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
	#print(df2.head())
	df = df2.transpose()
	print(df)
	df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
	train, test = df[df['is_train']==True], df[df['is_train']==False]
	features = df.columns[3:60483]
	y = pd.factorize(train['SURVIVAL'])[0]
	print(y)
	classif = RandomForestClassifier(n_jobs=2, random_state=0)
	classif.fit(train[features], y)
	print(classif.predict(df[features]))
	df3 = df
	df3['SURVIVAL'] = classif.predict(df[features])
	df4 = df3.sort_values(by=['SURVIVAL'])
	print(df4)
	df5_2 = df4.transpose()
	df5 = df5_2.sort_values(by=2)
	print(df5)
	foobar = df5.loc['SURVIVAL',:]
	colors = [colordict[i] for i in foobar]
	print(colors)
	df6 = df5.iloc[5:500,1:].apply(pd.to_numeric)
	sns.clustermap(df6,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
	#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
	#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
	#plt.yticks(range(data_present_length), geneNames, size='small')
	#plt.xticks(range(len(patients)), patients, size='small')
	plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
	#plt.switch_backend('SVG') 
	plt.savefig("/home/quirin/testproject/polls/static/test.png")
	print(logstring)
	script, div = components(plot)
	plot1 = plt.gcf()
	survival_yrs = df5.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
	survival_1 = {0:1.0}
	survival_2 = {0:1.0}
	survivor_nbr = sum(foobar)
	death_nbr = len(foobar) - sum(foobar)
	count1 = 1.0
	count2 = 1.0
	for j in range(1,len(survival_yrs)):
		if(foobar[j] == 0):
			survival_1[j] = count1 - (1/death_nbr)
			count1 = count1 - (1/death_nbr)
		if(foobar[j] == 1):
			survival_2[j] = count2 - (1/survivor_nbr)
			count2 = count2 - (1/survivor_nbr)
	trace1 = go.Scatter(
	x=list(survival_1.keys()),
	y=list(survival_1.values()),
	mode='lines+markers',
	name="'linear'",
	hoverinfo='name',
	line=dict(
	shape='hv'))
	trace2 = go.Scatter(
	x=list(survival_2.keys()),
	y=list(survival_2.values()),
	mode='lines+markers',
	name="'linear'",
	hoverinfo='name',
	line=dict(
	shape='hv'))
	#for i in range(0,25):
	#	temp1 = 0
	#	temp2 = 0
	#	for j in survival_yrs:
	#		if(foobar[j] == 0 and j >= i):
	#			temp1 = temp1+1
	#		if(foobar[j] == 1 and j >= i):
	#			temp2 = temp2+1
	#	survival_1[i] = temp1 / death_nbr
	#	survival_2[i] = temp2 / survivor_nbr
	data99 = [trace1,trace2]
	layout = dict(
	legend=dict(
	y=0.5,
 	traceorder='reversed',
	font=dict(size=16)))
	fig = dict(data=data99, layout=layout)
	#plot_3 = go.Figure(data=data99,layout=layout)
	#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
	#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
	plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
	#pio.write_image(fig, 'fig1.png')
	bla1 = sorted(survival_1.items())
	#figure_to_close = plt.gcf()
	#figure_to_close.clear()		
	plt.clf()
	if(len(survival_1)>0):
		x,y=zip(*bla1)
		plt.plot(x,y)
	bla2 = sorted(survival_2.items())
	z,a=zip(*bla2)
	plt.plot(z,a)
	#plt.show()
	fig = plt.gcf()
	#plot_mpl(fig, image='png', filename='/home/quirin/testproject/polls/static/testimage.png')
	#plt.show()
	fig.savefig("/home/quirin/testproject/polls/static/test_2.png")
	return(script,div,plot1,plot_div)
