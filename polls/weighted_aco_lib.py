from multiprocessing import Pool
import time
import pandas as pd
import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import scipy.sparse as sparse
from scipy.sparse import csr_matrix
from IPython.display import Audio, display
from sklearn.cluster.bicluster import SpectralCoclustering
from collections import Counter
import collections
from sklearn import preprocessing
flatten = lambda l: [item for sublist in l for item in sublist]
from scipy import stats
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
import seaborn as sns; sns.set(color_codes=True)
from multiprocessing import Pool
from numpy import linalg as LA
from sklearn.cluster import KMeans
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor    
from io import StringIO

def jac(x,y):
    if len(x)>0 and len(y)>0:
        return len(set(x).intersection(set(y)))/len((set(x).union(set(y))))
    else:
        return(0)
    
def jac_matrix(true,pred):
    res = np.zeros((len(true),len(true)))
    for i in range(len(true)):
        for j in range(len(true)):
            res[i,j] = jac(true[i],pred[j])
    cand1 = (res[0][0],res[1][1])
    cand2 = (res[0][1],res[1][0])
    if sum(cand1)>sum(cand2):
        return(cand1)
    else:
        return(cand2)
    
def matches(true1,true2,pred1,pred2):
    cand1 = (round(len(set(pred1).intersection(set(true2)))*100/len(true2)),
             round(len(set(pred2).intersection(set(true1)))*100/len(true1)))
    cand2 = (round(len(set(pred1).intersection(set(true1)))*100/len(true1)),
             round(len(set(pred2).intersection(set(true2)))*100/len(true2)))
    if sum(cand1)>sum(cand2):
        ans = cand1
    else:
        ans = cand2
    print(str(ans[0])+"%                  " + str(ans[1])+"%")
    


def joined_net(B,G):
    #joined net
    A_b = nx.adjacency_matrix(B).todense()
    A_b = A_b *1 ## trick to switch from boolean
    A_g = nx.adjacency_matrix(G).todense()
    n = len(A_g)
    A = A_b
    A[:n,:n] = A_g  
    return(A)





def hi(A_j,n,m):
    H = np.zeros((n+m,n+m))
    P = LA.matrix_power(A_j, 2)
    for node1 in range(n+m):
            for node2 in range(node1+1,n+m):
                if P[node1,node2] >0:
                    node_intersec = np.sum(np.multiply(A_j[node1,:],A_j[node2,:]))+A_j[node1,node2]
                    deg_node1 = np.sum(A_j[node1,:])
                    deg_node2 = np.sum(A_j[node2,:])
                    if node1<n and node2<n: #both beling to X
                        H[node1,node2] = node_intersec/(deg_node1+deg_node2+2)
                    if node1>=n and node2>=n: #both belong to Y
                        H[node1,node2] = 0.5*node_intersec/(deg_node1+deg_node2)
                    if (node2>=n) and (node1 <n): #node1 in X, node2 in Y or the other way around
                        H[node1,node2] = 4*node_intersec/(deg_node1+deg_node2)
                    H[node2,node1] = H[node1,node2]
    return(H*10)



def ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = True,show_plot = True, print_runs = True, save = None, show_nets = True):
    #fig = plt.figure(figsize=(10,8))
    #fig.suptitle('Your request is being processed...', fontsize=14, fontweight='bold')
    #plt.savefig("/home/quirin/testproject/polls/static/progress.png")
    #plt.close(fig)
    ge = GE.values
    H =H.astype(np.short)
    N = neigborhood(H,n,th)
    patients = np.arange(n,n+m)  
    
    cost = H/10
    cost = np.max(cost)-cost
    scores = []
    avs = []
    count_big = 0
    max_total_score = 0
    max_round_score = -100
    av_score = 0
    st = time.time()
    t0 = np.ones((n+m,n+m))*5
    t0 = t0.astype(np.short)
    probs= prob_upd(H,t0,a,b,n,th,N)
    end = time.time()
    flag = False
    score_change = []
    print ("Running time statistics:")
    print ("###############################################################")
    print("the joint graph has "+ str(n+m) + " nodes")
    print("probability update takes "+str(round(end-st,3)))
    W = 0
    while np.abs(max_round_score-av_score)>eps and count_big<times and (W<m/3):
        av_score = 0
        W = 0
        max_round_score = 0
        scores_per_round = []

        for i in range(K):
            #for each ant
            st = time.time()
            tot_score,gene_groups,patients_groups,new_scores,wars,no_int = ant_job(GE,N,H,th,clusters,probs,a,b,cost,m,n,patients,count_big,i,cost_limit,L_g_min,L_g_max,G,ge,print_runs)
            end = time.time()
            W = W+wars
            if count_big ==0 and i ==0:
                print("one ant run takes "+str(round(end-st,3)))
            scores_per_round.append(tot_score)
            av_score = av_score + tot_score
            if tot_score > max_round_score:
                max_round_score = tot_score
                solution = (gene_groups,patients_groups)
                full_scores = new_scores
                solution_big = (no_int,patients_groups)
            if count_big ==0 and i ==K-1:
                gs = 1.5*max_round_score

                t_max = (1/evaporation)*gs
                t_min = 0
   
                t0 = np.ones((n+m,n+m))*t_max
        #after all ants have finished:
        scores.append(scores_per_round)
        
        #saving rhe best overall solution
        if max_round_score>max_total_score:
            max_total_score = max_round_score
            best_solution = solution
            max_full_scores = full_scores 
            solution_big_best = solution_big
        score_change.append(round(max_round_score,3))
        print("Iteration # "+ str(count_big+1))
        print("best round score: " + str(round(max_round_score,3)))
        print("average score: " + str(round(av_score/K,3)))
        with open("/home/quirin/testproject/polls/static/output_console.txt", "w") as text_file:
        	#print("foobar")
        	text_file.write("Iteration # "+ str(count_big+1))
        	text_file.close()
        av_score = av_score/K
        avs.append(round(av_score,2))
        #print(scores)
        print(avs)
        #pher. and prob. updates
        t = pher_upd(t0,t_min,evaporation,max_full_scores,solution_big_best,flag)
        t0 = np.copy(t)
        
        probs= prob_upd(H,t,a,b,n,th,N)
        
        #visualization options:
        
        if show_pher:
            fig = plt.figure(figsize=(18,12))
            ax = fig.add_subplot(111)
            t_max = np.max(t)   
            cax = ax.matshow(t, interpolation='nearest',cmap=plt.cm.RdPu,vmin = t_min,vmax = t_max)
            plt.colorbar(cax)
            plt.title("Pheramones")
            plt.show(block=False)
            plt.close(fig)

        count_big = count_big +1
        if show_nets:
            features(solution, GE,G)    
        if show_plot:
            fig = plt.figure(figsize=(10,8))
            plt.boxplot(np.asarray(scores).T,manage_xticks = True, patch_artist=True)
            if opt!=None:
                plt.axhline(y=opt,label = "optimal solution score", c = "r")
            #plt.ylim((0,1))
            #plt.legend()
            #this was not commented before #plt.show(block=False)
            plt.savefig("/home/quirin/testproject/polls/static/progress.png")
            plt.close(fig)
        if len(set(score_change[:3])) ==1 and len(score_change)>3:
            flag = True
    if save != None:
        fig = plt.figure(figsize=(10,8))
        plt.boxplot(np.asarray(scores).T,manage_xticks = True)
        if opt!=None:
            plt.axhline(y=opt,label = "optimal solution score", c = "r")
        #plt.legend()
        plt.savefig(save+".png")
        plt.close(fig)
        
    #after the solutution is found we make sure to cluster patients the last time with that exact solution:
    data_new = ge[solution[0][0]+solution[0][1],:]
    kmeans = KMeans(n_clusters=2, random_state=0).fit(data_new.T)
    labels = kmeans.labels_
    patients_groups =[]
    for clust in range(clusters):
        wh = np.where(labels == clust)[0]
        group_p = [patients[i] for i in wh]
        patients_groups.append(group_p)
    if np.mean(ge[best_solution[0][0],:][:,(np.asarray(patients_groups[0])-n)])<np.mean(ge[best_solution[0][1],:][:,(np.asarray(patients_groups[0])-n)]):
        patients_groups = patients_groups[::-1]
    best_solution = [best_solution[0],patients_groups]
    
    print("best total score: "+str(max_total_score))
    #print_clusters(GE,best_solution)
    #features(best_solution, GE,G)
    return(best_solution,t,max_total_score,np.asarray(scores).T)
    
def neigborhood(H,n,th):
    N_per_patient = []
    dim = len(H)
    for i in range(n,dim):
        if th<0:
            N = np.where(H[i,:]>0.001)[0]
        else:
            rad = np.mean(H[i,:]) + th*np.std(H[i,:])
            N = np.where(H[i,:]>rad)[0]
        #N = np.where(H[i,:]>0)[0]
        N_per_patient.append(N)
    return N_per_patient
    
def prob_upd(H,t,a,b,n,th,N_per_patient):
    P_per_patient = []
    dim = len(H)
    temp_t = np.power(t,a)
    temp_H = np.power(H,b)
    temp = temp_t*temp_H 
    for i in range(n,dim):
        N_temp = N_per_patient[i-n]
        P = temp[:,N_temp]
        s = np.sum(P,axis = 1)
        s[s <1.e-4] = 1
        sum_p = 1/s
        sum_p = sum_p[:,None]
        P_new = P*sum_p[:np.newaxis]
        P_per_patient.append(P_new)

    return(P_per_patient)
        



def ant_job(GE,N,H,th,clusters,probs,a,b,cost,m,n,patients,count_big,count_small,cost_limit,L_g_min,L_g_max,G,ge,print_runs):
    if print_runs:
        print(str(count_big)+"."+str(count_small) + " run")
    paths = []
    wars = 0
    #set an ant on every patient
    
    for walk in range(m):
        k = cost_limit
        path = []
        start = patients[walk]
        Nn = N[walk] #neigbohood
        path.append(start)
        go = True
        P = probs[walk]            
        while go == True:
            P_new = P[start,:]
            #if there is any node inside the radious - keep mooving
            if np.sum(P_new)> 0.5:
                #transition:
                tr = np.random.choice(Nn,1,False,p = P_new)[0]
                c = cost[patients[walk],tr]
                #if there is any cost left we keep going
                if k-c >0:
                    path.append(tr)
                    start = tr
                    k = k - c
                #if not we are done and we save only genes from the path
                else:
                    go = False
            #no node to go - we are done and we save only genes from the path
            else:
                go = False
        path = np.asarray(path)
        #we are saving only genes
        path = path[path<n]
        paths.append(path)
        if len(path) == 0:
            wars = wars+1
            print("WARNING: emply path found")
    data_new = ge[list(set(flatten(paths))),:]
    kmeans = KMeans(n_clusters=2, random_state=0).fit(data_new.T)
    labels = kmeans.labels_
    gene_groups_set =[]
    patients_groups =[]
    for clust in range(clusters):
        wh = np.where(labels == clust)[0]
        group_g = [paths[i] for i in wh]
        group_g = flatten(group_g)
        gene_groups_set.append(set(group_g))
        #save only most common genes for a group
        group_p = [patients[i] for i in wh]
        patients_groups.append(group_p)
        
    #delete intersecting genes between groups
    
    I = set.intersection(*gene_groups_set)
    no_int =[list(gene_groups_set[i].difference(I)) for i in range(clusters)]
    gene_groups = no_int
    
    # make sure that gene clusters correspond to patients clusters:
    if np.mean(ge[gene_groups[0],:][:,(np.asarray(patients_groups[0])-n)])<np.mean(ge[gene_groups[1],:][:,(np.asarray(patients_groups[0])-n)]):
        patients_groups = patients_groups[::-1]
     
    gene_groups,sizes= clean_net(gene_groups,patients_groups, clusters,L_g_min,G,GE)
    

        
    new_scores = score(G,patients_groups,gene_groups,n,m,ge,sizes,L_g_min,L_g_max)
    tot_score = new_scores[0][0]*new_scores[0][1]+new_scores[1][0]*new_scores[1][1]   
    return(tot_score,gene_groups,patients_groups,new_scores,wars,no_int)
    
def pher_upd(t,t_min,p,scores,solution,flag):
    t = t*(1-p)
    t_new = np.copy(t)
    score = scores[0][0]*scores[0][1]+scores[1][0]*scores[1][1]
    for i in range(len(solution[0])):
        group_g = solution[0][i]
        group_p = solution[1][i]
        #score = scores[i][0]*scores[i][1]
        #ge_score = new_scores[i][0]*10
        #ppi_score = new_scores[i][1]*10
        for g1 in group_g:
            for p1 in group_p:
                t_new[g1,p1] = t[g1,p1]+ score
                t_new[p1,g1] = t[p1,g1]+ score
            for g2 in group_g:
                t_new[g1,g2] = t[g1,g2]+ score


    t_new[t_new < t_min] = t_min
            
    
    return(t_new)

def h_upd(H,scores,solution,p):
    H_new = np.copy(H)
    H_new = H*(1-p)
    score = (scores[0][0]*scores[0][1]+scores[1][0]*scores[1][1])*2
    for i in range(len(solution[0])):
        group_g = solution[0][i]
        group_p = solution[1][i]
        #score = scores[i][0]*scores[i][1]
        #ge_score = new_scores[i][0]*10
        #ppi_score = new_scores[i][1]*10
        for g1 in group_g:
            for p1 in group_p:
                H_new[g1,p1] = H[g1,p1]+ score
                H_new[p1,g1] = H[p1,g1]+ score
            for g2 in group_g:
                H_new[g1,g2] = H[g1,g2]+ score


    return(H_new)
    
    
    
    
def score(G,patients_groups,gene_groups,n,m,ge,sizes,L_g_min,L_g_max):
    clusters = len(patients_groups)
    conf_matrix = np.zeros((clusters,clusters))
    conect_ppi = []
    for i in range(clusters): #over genes
        group_g = np.asarray(gene_groups[i])
        s = sizes[i]
        if len(group_g)>0:
            for j in range(clusters): #over patients
                group_p = np.asarray(patients_groups[j])
                if len(group_p)>0:
                # gene epression inside the group
                    conf_matrix[i,j] = np.mean(ge[group_g,:][:,(group_p-n)])
            #ppi score    
            con_ppi = 1
            if s<L_g_min:
                con_ppi = s/L_g_min
            elif s>L_g_max:
                con_ppi = L_g_max/s
            conect_ppi.append(con_ppi)
        else:
            conect_ppi.append(0)           
    ans = []
    for i in range(clusters):
        all_ge = np.sum(conf_matrix[i,:])
        in_group = conf_matrix[i,i]
        out_group = all_ge - in_group
        ge_con = in_group-out_group
        #scaled = scaleBetween(num,0,0.5,0,1)
        ans.append((ge_con ,conect_ppi[i]))
        
    return(ans)
    

def aco_preprocessing(path_expr, path_ppi, col,log2, gene_list = None, size = None, sample= None):
    # path_expr - path for gene expression
    # path_ppi - path for ppi
    # col - split variable name (ONLY TWO CLASSES)
    # log2 - log2 transform
    #gene_list - preselected genes (if any)
    #size -  if genes are not preselected specify size of the gene set  for standard deviation selection
    # sample = None - all patients, otherwise specify fraction of patients taken
    expr = pd.read_csv(path_expr,sep = "\t") 
    expr = expr.set_index("Unnamed: 0")
    val1,val2 = list(set(expr[col]))
    group1_true = list(expr[expr[col]==val1].index)
    group2_true = list(expr[expr[col]==val2].index)
    patients_new = group1_true+group2_true
    if sample!=None:
        idx = list(expr.index)
        new_idx = np.random.choice(idx,int(sample*len(idx)),False)
        expr = expr.loc[new_idx]
        group1_true = list(expr[expr[col]==val1].index)
        group2_true = list(expr[expr[col]==val2].index)
        patients_new = group1_true+group2_true

    expr = expr.loc[patients_new]    
    net = pd.read_csv(path_ppi,sep = "\t", header= None)
    nodes_ppi = set(net[0]).union(set(net[1]))
    genes_ge = list(set(expr.columns) - set([col]))
    new_genes = [int(x) for x in genes_ge]
    intersec_genes = set.intersection(set(new_genes), set(nodes_ppi))
    genes_for_expr = [str(x) for x in list(intersec_genes)]
    expr = expr[genes_for_expr]
    #20188 genes
    if log2:
        expr = np.log2(expr)
    z_scores = stats.zscore(expr) 
    z_scores = pd.DataFrame(z_scores,columns = expr.columns, index = expr.index)
    if gene_list !=None and size == None:# gene list is given
        new_genes = [str(gene) for gene in gene_list] 
        
    elif gene_list == None and size!= None: #std selection
        std_genes = expr[genes_for_expr].std()
        std_genes, genes_for_expr = zip(*sorted(zip(std_genes, genes_for_expr)))
        genes_for_expr = genes_for_expr[len(std_genes)-size:]
        new_genes = list(genes_for_expr)
    elif gene_list == None and size == None: #all genes
        new_genes = genes_for_expr
    else:
        print("please specify gene selection method: predifined list, standart deviation filtering or none of them")
        return()

    expr = expr[new_genes]
    z_scores = z_scores[new_genes].values
    
    labels_B = dict()
    rev_labels_B = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in new_genes:
       labels_B[node] = g
       rev_labels_B[g] = node
       node = node+1
    for p in patients_new:
       labels_B[node] = p
       rev_labels_B[p] = node
       node = node+1
    

    #scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    #sim = scaler.fit_transform(expr)
    data_aco = pd.DataFrame(z_scores,columns= new_genes, index= patients_new)
    data_aco = data_aco.T
    n,m = data_aco.shape
    
    GE = pd.DataFrame(data_aco.values,index = np.arange(n), columns=np.arange(n,n+m))
    t = 2
    b = np.matrix(data_aco>t)
    b_sp = csr_matrix(b)
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(new_genes) and node2 in set(new_genes):    
            G.add_edge(rev_labels_B[node1],rev_labels_B[node2])
    A_new= nx.adj_matrix(G).todense()

    H = HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1)
    
    group1_true_ids= [rev_labels_B[x] for x in group1_true]
    group2_true_ids= [rev_labels_B[x] for x in group2_true]
    
    return B,G,H,n,m,GE,A_new,group1_true_ids,group2_true_ids,labels_B,rev_labels_B,val1,val2


def aco_preprocessing_ownfile(fh, fh_ppi, col,log2, gene_list = None, size = None, sample= None):
    # path_expr - path for gene expression
    # path_ppi - path for ppi
    # col - split variable name (ONLY TWO CLASSES)
    # log2 - log2 transform
    #gene_list - preselected genes (if any)
    #size -  if genes are not preselected specify size of the gene set  for standard deviation selection
    # sample = None - all patients, otherwise specify fraction of patients taken
    expr = pd.read_csv(fh,sep = "\t") 
    expr = expr.set_index("Unnamed: 0")
	#TODO: check if column 'prognosis' or 'cancer type' exists, set column based on this info
    if('cancer_type' in list(expr)):
    	col = 'cancer_type'
    else:
    	col = 'prognosis'
    val1,val2 = list(set(expr[col]))
    group1_true = list(expr[expr[col]==val1].index)
    group2_true = list(expr[expr[col]==val2].index)
    patients_new = group1_true+group2_true
    if sample!=None:
        idx = list(expr.index)
        new_idx = np.random.choice(idx,int(sample*len(idx)),False)
        expr = expr.loc[new_idx]
        group1_true = list(expr[expr[col]==val1].index)
        group2_true = list(expr[expr[col]==val2].index)
        patients_new = group1_true+group2_true

    expr = expr.loc[patients_new]    
    net = pd.read_csv(fh_ppi,sep = "\t", header= None)
    nodes_ppi = set(net[0]).union(set(net[1]))
    genes_ge = list(set(expr.columns) - set([col]))
    new_genes = [int(x) for x in genes_ge]
    intersec_genes = set.intersection(set(new_genes), set(nodes_ppi))
    genes_for_expr = [str(x) for x in list(intersec_genes)]
    expr = expr[genes_for_expr]
    #20188 genes
    if log2:
        expr = np.log2(expr)
    z_scores = stats.zscore(expr) 
    z_scores = pd.DataFrame(z_scores,columns = expr.columns, index = expr.index)
    if gene_list !=None and size == None:# gene list is given
        new_genes = [str(gene) for gene in gene_list] 
        
    elif gene_list == None and size!= None: #std selection
        std_genes = expr[genes_for_expr].std()
        std_genes, genes_for_expr = zip(*sorted(zip(std_genes, genes_for_expr)))
        genes_for_expr = genes_for_expr[len(std_genes)-size:]
        new_genes = list(genes_for_expr)
    elif gene_list == None and size == None: #all genes
        new_genes = genes_for_expr
    else:
        print("please specify gene selection method: predifined list, standart deviation filtering or none of them")
        return()

    expr = expr[new_genes]
    z_scores = z_scores[new_genes].values
    
    labels_B = dict()
    rev_labels_B = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in new_genes:
       labels_B[node] = g
       rev_labels_B[g] = node
       node = node+1
    for p in patients_new:
       labels_B[node] = p
       rev_labels_B[p] = node
       node = node+1
    

    #scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    #sim = scaler.fit_transform(expr)
    data_aco = pd.DataFrame(z_scores,columns= new_genes, index= patients_new)
    data_aco = data_aco.T
    n,m = data_aco.shape
    
    GE = pd.DataFrame(data_aco.values,index = np.arange(n), columns=np.arange(n,n+m))
    t = 2
    b = np.matrix(data_aco>t)
    b_sp = csr_matrix(b)
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(new_genes) and node2 in set(new_genes):    
            G.add_edge(rev_labels_B[node1],rev_labels_B[node2])
    A_new= nx.adj_matrix(G).todense()

    H = HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1)
    
    group1_true_ids= [rev_labels_B[x] for x in group1_true]
    group2_true_ids= [rev_labels_B[x] for x in group2_true]
    
    return B,G,H,n,m,GE,A_new,group1_true_ids,group2_true_ids,labels_B,rev_labels_B,val1,val2


def aco_preprocessing_strings(expr_str, ppi_str, col,log2, gene_list = None, size = None, sample= None):
    # path_expr - path for gene expression
    # path_ppi - path for ppi
    # col - split variable name (ONLY TWO CLASSES)
    # log2 - log2 transform
    #gene_list - preselected genes (if any)
    #size -  if genes are not preselected specify size of the gene set  for standard deviation selection
    # sample = None - all patients, otherwise specify fraction of patients taken
    EXPRDATA = StringIO(expr_str)
    expr = pd.read_csv(EXPRDATA,sep = "\t") 
    expr = expr.set_index("Unnamed: 0")
	#TODO: check if column 'prognosis' or 'cancer type' exists, set column based on this info
    if('cancer_type' in list(expr)):
    	col = 'cancer_type'
    else:
    	col = 'prognosis'
    val1,val2 = list(set(expr[col]))
    group1_true = list(expr[expr[col]==val1].index)
    group2_true = list(expr[expr[col]==val2].index)
    patients_new = group1_true+group2_true
    if sample!=None:
        idx = list(expr.index)
        new_idx = np.random.choice(idx,int(sample*len(idx)),False)
        expr = expr.loc[new_idx]
        group1_true = list(expr[expr[col]==val1].index)
        group2_true = list(expr[expr[col]==val2].index)
        patients_new = group1_true+group2_true

    expr = expr.loc[patients_new]
    PPIDATA = StringIO(ppi_str)    
    net = pd.read_csv(PPIDATA,sep = "\t", header= None)
    nodes_ppi = set(net[0]).union(set(net[1]))
    genes_ge = list(set(expr.columns) - set([col]))
    new_genes = [int(x) for x in genes_ge]
    intersec_genes = set.intersection(set(new_genes), set(nodes_ppi))
    genes_for_expr = [str(x) for x in list(intersec_genes)]
    expr = expr[genes_for_expr]
    #20188 genes
    if log2:
        expr = np.log2(expr)
    z_scores = stats.zscore(expr) 
    z_scores = pd.DataFrame(z_scores,columns = expr.columns, index = expr.index)
    if gene_list !=None and size == None:# gene list is given
        new_genes = [str(gene) for gene in gene_list] 
        
    elif gene_list == None and size!= None: #std selection
        std_genes = expr[genes_for_expr].std()
        std_genes, genes_for_expr = zip(*sorted(zip(std_genes, genes_for_expr)))
        genes_for_expr = genes_for_expr[len(std_genes)-size:]
        new_genes = list(genes_for_expr)
    elif gene_list == None and size == None: #all genes
        new_genes = genes_for_expr
    else:
        print("please specify gene selection method: predifined list, standart deviation filtering or none of them")
        return()

    expr = expr[new_genes]
    z_scores = z_scores[new_genes].values
    
    labels_B = dict()
    rev_labels_B = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in new_genes:
       labels_B[node] = g
       rev_labels_B[g] = node
       node = node+1
    for p in patients_new:
       labels_B[node] = p
       rev_labels_B[p] = node
       node = node+1
    

    #scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    #sim = scaler.fit_transform(expr)
    data_aco = pd.DataFrame(z_scores,columns= new_genes, index= patients_new)
    data_aco = data_aco.T
    n,m = data_aco.shape
    
    GE = pd.DataFrame(data_aco.values,index = np.arange(n), columns=np.arange(n,n+m))
    t = 2
    b = np.matrix(data_aco>t)
    b_sp = csr_matrix(b)
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(new_genes) and node2 in set(new_genes):    
            G.add_edge(rev_labels_B[node1],rev_labels_B[node2])
    A_new= nx.adj_matrix(G).todense()

    H = HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1)
    
    group1_true_ids= [rev_labels_B[x] for x in group1_true]
    group2_true_ids= [rev_labels_B[x] for x in group2_true]
    #print(group1_true + "babaaba")
    return B,G,H,n,m,GE,A_new,group1_true_ids,group2_true_ids,labels_B,rev_labels_B,val1,val2,group1_true,group2_true


def HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1):
    scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))#
    H_g_to_g = (data_aco.T.corr())*gtg_weight
    H_p_to_p = data_aco.corr()*ptp_weight
    H_g_to_g = scaler.fit_transform(H_g_to_g)
    H_p_to_p = scaler.fit_transform(H_p_to_p)
    H_g_to_p = scaler.fit_transform(data_aco)
    H_full_up = np.concatenate([H_g_to_g,H_g_to_p*gtp_weight], axis = 1)
    H_full_down = np.concatenate([H_g_to_p.T*gtp_weight,H_p_to_p], axis = 1)
    H_full =  np.concatenate([H_full_up,H_full_down], axis = 0)*10
#    H_full[H_full < 1] = 1
#    np.fill_diagonal(H_full, 1)
    np.fill_diagonal(H_full, 0)
    return(H_full)
    

def most_common(lst,top,L_g):
    data = Counter(lst)
    l = len(data)
    take = int(top*l)
    if take == 0:
        take = 1
    if take <L_g:
        take = L_g
    count = data.most_common(take)
    genes = [x[0] for x in count]
    return genes



def scaleBetween(unscaledNum, minAllowed, maxAllowed, min_cur, max_cur):
  return (maxAllowed - minAllowed) * (unscaledNum - min_cur) / (max_cur
         - min_cur) + minAllowed
def print_clusters(GE,solution):
    grouping_p = []
    p_num = list(GE.columns)
    for p in p_num:
        if p in solution[1][0]:
            grouping_p.append(1)
        else:
            grouping_p.append(2)
    grouping_p = pd.DataFrame(grouping_p,index = p_num)
    grouping_g = []
    g_num = list(GE.index)
    for g in g_num:
        if g in solution[0][0]:
            grouping_g.append(1)
        elif  g in solution[0][1]:
            grouping_g.append(2)
        else:
            grouping_g.append(3)
            
    grouping_g = pd.DataFrame(grouping_g,index = g_num)
    species = grouping_p[0]
    lut = {1: '#A52A2A', 2: '#7FFFD4'}
    row_colors = species.map(lut)
    species = grouping_g[0]
    lut = {1: '#A52A2A', 2: '#7FFFD4', 3:'#FAEBD7'}
    col_colors = species.map(lut)
    sns.clustermap(GE.T, row_colors=row_colors, col_colors = col_colors,figsize=(15, 10))
    
def features(solution, GE,G,pos = None):
    genes1,genes2 = solution[0]
    patients1, patients2 = solution[1]
    
    means1 = list(np.mean(GE[patients1].loc[genes1],axis = 1)-np.mean(GE[patients2].loc[genes1],axis = 1).values)
    means2 = list(np.mean(GE[patients1].loc[genes2],axis = 1)-np.mean(GE[patients2].loc[genes2],axis = 1).values)
    G_small = nx.subgraph(G,genes1+genes2)
    
    fig = plt.figure(figsize=(15,10))
    vmin = min(means1+means2)
    vmax = max(means1+means2)
    if pos == None:
        pos = nx.spring_layout(G_small)
    ec = nx.draw_networkx_edges(G_small,pos)
    nc1 = nx.draw_networkx_nodes(G_small,nodelist =genes1, pos = pos,node_color=means1, node_size=200,alpha=1.0,
                                 vmin=vmin, vmax=vmax,node_shape = "^",cmap =plt.cm.viridis)
    nc2 = nx.draw_networkx_nodes(G_small,nodelist =genes2, pos = pos,node_color=means2, node_size=200,
                                 alpha=1.0,
                                 vmin=vmin, vmax=vmax,node_shape = "o",cmap =plt.cm.viridis)
    nx.draw_networkx_labels(G_small,pos)
    plt.colorbar(nc1)
    plt.axis('off')
    
    plt.show(block=False)
    plt.close(fig)
    
def stability_plot(data, labels, name = None):
    jaccards = []
    for categ in data:       
        jk = []
        for i in range(len(categ)):
            for j in range(i+1,len(categ)):
                jk.append(jac(categ[i],categ[j]))
        jaccards.append(jk)
    fig, ax = plt.subplots(figsize=(15, 10))
    bplot1 = ax.boxplot(jaccards,
                             vert=True,  # vertical box alignment
                             patch_artist=True,  # fill with color
                             labels=labels)  # will be used to label x-ticks
    plt.ylabel("jaccard index")
    plt.xlabel("required gene module")
    plt.ylim(0,1)
    if name!=None:
        plt.savefig(name+".png")
    plt.show(block=False)
    plt.close(fig)


    
    
def clean_net(gene_groups,patients_groups, clusters,L_g,G,GE):    
    genes_components = []
    sizes = []
    for clust in range(clusters):
        group_g = gene_groups[clust]
        if clust == 0:
            not_clust = 1
        else:
            not_clust = 0
        if len(group_g)>=L_g:
            g = nx.subgraph(G,group_g)
            comp_big = max(nx.connected_component_subgraphs(g), key=len)
            dg = dict(nx.degree(comp_big))
            ones = [x for x in dg if dg[x]==1]
            nodes = list(comp_big.nodes)
            size_comp = len(nodes)
            max_out = len(nodes)- L_g
            while max_out >0:
                dif = np.mean(GE[patients_groups[clust]].loc[ones],axis = 1)-np.mean(GE[patients_groups[not_clust]].loc[ones],axis = 1)
                dif = dif.sort_values()
                ones = list(dif[dif<2].index)
                if len(ones)>0:
                    if len(ones)<=max_out:
                        outsiders = list(ones)
                    if len(ones) > max_out:
                        outsiders = list(ones)[:max_out]
     
                    nodes  = list(set(nodes) - set(outsiders))
                    g = nx.subgraph(G,nodes)
                    comp_big = g
                    dg = dict(nx.degree(comp_big))
                    ones = [x for x in dg if dg[x]==1]
                    nodes = list(comp_big.nodes)
                    size_comp = len(nodes)
                    max_out = len(nodes)- L_g
                else:
                    max_out = 0
                    
            group_g = nodes
        elif len(group_g)>0:
            g = nx.subgraph(G,group_g)
            comp_big = max(nx.connected_component_subgraphs(g), key=len)
            nodes = list(comp_big.nodes)
            size_comp = len(nodes)
        else:
            size_comp = 0
            
        genes_components.append(group_g)
        sizes.append(size_comp)
    return genes_components,sizes

def sim_data(genes1,genes2,background,patients1,patients2,dens):
    n = genes1+genes2+background
    m = patients1 +patients2
    
    
    genes = np.arange(n)
    groups_genes = list(np.ones(genes1))+list(np.ones(genes2)*2)+list(np.ones(background)*3)
    groups_p = [1 if node<patients1 else 2 for node in range(m)]
    
    to_sparce = 0.3 #to sparcify bipartite
    to_mix = 0.99 # to mix edges berween groups
    b = np.zeros((n,m))
    ge = np.random.normal(0,1,n*m).reshape(n,m)

    for patient in range(m):
        for gene in range(n):
            p_gr = groups_p[patient]
            g_gr = groups_genes[gene]
            if p_gr ==1 and g_gr == 1: #all up
                ge[gene,patient] = np.random.normal(1,0.35,1)
            elif p_gr ==2 and g_gr == 2:
                ge[gene,patient] = np.random.normal(1,0.35,1) #also up
            elif p_gr ==1 and g_gr == 2:
                ge[gene,patient] = np.random.normal(-1,0.35,1) #down
            elif p_gr ==2 and g_gr == 1:
                ge[gene,patient] = np.random.normal(-1,0.35,1) #down
    for patient in range(m):
        for gene in range(genes1+genes2):
            prob = np.random.uniform(0,1)
            if prob>0.9:
                ge[gene,patient] = np.random.normal(0,1,1)
                
    for gene in range(genes1+genes2,n):
        prob = np.random.uniform(0,1)
        if prob<0.05:


            for patient in range(m):
                if groups_p[patient] ==1: #all up
                    ge[gene,patient] = np.random.normal(0.3,0.35,1)
                else:
                    ge[gene,patient] = np.random.normal(-0.3,0.35,1)
        if prob>0.05 and prob<0.1:
            for patient in range(m):
                if groups_p[patient] ==1: #all up
                    ge[gene,patient] = np.random.normal(-0.3,0.35,1)
                else:
                    ge[gene,patient] = np.random.normal(0.3,0.35,1)
                    
                
    g1 = nx.barabasi_albert_graph(genes1,1)
    g2 = nx.barabasi_albert_graph(genes2,1)
    g3 = nx.barabasi_albert_graph(background, 1)
    G = nx.disjoint_union(g1,g2)
    G = nx.disjoint_union(G,g3)
    for _ in range( int(dens*n)):
        node1 = np.random.randint(0,genes1)
        node2 = np.random.randint(genes1,genes1+genes2)
        node3_1 = np.random.randint(genes1+genes2,n)
        node3_2 = np.random.randint(genes1+genes2,n)
        G.add_edges_from([(node1,node3_1),
                          (node2,node3_2)])
     
    d =  nx.density(G)
    count = 0 
    while d>0.002 and count<10:
        
        node3_1 = np.random.randint(genes1+genes2,n)
        node3_2 = np.random.randint(genes1+genes2,n)
        count = count+1
        if G.has_edge(node3_1,node3_2):
            G.remove_edge(node3_1,node3_2)
            d =  nx.density(G) 
        
    #A_g = nx.adj_matrix(G).todense() *1
    b_sp = csr_matrix(b) #sparse matrix for making bipartite graph
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    GE = pd.DataFrame(ge,index = np.arange(n),columns = np.arange(n,n+m))
    H = HI_big(GE,1,1,1)
    return(B,GE,G,H,d,n,m)
