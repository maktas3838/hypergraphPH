# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:20:47 2022

@author: maktas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:08:33 2021

@author: mehmet
"""

#import dionysus as d
import numpy as np
import networkx as nx
import ast
from itertools import combinations 
import time

file_name = "text_bbc" # change file name here

file = open("Datasets/" + file_name + "/" + file_name + "_nverts.txt", "r")
file_simplices = open("Datasets/" + file_name + "/" + file_name + "_simplices.txt", "r")
file_songs = open("Datasets/" + file_name + "/number_of_songs.txt", "r")

songs = []
for x in file_songs:
    songs.append(int(x))
print(len(songs))
file_songs.close()

max_ss=25 # max simplex size
comb={}

for ss in range(2, max_ss):
    cc=[]
    for si in range(1,ss):
        comb_list=list(combinations(list(range(ss)), si))
        cc.append(comb_list)
    comb[ss]=cc

first = 1
ind = 0

for n in range(len(songs)):
    if (songs[n] == 0): # skip songs that have no simplices
        continue
    print(n)
    
    start_time = time.time()
    ind=ind+1
    simplices = []
    time_counter = 0
    line_counter = 1
    simp_notime=[]
    G = nx.Graph()

    for x in file:
        nline = int(x)
        if nline == 0: # skip non-simplex
            continue
        count = 0
        array = []
        for y in file_simplices:  # loop will run n lines in simplices file corresponding to the number in nverts
            if int(y) not in array:
                array.append(int(y))
            count = count + 1
            if count == nline:
                time_counter = time_counter + len(array)
                array.sort()
                for v in array:
                    if [v] not in [item[0] for item in simplices]:
                        simp_notime.append([v])
                        simplices.append(([v], time_counter))
                if array not in [item[0] for item in simplices]:
                    simp_notime.append(array)
                    simplices.append((array, time_counter))
                break
        if line_counter == songs[n]:
            break
        line_counter += 1

    # # skip (exclude) hypergraphs have simplices' size greater than 25
    # list_len = [len(i[0]) for i in simplices]
    # if max(list_len) > 25:
    #     continue

    simplices.sort(key = lambda x: x[1]) # sort by time
    simplexClosure=simplices.copy()
    missing=[]
    
    sl=[len(x) for x in simp_notime]
    max_value = None
    simp_notime2=simp_notime.copy()
    for num in sl:
        if (max_value is None or num > max_value):
            max_value = num
            
    print(ind, max_value)
    if max_value<=25:
        for s in simplices:
            simp=s[0]
            ls=len(simp)
            if ls>1:
                aa=comb[ls]
                for i in range(len(aa)):
                    bb=aa[i]
                    for j in range(len(bb)):
                        cl=list(bb[j])
                        a=[]
                        for k in cl:
                            a.append(simp[k])
                        out = [t for t in a]
                        if out not in simp_notime:
                            simplexClosure.append((out,s[1]))    # we add the missing hyperedges at the minimum time value observed in other hyperedges
                            simp_notime.append(out)
                            missing.append(out)
    else:
        print('The song with max size is ...', max_value)
        print('This happens on ...', ind)
    
    arr = []
    dic = {}
    value = 1 
    for s in simplices:
        each = []
        if len(s[0]) >= 2:
            string_ints = [str(inl) for inl in s[0]]
            stri = "".join(string_ints)
            each = [int(stri)]
        else:
            each = [s[0][0]] #to get 1 instead of [1]
        dic[each[0]] = value # make simplices as key in dictionary 
        G.add_node(value, time = s[1])
        arr.append(([value], s[1]))
        value +=1
    
    for s in simplices:
        each = []
        if len(s[0]) >= 2:
            string_ints = [str(inl) for inl in s[0]]
            stri = "".join(string_ints)
            each = [ dic[int(stri)] ] # get value of that simplex(key) in dictionary 
        else:
            each = [ dic[s[0][0]] ] 
            
        for s_check in simplices:
            org = each.copy()
            if len(s[0]) >= len(s_check[0]):
                continue
            else:
                if set(s[0]).issubset(set(s_check[0])):
                    string_int = [str(inl) for inl in s_check[0]]
                    strl = "".join(string_int)
                    org.append( dic[int(strl)] )
                    if s[1] > s_check[1]:
    
                        arr.append((org,s[1]))  
                        G.add_edge(org[0],org[1], time = s[1])
                    else:
                        arr.append((org,s_check[1]))  
                        G.add_edge(org[0],org[1], time = s_check[1])
                        
    ns=len(simplices) # number of simplices in the barycentric subdivision
    nsc=len(simplexClosure)
    edge_node=[] # find the nodes in the original graph where we need to add an edge to the qutient node
    for i in range(ns,nsc):
        scs=simplexClosure[i][0]
        if len(scs)>1:
            aa=comb[len(scs)]
            for i in range(len(aa)):
                bb=aa[i]
                for j in range(len(bb)):
                    cl=list(bb[j])
                    a=[]
                    for k in cl:
                        a.append(scs[k])
                    out = [t for t in a]
                    if out in simp_notime2:
                        ind_node=simp_notime2.index(out)
                        if ind_node not in edge_node:
                            edge_node.append(ind_node)
        else:
            for simm in scs:
                if simm in simp_notime2:
                    ind_node=simp_notime2.index(simm)
                    if ind_node not in edge_node:
                        edge_node.append(ind_node)
                
        for simp in simp_notime2:
            if len(simp)>len(scs):
                if set(scs)<set(simp):
                    ind_node=simp_notime2.index(simp)
                    if ind_node not in edge_node:
                        edge_node.append(ind_node)
                        continue
                
        
    #G = M.to_undirected()
    cliques = [c for c in nx.enumerate_all_cliques(G)]
    
    # find clique time
    arr2=[]
    for c in cliques:
        if len(c)==1:
            for k in range(len(arr)):
                if c==arr[k][0]:
                    arr2.append((c,arr[k][1]))
                    break
        else:
            com = [(u,v) for u,v in combinations(c,2)] # find each pair of clique
            max_time = 0
            for e in com:
                if G.edges[e]['time'] > max_time:
                    max_time = G.edges[e]['time']
            arr2.append((c,max_time))
    '''
    f = d.Filtration()
    
    for vertices, time1 in arr2:
        f.append(d.Simplex(vertices, time1))
    
    f.sort()
    m = d.homology_persistence(f)
    
    # init_diagrams() method:
    dgms = d.init_diagrams(m,f)
    print(dgms)
    

    if first == 1:
        out = open("Datasets/" + file_name + "/relative_" + file_name + "_output.txt","w")
        first += 1
    else:
        out = open("Datasets/" + file_name + "/relative_" + file_name + "_output.txt","a")
    for i, dgm in enumerate(dgms):
        for pt in dgm:
            out.write(str(i) + " " + str(pt.birth) + " " + str(pt.death) + "\n")
    out.write("\n")
    out.close()
    
    '''
    #d.plot.plot_bars(dgms[0], show = True)
    #d.plot.plot_bars(dgms[1], show = True)

    # d.plot.plot_bars(dgms[0], show = True)
    # # d.plot.plot_bars(dgms[2], show = True)

    # # wasserstein_distance() computes q-th Wasserstein distance between a pair of persistence diagrams. 
    # # bottleneck_distance() computes the bottleneck distance.

    # f1 = d.fill_rips(np.random.random((20, 2)), 2, 1)
    # m1 = d.homology_persistence(f1)
    # dgms1 = d.init_diagrams(m1, f1)
    # f2 = d.fill_rips(np.random.random((20, 2)), 2, 1)
    # m2 = d.homology_persistence(f2)
    # dgms2 = d.init_diagrams(m2, f2)
    # wdist = d.wasserstein_distance(dgms1[1], dgms2[1], q=2)
    # print("2-Wasserstein distance between 1-dimensional persistence diagrams:", wdist)
    # bdist = d.bottleneck_distance(dgms1[1], dgms2[1])
    # print("Bottleneck distance between 1-dimensional persistence diagrams:", bdist)