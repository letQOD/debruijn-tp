#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 19:11:22 2019

@author: lfabre
"""

###Requirement###

import os
import argparse
import networkx as nx
from networkx import algorithms
import statistics

############################################################
############# Identification des -mers uniques #############
############################################################

### parsing du fastq ###
def read_fastq(fqfile):
    """returns a sequence iterator, skips comment and score lines"""    
    fq=open(fqfile, 'r')
    for line in fq:
        yield next(fq).strip()
        next(fq)
        next(fq)
        
### identification des k-mers uniques ###
def cut_kmer(seq, k):   
    """returns a k-mer iterator"""                
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

### creation du dictionnaire de k-mers ###
def build_kmer_dict(fqfile, k): 
    """returns a kmer dictionary as key=k-mer, value= n k-mer"""
    kmer_dic = {}
    for i in read_fastq(fqfile):
        for kmer in cut_kmer(i, k):
            if not kmer in kmer_dic:
                kmer_dic[kmer] = 1
            else:
                kmer_dic[kmer] += 1
    return kmer_dic


############################################################
############ Construction de l arbre de Bruijn #############
############################################################

def build_graph(kmer_dic):
    """returns a graph of k-mers"""
    G = nx.DiGraph()    
    for kmer in kmer_dic:
        G.add_edge(kmer[0:len(kmer)-1], kmer[1:len(kmer)], weight=kmer_dic[kmer])
    return G


############################################################
############ Parcours du graphe de Bruijn ##################
############################################################

def get_starting_nodes(G):
    """returns input node of a graph"""
    start_node_list = []
    for node in G:
        start_node = list(G.predecessors(node))
        if not start_node:
            start_node_list.append(node)
    return start_node_list


def get_sink_nodes(G):
    """returns output nodes of a graph"""
    end_node_list = []
    for node in G:
        end_node = list(G. successors(node))
        if not end_node:
            end_node_list.append(node)
    return end_node_list

def get_contigs(G, start_node_list, end_node_list):
    """Get all paths from input to output nodes : returns a list of contigs and their size"""
    contigs = []
    for source in start_node_list:
        for target in end_node_list:
            if algorithms.has_path(G ,source, target)==True:
                path = algorithms.shortest_path(G, source, target)
                contig=path[0]
                for i in range(len(path)-1):
                    contig +=path[i+1][-1]
                contigs.append((contig, len(contig)))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(tupple, fqfile):
    """returns a fasta formatted file"""
    file = open(fqfile, 'w+')
    for i in range(len(tupple)):
        file.write('>contig_' + str(i) + ' len=' + str(tupple[i][1]) + '\n' 
                   + str(fill(tupple[i][0])) + '\n')
    file.close()


############################################################
########## Simplification du graphe de Bruijn ##############
############################################################
 
### resolution des bulles ###

def std(val_list):
    """returns the standard deviation"""
    return statistics.pstdev(val_list) 


def path_average_weight(G, path):
    """returns the mean weight of the path"""
    weight_list = []
    for i,j,dic_weight in G.subgraph(path).edges(data=True):
        weight_list.append(dic_weight["weight"])
    return statistics.mean(weight_list)
    

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    #return(clean_graph)
    pass

def select_best_path(graph, path_list,pathlen_list, meanweight_list, delete_entry_node=False, delete_sink_node=False):
    #return(clean_graph)
    pass

def solve_bubble(graph, ancestor, descendant):
    #return(clean_graph)
    pass

def simplify_bubbles(graph):
    #return(clean_graph)
    pass


### detection des tips ###

def solve_entry_tips(graph, entry_node_list):
    #return(cleangraph)
    pass

def solve_out_tips(graph, end_node_list):
    #return(cleangraph)
    pass




############################################################
####################### Main ###############################
############################################################

    
def main():
    """main"""

#arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", metavar="--fastq", help="inputfile fastq single end",type=str, required=True)
    parser.add_argument("-k", metavar="--kmer_size", help="k-mer size", type=int, default=21)
    parser.add_argument("-o", metavar="--output_file", help="output file contigs", type=str, required=True)
    args = parser.parse_args()

#Lancement des fonctions    
    
    kmer_dic = build_kmer_dict(args.i, int(args.k))
    G = build_graph(kmer_dic)
    starting_nodes = get_starting_nodes(G)
    sink_nodes = get_sink_nodes(G)
    contigs_tuple = get_contigs(G, starting_nodes, sink_nodes)
    save_contigs(contigs_tuple, args.o)
    
    
if __name__ == "__main__":
    main()





