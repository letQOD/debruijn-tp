#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 19:11:22 2019

@author: lfabre
"""

###Requirement###
import sys
import os
import argparse
import networkx as nx


############################################################
############# Identification des -mers uniques #############
############################################################

### parsing du fastq ###
def read_fastq(fqfile):         ###fqfile=arg.i
    fq = open(fqfile, 'r')
    for line in fq:
        yield next(fq).strip()
        next(fq)
        next(fq)
        
### identification des k-mers uniques ###
def cut_kmer(seq, k):                       ###k=arg.k , seq= read
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

### creation du dictionnaire de k-mers ###
def build_kmer_dict(fqfile, k):            ###fqfile=arg.i, k=arg.k
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
    G = nx.DiGraph()    
    for kmer, val in kmer_dic.items():
        G.add_edge(kmer[0:len(kmer)-1], kmer[1:len(kmer)], weight=val)  ####weight renseigne quand ???
    return G












############################################################
############ Main #############
############################################################

def main():

#arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="inputfile fastq single end")
    parser.add_argument("-k", help="k-mer size", default=21)
    parser.add_argument("-o", help="output file contigs")
    args = parser.parse_args()

#Lancement des fonctions    
    
    kmer_dic = build_kmer_dict(args.i, int(args.k))
    G = build_graph(kmer_dic)
    #starting_nodes = get_starting_nodes(G)
    #sink_nodes = get_sink_nodes(G)
    #contigs_tuple = get_contigs(G, starting_nodes, sink_nodes)
    #save_contigs(contigs_tuple, 'exit.txt')
    
    
if __name__ == "__main__":
    main()





