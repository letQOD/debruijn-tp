#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 19:11:22 2019

@author: lfabre
"""

###Requirement###
import sys
import sys
import os
import argparse
import networkx as nx


############################################################
############# Identification des -mers uniques #############
############################################################

### parsing du fastq ###
def read_fastq(fqfile):         ###fqfile=arg.i
	fq = open (fqfile,'r')
	for line in fq:
		yield next(fq).strip()
		next(fq)
		next(fq)
        
### identification des k-mers uniques ###
def cut_kmer(seq, k):                       ###k=arg.k , seq= read
    for i in range(len(seq)-kmer_size+1):
        yield seq[i:i+kmer_size]

### creation du dictionnaire de k-mers ###
def build_kmer_dict (fqfile, k):            ###fqfile=arg.i, k=arg.k
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











##### MAIN #####
def main():

#arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="inputfile fastq single end")
	parser.add_argument("-k", help="k-mer size", default=21)
	parser.add_argument("-o", help= "output file contigs")
	args = parser.parse_args()
