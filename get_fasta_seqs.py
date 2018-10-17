#!/usr/bin/env python

############################
# Author: Anna Mathioudaki #
############################

# <><><> Required Packages <><><>
import pandas as pd
import argparse
import re
import os

# <><><> Create Parser Argument <><><>
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", required = True, help="Input file of protein domains")
parser.add_argument("-fa", "--fasta", required = True, help="Fasta input file")# type=argparse.FileType('r') if i contain this argument and then i try to open the file-> error because alr$
parser.add_argument("-ev", "--evalue", help="e-value required cutoff", required = False ) #optional
parser.add_argument("-o", "--output", help="Output the result to a file")

args = parser.parse_args()

final_df = pd.DataFrame()

def file_validation(parser, arg):
    """Validate if a file exists"""
    if not os.path.exists(str(arg)):
        parser.error("{}{}{}".format("File Error.", arg, " does not exist."))
    else:
        pass

check_list = list() #list of domains extracted from pipeline

def get_coordinates(inputfile):
    file_validation(parser, inputfile)
    with open(str(inputfile), "r") as f:
        if args.evalue:
            print '{}{}'.format("E-value threshold is equal to:", args.evalue)
            flag = 1
        else:
            print "No evalue threshold was set"
            flag = 0 
        i = 0
        for line in f:
            columns = line.split("\t")
            all_domains = (columns[1]).split('~')
            for domain in all_domains:
                if 'NB-ARC' in domain:
                    nbarc = domain.replace("(", ",").replace("=", ",").replace(")", ",")
                    sections = nbarc.split(",")
                    start = sections[2]
                    stop = sections[4]
                    evalue = sections[6]
                    if flag == 1:
                        if evalue >= args.evalue:
                            final_df.loc[i, 0] = columns[0]
                            final_df.loc[i, 1] = start
                            final_df.loc[i, 2] = stop
                            final_df.loc[i, 3] = evalue
                    else:
                        check_list.append(columns[0])
                        final_df.loc[i, 0] = columns[0]
                        final_df.loc[i, 1] = start
                        final_df.loc[i, 2] = stop
                        final_df.loc[i, 3] = evalue
                    i = i + 1
        final_df.columns = ["proteinID", "Start", "End", "Evalue"]

def fasta_iteration(fasta_file, output_file):
    """Yield headers and sequences from a fasta file, to a new fasta file
    The new file will contain IDs of proteins included in a given matrix"""

    file_validation(parser, fasta_file)
    seq_dict = {}

    with open(fasta_file, "r") as fh:
        for line in fh:         
            if line.startswith('>') and line[0:line.find(' ')].replace('>', '') in check_list:
                seq = ""
                id =  line[0:line.find(' ')].replace('>', '')
                flag = 1
            elif line.startswith('>'):
                flag = 0
            else:
                if flag == 1:
                     seq = seq + str(line)
                     seq = seq.replace('\n', '')
                     seq_dict[id] = seq

    with open(output_file, "w") as wf:
        for id in seq_dict.keys():
            if id in check_list:
                tmp_df = final_df.loc[final_df["proteinID"] == id]
                for j in range(0, tmp_df.shape[0]):
                    begin = int(tmp_df.iloc[j, 1])
                    end = int(tmp_df.iloc[j, 2])
                    domain = seq_dict[id][begin:end]
                    wf.write('>' + id +'_'+str(j)+ '\n' + domain + '\n')
            else:
                pass

if __name__ == '__main__':
     get_coordinates(args.input_file)
     check_list = set(check_list)
     fasta_iteration(args.fasta, args.output)
