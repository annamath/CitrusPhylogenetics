#!/usr/bin/env python

# <><><> Required Packages <><><>
from Bio import AlignIO
import argparse

# <><><> Create Parser Argument <><><>
parser = argparse.ArgumentParser()

parser.add_argument("-fa", required=True, help="Fasta file to convert to phylip format")
parser.add_argument("-o", required=True, help="Output file in phylip format")

args = parser.parse_args()

def file_validation(parser, arg):
     """Validate if a file exists"""
     if not os.path.exists(str(arg)):
        parser.error("{}{}{}".format("File Error.", arg, " does not exist."))
     else:
        pass

if __name__ == "__main__":
     file_validation(parser, args.fa)
     AlignIO.convert(args.fa, "fasta", args.o, "phylip-relaxed")
