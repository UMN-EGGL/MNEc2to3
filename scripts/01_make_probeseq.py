#!/usr/bin/python3

import argparse
from locuspocus import Fasta
from minus80.Tools import available
import sys
import pysam

def print_flank_sequence(vcf,fasta,window=35):                                     
    vcf = pysam.VariantFile(vcf)                                                   
    if not available("Fasta",'temp'): 
        fasta = Fasta.from_file('temp',fasta)                                                 
    else:
        fasta = Fasta('temp')

    for i,var in enumerate(vcf):                                                   
        flank5 = ''.join(fasta[var.chrom][var.pos-window:var.pos-1]).upper()                
        flank3 = ''.join(fasta[var.chrom][var.pos+1:var.pos+window]).upper()                
        kmer = flank5 + var.ref + flank3                                           
        print(">{}.{}\n{}".format(var.chrom, var.pos, kmer))                       
                 

if __name__ == '__main__':
    # Drop into ipdb if there is an error	
    #from IPython.core import ultratb
    #sys.excepthook = ultratb.FormattedTB(
    #    mode='Verbose', color_scheme='Linux', call_pdb=1
    #)      
 
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--vcf'
    )
    parser.add_argument(
        '--fasta'
    )
    args = parser.parse_args()
    print_flank_sequence(args.vcf,args.fasta)
