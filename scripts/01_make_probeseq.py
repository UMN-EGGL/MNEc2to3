#!/usr/bin/python3

import argparse
from locuspocus import Fasta
from minus80.Tools import available
import sys
import pysam

def print_flank_sequence(vcf,fasta,out,window=30):                                     
    vcf = pysam.VariantFile(vcf)                                                   
    if not available("Fasta",'temp'): 
        fasta = Fasta.from_file('temp',fasta)                                                 
    else:
        fasta = Fasta('temp')
    

    with open(out,'w') as OUT:
        for i,var in enumerate(vcf):                                                   
            # Grab the FASTA flank sequence
            if fasta[var.chrom][var.pos] != var.ref:
                print(f"{var.id} REF does not match the FASTA. oops.")
            # Grab the probe sequence anyways
            kmer = ''.join(fasta[var.chrom][var.pos+1:var.pos+window]).upper()
            print(">{}.{}_{}\n{}".format(var.chrom, var.pos, "1", kmer),file=OUT)
            kmer = ''.join(fasta[var.chrom][var.pos-window:var.pos-1]).upper()
            print(">{}.{}_{}\n{}".format(var.chrom, var.pos, "2", kmer),file=OUT)                       
                 

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
    parser.add_argument(
        '--out',
        default='probeseq.fa'
    )
    args = parser.parse_args()
    print_flank_sequence(args.vcf,args.fasta,args.out)
