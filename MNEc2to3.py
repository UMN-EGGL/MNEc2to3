#!/usr/bin/python3

import argparse,datetime,os,glob
import pandas as pd
import numpy as np

def fixGenos(x):
    # update genotypes in case of REF/ALT allele flips
    if x == '0/0':
        return '1/1'
    elif x == '1/1':
        return '0/0'
    else:
        return x

def sortChrom(chrom,outfile):
    c = pd.read_table('.'+chrom+'.vcf')
    c.sort_values(by=['POS'],inplace=True)
    c.to_csv(outfile,index=False,header=False,sep='\t')

def writeVCF(remap,infile,outfile,chrUn=False):
    # basing re-mapping on EquCab2 chromosome & positions instead of SNP IDs
    remap['id'] = remap.EC2_chrom.map(str)+'.'+remap.EC2_pos.map(str)
    remap_index = pd.Series(remap.index,index=remap.id).to_dict()
    # write new VCF header
    out = open(outfile,'w')
    now = datetime.datetime.now()
    header = ['##fileformat=VCFv4.2\n',
            '##filedate='+str(datetime.datetime.now()).split(' ')[0]+'\n',
            '##source=MNEc2to3.py\n',
            '##INFO=<ID=Name,Number=A,Type=String,Description="SNP Name">\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',]
    if chrUn:
        contigs = [line for line in open('EquCab3_contigs.txt','r')]
    else:
        contigs = [line for line in open('EquCab3_contigs.txt','r') if not 'chrUn' in line]
    header = header+contigs
    for line in header:
        out.write(line)
    # create temporary files for each chromosome to pre-sort
    chroms = [contig.split(',')[0].split('=')[2] for contig in contigs]
    tmp_vcfs = {}
    for chrom in chroms:
        tmp_vcfs[chrom] = open('.'+chrom+'.vcf','a')
    # process input VCF by line
    vcf = open(infile,'r')
    pairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    print('Reading VCF...')
    ctr = 0
    for line in vcf:
        ctr += 1
        if ctr % 100000 == 0:
            print('\tProcessed',ctr,'lines...')
        # skip existing header. Contigs are going to be wrong now.
        if line.startswith('##'):
            if line.startswith('##FORMAT'):
                out.write(line)
            continue
        elif line.startswith('#CHROM'):
            header = line.strip().split()
            # make sure we know where samples start
            if header[8] != 'FORMAT':
                raise Exception('Check your VCF format!')
            out.write(line)
            for chrom in chroms:
                tmp_vcfs[chrom].write(line)
        else:
            l = line.strip().split()
            try:
                # find the re-mapping info for this SNP
                index = remap_index[l[0]+'.'+l[1]]
                # update chromosome & position
                l[0] = remap.loc[index,'EC3_chrom']
                l[1] = str(remap.loc[index,'EC3_pos'])
                ref = l[3]
                alt = l[4]
                # update info field; existing info is likely invalid now
                l[7] = 'Name='+remap.loc[index,'Name']
                # check to see if there is a REF/ALT allele flip +/- strand flip
                if (ref == remap.loc[index,'EC3_ALT'] and alt == remap.loc[index,'EC3_REF']) or (ref == pairs[remap.loc[index,'EC3_ALT']]):
                    l[3] = remap.loc[index,'EC3_REF']
                    l[4] = remap.loc[index,'EC3_ALT']
                    l[8] = 'GT'
                    for i in range(9,len(l)):
                        geno = l[i].split(':')[0].replace('|','/')
                        l[i] = fixGenos(geno)
                # or a strand flip without an allele flip
                elif ref == pairs[remap.loc[index,'EC3_REF']]:
                    l[3] = remap.loc[index,'EC3_REF']
                    l[4] = remap.loc[index,'EC3_ALT']
                # or if there is no ACTG ALT allele (non-polymorphic)
                elif ref == remap.loc[index,'EC3_REF'] and alt not in pairs:
                    l[4] = remap.loc[index,'EC3_ALT']
                # skip if there is no ACTG REF allele
                elif (ref not in pairs) or (',' in alt):
                    continue
                # re-write line and print to new VCF
                st = ''
                for i in range(len(l)):
                    st += l[i]+'\t'
                tmp_vcfs[l[0]].write(st.strip()+'\n')
            # skip if not in re-map file
            except KeyError:
                continue
    print('Printing VCF...')
    for chrom in chroms:
        print('\t',chrom,sep='')
        tmp_vcfs[chrom].close()
        sortChrom(chrom,out)
        os.remove('.'+chrom+'.vcf')



if __name__ == '__main__':
 
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--vcf', '-v', help="VCF file name", required=True
    )
    parser.add_argument(
        '--chip', '-c', help="SNP chip name. Choose from MNEc2M, MNEc670k, SNP50, or SNP70",required=True
    )
    parser.add_argument(
        '--out','-o', help="Output file name", required=True
    )
    parser.usage=parser.format_usage().replace('usage:','')+'\nYou must include a remap table (e.g., MNEc2M.unique remap.FINAL.csv.gz) in\nthe same directory as this script. All final remap tables can be downloaded\nfrom https://www.animalgenome.org/repository/pub/UMN2018.1003/.'
    args = parser.parse_args()
    if args.chip not in ['MNEc2M', 'MNEc670k', 'SNP50', 'SNP70']:
        raise Exception('Invalid SNP chip name. Choose from MNEc2M, MNEc670k, SNP50, or SNP70')
    chip = args.chip
    try:
        writeVCF(pd.read_csv(glob.glob(chip+'.unique_remap.FINAL.csv*')[0]),args.vcf,args.out)
    except FileNotFoundError:
        print('You must include a remap table (e.g., MNEc2M.unique remap.FINAL.csv.gz) in\nthe same directory as this script. All final remap tables can be downloaded\nfrom https://www.animalgenome.org/repository/pub/UMN2018.1003/.')
