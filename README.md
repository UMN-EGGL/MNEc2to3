# MNEc2to3
## Files and script for mapping probe coordinates from EquCab2 to EquCab3

This is an open project that aims to provide updated and quality controlled genome coordinates for the 
[MNEc2M chip](https://doi.org/10.1186/s12864-017-3943-8). It contains scripts and data files related to the 
remapping of probe sequences. 

This project is designed to be open and welcomes contributions and discussion from the community in order to develop
the best product possible. The data contained herein is available under the [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) 
license and the code is available under the MIT license. People who contribute to the project will be noted on this page
and will be included on any scientific publications that directly result from this work.

## Background
As described in a [previous](https://doi.org/10.1186/s12864-017-3943-8) publication, the equine community banded together to develop
a high density SNP chip based on variants discovered in >20 different horse breeds. The initial discovery effort found over 
22 million sites that were likely polymorphic and went through several filtering steps to identify a subset of 2 million and 
670K SNPs (a proper subset of the 2M) that were put on a commercially available array.

### Data Definition: 
- **2M SNPs**: the SNPs assayed by the limited supply MNEc2M test chip
- **670K SNPs**: a subset of the 2M SNPs available on the commercial array

To be as exhaustive as possible, the initial SNP discovery effort included all the autosomes cataloged by EquCab2 as well as
ChrX and two artificially build chromosomes we reffered to as chrUn1 and chrUn2. These contained contigs that were yet unmapped
as of EquCab2, many of which were successfully mapped in EquCab3. Therefore in this remapping process we are dealing with three
different reference genomes:

### Data Definition:
- **MNEc_Ec2**: The custom version of EquCab2 that was used in the 2M SNP discovery (contains EquCab2+ChrUn1+ChrUn2)
- **EquCab2**: The reference genome [available from NCBI](ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/ARCHIVE/BUILD.2.2/)
- **EquCab2**: The newset reference genome also [available from NCBI](ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus)

## Description of Workflow
Due to the fact that we are, for now, stuck with the design of the MNEc 2M and 670K SNPs, the unit of data considered here is
the probe sequence that is actually being hybridized by the array. In other words, what we are mapping are the probe sequences 
that were put on the actual array in relation to their coordinates from EquCab2 to EquCab3. Each SNP on the chip is interrogated
by the array by hybridizing a probe to genomic DNA. This probe was designed using the MNEc_Ec2 reference genome.

The simplest use case here is to extract the probe sequence used by the MNEc chip and BLAST it to the new genome build. The first
script in the `scripts/` directory called `01_make_probeseq.py` does just that. The script takes an input a Fasta file as well as
a VCF file containing the MNEc SNP coordinates and outputs a FASTA file containing probe sequnces. The Fasta file and VCF are not 
included in this project due to limitations of GitHub hosting. But the probesequences generated from the **MNEc_Ec2** reference
genome used in the array design are available in `data/MNEc2M_probeseq.fa.gz`.

From here, the probe sequences were BLASTED against three different reference genomes: MNEc_Ec2 (as a sanity check), EquCab2 (as 
another sanity check), and EquCab3. The code used to do this is in `scripts/02_do_blast.sh`. Again, some of the input files are
not available because of their size, but each of the blast output files can be found in `data/`:
`data/MNEc2M.blast.EquCab2.txt.gz`,`data/MNEc2M.blast.EquCab3.txt.gz`, and `data/MNEc2M.blast.MNEc_Ec2.txt.gz`. These will be the 
most useful files for doing the remapping analysis and the updated probe quality control assessment. 
