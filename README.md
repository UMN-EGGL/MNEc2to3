# MNEc2to3
## Files and scripts for mapping probe coordinates from EquCab2 to EquCab3

This is an open project that aims to provide updated and quality controlled genome coordinates for the 
[MNEc2M chip](https://doi.org/10.1186/s12864-017-3943-8). It contains scripts and data files related to the 
remapping of probe sequences. 

This project is designed to be open and welcomes contributions and discussion from the community in order to develop
the best product possible. The data contained herein is available under the [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) 
license and the code is available under the MIT license. People who contribute to the project will be noted on this page
and will be included on any scientific publications that directly result from this work.

**To download the final EquCab3 coordinate lists for all equine SNP arrays, click [here](https://www.animalgenome.org/repository/pub/UMN2018.1003/).**

We are currently in the process of developing software that will remap your SNP array coordinates for you. In the meantime, the Python script `MNEc2to3.py` is included here. It uses the remap tables from the link above to update SNP coordinates and handles strand flips and translocations for VCF files only. For detailed usage instructions, run the following command:

``python3 MNEc2to3.py --help``

## Background
As described in a [previous](https://doi.org/10.1186/s12864-017-3943-8) publication, the equine community banded together to develop
a high density SNP chip based on variants discovered in >20 different horse breeds. The initial discovery effort found over 
22 million sites that were likely polymorphic and went through several filtering steps to identify a subset of 2 million and 
670K SNPs (a proper subset of the 2M) that were put on a commercially available array.

### Data Definition: 
- **2M SNPs**: the SNPs assayed by the limited supply MNEc2M test chip
- **670K SNPs**: a subset of the 2M SNPs available on the commercial array

To be as exhaustive as possible, the initial SNP discovery effort included all the autosomes cataloged by EquCab2 as well as
ChrX and two artificially built chromosomes we referred to as chrUn1 and chrUn2. These contained contigs that were yet unmapped
as of EquCab2, many of which were successfully mapped in EquCab3. Therefore in this remapping process we are dealing with three
different reference genomes:

### Data Definition:
- **MNEc_Ec2**: The custom version of EquCab2 that was used in the 2M SNP discovery (contains EquCab2+ChrUn1+ChrUn2)
- **EquCab2**: The reference genome [available from NCBI](ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/ARCHIVE/BUILD.2.2/)
- **EquCab3**: The newest reference genome also [available from NCBI](ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus)

## Description of Workflow
The main focus and the primary unit of data considered here is the probe sequence that is actually being hybridized by the array.
In other words, what we are mapping are the probe sequences 
that were put on the actual array in relation to their coordinates from EquCab2 to EquCab3. Each SNP on the chip is interrogated
by the array by hybridizing a probe to genomic DNA. This probe was designed using the MNEc_Ec2 reference genome.

The simplest use case here is to extract the probe sequence used by the MNEc chip and BLAST it to the new genome build. The first
script in the `scripts/` directory called `01_make_probeseq.py` does just that. The script takes an input a Fasta file as well as
a VCF file containing the MNEc SNP coordinates and outputs a FASTA file containing probe sequences. The Fasta file and VCF are not 
included in this project due to limitations of GitHub hosting. But the probe sequences generated from the **MNEc_Ec2** reference
genome used in the array design are available in `data/MNEc2M_probeseq.fa.gz`.

From here, the probe sequences were BLASTED against three different reference genomes: MNEc_Ec2 (as a sanity check), EquCab2 (as 
another sanity check), and EquCab3. The code used to do this is in `scripts/02_do_blast.sh`. Briefly, the script builds some
custom blast databases with megablast indexes and then performs a BLAST on each of the MNEc 2M probe sequences and only retains
100%matches (of the 71mer).

Again, some of the input files are
not available because of their size, but each of the blast output files can be found in `data/`:
`data/MNEc2M.blast.EquCab2.txt.gz`,`data/MNEc2M.blast.EquCab3.txt.gz`, and `data/MNEc2M.blast.MNEc_Ec2.txt.gz`. These will be the 
most useful files for doing the remapping analysis and the updated probe quality control assessment. 

**Note:** The output of the BLAST program (and the headers of the blast files above) are: 
```
MNEc Coordinate | BLASTED Chromosome | % Identity | Blasted start | Blasted end
```

The final step in the workflow is performed by `scripts/03_process_BLAST.py`. This script borrows some annotations about the 
SNPs from [PonyTools](https://github.com/schae234/PonyTools) and iterates over the results produced by BLAST in order to add
some information about the mapping. The results from this script are in `data/MNEc2M.probe_blast_counts.csv.gz`. The table looks
like:
```
In [135]: info.head(n=20)
Out[135]: 
                               MNEcID       AffySNPID    ProbeSetID  chrom       pos AlleleA AlleleB REF ALT  in670  legacy              id  MNEc_Ec2  EquCab3  EquCab2
0               MNEc.2.10.9554699.VIP  Affx-101861294  AX-104415554  chr10   9554699       C       G   C   G  False    True   chr10.9554699         1        1        1
1               MNEc.2.10.9554699.VIP  Affx-101861294  AX-103524731  chr10   9554699       C       G   C   G  False    True   chr10.9554699         1        1        1
2               MNEc.2.2.13074277.VIP  Affx-103007064  AX-103861101   chr2  13074277       T       C   C   T  False    True   chr2.13074277         1        1        1
3               MNEc.2.2.13074277.VIP  Affx-103007064  AX-102948897   chr2  13074277       T       C   C   T  False    True   chr2.13074277         1        1        1
4              MNEc.2.21.30666626.VIP  Affx-102471109  AX-104189538  chr21  30666626       A       G   G   A  False    True  chr21.30666626         1        1        1
5              MNEc.2.21.30666626.VIP  Affx-102471109  AX-104916861  chr21  30666626       A       G   G   A  False    True  chr21.30666626         1        1        1
6              MNEc.2.23.22999655.VIP  Affx-101668032  AX-104553122  chr23  22999655       A       C   C   A  False    True  chr23.22999655         1        1        1
7              MNEc.2.23.22999655.VIP  Affx-101668032  AX-104971885  chr23  22999655       A       C   C   A  False    True  chr23.22999655         1        1        1
8              MNEc.2.18.66493737.VIP  Affx-102760722  AX-104041705  chr18  66493737       T       C   T   C   True    True  chr18.66493737         1        1        1
9              MNEc.2.18.66493737.VIP  Affx-102760722  AX-103119925  chr18  66493737       T       C   T   C   True    True  chr18.66493737         1        1        1
10      MNEc.2.6.6431969.BIEC2-985609  Affx-101987248  AX-104524509   chr6   6431969       A       G   A   G   True    True    chr6.6431969         1        1        1
11      MNEc.2.6.6431969.BIEC2-985609  Affx-101987248  AX-103467496   chr6   6431969       A       G   A   G   True    True    chr6.6431969         1        1        1
12      MNEc.2.6.6437698.BIEC2-985637  Affx-102389295  AX-104231007   chr6   6437698       A       G   G   A   True    True    chr6.6437698         1        1        1
13      MNEc.2.6.6562393.BIEC2-985785  Affx-101867915  AX-104411696   chr6   6562393       A       G   A   G   True    True    chr6.6562393         1        1        1
14      MNEc.2.6.6562393.BIEC2-985785  Affx-101867915  AX-103521326   chr6   6562393       A       G   A   G   True    True    chr6.6562393         1        1        1
15      MNEc.2.6.6671817.BIEC2-985848  Affx-101434215  AX-104672882   chr6   6671817       T       C   C   T   True    True    chr6.6671817         1        1        1
16      MNEc.2.6.6671817.BIEC2-985848  Affx-101434215  AX-103715826   chr6   6671817       T       C   C   T   True    True    chr6.6671817         1        1        1
17      MNEc.2.6.6698388.BIEC2-985859  Affx-102562234  AX-102951995   chr6   6698388       A       G   G   A   True    True    chr6.6698388         1        1        1
18      MNEc.2.6.6698388.BIEC2-985859  Affx-102562234  AX-102951996   chr6   6698388       A       G   G   A   True    True    chr6.6698388         1        1        1
19  MNEc.2.20.33196687.MHCClIIDRB2-38  Affx-102831127  AX-104005803  chr20  33196687       T       C   T   C  False    True  chr20.33196687         2        2        2
```

The first three columns (`MNEcID,AffySNPID,ProbeSetID`) are IDs from the Chip and the corresponding Affymetrix probe and SNP IDs. The next six columns `chrom,pos,AlleleA,AlleleB,REF,ALT` are the original MNEc_Ec2 coordinates and allele information. The next three columns `in670,legacy,id` contain some very useful information:
- in670: indicates membership on the 670k commercial array
- legacy: indicates if the SNP was on the Illumina 54K and 65K arrays (all SNPs from legacy chips were included on the MNEc array)
- id: this is an ID indicating the original MNEc_Ec2 SNP position

The last three columns: `MNEc_Ec2,EquCab3,EquCab2` contain the number of times the probe sequence for each SNP BLASTED to each
reference genome. You can see this is mostly 1 here which is a good indication. You could confidently find the updated coordinates 
in EquCab3 by searching for the ID in the BLAST output file. 

However, for the last SNP we see that, in each genome, the probe sequence blasted to two places. However we also see that the SNP
was a legacy SNP. The SNP was included on the chip regardless of its initial quality control status due to requirements for 
backwards compatibility with the previous generation chips. These are the columns that will indicate successful remapping of probe
sequence from EquCab2 to EquCab3.

## Analysis
### Mean and Median number of 100% probe mappings from EquCab2 to EquCab3:
```
In [148]: info.EquCab3.mean()
Out[148]: 1.0436152098740357

In [149]: info.EquCab3.median()
Out[149]: 1.0
```
On average probes had a unique mapping from EquCab2 to EquCab3.

### The distribution 100% probe mappings from EquCab2 to EquCab3

| # probe BLAST matches | Frequency |                                                 
| --- | --- |                                                                      
|1    | 1690599|                                                                   
|2    |   22125|                                                                   
|3    |   13865|                                                                   
|4    |    5866|                                                                   
|5    |     373|                                                                   
|6    |     116|                                                                   
|7    |      49|                                                                   
|8    |      30|                                                                   
|9    |      12|                                                                   
|10   |      12|                                                                   
|11   |       7|                                                                   
|12   |       5|                                                                   
|13   |       7|                                                                   
|14   |       5|                                                                   
|16   |       1|                                                                   
|17   |       5|                                                                   
|19   |       3|                                                                   
|20   |       2|                                                                   
|21   |       3|                                                                   
|27   |       4|                                                                   
|35   |       2|                                                                   
|112  |       2|                                                                   
|125  |       2|                                                                   
|152  |       2|                                                                   
|181  |       2|                                                                   
|224  |       2|                                                                   
|254  |       2|                                                                   
|266  |       2|                                                                   
|294  |       2|                                                                   
|319  |       2|                                                                   
|415  |       2|    


### Discussion
Discussion can happen by opening a GitHub issue [here](https://github.com/UMN-EGGL/MNEc2to3/issues). This lets us track progress
and discuss issues using tools build into Github.
