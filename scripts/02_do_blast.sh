

# Assumes that FASTA files are in /project/fasta and Blast databases are in /project/blast 
cd /project/blast

# There are three FASTA files of interest...
# 1) The one we used the build the MNEC2M chip, which was a hybrid containing a 
# 2) Ensemble EquCab2
# 3) Ensemble EquCab3 

makeblastdb -in /project/Data/Fasta/EquCab2_wChrun1_2/Equus_cab_nucl_wChrUn1_2.fasta -dbtype nucl -parse_seqids -out MNEc_Ec2
makeblastdb -in /project/Data/Fasta/EquCab2/EquCab2.fa -dbtype nucl -parse_seqids -out EquCab2
makeblastdb -in /project/Data/Fasta/EquCab3/EquCab3.fa -dbtype nucl -parse_seqids -out EquCab3

# Make an index with the largest word size possible (15) in the /project/blast directory
makembindex -input MNEc_Ec2 -iformat blastdb -nmer 15
makembindex -input EquCab2 -iformat blastdb -nmer 15
makembindex -input EquCab3 -iformat blastdb -nmer 15

# Perform a blastn with each blast database using the index looking for
# perfect 70mer identity
parallel  'blastn -db /project/blast/{} \
           -word_size 70 -perc_identity 100 \
           -outfmt "6 qseqid sseqid pident sstart send" \
           -query MNEc2M_probeseq.fa -use_index True \
           -index_name /project/blast/{} -num_threads 4 \
           -out MNEc2M.blast.{}.txt' ::: EquCab2 EquCab3 MNEc_Ec2
