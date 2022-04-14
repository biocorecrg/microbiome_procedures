# This procedure is for the analysis of 16S rRNA amplicons sequenced on MiSeq
# following the protocol at https://www.mothur.org/wiki/MiSeq_SOP
#
# Last time tested using mothur v.1.44.1 
# running it from the Singularity image mothur_1.44.1--h59c8b9c_1.sif
# made with the command from docker pull quay.io/biocontainers/mothur:<tag>
#   singularity build mothur_1.44.1--h59c8b9c_1.sif docker://biocontainers/mothur:1.44.1--h59c8b9c_1
#
# using the latest RDP reference trainset (see below the link) and SILVA v.132
######################
# For the initial commit, it was tested on mothur version 1.39.5
# Run time on 44 CPUs for three mock samples (*V3V4*) in /data (~450K initial reads) was ~20 min.
#
# For runs with real samples, time depends on reachness of samples. 
# Thus a run of ~5.5M reads (44 samples) took ~6 hrs; of 20.7M reads (96 samples), ~13 hrs; 
# and yet another run of ~19.8M reads(96 samples) took ~7 hrs.
#
#########################################################################################
# Run everything on control positive samples first (MOCK), preferably step-by-step to adjust filtering parameters for contigs and alignments
# and to estimate the error.
# Then run everything at once on all samples.
#########################################################################################
# 
#	To run this procedure:
#
# 1. Create a stability file (a text file with "ID PAIR1 PAIR2", separated by space)
#    using in bash, for a example:
#      $ for i in ../data/*_R1_*.gz; do file=`basename $i`; echo ${file%-*} $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#    or
#      $ for i in ../data/*_R1_*.gz; do file=`basename $i`; echo ${file:0:5} $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#    or, for files like this 30427-90_S46_L001_R1_001.fastq.gz 
#	 or 30745-LAL-0025-A_S29_L001_R1_001.fastq.gz, with the same first 6 characters
#      $ for i in ../../data/*_R1_*.gz; do file=`basename $i`; file=${file:6}; file=${file%S*}; file=`echo $file | sed s/[_-]//g`; echo $file  $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#	 and remove from this file unwanted samples (e.g., those that didn't pass QC).
#
#	 Make sure to keep sample names as short as possible (2-4 characters) and remove from them 
#      any non-letter characters, otherwise it will be impossible to filter them out inside mothur.
#
# An example of stability.file
# A1 ../data/A1-V3V4_S3_L001_R1_001.fastq.gz ../data/A1-V3V4_S3_L001_R2_001.fastq.gz
# A2 ../data/A2-V3V4_S11_L001_R1_001.fastq.gz ../data/A2-V3V4_S11_L001_R2_001.fastq.gz
# A3 ../data/A3-V3V4_S19_L001_R1_001.fastq.gz ../data/A3-V3V4_S19_L001_R2_001.fastq.gz
# B1 ../data/B1-V3V4_S27_L001_R1_001.fastq.gz ../data/B1-V3V4_S27_L001_R2_001.fastq.gz
# B2 ../data/B2-V3V4_S35_L001_R1_001.fastq.gz ../data/B2-V3V4_S35_L001_R2_001.fastq.gz
#
# 2. If you already worked with the same 16S primers and same SILVA version before, 
#	 provide the correct path (inputdir) to the directory where these two files are
#	 in the command align.seqs() below
#	 or copy these files in the current directory: 
#		silva.nr_v132.pcr.good.align 
#    	silva.nr_v132.pcr.good.8mer
#    
#	Otherwise, make a reference database for the alignment (this takes a while to match primers).
#   You will need files silva.nr_v132.align and pcrTest.oligos
#   pcrTest.oligos for V3-V4 primers looks like that:
# forward CCTACGGGNGGCWGCAG
# reverse GACTACHVGGGTATCTAATCC
#
#		mothur > pcr.seqs(fasta=silva.nr_v132.align, oligos=./pcrTest.oligos, keepdots=F, processors=44)
#	    mothur > summary.seqs(fasta=silva.nr_v132.pcr.align) 	 
#		mothur > screen.seqs(fasta=silva.nr_v132.pcr.align, start=4965, end=21977)
#		mothur > summary.seqs(fasta=silva.nr_v132.pcr.good.align)
#		
# 3. Download the latest reference from https://mothur.org/wiki/RDP_reference_files
#		wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.pds.tgz
#		tar -xvzf trainset18_062020.pds.tgz 
#	 Make sure to provide correct path and versions in the command classify.seqs() below.
#
# 4. Put rename_oturep.py in the current folder or modify the run_mothur.sh to provide the path.
#
# 5. Make the following 4 files (as specified below: 
#     (1) batch_start.txt that contains the mothur commands until the command screen.seqs() after contig assembly
#		  (2) batch_align.txt 

#
#    and batch_tre.txt that contains the last two commands to generate the OTU tree:
# 		pairwise.seqs(processors=44, fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)
# 		clearcut(phylip=clean_repFasta.phylip.dist)
#
# 6. Make the file run_mothur.sh to run mothur in a batch mode on your cluster,
#    for example, using SGE (modify as needed, also change the paramater "processors" in commands below if needed)
#      #!/bin/bash
#      #$ -cwd
#      #$ -q [your queue]
#      #$ -pe smp 44
#      #$ -j y
#      #$ -l virtual_free=450G
#      #$ -l h_rt=100:00:00
#	 mothur batch_start.txt
#	#python rename_oturep.py stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.rep.fasta
#	#mothur batch_tre.txt
#
# 7. Launch the pipeline:
#      $ qsub -N mothur run_mothur.sh&
#
# 8. Look in the log file at the output of the first two commands
#		    make.contigs(processors=44, file=stability.file)
#		    summary.seqs(fasta=stability.trim.contigs.fasta)
#   to decide on the parameters in the command
#     screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=500)
#
# 9. Comment the first two lines and run mothur.sh with batch.txt
#
#
# 10. Evaluate the number of assembled contigs and remove manually samples with low number of contigs
#       $mothur > remove.groups(group=stability.contigs.groups, groups=A-B-C-45-66A)
#
# 11. Replace content of batch.txt with this file.
#    Uncomment two running lines in run_mothur.sh.
#    Launch the pipeline again.
#
#########################################################################################
#
#   Estimation of the sequencing error rate and sequencing depth can be performed on MOCK samples 
#   after running the pipeline. The procedure is described in the end of this file.
#	
#########################################################################################
########### This is batch_start.txt #########

# Check the SILVA alignment 
summary.seqs(fasta=silva.nr_v132.pcr.good.align, inputdir=/db/silva/mothur/silva_v132/silva_CRG_V3_V4_primers/, outputdir=.)

# Assemble contigs (even though some might not overlap)
make.contigs(processors=44, file=stability.file)
summary.seqs(fasta=stability.trim.contigs.fasta)
# Output File Names:
# stability.trim.contigs.summary

########### this is batch_align.txt  #########

# Decide on the parameters in screen.seqs() looking at the summary.seqs() output for assembled contigs in the end of log files 
# Remove contigs with ambiguous bases and longer than maxlength (in bp) and with maxhomop

screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=500, maxhomop=7)
summary.seqs(fasta=stability.trim.contigs.good.fasta)

# Leave unique contigs only by removing duplicates (to speed up the further mapping)
# This step requires a lot of memory and cpu, so do it through batch submission
unique.seqs(fasta=stability.trim.contigs.good.fasta)
# Output File Names:
# stability.trim.contigs.good.names
# stability.trim.contigs.good.unique.fasta

# Generate a table where the rows are the names of unique sequences and # the columns are the names of groups. 
# For each unique sequence, the table contains counts in each sample (group).
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
# Output File Names:
# stability.trim.contigs.good.count_table

# Align contigs to the reference alignment
align.seqs(processors=44, fasta=stability.trim.contigs.good.unique.fasta, reference=/db/silva/mothur/silva_v132/silva_CRG_V3_V4_primers/silva.nr_v132.pcr.good.align, flip=t)
#align.seqs(processors=44, fasta=stability.trim.contigs.good.unique.fasta, reference=silva.nr_v132.pcr.good.align, flip=t)
# Output File Names:
# stability.trim.contigs.good.unique.align
# stability.trim.contigs.good.unique.align.report
# stability.trim.contigs.good.unique.flip.accnos

# this step cannot be skipped because the output *summary file is used in the next step
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)

########### this is batch_classify.txt  (including get.oturep() to make the tree) #########

# Clean the alignment: (change the start and the end) and remove alignments with more than 8 homopolymers 
# (the maximum number found in the reference alignment) and minlength - 
# to choose these parameters, see the summary.seqs() output for alignment in the end of log files

screen.seqs(processors=44, fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=4965, end=21977, maxhomop=8, minlength=405)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)

# Filter the alignment to remove the overhangs at both ends and remove dots (gaps) in the alginment 
filter.seqs(processors=44, fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

# Pre-cluster sequences allowing for up to 4 (in mothur SOP, 2 are used) differences between sequences. 
# This command will split the sequences by groups and then sort them by abundance and 
# go from the most to least abundant, identify sequences that are within 4nt from each other, and merge them.
pre.cluster(processors=44, diffs=2, fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

# Remove clusters with one read (not in mothur SOP, so be careful further on with file names!)
split.abund(cutoff=1, accnos=true, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

# Remove chimeras. 
# Chimeras are sequences formed from two or more biological sequences joined together. 
# Amplicons with chimeric sequences can form during PCR when closely related sequences are amplified.
# The majority of chimeras are believed to arise from incomplete extension: During subsequent 
# cycles of PCR, a partially extended strand can bind to a template derived from a different 
# but similar sequence. This then acts as a primer that is extended to form a chimeric sequence.
# Remove chimera sequences from count file (mothur: "if a sequence is flagged as chimeric in one sample, 
# the default (dereplicate=F) is to remove it from all samples. Our experience suggests 
# that this is a bit aggressive since we've seen rare sequences get flagged as chimeric 
# when they're the most abundant sequence in another sample.")
chimera.vsearch(processors=44, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.count_table, dereplicate=t)

# Remove chimera sequences from fasta file
remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.accnos)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table)

# Classify sequences using RDP trainset to remove non bacteria taxa
classify.seqs(cutoff=80, processors=44, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, reference=./trainset18_062020.pds/trainset18_062020.pds.fasta, taxonomy=./trainset18_062020.pds/trainset18_062020.pds.tax)
# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.tax.summary
# The summary file contains read counts by taxa levels and samples

remove.lineage(taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy)
summary.tax(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table)

###### this is batch_otu.txt #######
#
## NOTE: if the log file from running the batch_classify.txt says *** No contaminants to remove ***
# files *pick.pick.* won't be generated; therefore, to proceed with the code below, you need to copy a few files:
cp -R -u -p stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta
cp -R -u -p stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table


# Calculate pairwise distances among unique sequences 
dist.seqs(processors=44, output=lt, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta, cutoff=0.03)

# Cluster sequences into OTUs based on distances
cluster(phylip=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table)

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, label=0.03)
# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.shared

# Get the consensus taxonomy for each OTU 
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.taxonomy, label=0.03)
# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.0.03.cons.taxonomy
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.0.03.cons.tax.summary

# Bin sequences into phylotypes according to their taxonomic classification
phylotype(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.taxonomy)

# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.sabund
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.rabund
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.list

# Get the shared file at the genus level (label=1; up to corresponding to up Kingdom level)
make.shared(label=1, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table)

# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.shared

classify.otu(label=1, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.taxonomy)

# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.cons.taxonomy
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.cons.tax.summary

# The first file shows OTUs at the genus level. If exploring this file, undesired taxa are found, 
# you can go back to remove.lineage() and remove them, then repeat all commands from that point down.
# The summary file is the table of OTU counts by samples and taxa levels

# Get the final number of sequences in each sample
count.groups(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.shared)

# Calculate alpha diversity (all provided measures; subsampling will be done by the number of reads in the smallest file; 
# alternatively, use the parameter subsample=10000, then the samples with smaller amount of reads will be eliminated;
# in this case, make sure to remove.groups to make all files that will be used by PhyloSeq to correspond)
summary.single(label=0.03, subsample=T, shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.shared)

# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.groups.ave-std.summary
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.groups.summary


#######################################################################

### Make a tree of OTU representative sequences ###

# 1. Get representative sequences for each OTU
get.oturep(sorted=bin, label=1, method=abundance, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.list, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta, name=stability.trim.contigs.good.names)
# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.rep.names
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.rep.fasta

# 2. Rename sequence names into OTU names (needs Python script in this repository)
# python rename_oturep.py stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.rep.fasta

# Place the code below (uncomment pairwise.seqs and clearcut) in batch_tre.txt to run via queue

# 3. # calculate distances between sequences in clean_repFasta.fasta
# pairwise.seqs(processors=44, fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)

# 4. Build a tree for representative sequences (i.e., final OTUs)
# clearcut(phylip=clean_repFasta.phylip.dist)

#######################################################################

# File to see read counts by taxa levels and samples is
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.tax.summary

# for a specific sample and taxa level it can be explored, e.g., using the command
#$ cut -f1-3,99 stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.tax.summary | awk '{if ( ($1 == 2 && $4!=0) || ($1 == 3 && $4!=0 ) || ($1 == 4 && $4!=0) || ($1 == 5 && $4 !=0) || ($1==6 && $4!=0) ) print $0;}'

#######################################################################

# Files to pass to phyloseq for downstream analysis:

# metadata --> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.groups.ave-std.summary
# OTU_taxonomy --> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.1.cons.taxonomy
# OTU_counts --> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.tx.shared
# OTU_tree --> clean_repFasta.phylip.tre

# For 3 mock samples in /data, 51 final OTUs was obtained

#######################################################################

# Assessing sequencing error rates using MOCK communities
# Run the commands below from within mothur

# mothur > get.groups(groups=HM782-HM783, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta)

# Selected 249 sequences from your fasta file.
# Selected 287591 sequences from your count file.

# Output File names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table

# Copy from the previous project or donwload from the BEI resources the fast file with MOCK species sequences
#  cp /users/project/sacalalengua/sacalall2/anno/HM782_3_MOCK.fasta .

# mothur > seq.error(reference=HM782_3_MOCK.fasta, aligned=F, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table)

# Output:
# Multiply error rate by 100 to obtain the percent sequencing errors.
# Overall error rate:	5.95575e-05
# Errors	Sequences
# 0	285321
# 1	887
# 2	439
# ...

# The error rate is below 0.006%

# Cluster sequences in OTUs
# mothur > dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta, cutoff=0.03)
# mothur > cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table)
# mothur > make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, label=0.03)
# mothur > rarefaction.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.shared)

# How many OTUs is MOCK samples
# $ head *HM*rabund | cut -f2
# ==> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.HM782.rabund <==
# 30
# ==> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.HM783.rabund <==
# 77

# Make rarefaction plots in R 
#R
#setwd(" [current directory] ")
#data <- read.table(file="stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.groups.rarefaction", header=T)
#pdf("rarefaction_plot_mock.pdf")
#plot(x=data$numsampled, y=data$X0.03.HM783, xlab="Number of Tags Sampled",ylab="OTUs", type="l", col="red", font.lab=3, main="Rarefaction plot mock community run 1")
#points(x=data$numsampled, y=data$X0.03.HM782, type="l", col="black")
#legend("bottomright", c("HM783D (staggered)", "HM782D (equimolar)"), fill=c("red", "black"))
#dev.off()
#quit()

# Explore the plots
# $ gnome-open rarefaction_plot_mock.pdf

####### THE END ###############################################################################################
