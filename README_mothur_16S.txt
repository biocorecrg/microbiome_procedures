# This procedure is for the analysis of 16S rRNA amplicons sequenced on MiSeq
# following the protocol at https://www.mothur.org/wiki/MiSeq_SOP
#
# Tested on mothur version 1.39.5
# Run time on 44 CPUs for three mock samples (*V3V4*) in /data (~450K initial reads) was ~20 min,
# for ~5.5M reads, ~6 hrs.
#
#########################################################################################
# 
#	To run this procedure:
#
# 1. Create a stability file (a text file with "ID PAIR1 PAIR2")
#    using in bash, for a example:
#      $ for i in ../data/*_R1_*.gz; do file=`basename $i`; echo ${file%-*} $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#    or
#      $ for i in ../data/*_R1_*.gz; do file=`basename $i`; echo ${file:0:5} $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#    and remove from this file unwanted samples (e.g., those that didn't pass QC).
#
#	 Make sure to keep sample names as short as possible (2-4 characters) and remove from them 
#      any non-letter characters, otherwise it will be impossible to filter them out inside mothur.
#
# 2. If you already worked with the same 16S primers before, 
#	 provide the correct path (inputdir) to the directory where these two files are
#	 in the command align.seqs() below
#	 or copy these files in the current directory: 
#		silva.nr_v132.pcr.good.align 
#    	silva.nr_v132.pcr.good.8mer
#    
#	Otherwise, make a reference database for the alignment (this takes a while to match primers).
#   You will need files silva.nr_v132.align and pcrTest.oligos
#		mothur > pcr.seqs(fasta=silva.nr_v132.align, oligos=./pcrTest.oligos, keepdots=F, processors=44)
#	    mothur > summary.seqs(fasta=silva.nr_v132.pcr.align) 	 
#		mothur > screen.seqs(fasta=silva.nr_v132.pcr.align, start=4965, end=21977)
#		mothur > summary.seqs(fasta=silva.nr_v132.pcr.good.align)
#		
# 3. Download the latest reference from https://mothur.org/wiki/RDP_reference_files
#		wget https://mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz
#		tar -xvzf Trainset16_022016.rdp.tgz -C Trainset/.
#	 Make sure to provide correct path and versions in the command classify.seqs() below.
#
# 4. Put rename_oturep.py in the current folder or modify the run_mothur.sh to provide the path.
#
# 5. Make two files: batch.txt that contains everything in this file and batch_tre.txt that contains the last two commands to generate the OTU tree.
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
#	mothur batch.txt
#	python rename_oturep.py stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.rep.fasta
#	mothur batch_tre.txt
#
# 7. Launch the pipeline:
#      $ qsub -N mothur run_mothur.sh&
#
#
#   Estimation of the sequencing error rate and sequencing depth can be performed on MOCK samples 
#   after running the pipeline. The procedure is described in the end of this file.
#	
#########################################################################################

# Check the SILVA alignment 
summary.seqs(fasta=silva.nr_v132.pcr.good.align, inputdir=/db/silva/mothur/silva_v132/silva_CRG_V3_V4_primers/, outputdir=.)


# Assemble contigs (even though some might not overlap)
make.contigs(processors=44, file=stability.file)
summary.seqs(fasta=stability.trim.contigs.fasta)
# Output File Names:
# stability.trim.contigs.summary

# Remove contigs with ambiguous bases and longer than 500bp
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=500)
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
#align.seqs(processors=44, fasta=stability.trim.contigs.good.unique.fasta, reference=/db/silva/mothur/silva_v132/silva_CRG_V3_V4_primers/silva.nr_v132.pcr.good.align, flip=t)
align.seqs(processors=44, fasta=stability.trim.contigs.good.unique.fasta, reference=silva.nr_v132.pcr.good.align, flip=t)
# Output File Names:
# stability.trim.contigs.good.unique.align
# stability.trim.contigs.good.unique.align.report
# stability.trim.contigs.good.unique.flip.accnos

# this step cannot be skipped because output *summary file is used in the next step
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)

# Clean the alignment (change the start and the end) and remove alignments with more than 8 homopolymers (the maximum number found in the reference alignment)
screen.seqs(processors=44, fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=4965, end=21977, maxhomop=8)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)

# Filter the alignment to remove the overhangs at both ends and remove dots in the alginment to make it shorter
filter.seqs(processors=44, fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

# Pre-cluster sequences allowing for up to 4 (in mothur SOP, 2 are used) differences between sequences. 
# This command will split the sequences by groups and then sort them by abundance and 
# go from most abundant to least and identify sequences that are within 4nt from each other.
# If they are, then they get merged. 
pre.cluster(processors=44, fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=4)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

# Remove clusters with one read (not in mothur SOP, so be careful further on with file names!)
split.abund(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, cutoff=1, accnos=true)

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

# Classify sequences
classify.seqs(cutoff=80, processors=44, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, reference=./Trainset/trainset16_022016.pds/trainset16_022016.pds.fasta, taxonomy=./Trainset/trainset16_022016.pds/trainset16_022016.pds.tax)
# Output File Names:
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy
# stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.tax.summary
# The summary file contains read counts by taxa levels and samples

remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table)

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
# alternatively, use the parameter subsample=10000, then the samples with smaller amount of reads will be eliminated)
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

# Place the code below in another batch_tre.txt and run via queue

# 3. # calculate distances between sequences in clean_repFasta.fasta
# pairwise.seqs(processors=44, fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)

# 4. Build a tree for representative sequences (i.e., final OTUs)
# clearcut(phylip=clean_repFasta.phylip.dist)

#######################################################################

# Files to pass to phyloseq for downstream analysis:

# metadata --> stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.phylip.opti_mcc.groups.summary
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