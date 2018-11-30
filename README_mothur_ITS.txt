# This procedure is for the analysis of fungi ITS amplicons sequenced on MiSeq
# Tested on mothur version 1.39.5
# Running time on 44 CPUs was ~30 min for initial 1M reads and ~3hrs for 5,5M reads.
#
#########################################################################################
# To run this procedure:
#
# 1. Create a stability file (a text file with "ID PAIR1 PAIR2")
#    using in bash, for a example:
#      $ for i in ../data/*_R1_*.gz; do file=`basename $i`; echo ${file%-*} $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file
#    and remove from this file unwanted samples (e.g., those that didn't pass QC)
#
# 2. Download the mothur release of the UNITE ITS database from https://unite.ut.ee/repository.php and
#    place files UNITEv6_sh_dynamic_s.fasta and UNITEv6_sh_dynamic_s.tax in the current folder.
#
# 3. Put rename_oturep.py in the current folder or modify the run_mothur.sh to provide the path.
#
# 4. Make two files: batch.txt that contains everything in this file and batch_tre.txt that contains the last two commands to generate the OTU tree.
#
# 5. Make the file run_mothur.sh to run mothur in a batch mode on your cluster,
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
# 6. Launch the pipeline:
#      $ qsub -N mothur run_mothur.sh&
#########################################################################################

# This procedure avoids pre-clustering and removing rare clusters.
# It produces a lot of OTUs, many of which are singletons of the same species and can be 
# collapsed and filtered further in phyloseq (it is inconvenient to do in mothur)
# Thus 4881 OTUs produced in the end of this procedure, after removing (*) in taxa and collapsing by same species
# became 846 OTUs.
#
####################################################################################################################

# Assemble contigs (even though some might not overlap)
make.contigs(processors=44, file=stability.file)
summary.seqs(fasta=stability.trim.contigs.fasta)

# Remove low quality contigs, contigs longer than maxlength, contigs with ambiguities and containing homopolimers longer 13nt
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxlength=554, summary=stability.trim.contigs.summary,  maxhomop=13, maxambig=0)
summary.seqs(fasta=stability.trim.contigs.good.fasta)

# Remove dublicates. Leave only unique contigs to reduce computation
unique.seqs(fasta=stability.trim.contigs.good.fasta)
# Output File Names:
# stability.trim.contigs.good.names
# stability.trim.contigs.good.unique.fasta

count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
# Output File Names:
# stability.trim.contigs.good.count_table

summary.seqs(count=stability.trim.contigs.good.count_table, fasta=stability.trim.contigs.good.unique.fasta)
# Output File Names:
# stability.trim.contigs.good.unique.summary

# Remove chimeras (sequences formed from two or more biological sequences joined together)
chimera.vsearch(processors=44, fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)
remove.seqs(fasta=stability.trim.contigs.good.unique.fasta, accnos=stability.trim.contigs.good.unique.denovo.vsearch.accnos, count=stability.trim.contigs.good.count_table)

# Output File Names:
# stability.trim.contigs.good.unique.pick.fasta
# stability.trim.contigs.good.pick.count_table

summary.seqs(fasta=stability.trim.contigs.good.unique.pick.fasta, count=stability.trim.contigs.good.pick.count_table)

# Classify contigs mapping them to the UNITE ITS database (download the mothur release of the database from https://unite.ut.ee/repository.php)
# using the manually curated dataset of ITS sequences and taxonomy (dynamic files)
# Using default cutoff of 80%
classify.seqs(processors=44, fasta=stability.trim.contigs.good.unique.pick.fasta, count=stability.trim.contigs.good.pick.count_table, reference=UNITEv6_sh_dynamic_s.fasta, taxonomy=UNITEv6_sh_dynamic_s.tax)
# Output File Names:
# stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.taxonomy
# stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.tax.summary
# stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.flip.accnos

# Remove unwanted lineages (Bacteria-Animalia-Plantae-unclassified-Plantae_unclassified-unknown-Protista)
remove.lineage(fasta=stability.trim.contigs.good.unique.pick.fasta, count=stability.trim.contigs.good.pick.count_table, taxonomy=stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.taxonomy, taxon=Protozoa-Chromista-Eukaryota_kgd_Incertae_sedis-Bacteria-Animalia-Plantae-Plantae_unclassified-unknown-Protista-Fungi_unclassified-unclassified_Fungi)
# Output File Names:
# stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.pick.taxonomy
# stability.trim.contigs.good.unique.pick.pick.fasta
# stability.trim.contigs.good.pick.pick.count_table

summary.seqs(fasta=stability.trim.contigs.good.unique.pick.pick.fasta, count=stability.trim.contigs.good.pick.pick.count_table)

# Cluster contigs into OTUs using cutoff=0.05 and method=agc (VSEARCH clustering method because the distance matrix cannot be calculated due to different length of contigs)
cluster(fasta=stability.trim.contigs.good.unique.pick.pick.fasta, count=stability.trim.contigs.good.pick.pick.count_table, method=agc, cutoff=0.05)
# Output File Names:
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.list

# Make OTU table (the output file stability.*.unique_list.shared contains read counts by OTUs for each sample)
make.shared(list=stability.trim.contigs.good.unique.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.pick.pick.count_table)
# Output File Names:
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.shared

# Classify OTUs based on UNITE ITS taxonomy
classify.otu(list=stability.trim.contigs.good.unique.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.pick.UNITEv6_sh_dynamic_s.wang.pick.taxonomy)
# Output File Names:
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.cons.taxonomy
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.cons.tax.summary

# Get final number of sequences in each sample
count.groups(shared=stability.trim.contigs.good.unique.pick.pick.agc.unique_list.shared)

# Calculate alpha diversity (all provided measures; subsampling will be done by the number of reads in the smallest file; alternatively, use the parameter subsample=10000, then the samples with smaller amount of reads will be eliminated)
summary.single(shared=stability.trim.contigs.good.unique.pick.pick.agc.unique_list.shared, subsample=T)
# Output File Names:
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.groups.summary


### Make a tree of OTU representative sequences ###
# 1. Get representative sequences for each OTU
get.oturep(list=stability.trim.contigs.good.unique.pick.pick.agc.unique_list.list, fasta=stability.trim.contigs.good.unique.pick.fasta, name=stability.trim.contigs.good.names, sorted=bin, method=abundance)
# Output File Names:
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.rep.names
# stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.rep.fasta

# 2. Rename sequence names into OTU names (needs Python script in this repository)
# python rename_oturep.py stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.rep.fasta

# Place the code below in another batch_tre.txt and run via queue

# 3. # calculate distances between sequences in clean_repFasta.fasta
# pairwise.seqs(processors=44, fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)

# 4. Build a tree for representative sequences (i.e., final OTUs)
# clearcut(phylip=clean_repFasta.phylip.dist)
######

# Files to pass to phyloseq for downstream analysis:
# metadata --> stability.trim.contigs.good.unique.pick.pick.agc.unique_list.groups.summary
# OTU_taxonomy --> stability.trim.contigs.good.unique.pick.pick.agc.unique_list.0.05.cons.taxonomy
# OTU_counts --> stability.trim.contigs.good.unique.pick.pick.agc.unique_list.shared
# OTU_tree --> clean_repFasta.phylip.tre
