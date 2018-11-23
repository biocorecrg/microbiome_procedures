# This procedure for analysis of fungi ITS amplicons is based on mothur version 1.39.5
# Adapted with changes from https://github.com/krmaas/bioinformatics/blob/master/mothur.fungal.batch

# Create a stability file (a text file with "ID PAIR1 PAIR2")
# for a example:
for i in *_R1_*.gz; do echo `basename $i _R1_fastq.gz` $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file

# Assemble contigs (even though some might not overlap)
make.contigs(processors=12, file=stability.file)
summary.seqs(fasta=stability.trim.contigs.fasta)

# Remove low quality contigs, contigs longer than maxlength, contigs with ambiguities and containing homopolimers longer 13nt
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxlength=554, summary=stability.trim.contigs.summary,  maxhomop=13, maxambig=0)

# Remove dublicates. Leave only unique contigs to reduce computation
unique.seqs(fasta=stability.trim.contigs.good.fasta)
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
summary.seqs(count=stability.trim.contigs.good.count_table, fasta=stability.trim.contigs.good.unique.fasta)

# Pre-cluster contigs at ~1% difference to reduce computation time (4 different nucleotides in an amplicon of 300-500nt). 
# (this step is better to run in bach mode since it is slow)
pre.cluster(fasta=stability.trim.contigs.good.unique.fasta, diffs=4, count=stability.trim.contigs.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.fasta, count=stability.trim.contigs.good.unique.precluster.count_table)

# Remove clusters consisting of only one contig
split.abund(fasta=stability.trim.contigs.good.unique.precluster.fasta, count=stability.trim.contigs.good.unique.precluster.count_table, cutoff=1, accnos=true)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.precluster.abund.count_table) 

# Remove chimeras (sequences formed from two or more biological sequences joined together)
chimera.vsearch(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.precluster.abund.count_table)
remove.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, accnos=stability.trim.contigs.good.unique.precluster.abund.denovo.vsearch.accnos, count=stability.trim.contigs.good.unique.precluster.abund.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table)

# Classify contigs mapping them to the UNITE ITS database (download the mothur release of the databae from https://unite.ut.ee/repository.php)
# using the manually curated dataset of ITS sequences and taxonomy (dynamic files)
classify.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table, reference=UNITEv6_sh_dynamic_s.fasta, taxonomy=UNITEv6_sh_dynamic.tax, cutoff=60, processors=44)

# Remove unwanted lineages (Bacteria-Animalia-Plantae-unclassified-Plantae_unclassified-unknown-Protista)
remove.lineage(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table, taxonomy=stability.trim.contigs.good.unique.precluster.abund.pick.UNITEv6_sh_dynamic_s.wang.taxonomy, taxon=Bacteria-Animalia-Plantae-unclassified-Plantae_unclassified-unknown-Protista)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table)

# Cluster contigs into OTUs using cutoff=0.05 and method=agc (VSEARCH clustering method because the distance matrix cannot be calculated due to different length of contigs)
cluster(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table, method=agc, cutoff=0.05)

# Make OTU table (the output file stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared contains read counts by OTUs for each sample)
make.shared(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table)

# Classify OTUs based on UNITE ITS taxonomy
classify.otu(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.precluster.abund.pick.UNITEv6_sh_dynamic_s.wang.pick.taxonomy)

# Get final number of sequences in each sample
count.groups(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared)

# Calculate alpha diversity (all provided measures; samples with less than subsample=10000 reads are eliminated)
summary.single(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, subsample=10000)

# Calculate beta diversity
dist.shared(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, calc=braycurtis-jest-thetayc, subsample=10000)

# Make a rarefied OTU table for heatmaps, indicator species, etc.
sub.sample(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, size=10000)
 
summary.tax(taxonomy=current, count=current)

### Make a tree of OTU representative sequences ###
# 1. Get representative sequences for each OTU
get.oturep(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, name=stability.trim.contigs.good.names, sorted=bin, method=abundance)

# 2. Rename sequence names into OTU names (needs Python script in this repository)
python rename_oturep.py stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.0.05.rep.fast

# 3. # calculate distances between sequences in clean_repFasta.fasta
pairwise.seqs(fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)

# 4. Build a tree for representative sequences (i.e., final OTUs)
clearcut(phylip=clean_repFasta.phylip.dist)
######

get.current()

###### Files needed for the downstream analysis in R package phyloseq: ######
# 1. OTU counts: stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared
# 2. OTU taxonomy: stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.0.05.cons.taxonomy
# 3. OTU tree: clean_repFasta.phylip.tre
# 4. Samples metadata with alpha-diversity indexes: stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.groups.summary
######
