# This procedure is based on mothur version 1.39.5
# First you need to create a stability file (a text file with "ID PAIR1 PAIR2")
# This comes from https://github.com/krmaas/bioinformatics/blob/master/mothur.fungal.batch
#example:
for i in *_R1_*.gz; do echo `basename $i _R1_fastq.gz` $i `echo $i | sed s/_R1_/_R2_/g`; done >> stability.file

# Assemble contigs (even though some might not overlap)
make.contigs(processors=12, file=stability.file)
summary.seqs(fasta=stability.trim.contigs.fasta)

# Removing low quality contigs. We try to be strict asking for no ambiguities and max 13 homopolimer = stretch of 13 or more of the same base in a row)
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxlength=554, summary=stability.trim.contigs.summary,  maxhomop=13, maxambig=0)

# Make unique sequences to reduce the computation
unique.seqs(fasta=stability.trim.contigs.good.fasta)
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
summary.seqs(count=stability.trim.contigs.good.count_table, fasta=stability.trim.contigs.good.unique.fasta)

#Pre-cluster. We used 1% difference to reduce computation time (4 differences in an amplicon ~300/500). You might want to run this step in bach mode since is slow.
pre.cluster(fasta=stability.trim.contigs.good.unique.fasta, diffs=4, count=stability.trim.contigs.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.fasta, count=stability.trim.contigs.good.unique.precluster.count_table)

#Remove rare clusters (i.e. cluster composed by only one contig, likley to be porduced by sequencing errors)
split.abund(fasta=stability.trim.contigs.good.unique.precluster.fasta, count=stability.trim.contigs.good.unique.precluster.count_table, cutoff=1, accnos=true)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.precluster.abund.count_table) 


#Scan for chimeras
chimera.vsearch(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.precluster.abund.count_table)
remove.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.fasta, accnos=stability.trim.contigs.good.unique.precluster.abund.denovo.vsearch.accnos, count=stability.trim.contigs.good.unique.precluster.abund.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table)

#Classify your reads according to silva
#RDP classifier using silva as the reference (similar results as RDP reference but some are classified to species)
classify.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table, reference=UNITEv6_sh_dynamic_s.fasta, taxonomy=UNITEv6_sh_dynamic.tax, cutoff=60, processors=44)

#Remove unwanted lineages (Bacteria-Animalia-Plantae-unclassified-Plantae_unclassified-unknown-Protista)
remove.lineage(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.count_table, taxonomy=stability.trim.contigs.good.unique.precluster.abund.pick.UNITEv6_sh_dynamic_s.wang.taxonomy, taxon=Bacteria-Animalia-Plantae-unclassified-Plantae_unclassified-unknown-Protista)
summary.seqs(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table)

#cluster using vsearch, default level is cutoff=0.03, since these are ITS we'll use 5% instead
# We cannot calculate the distance like with bacteria because the sequence lengths are different
cluster(fasta=stability.trim.contigs.good.unique.precluster.abund.pick.pick.fasta, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table, method=agc, cutoff=0.05)

#make OTU matrix 
make.shared(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table)

#classify each OTU, used the RDP classification 100% means all seqs in that OTU match at that classification level
classify.otu(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, count=stability.trim.contigs.good.unique.precluster.abund.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.precluster.abund.pick.UNITEv6_sh_dynamic_s.wang.pick.taxonomy)

#check number of sequences in each sample
count.groups(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared)

#alpha diversity (all provided measures; samples with less than 10000 reads are eliminated)
summary.single(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, subsample=10000)

#beta diversity
dist.shared(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, calc=braycurtis-jest-thetayc, subsample=10000)

#make a rarefied OTU table for heatmaps, indicator species, etc
sub.sample(shared=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.shared, size=10000)
 
summary.tax(taxonomy=current, count=current)

### Make a tree of OTU representative sequences. You need rename_oturep.py ###
# 1. Get representative sequences for each OTU
get.oturep(list=stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.list, fasta=stability.trim.contigs.good.unique.precluster.abund.pick.fasta, name=stability.trim.contigs.good.names, sorted=bin, method=abundance)

# 2. Rename sequence names into OTU names
python rename_oturep.py stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.0.05.rep.fast

# 3. # calculate distances between sequences in clean_repFasta.fasta
pairwise.seqs(fasta=clean_repFasta.fasta, cutoff=0.05, output=lt)

# 4. Build a tree for representative sequences (i.e., final OTUs)
clearcut(phylip=clean_repFasta.phylip.dist)
######

get.current()

### Files needed for downstream analysis in R package phyloseq: ###
# 1. OTU counts: 
# 2. OTU taxonomy: stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.0.05.cons.taxonomy
# 3. OTU tree: clean_repFasta.phylip.tre
# 4. Samples metadata with alpha-diversity indexes: stability.trim.contigs.good.unique.precluster.abund.pick.pick.agc.unique_list.groups.summary
######
