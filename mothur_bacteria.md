## This procedures has been used in SacaLaLengua 2 with mothur v.1.39.5
# March 2018

# First you need to create a stability file (a text file with "ID PAIR1 PAIR2")
for i in *R1_001.fastq; do echo `basename $i| sed s/_L001_R1_001.fastq//g | sed s/-/_/g` $i `basename $i| sed s/_R1_/_R2_/g`; done > stability.files


# Make contigs: create contigs from forward and reverse fastq files.
make.contigs(file=stability.files, processors=42)

# Summarize the quality of sequences
summary.seqs(fasta=stability.trim.contigs.fasta) # Stats #1 in stats_mothur_run6.txt

---------------
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	297	297	0	3	1
2.5%-tile:	1	440	440	0	4	373665
25%-tile:	1	460	460	0	4	3736648
Median: 	1	465	465	1	5	7473296
75%-tile:	1	465	465	9	6	11209944
97.5%-tile:	1	466	466	32	7	14572927
Maximum:	1	602	602	337	301	14946591
Mean:	1	462.394	462.394	6.49305	5.39431
# of Seqs:	14946591

Output File Names: 
stability.trim.contigs.summary
---------------


# Keep sequences that fulfill certain user defined criteria. 	
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=4,minlength=440, maxlength=466, processors=42)

# Returns only the unique sequences found in a fasta-formatted sequence file and a file that indicates those sequences that are identical to the reference sequence.
unique.seqs(fasta=stability.trim.contigs.good.fasta)

# Counts the number of sequences represented by the representative sequence in a name file.
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)

summary.seqs(count=stability.trim.contigs.good.count_table) # Stats #2 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	440	440	0	3	1
2.5%-tile:	1	442	442	0	4	246544
25%-tile:	1	460	460	0	4	2465439
Median: 	1	465	465	0	5	4930877
75%-tile:	1	465	465	1	6	7396315
97.5%-tile:	1	465	465	4	7	9615210
Maximum:	1	466	466	4	154	9861753
Mean:	1	460.825	460.825	0.591156	4.96562
# of unique seqs:	5291858
total # of seqs:	9861753

Output File Names: 
stability.trim.contigs.good.unique.summary


# SILVA SB prep: this was done in round 2 - run 1 - no need to repeat here; results linked to current directory
	# pcr.seqs(fasta=silva.nr_v132.align, keepdots=F, processors=40, oligos=pcrTest.oligos)
	# summary.seqs(fasta=silva.nr_v132.pcr.align)
	# screen.seqs(fasta=silva.nr_v132.pcr.align, start=4965, end=21977)
	# summary.seqs(fasta=silva.nr_v132.pcr.good.align)
	# rename.file(input=silva.nr_v132.pcr.good.align, new=silva.v34.fasta)
###############################

# Aligns a user-supplied fasta-formatted candidate sequence file to a user-supplied fasta-formatted template alignment (here, SILVA database).
	# find the closest template for each candidate using kmer (default) searching, blastn, or suffix tree searching
	# make a pairwise alignment between the candidate and de-gapped template sequences using the Needleman-Wunsch (default), Gotoh, or blastn algorithms
	# re-insert gaps to the candidate and template pairwise alignments using the NAST algorithm so that the candidate sequence alignment is compatible with the original template alignment.

align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v34.fasta, flip=true, processors=42)

summary.seqs(fasta=stability.trim.contigs.good.unique.align) # Stats #3 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	4964	4965	1	0	1	1
2.5%-tile:	4965	21977	403	0	4	132297
25%-tile:	4965	21977	422	0	4	1322965
Median: 	4965	21977	427	1	4	2645930
75%-tile:	4965	21977	427	2	5	3968894
97.5%-tile:	4965	21977	427	4	7	5159562
Maximum:	21977	21980	436	4	137	5291858
Mean:	4966.09	21976.6	422.48	1.06441	4.75694
# of Seqs:	5291858

Output File Names: 
stability.trim.contigs.good.unique.summary

screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=21977, end=4965, maxhomop=12, minlength=403, processors=42)

summary.seqs(fasta=current, count=current) # Stats #4 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	4964	21122	403	0	3	1
2.5%-tile:	4965	21977	405	0	4	241980
25%-tile:	4965	21977	422	0	4	2419792
Median: 	4965	21977	427	0	4	4839584
75%-tile:	4965	21977	427	1	5	7259376
97.5%-tile:	4965	21977	427	4	7	9437188
Maximum:	6234	21980	436	4	12	9679167
Mean:	4965.02	21977	423.197	0.589182	4.57573
# of unique seqs:	5179437
total # of seqs:	9679167

Output File Names: 
stability.trim.contigs.good.unique.good.summary

filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=42)

unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

summary.seqs(fasta=current, count=current) # Stats #5 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	1069	381	0	3	1
2.5%-tile:	1	1071	384	0	4	241980
25%-tile:	1	1071	401	0	4	2419792
Median: 	1	1071	406	0	4	4839584
75%-tile:	1	1071	406	1	5	7259376
97.5%-tile:	1	1071	406	4	7	9437188
Maximum:	5	1071	415	4	12	9679167
Mean:	1.00009	1071	402.201	0.589182	4.57298
# of unique seqs:	3742349
total # of seqs:	9679167

Output File Names: 
stability.trim.contigs.good.unique.good.filter.unique.summary

pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=4, processors=40)

#Remove rare clusters
split.abund(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, cutoff=1, accnos=true)

summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.count_table) # Stats #6 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	1071	381	0	3	1
2.5%-tile:	1	1071	384	0	4	193245
25%-tile:	1	1071	401	0	4	1932444
Median: 	1	1071	406	0	4	3864887
75%-tile:	1	1071	406	0	5	5797330
97.5%-tile:	1	1071	406	0	7	7536528
Maximum:	3	1071	411	4	10	7729772
Mean:	1.00001	1071	402.166	0.0134308	4.52566
# of unique seqs:	66151
total # of seqs:	7729772

Output File Names: 
stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.summary

chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.count_table, dereplicate=t, processors=42)

remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.accnos)

summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table) # Stats #7 in stats_mothur.txt
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	1071	381	0	3	1
2.5%-tile:	1	1071	384	0	4	187166
25%-tile:	1	1071	401	0	4	1871655
Median: 	1	1071	406	0	4	3743309
75%-tile:	1	1071	406	0	5	5614963
97.5%-tile:	1	1071	406	0	7	7299451
Maximum:	1	1071	411	4	10	7486616
Mean:	1	1071	402.125	0.0118596	4.52605
# of unique seqs:	36424
total # of seqs:	7486616

Output File Names: 
stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.summary

# assign sequences to taxonomy (default method "wang")
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, reference=trainset16_022016.pds.fasta, taxonomy=trainset16_022016.pds.tax, cutoff=80, processors=40)

# reads a taxonomy file and a taxon and generates a new file that contains only the sequences not containing that taxon.
remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#### mock community ####

# selects sequences from a specific group or set of groups (here select mock community samples)
get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta, groups=HM_782D_S1-HM_783D_S191)

# reads a fasta file and searches for errors in sequence compared to a reference file; assesses error rate
seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, reference=HM782_3_MOCK.fasta, aligned=F)

# grep "Overall error rate" *.logfile
# Overall error rate 6th round: 0.000177197
# first round: 5.05289e-05

# generate an uncorrected pairwise distances between aligned DNA sequences.
dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta, cutoff=0.03, processors=42)

# "opticlust" is the default clustering method for version v.1.39.0 of mothur.
	# opti: OTUs are assembled using metrics to determine the quality of clustering.
cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table)

# reads a list and group file or biom file and creates a .shared file as well as a rabund file for each group.

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, label=0.03)

# generate intra-sample rarefaction curves on multiple samples ("shared"): provide a way of comparing the richness observed in different samples.
rarefaction.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.shared)

######## in R ########
R
setwd("/nfs/users/project/sacalalengua/sacalall2/analysis_6")
data<-read.table(file="stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.groups.rarefaction", header=T)
pdf("rarefaction_plot_mock_run6.pdf")
plot(x=data$numsampled, y=data$X0.03.HM_783D_S191, xlab="Number of Tags Sampled",ylab="OTUs", type="l", col="red", font.lab=3, main="Rarefaction plot mock community run 6")
points(x=data$numsampled, y=data$X0.03.HM_782D_S1, type="l", col="black")
legend("bottomright", c("HM_783D (staggered)", "HM_782D (equimolar)"), fill=c("red", "black"))
dev.off()



### back to mothur ###

# check number of OTUS for mock community samples
# head *HM*rabund | cut -f2
# ==> HM_782D_S1.rabund <==
# 28
# ==> HM_783D_S191.rabund <==
# 44


#### again our samples ####

# removes sequences from a specific group or set of groups: here we remove the mock samples from the others
remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.taxonomy, groups=HM_782D_S1-HM_783D_S191)

# generate an uncorrected pairwise distances between aligned DNA sequences.
dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.fasta, cutoff=0.03, processors=40)

# "opticlust" is the default clustering method for version v.1.39.0 of mothur.
	# opti: OTUs are assembled using metrics to determine the quality of clustering.
cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table)

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, label=0.03)

# get a consensus taxonomy for an otu. 
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.pick.taxonomy, label=0.03)

# assign sequences to OTUs based on their taxonomy (taxonomy was obtained with "classify.seqs"
phylotype(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.pick.taxonomy)

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.pick.tx.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, label=1)

classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.pick.tx.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.denovo.vsearch.pick.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pds.wang.pick.pick.taxonomy, label=1)


