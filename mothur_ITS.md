# This procedure is based on mothur version 1.39.5
# First you need to create a stability file (a text file with "ID PAIR1 PAIR2")
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

#Pre-cluster. We used 1% difference to reduce computation time (4 differences in an amplicon ~300/500). You might want to run this step in bach mode since is slow.
pre.cluster(fasta=stability.trim.contigs.good.unique.fasta, diffs=4, count=stability.trim.contigs.good.count_table)
summary.seqs(fasta=current, count=current)

#Scan for chimeras
chimera.vsearch(fasta=stability.trim.contigs.good.good.unique.precluster.fasta, count=stability.trim.contigs.good.good.unique.precluster.count_table)
remove.seqs(fasta=current, accnos=current, count=current)
summary.seqs(fasta=current, count=current)

#Classify your reads according to silva
#RDP classifier using silva as the reference (similar results as RDP reference but some are classified to species)
classify.seqs(fasta=current, count=current, reference=UNITEv6_sh_dynamic.fasta, taxonomy=UNITEv6_sh_dynamic.tax, cutoff=60)



get.current()
