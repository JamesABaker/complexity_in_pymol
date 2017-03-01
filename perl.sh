#Clean previous run results
rm inputseq.fasta
rm TMsegments.txt
rm inputseq.fasta

touch inputseq.fasta
touch TMsegments.txt
touch inputseq.fasta

#run TMSOC
perl TMSOC.pl inputseq.fasta TMsegments.txt >> TMSOCoutput.txt
