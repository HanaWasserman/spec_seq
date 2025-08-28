# spec_seq
## Scripts for processing spec-seq library sequencing data
1. Pre-process fastq file to get reverse complement and extract sequence to ```<input>_rc_seqs.txt``` output file:
```
bash preprocess_specseq_library.bs <input>.fastq
```
2. Do preliminary alignment of reads to primer sequence. This is helpful for some basic QC. I split the read file and then batch alignment jobs and concatenate them when done:
```
split -l 200 -d <input>_rc_seqs.txt <input>_chunk_
sbatch run_align_primers_job.slurm
...
cat <input>_chunk_* > <input>_rc_primers_aligned.csv
```
3. 
