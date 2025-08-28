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
3. Separate out reads that have a perfect match to reference library design and imperfect matches. See ```ETS_library_30_18_1_24_1_18_30_11152024.txt``` for example library design file. Perfect matches determined by ROI+Barcode1 match (ROI is the binding site). Perfect matches are aggregated for read coverage:
```
python delineate_matches.py <input>_rc_primers_aligned.csv <library_design>.txt
```
**Optional**
4. Rescue imperfect reads with Smith-Waterman alignment. Again, split file and batch alignment jobs. Customize ```align_imperfect.py``` file for desired alignment scoring parameters:
```
split -l 200 -d <input>_imperfect_matches.csv <input>_imperfect_matches_chunk_
sbatch run_align_imperfect_job.slurm 
```
