import swalign
import pandas as pd
import sys
import os

# Script for initial, easy alignment of spec-seq reads to primer sequences (useful for rough QC)

# Customize Smith-Waterman alignment parameters as desired
match = 1
mismatch = -1
gap_penalty = -4
gap_extension_penalty = 0
sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(match, mismatch), gap_penalty, gap_extension_penalty)

# Replace with desired primer sequence 
primers = 'TACGCATACACGCACACGCACATACACATA'

aligned_seqs = []
def align_chunk(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            seq = line.strip()
            alignment = sw.align(seq, primers)
            aligned_seqs.append((seq, alignment.score, alignment.r_pos, alignment.matches))
    aligned_seqs_df = pd.DataFrame(aligned_seqs, columns=['seq','score','r_pos','matches'])
    aligned_seqs_df.to_csv(f'{os.path.splitext(file_path)[0]}.csv', index=False, header=False)

if __name__== "__main__":
    file_path = sys.argv[1]
    align_chunk(file_path)
