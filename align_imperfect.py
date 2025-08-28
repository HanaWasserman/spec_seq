import pandas as pd
import sys
import os
import Levenshtein
import swalign
from tqdm import tqdm
tqdm.pandas()
pd.options.mode.chained_assignment = None

match = 1
mismatch = -1
gap_penalty = -1
sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(match, mismatch), gap_penalty)
def align_subseq(df_row, ref_subseq_lst, meth):
    subseq = df_row['subseq']
    aligns = []
    for ref_subseq in ref_subseq_lst['ref_subseq']:
        if meth=='levenshtein':
            score = Levenshtein.distance(ref_subseq, subseq)
        elif meth=='smith_waterman':
            alignment = sw.align(ref_subseq, subseq)
        aligns.append((ref_subseq, alignment.score, alignment.r_pos))
        if alignment.score > 57:
            break
    sorted_aligns = sorted(aligns, key=lambda x: x[1], reverse=True)
    try:
        nearest_lev_dist = sorted_aligns[0][1]-sorted_aligns[1][1]
    except:
        nearest_lev_dist = -1
    return pd.Series(sorted_aligns[0] + (nearest_lev_dist,), index=['ref_subseq','dist_score','align_pos','nearest_dist_score'])

def rescue_seqs(lib, rescue_sqs_path):
    reads = pd.read_csv(rescue_sqs_path, names=['l_spacer','roi','r_spacer','r_barcode','subseq','count'])
    reads[['ref_subseq', 'dist_score', 'align_pos', 'nearest_dist_score']] = reads.progress_apply(align_subseq, ref_subseq_lst=lib, meth='smith_waterman', axis=1)
    reads_merge = pd.merge(reads, lib[['id', 'ref_subseq']], on='ref_subseq')
    reads_merge.to_csv(f'{os.path.splitext(rescue_sqs_path)[0]}_mapped.csv', index=False, header=False)

if __name__== "__main__":
    libr_path = f'{os.getcwd()}/ETS_library_30_18_1_24_1_18_30_11152024.txt'
    libr = pd.read_table(libr_path, names=['id', 'seq'])
    libr[['left_primer', 'l_barcode', 'sp1', 'roi', 'sp2', 'r_barcode', 'right_primer']] = libr['seq'].str.split('-',expand=True)
    libr['ref_subseq'] = libr['l_barcode'].str[2:] + libr['sp1'] + libr['roi'] + libr['sp2'] + libr['r_barcode']
    rescue_seqs_path = sys.argv[1]
    rescue_seqs(libr, rescue_seqs_path)


