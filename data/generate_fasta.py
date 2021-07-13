import pybedtools
from pybedtools import BedTool
import pandas as pd
import sys
import numpy as np


class GetFasta():
    def __init__(self, peak_bed, snp_bed, ref_fasta):
        self.peak_bed = peak_bed
        self.snp_bed = snp_bed
        self.ref_fasta = ref_fasta

    # read in the bed
    def read_bed(self):
        peak_bed_df = pd.read_csv(self.peak_bed, sep="\t", header=None)
        snp_bed_df = pd.read_csv(self.snp_bed, sep="\t", header=None)
        combined = pd.concat([peak_bed_df, snp_bed_df], axis=1)
        combined.columns = ['peak_chr', 'peak_start', 'peak_end',
                            'snp_chr', 'snp_start', 'snp_end',
                            'snp_id', 'wt', 'mutant']

        # check if the the snp is within peak
        if combined['peak_chr'].equals(combined['snp_chr']) and \
                all(combined['peak_start'] <= combined['snp_start'])and \
                all(combined['peak_end'] >= combined['snp_end']):
            pass
        else:
            raise ValueError('Given WT nucleic acid does not match with the one retrieve in fasta ')

        return peak_bed_df, snp_bed_df, combined

    # compute the location index of snp within the peak
    def compute_pos_index(self, combined):
        return list(combined['snp_start'] - combined['peak_start'])

    # generate wt and mutant fasta
    # TODO check why the snp info and the returned nucliec does not match - returned one pos further
    def get_fasta(self, snp_pos, combined):
        peak_bed_bt = BedTool(self.peak_bed)
        ref_bt = BedTool(self.ref_fasta)
        peak_fa = peak_bed_bt.sequence(fi=ref_bt)
        peak_fa = open(peak_fa.seqfn).read()
        peak_fa_list = self._process_seq(peak_fa, combined)
        mutant_peak_fa_list = self._generate_mutant_fa(combined, peak_fa_list)

        return np.array(peak_fa_list), np.array(mutant_peak_fa_list)

    def _process_seq(self, peak_fa, combined):
        n_peak = combined.shape[0]
        index_list = [i*2+1 for i in range(n_peak)]
        peak = peak_fa.split("\n")
        peak_list = [peak[i].upper() for i in index_list] # TODO: change to upper
        print(peak_list)
        return peak_list

    def _generate_mutant_fa(self, combined, peak_fa_list):

        snp_index = combined['snp_start'] - combined['peak_start'] - 1
        snp_index_list = []
        for i in snp_index:
            snp_index_list.append(i)

        mutant_peak_fa_list = []
        for i in range(len(peak_fa_list)):
            fa = peak_fa_list[i]
            snp_index = snp_index_list[i]

            if fa[snp_index] != combined.iloc[i, 7]:
                raise ValueError('Given WT nucleic acid does not match with the one retrieve in fasta')
            else:
                fa = list(fa)
                fa[snp_index] = combined.iloc[i, 8].upper()
                mutant_fa = ''.join(fa)
                mutant_peak_fa_list.append(mutant_fa)

        return mutant_peak_fa_list