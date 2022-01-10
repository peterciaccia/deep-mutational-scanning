"""
Created by Peter Ciaccia
"""

# built-ins
import os
import math

# dependencies
import numpy as np
import pylab as plt
import matplotlib.patches as patches
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import pairwise2
from dotenv import load_dotenv

# loads environment variables
load_dotenv()


class Variant:

    def __init__(self, record, alignment):
        self.record = record
        self.alignment = alignment

        print(self.record.id)


class Analyzer:

    def __init__(self, filename, outdir, reanalyze=False):
        self.filename = filename
        self.outdir = outdir
        self.reference_records = [record for record in
                                  SeqIO.parse(os.getenv("REFERENCE_SEQ_PATH"), 'fasta')]
        self.concat_reference_seq = ''.join([str(record.seq) for record in self.reference_records])

    @property
    def parser(self):
        return SeqIO.parse(self.filename, 'fastq')

    def _make_read_df(self, size=None):

        records = []
        i = 0
        for fastq_record in self.parser:
            if size is not None:
                i += 1
                if i > size:
                    break
            records.append([fastq_record.id, str(fastq_record.seq)])
        df = pd.DataFrame(records, columns=['id', 'seq'])
        return df

    def plot_fastq_gc_content(self, limit=None):

        # returns normal distriubtion around mean
        def norm_pdf(x, mean, sd):
            variance = float(sd) ** 2
            denominator = (2 * math.pi * variance) ** .5
            numerator = math.exp(-(float(x)-float(mean))**2/(2*variance))
            return numerator/denominator

        # configures plots
        f, ax = plt.subplots(figsize=(12, 5))
        df = self._make_read_df(size=limit)
        gc = df.seq.apply(lambda x: GC(x))
        gc.hist(ax=ax, bins=50, color='black', grid=False, histtype='step', lw=2)
        ax.set_xlim((0, 100))
        x_ = np.arange(1, 100, 0.1)
        f = [norm_pdf(i, gc.mean(), gc.std()) for i in x_]
        ax2 = ax.twinx()
        ax2.plot(x_, f)
        ax2.set_ylim(0, max(f))
        ax.set_title('GC content', size=15)

        # saves figure
        gc_plot_path = os.path.join(self.outdir, "gc_FWD_and_REV.png")
        plt.savefig(gc_plot_path, transparent=True)
        plt.clf()
        return

    def make_quality_score_summary(self):

        scores = [record.letter_annotations['phred_quality'] for record in self.parser]
        df = pd.DataFrame(scores)
        length = len(df.T) + 1

        f, ax = plt.subplots(figsize=(12, 5))
        rect = patches.Rectangle((0, 0), length, 20, linewidth=0, facecolor='r', alpha=.4)
        ax.add_patch(rect)
        rect = patches.Rectangle((0, 20), length, 8, linewidth=0, facecolor='yellow', alpha=.4)
        ax.add_patch(rect)
        rect = patches.Rectangle((0, 28), length, 12, linewidth=0, facecolor='g', alpha=.4)
        ax.add_patch(rect)
        df.mean().plot(ax=ax, c='black')
        boxprops = dict(linestyle='-', linewidth=1, color='black')
        df.plot(kind='box', ax=ax, grid=False, showfliers=False,
                color=dict(boxes='black', whiskers='black'))
        ax.set_xticks(np.arange(0, length, 5))
        ax.set_xticklabels(np.arange(0, length, 5))
        ax.set_xlabel('position(bp)')
        ax.set_xlim((0, length))
        ax.set_ylim((0, 40))
        ax.set_title('per base sequence quality')

        quality_score_path = os.path.join(self.outdir, "quality_hist_FWD_and_REV.png")
        plt.savefig(quality_score_path, transparent=True)
        plt.clf()
        return

    def sort_records(self):

        def get_max_alignment():
            alignments = pairwise2.align.globalxx(sequenceA=reference.seq, sequenceB=record.seq)
            maximum = max([alignment for alignment in alignments], key=lambda x: x.score)
            return maximum

        gem_records = {}
        zpm_records = {}
        for record in self.parser:
            score_dict = {}
            for reference in self.reference_records:
                score_dict[reference.id] = get_max_alignment()

            filtered_alignment = max(score_dict.items(), key=lambda x: x[1].score)
            if filtered_alignment[0] == 'ZPM_amplicon':
                zpm_records[record.id] = Variant(record, filtered_alignment[1])
            elif filtered_alignment[0] == 'GEM_amplicon':
                gem_records[record.id] = Variant(record, filtered_alignment[1])
            # sorted_alignments = sorted([alignment for alignment in alignments],
            #                            key=lambda x: x.score)


def get_data_paths(get_paths=False):
    data_dir = os.getenv("DATADIR")
    file_list = sorted(
        [file
         for file in os.listdir(data_dir)
         if file.endswith(".fastq")],
        key=lambda x: int(x.split('/')[-1].split('_')[0])
    )

    if get_paths:
        path_list = [os.path.join(data_dir, file)
                     for file in file_list]
        return path_list
    else:
        return file_list


def get_output_paths(outfile=None):
    if not os.path.isdir(os.getenv("OUTDIR")):
        os.mkdir(os.getenv("OUTDIR"))
    outpath_list = []
    for file in get_data_paths():
        filename_as_dirname = os.path.splitext(file)[0]
        output_path = os.path.join(os.getenv("OUTDIR"), filename_as_dirname)
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        if outfile is not None:
            output_path = os.path.join(output_path, outfile)
        outpath_list.append(output_path)
    return outpath_list


def run():

    file_list = get_data_paths()

    # makes output directories if not exists


    # runs analysis for each fastq file
    # analyses = {}
    # for path, outdir_ in zip(path_list, outdir_list):
    #     filename = os.path.basename(path)
    #     analyses[filename] = Analyzer(path, outdir_, reanalyze=True)
    #
    # for k, v in analyses.items():
    #     v.sort_records()
    #     v.make_quality_score_summary()
    #     v.plot_fastq_gc_content()

    for k, v in analyses.items():
        v.sort_records()
        v.make_quality_score_summary()
        v.plot_fastq_gc_content()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()
