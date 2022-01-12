"""
Created by Peter Ciaccia

Acknowledgments
Credit to dmnfarrell (ORDID 0000-0003-3020-7945) for implementations of gc_content, quality score visualizations from
which Analyzer.plot_fastq_gc_content and Analyzer.make_quality_score_summary were derived:
https://dmnfarrell.github.io/python/fastq-quality-python
"""

# built-ins
import os
import math
import subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from pathlib import Path
import shutil

# dependencies
import numpy as np
import pylab as plt
import matplotlib.patches as patches
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils import GC
from dotenv import load_dotenv

# loads environment variables
load_dotenv()


class Variant:

    def __init__(self):
        pass


class Analyzer:

    def __init__(self, analysis_dir, reanalyze=False):
        self.analysis_dir = analysis_dir
        self.merged_reads_path = os.path.join(self.analysis_dir, 'merged_reads.fastq')
        self.alignment_chunks_dir = os.path.join(self.analysis_dir, 'alignment_chunks')
        self.reference_records = [record for record in
                                  SeqIO.parse(os.getenv("REFERENCE_SEQ_PATH"), 'fasta')]
        self.concat_reference_seq = ''.join([str(record.seq) for record in self.reference_records])

    def get_abs_path(self, filename):
        return os.path.join(self.analysis_dir, filename)

    def get_unaligned_chunk_path(self, filename):
        return os.path.join(self.alignment_chunks_dir, 'input', filename)

    def get_aligned_chunk_path(self, filename):
        return os.path.join(self.alignment_chunks_dir, 'output', filename)

    def fastq_parser(self, filename='merged_reads.fastq'):
        file_path = self.get_abs_path(filename)
        return SeqIO.parse(file_path, 'fastq')

    def _make_read_df(self, size=None):

        records = []
        i = 0
        for fastq_record in self.fastq_parser():
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
        gc_plot_path = os.path.join(self.analysis_dir, "gc_FWD_and_REV.png")
        plt.savefig(gc_plot_path, transparent=True)
        plt.close()
        return

    def make_quality_score_summary(self):

        scores = [record.letter_annotations['phred_quality'] for record in self.fastq_parser()]
        df = pd.DataFrame(scores)
        length = len(df.T) + 1

        f, ax = plt.subplots(figsize=(24, 5))
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
        ax.set_xticks(np.arange(0, length, 25))
        ax.set_xticklabels(np.arange(0, length, 25))
        ax.set_xlabel('position(bp)')
        ax.set_xlim((0, length))
        ax.set_ylabel('quality score')
        ax.set_ylim((0, 40))
        ax.set_title('per base sequence quality')

        quality_score_path = os.path.join(self.analysis_dir, "quality_hist_FWD_and_REV.png")
        plt.savefig(quality_score_path, transparent=True)
        plt.close()
        return

    @staticmethod
    def _wrap_fasta(sequence, fasta_line_length=60):

        return '\n'.join([sequence[i:i + fasta_line_length]
                          for i in range(0, len(sequence), fasta_line_length)]
                         )

    def _get_reference_fasta(self):

        with open(os.getenv("REFERENCE_CONCAT_AMPLICON_PATH"), 'r') as in_handle:
            record = SeqIO.read(in_handle, "fasta")
            (name, seq) = (record.id, str(record.seq))
            return f">{name}\n{self._wrap_fasta(seq)}\n"

    def chunk_reads(self, reanalyze=True, chunk_size=180):

        unaligned_fasta_dir = os.path.join(self.alignment_chunks_dir, 'input')
        aligned_fasta_dir = os.path.join(self.alignment_chunks_dir, 'output')
        if reanalyze:
            try:
                shutil.rmtree(self.alignment_chunks_dir)
                print('Directories emptied')
            except FileNotFoundError:
                pass
        if not os.path.isdir(self.alignment_chunks_dir):
            os.makedirs(unaligned_fasta_dir)
            os.makedirs(aligned_fasta_dir)

        count = 0
        chunk = 0
        with open(self.merged_reads_path, 'r') as in_handle:
            # parser = SeqIO.parse(in_handle, 'fastq')
            chunk_str = self._get_reference_fasta()
            for name, seq, qual in FastqGeneralIterator(in_handle):
                count += 1
                seq = self._wrap_fasta(seq)
                chunk_str += f'>{name}\n{seq}\n'
                if count == chunk_size:
                    chunk_name = f'chunk_{str(chunk).zfill(8)}'
                    unaligned_chunk_path = os.path.join(unaligned_fasta_dir,
                                                        f'{chunk_name}.fa')
                    with open(unaligned_chunk_path, 'w') as f:
                        f.write(chunk_str)
                    count = 0
                    chunk += 1
                    chunk_str = self._get_reference_fasta()

    def align_chunks(self):

        def worker(unaligned_path, aligned_path):
            print(aligned_path)
            out_handle = open(aligned_path, 'w')
            process = subprocess.Popen(mafft_args,
                                       stdout=out_handle,
                                       stderr=subprocess.PIPE
                                       )
            stuff = process.communicate()
            out_handle.close()

            out = stuff[0]
            err = stuff[1]
            print(f"aligned chunk {os.path.basename(unaligned_path)}")

        # pool = mp.Pool()
        pool = ThreadPool()

        for file in sorted(os.listdir(os.path.join(self.alignment_chunks_dir, 'input'))):
            unaligned_chunk_path = self.get_unaligned_chunk_path(file)
            aligned_chunk_path = self.get_aligned_chunk_path(file)
            if os.path.isfile(aligned_chunk_path):
                os.remove(aligned_chunk_path)
            Path(aligned_chunk_path).touch()

            mafft_args = [
                'mafft',
                '--adjustdirection',
                # '--localpair', '--lop', '-10.00',
                '--quiet',
                '--clustalout',
                unaligned_chunk_path,
                # '>',
                # aligned_chunk_path
            ]
            pool.apply_async(worker, (unaligned_chunk_path, aligned_chunk_path))

        pool.close()
        pool.join()


def get_raw_data_paths(get_paths=False):
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


def get_data_analysis_paths(filename=None):
    if not os.path.isdir(os.getenv("OUTDIR")):
        os.mkdir(os.getenv("OUTDIR"))
    outpath_list = []
    for file in get_raw_data_paths():
        filename_as_dirname = os.path.splitext(file)[0]
        output_path = os.path.join(os.getenv("OUTDIR"), filename_as_dirname)
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        if filename is not None:
            output_path = os.path.join(output_path, filename)
        outpath_list.append(output_path)
    return outpath_list


def generate_statistics():

    input_path_list = get_data_analysis_paths(filename='merged_reads.fastq')
    outdir_list = get_data_analysis_paths()

    analyses = {}
    for outdir_ in outdir_list:
        analyses[outdir_] = Analyzer(outdir_)

    for k, v in analyses.items():
        v.make_quality_score_summary()
        v.plot_fastq_gc_content()


def merge_reads(check_first=True):

    for input_path, output_path, outu_path in zip(get_raw_data_paths(get_paths=True),
                                                  get_data_analysis_paths(filename='merged_reads.fastq'),
                                                  get_data_analysis_paths(filename='unmerged_reads.fastq')):
        if check_first:
            if os.path.isfile(output_path):
                continue

        bbmerge_args = [
            'bbmerge-auto.sh',
            f'in={input_path}',
            f'out={output_path}',
            f'outu={outu_path}'
        ]
        additional_args = [
            'rem',
            'k=62',
            'extend2=50',
            'ecct'
        ]
        bbmerge_args.extend(additional_args)
        process = subprocess.Popen(bbmerge_args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()


def get_variants():

    for analysis_path in get_data_analysis_paths():
        align_path = os.path.join(os.getenv('PROJECTPATH'), 'align_merged_reads.sh')
        process = subprocess.Popen([align_path, analysis_path])
        stdout, stderr = process.communicate()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # merge_reads()
    # generate_statistics()
    get_variants()
