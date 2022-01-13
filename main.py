"""
Created by Peter Ciaccia

Acknowledgments
Credit to dmnfarrell (ORDID 0000-0003-3020-7945) for implementations of gc_content, quality score visualizations from
which Analyzer.plot_fastq_gc_content and Analyzer.make_quality_score_summary were derived:
https://dmnfarrell.github.io/python/fastq-quality-python
"""

# built-ins
import logging
import os
import math
import subprocess
import multiprocessing as mp
from functools import partial
from multiprocessing.pool import ThreadPool
from pathlib import Path
import shutil
from io import StringIO
import re
from collections import Counter

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

# set up logging formatting
formatter = logging.Formatter(fmt="%(asctime)s %(module)s::%(lineno)d %(levelname)s %(message)s")


def setup_logger(log_path, level=logging.INFO):
    name = os.path.basename(log_path)
    log_file = os.path.join(log_path, 'log.out')

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.addHandler(stream_handler)

    return logger


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

        with open(os.getenv("REFERENCE_SEQ"), 'r') as in_handle:
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
    # makes directories if they do not exist
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


def call_variant_freqs():

    align_merged_reads_path = os.path.join(os.getenv('PROJECTPATH'), 'align_merged_reads.sh')
    for analysis_path in get_data_analysis_paths():
        variants_path = os.path.join(analysis_path, 'var.vcf')
        process = subprocess.Popen([align_merged_reads_path, analysis_path],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        with open(variants_path, 'r') as f:
            """bcftools mpileup uses "##" for comment lines, which read_csv does not support.
               Pre-processes var.vcf to replace "##" with "%" """
            file_obj = StringIO(re.sub(r"^#{2,}", "%", f.read(), flags=re.M))
        df = pd.read_csv(file_obj, comment="%", sep="\t")
        # TODO: get frequencies of each position


def make_variant_sequence(mismatch_list_, position, read_, reference, verbose=False, logger=None):

    if logger is None:
        raise NotImplementedError

    nums = [int(x) for x in mismatch_list_[0:len(mismatch_list_):2]]
    nucs = mismatch_list_[1:len(mismatch_list_):2]

    # discards reads with indels
    if "^" in ''.join(nucs):
        logger.debug(f"indel in {mismatch_list_}")
        return None
    # returns reference if no mutations identified
    if len(mismatch_list_) == 1:
        return reference
    variant_seq = reference[:position-1] + read_ + reference[position-1+len(read_):]
    if len(variant_seq) != len(reference):
        print(len(variant_seq))

    if verbose:
        print(nums, nucs, position)
        print(f"read:\t{read_[:nums[0]-1]} "
              f"{read_[nums[0]-1:nums[0]]} "
              f"{read_[nums[0]:nums[0]+4]}")
        print(f"ref:\t{reference.seq[position:position+nums[0]-1]} "
              f"{reference.seq[position+nums[0]-1:position+nums[0]]} "
              f"{reference.seq[position+nums[0]:position+nums[0]+4]}")
        print(variant_seq.seq[position:position+30])
        print(reference.seq[position:position+30])
        print('')

    # print(len(reference), len(variant_seq))
    if len(reference) != len(variant_seq):
        logger.warning("Variant length is not the same as reference length:\n"
                       f"Mismatches:\t{mismatch_list_}\n"
                       f"Variant:\t{variant_seq.seq}\n"
                       f"Reference:\t{reference.seq}"
                       )
        # print(nums, nucs, position)
        # print(variant_seq.translate().seq)
    return variant_seq


def worker(reference_nuc, analysis_path):

    def map_aa_pos(i):
        """maps pseudoORF aa indices to aa indices in biosensor"""
        if i <= 373:
            return i + 139
        else:
            return i - 373

    logger = setup_logger(analysis_path)
    sam_path = os.path.join(analysis_path, "aligned_reads_sorted.sam")
    header_names = [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
        "mismatch_num",
        "mismatch"
    ]
    dtypes = [
        str,
        int,
        str,
        int,
        int,
        str,
        str,
        str,
        int,
        str,
        str,
        str,
        str
    ]
    dtype_dict = dict(zip(header_names, dtypes))
    df = pd.read_csv(sam_path, comment="@", sep="\t", usecols=range(13),
                     names=header_names, dtype=dtype_dict)

    logger.info(f"Analyzing {analysis_path}")
    mutations = []
    for name, pos, seq, map_col in zip(df['QNAME'],
                                       df['POS'],
                                       df['SEQ'],
                                       df['mismatch']):

        if str(map_col).startswith('MD'):
            sequence = map_col.split(":")[2]
            mismatch_list = re.split(r"(\D+)", sequence)
            variant_nuc = make_variant_sequence(mismatch_list, pos, seq,
                                                reference_nuc, logger=logger)
            if variant_nuc is None:
                continue
            variant_pept = variant_nuc.translate()

            reference_pept = reference_nuc.translate()
            aa_pos = map(map_aa_pos, range(len(reference_pept.seq)))
            for i, ref_aa, var_aa in zip(aa_pos,
                                         list(reference_pept.seq),
                                         list(variant_pept.seq)):
                if ref_aa != var_aa:
                    mutations.append(''.join([ref_aa, str(i+1), var_aa]))

    counts = Counter(mutations).most_common()

    with open(os.path.join(analysis_path, "substitution_counts.txt"), "w") as f:
        for k, v in counts:
            f.write(f"{k}\t{v}\n")


def extract_variants_from_sam():

    with open(os.getenv("REFERENCE_SEQ"), "r") as f:
        reference_nuc_ = SeqIO.read(f, format="fasta")

    pool = mp.Pool()
    analysis_paths_ = [path_ for path_ in get_data_analysis_paths()]
    func = partial(worker, reference_nuc_)
    pool.map(func, analysis_paths_)


def run(rerun_all=False, rerun_merge=False, rerun_stats=False,
        rerun_variant_call=False, rerun_extract_variants=False):

    if rerun_all:
        merge_reads(check_first=False)
        generate_statistics()
        call_variant_freqs()
        extract_variants_from_sam()
    else:
        if rerun_merge:
            merge_reads(check_first=False)
        if rerun_stats:
            generate_statistics()
        if rerun_variant_call:
            call_variant_freqs()
        if rerun_extract_variants:
            extract_variants_from_sam()


if __name__ == '__main__':
    run(
        rerun_extract_variants=True
    )
