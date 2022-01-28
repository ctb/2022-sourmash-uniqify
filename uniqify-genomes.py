#! /usr/bin/env python
"""
Peform an iterative, greedy clustering of a pile of sourmash signatures.

This is a practical alternative to a more principled clustering; see
https://github.com/ctb/2017-sourmash-cluster

Author: C. Titus Brown, github.com/ctb/, titus@idyll.org.

This code is under CC0.
"""
import sys
import argparse
import random
import csv
import screed
import os
import shutil
import gzip

import sourmash
from sourmash import load_file_as_signatures
from sourmash.logging import notify


def main():
    p = argparse.ArgumentParser()
    p.add_argument('signature_sources', nargs='+',
                   help='signature files, directories, and sourmash databases')
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--moltype', default='DNA')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--threshold', type=float, default=0.2)
    p.add_argument('--prefix', default='clust.',
                   help='output filename prefix (can be a path)')
    p.add_argument('--max-containment', action='store_true',
                   help="Use max containment instead of similarity to cluster.")
    p.add_argument('--merge-files', action='store_true',
                   help="Merge fasta files, starting with founder.")
    args = p.parse_args()

    # @CTB ovviously
    base_mh = sourmash.MinHash(n=0, ksize=31, scaled=1000)

    siglist = []
    for filename in args.signature_sources:
        notify(f'loading sequences from {filename}')
        with screed.open(filename) as screed_iter:
            mh = base_mh.copy_and_clear()
            for record in screed_iter:
                mh.add_sequence(record.sequence, force=True)

            siglist.append((filename, mh))

    notify(f'loaded {len(siglist)} files total.')

    notify(f'setting random number seed to {args.seed} and shuffling seqs')
    random.seed(args.seed)
    random.shuffle(siglist)

    cluster_summary = []

    pass_n = 0
    while len(siglist):
        notify(f'starting pass {pass_n+1}')
        (founder_from, founder) = siglist.pop()
        cluster_summary.append((founder_from, founder, pass_n, 'founder'))

        cluster = []
        leftover = []
        for (sig_from, sketch) in siglist:
            if args.max_containment:
                score = sketch.max_containment(founder)
            else:
                score = sketch.similarity(founder)

            if score >= args.threshold:
                cluster.append((sig_from, sketch))
                cluster_summary.append((sig_from, sketch, pass_n, 'member'))
            else:
                leftover.append((sig_from, sketch))

        if cluster:
            notify(f'clustered {len(cluster)} genomes with founder {founder_from}')

            # output individual founder/cluster
            if args.merge_files:
                with gzip.open(f"{args.prefix}cluster.{pass_n}.fa.gz", "wt") as fp:
                    for record in screed.open(founder_from):
                        fp.write(f">{record.name}\n{record.sequence}\n")
                    for ff, _ in cluster:
                        for record in screed.open(ff):
                            fp.write(f">{record.name}\n{record.sequence}\n")
                pass
            else:
                prefix = f'{args.prefix}cluster.{pass_n}'
                try:
                    os.mkdir(prefix)
                except FileExistsError:
                    print(f"note - '{prefix}' already exists.")

                # shutil copy instead
                shutil.copy(sig_from, prefix)
                for (ff, _) in cluster:
                    shutil.copy(ff, prefix)
                print(f'saved founder and {len(cluster)} files to {prefix}')

        else:
            notify(f'founder {founder_from} is a singleton.')

            if args.merge_files:
                with gzip.open(f"{args.prefix}cluster.{pass_n}.fa.gz", "wt") as fp:
                    for record in screed.open(founder_from):
                        fp.write(f">{record.name}\n{record.sequence}\n")
            else:
                prefix = f'{args.prefix}cluster.{pass_n}'
                try:
                    os.mkdir(prefix)
                except FileExistsError:
                    print(f"note - '{prefix}' already exists.")

                shutil.copy(founder_from, prefix)
                print(f'copied singleton file to {prefix}')

        siglist = leftover
        pass_n += 1

    if 0:
        # output summary spreadsheet
        headers = ['origin_path', 'name', 'filename', 'md5sum', 'cluster', 'member_type']
        csv_name = f'{args.prefix}summary.csv'

        with open(csv_name, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(headers)

            for (origin_path, sig, cluster_n, member_type) in cluster_summary:
                name = str(sig)
                filename = sig.filename
                md5sum = sig.md5sum()

                w.writerow([origin_path, name, filename, md5sum, cluster_n, member_type])

        notify(f"wrote {len(cluster_summary)} entries to clustering summary at '{csv_name}'")


if __name__ == '__main__':
    sys.exit(main())
