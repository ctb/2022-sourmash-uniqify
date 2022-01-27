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
    p.add_argument('--prefix', default='cluster',
                   help='output filename prefix (can include directories)')
    p.add_argument('--max-containment', action='store_true',
                   help="Use max containment instead of similarity to cluster.")
    args = p.parse_args()

    siglist = []
    for filename in args.signature_sources:
        notify(f'loading from {filename}')
        m = 0
        for sig in load_file_as_signatures(filename,
                                           select_moltype=args.moltype,
                                           ksize=args.ksize):
            m += 1
            siglist.append((filename, sig))
        notify(f'...got {m} signatures.')

    notify(f'loaded {len(siglist)} signatures total.')

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
        for (sig_from, sig) in siglist:
            if args.max_containment:
                score = sig.max_containment(founder)
            else:
                score = sig.similarity(founder)

            if score >= args.threshold:
                cluster.append((sig_from, sig))
                cluster_summary.append((sig_from, sig, pass_n, 'member'))
            else:
                leftover.append((sig_from, sig))

        if cluster:
            notify(f'clustered {len(cluster)} signature(s) with founder sig {str(founder)[:30]}...')

            prefix = f'{args.prefix}.cluster.{pass_n}'
            with open(f'{prefix}.founder.sig', 'wt') as fp:
                sourmash.save_signatures([founder], fp)
            with open(f'{prefix}.cluster.sig', 'wt') as fp:
                cluster_sigs = [ x[1] for x in cluster ]
                sourmash.save_signatures(cluster_sigs, fp)

            print(f'saved founder and {len(cluster)} signatures to {prefix}.*')
        else:
            notify(f'founder sig {str(founder)[:30]}... is a singleton.')

            prefix = f'{args.prefix}.cluster.{pass_n}'
            with open(f'{prefix}.founder.sig', 'wt') as fp:
                sourmash.save_signatures([founder], fp)
            print(f'saved singleton signature to {prefix}.*')

        siglist = leftover
        pass_n += 1

    # output summary spreadsheet
    headers = ['origin_path', 'name', 'filename', 'md5sum', 'cluster', 'member_type']
    csv_name = f'{args.prefix}.summary.csv'

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
