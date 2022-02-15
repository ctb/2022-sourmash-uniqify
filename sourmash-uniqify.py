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
from collections import namedtuple

import sourmash
from sourmash import load_file_as_signatures
from sourmash.logging import notify


ClusterMember = namedtuple('ClusterMember',
                           'founder_from, sigobj, cluster_num, member_type, cluster_type')


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
    p.add_argument('--merge-cluster-signatures', action='store_true',
                   help="Merge cluster signatures into founder.")
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
        print('new founder', founder_from)
        founder_info = ClusterMember(founder_from=founder_from,
                                     sigobj=founder,
                                     cluster_num=pass_n,
                                     member_type='founder',
                                     cluster_type='singleclust')

        cluster = []
        leftover = []
        for (sig_from, sig) in siglist:
            if args.max_containment:
                score = sig.max_containment(founder)
            else:
                score = sig.similarity(founder)

            if score >= args.threshold:
                cluster.append((sig_from, sig))
                member_info = ClusterMember(founder_from=sig_from,
                                            sigobj=sig,
                                            cluster_num=pass_n,
                                            member_type='member',
                                            cluster_type='multiclust')
                cluster_summary.append(member_info)
            else:
                leftover.append((sig_from, sig))

        if cluster:
            member_info = ClusterMember(founder_from=founder_from,
                                        sigobj=founder,
                                        cluster_num=pass_n,
                                        member_type='founder',
                                        cluster_type='multiclust')
            cluster_summary.append(member_info)

            notify(f'clustered {len(cluster)} signature(s) with founder sig {str(founder)[:30]}...')

            # output individual founder/cluster
            prefix = f'{args.prefix}cluster.{pass_n}'
            if not args.merge_cluster_signatures:
                with open(f'{prefix}.founder.sig', 'wt') as fp:
                    sourmash.save_signatures([founder], fp)
                with open(f'{prefix}.cluster.sig', 'wt') as fp:
                    cluster_sigs = [ x[1] for x in cluster ]
                    sourmash.save_signatures(cluster_sigs, fp)
                print(f'saved founder and {len(cluster)} signatures to {prefix}.*')
            # merge!
            else:
                mh = founder.minhash.to_mutable()
                for filename, ss in cluster:
                    mh += ss.minhash
                founder.minhash = mh
                print(f'merged {len(cluster) - 1} sketches into founder and saved to {prefix}.*')

        else:
            notify(f'founder sig {str(founder)[:30]}... is a singleton.')
            cluster_summary.append(founder_info)

            prefix = f'{args.prefix}cluster.{pass_n}'
            with open(f'{prefix}.founder.sig', 'wt') as fp:
                sourmash.save_signatures([founder], fp)
            print(f'saved singleton signature to {prefix}.*')

        siglist = leftover
        pass_n += 1

    # output summary spreadsheet
    headers = ['origin_path', 'name', 'filename', 'md5sum', 'cluster_num', 'member_type', 'cluster_type']
    csv_name = f'{args.prefix}summary.csv'

    with open(csv_name, 'w', newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=headers)
        w.writeheader()

        for member_info in cluster_summary:
            d = dict(member_info._asdict())
            del d['sigobj']
            del d['founder_from']
            d['origin_path'] = member_info.founder_from
            d['name'] = str(member_info.sigobj)
            d['filename'] = sig.filename
            d['md5sum'] = sig.md5sum()

            w.writerow(d)

    notify(f"wrote {len(cluster_summary)} entries to clustering summary at '{csv_name}'")


if __name__ == '__main__':
    sys.exit(main())
