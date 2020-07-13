# Thomas J. Savitsky
# July 10, 2020

import argparse
import sys
import kpolymatroid

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='kpolyext.py',
        description='Generate single-element extensions of a k-polymatroid')
    parser.add_argument('--label', help='Just canonically label the '\
        'polymatroids.  Do not generate extensions.', action='store_true',
        default=False)
    parser.add_argument('--all', help='Generate all extensions without '\
        'doing isomorphism checks. Defaults to False.', action='store_true',
        default=False)
    parser.add_argument('-c', help='The maximum rank of the new element. '\
            'Defaults to 1.', action='store',
            dest='max_rank', default=1, type=int)
    parser.add_argument('--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('infile', help='name of input file',
        nargs='?', default=None, action='store')
    parser.add_argument('outfile', help='name of output file',
        nargs='?', default=None, action='store')

    results = parser.parse_args()

    if results.infile is None:
        inf = sys.stdin
    else:
        inf = open(results.infile, 'r')

    if results.outfile is None:
        outf = sys.stdout
    else:
        outf = open(results.outfile, 'w')

    for line in inf.readlines():
        kpm = kpolymatroid.KPolyMatroid()
        kpm.read_flats_str(line.rstrip())
        if results.label:
            canon = kpm.canonical_label()
            outf.write(str(canon)+'\n')
            continue
        kpm.extend(results.all, results.max_rank, outf)

    inf.close()
    outf.close()
