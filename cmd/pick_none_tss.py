import random

from utils.feat_utils import bedIterator
import argparse


def randomPick(bed, output, genLen):
    readLen = 2000
    bedfileList = [bed]
    chain = []
    for bed in bedIterator(bedfileList):
        chain.extend(list(range(bed[1]-readLen, bed[2]+readLen)))
    chain = set(chain)
    i = 0
    with open(output, 'w+') as f:
        while i < genLen:
            start = random.randint(0, 249133423)
            if start in chain:
                continue
            else:
                i += 1
                f.write('1\t%d\t%d\n' % (start, start+readLen))
    print('write to %s' % output)


parser = argparse.ArgumentParser(description="Pick none tss zone")
parser.add_argument("bed", help="bed file")
parser.add_argument('-o', '--output', help='output file',
                    dest='output', default='build/none_tss.bed')
parser.add_argument('-l', '--gen-len', help='generate read length',
                    type=int, dest='genLen', default=6500)
args = parser.parse_args()
randomPick(args.bed, args.output, args.genLen)
