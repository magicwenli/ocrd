'''
Database_web_link = http://bioinf.xmu.edu.cn/PaGenBase/index.jsp
'''
import os
import glob
from argparse import ArgumentParser

parser = ArgumentParser(description='Extract gene infos from PaGenBase file.')
parser.add_argument('path', help='Path of files that need to be extracted.')
parser.add_argument('-o', '--output', dest='output',
                    help='Output dir, default to be *path*.')
parser.add_argument('-g', '--GRCh37', dest='grch37',
                    help='Fullname of GRCh37.gene.bed')
parser.add_argument('-sm', '--spm-min', dest='sm', type=float,
                    help='Min SPM, default 0.5', default=0.5)
parser.add_argument('-sn', '--spm-max', dest='sn', type=float,
                    help='Max SPM, default 1.0', default=1.0)
args = parser.parse_args()

raw_files = glob.glob(os.path.join(args.path, "*.txt"))
grch37 = ""
try:
    with open(args.grch37, "r") as f:
        grch37 = f.readlines()
except Exception as e:
    print(e)

grch37_genes = {}

for (i, x) in enumerate(grch37):
    gene = x.split("\t")[4].replace('\n', "")
    grch37_genes[gene] = i

for raw in raw_files:
    with open(raw, "r") as f:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        new_file_gene = os.path.join(args.output, os.path.splitext(
            os.path.basename(raw))[0]+"."+str(args.sm)+"_"+str(args.sn)+".gene.txt")
        new_file_TSS = os.path.join(args.output, os.path.splitext(
            os.path.basename(raw))[0]+"."+str(args.sm)+"_"+str(args.sn)+".TSS.bed")
        out_gene = open(new_file_gene, "w")
        out_TSS = open(new_file_TSS, "w")

        for line in f.readlines():
            if line.startswith(('^', '!', 'Gene')):
                continue
            line = line.strip().split("\t")
            try:
                spm = float(line[4])
            except ValueError:
                continue
            if spm < args.sn and spm > args.sm:
                continue
            if line[0] in grch37_genes:
                grch = grch37[grch37_genes[line[0]]]
                out_gene.write(grch)
                grch_list = grch.split("\t")
                if grch_list[3] == "+":
                    left = int(grch_list[1])-1000
                    out_TSS.write("\t".join([grch_list[0], str(left), str(
                        left+2000), grch_list[3], grch_list[4]]))
                else:
                    right = int(grch_list[2])-1000
                    out_TSS.write("\t".join([grch_list[0], str(right), str(
                        right+2000), grch_list[3], grch_list[4]]))
        out_gene.close()
        out_TSS.close()
        print("save to {}, {}".format(new_file_gene, new_file_TSS))
