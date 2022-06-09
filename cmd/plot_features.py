import argparse

import pandas as pd
from utils.draw_utils import plot_feature_line

parser = argparse.ArgumentParser(description="plot feature line")
parser.add_argument("csv", help="csv file")
parser.add_argument('-o', '--output', help='output file', dest='output')
args = parser.parse_args()

data = pd.read_csv(args.csv)

titles = list(data.head())
titles.remove('classType')
titles.remove('conting')
titles.remove('start')
titles.remove('end')

for title in titles:
    plot_feature_line(data[['classType', title]],
                      title=title, save_path=args.output)
