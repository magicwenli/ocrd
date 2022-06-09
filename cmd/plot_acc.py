import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams["font.sans-serif"] = ["simsun"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题
plt.rcParams.update({'font.size': 12})

source = "/Source/scripts/build/2022-06-05_16-48-22/result.csv"
workdir = os.path.dirname(source)


def plot_groupbar(series1, series2, save_path, ylabel, xlabel, labels):
    fig, ax = plt.subplots()
    width = 0.3  # the width of the bars
    x = np.arange(len(labels))
    rects1 = ax.bar(x - width/2-0.025, series1, width, label='普通随机森林分类器')
    rects2 = ax.bar(x + width/2+0.025, series2, width, label='置信学习预处理的随机森林分类器')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    # ax.set_title('Scores by group and gender')
    ax.set_xticks(x, labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()
    fig.set_size_inches(8, 4)
    fig.savefig(os.path.join(workdir, save_path+'.pdf'),
                bbox_inches='tight', dpi=300)


def plot_groupbar2(series1, series2, save_path, ylabel, xlabel, labels):
    fig, ax = plt.subplots()
    width = 0.3  # the width of the bars
    x = np.arange(len(labels))
    rects1 = ax.bar(x - width/2-0.025, series1, width, label='普通随机森林分类器')
    rects2 = ax.bar(x + width/2+0.025, series2, width, label='置信学习预处理的随机森林分类器')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 1.1)
    # ax.set_title('Scores by group and gender')
    ax.set_xticks(x, labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()
    fig.set_size_inches(8, 4)
    fig.savefig(os.path.join(workdir, save_path+'.pdf'),
                bbox_inches='tight', dpi=300)


data = pd.read_csv(source, header=None, names=[
                   "fraction", "isCL", "noiseRatio", "accscore", "aucscore", "recallsocre", "f1score"], index_col=False)

no_noise = data[data.noiseRatio == 0]
no_noise_cl = no_noise[no_noise.isCL == 1].round(2)
no_noise_no_cl = no_noise[no_noise.isCL == 0].round(2)
labels = ['1X', '2X', '5X', '15X', '30X', '100X']
plot_groupbar(no_noise_no_cl.aucscore.to_numpy(),
              no_noise_cl.aucscore.to_numpy(), "auc_no_noise", 'AUC Score', "覆盖度", labels)
plot_groupbar(no_noise_no_cl.accscore.to_numpy(),
              no_noise_cl.accscore.to_numpy(), "acc_no_noise", 'Accuracy Score', "覆盖度", labels)
no_noise[no_noise.isCL == 0].to_latex(os.path.join(workdir, 'no_noise.tex'))
no_noise[no_noise.isCL == 1].to_latex(os.path.join(workdir, 'no_noise_cl.tex'))


full_coverage = data[data.fraction == 1.0]
full_coverage_cl = full_coverage[full_coverage.isCL == 1].round(2)
full_coverage_no_cl = full_coverage[full_coverage.isCL == 0].round(2)
labels = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6']
plot_groupbar2(full_coverage_no_cl.aucscore.to_numpy(),
               full_coverage_cl.aucscore.to_numpy(), "auc_full_coverage", 'AUC Score', "标签噪声比例", labels)
full_coverage[full_coverage.isCL == 1].to_latex(
    os.path.join(workdir, 'full_coverage.tex'))

tri_coverage = data[data.fraction == 0.3]
tri_coverage_cl = tri_coverage[tri_coverage.isCL == 1].round(2)
tri_coverage_no_cl = tri_coverage[tri_coverage.isCL == 0].round(2)
labels = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6']
plot_groupbar2(tri_coverage_no_cl.aucscore.to_numpy(),
               tri_coverage_cl.aucscore.to_numpy(), "auc_tri_coverage", 'AUC Score', "标签噪声比例", labels)
tri_coverage[tri_coverage.isCL == 1].to_latex(
    os.path.join(workdir, 'tri_coverage.tex'))
