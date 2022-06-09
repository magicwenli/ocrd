import warnings
from os.path import join as pjoin
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
from scipy.fftpack import tilbert
import seaborn as sns
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)

plt.rcParams["font.sans-serif"] = ["simsun"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题
plt.rcParams.update({'font.size': 12})
# style.use('ggplot')


def plot_one(y, x=[], title="", xlabel="position", ylabel="value",  save_path="./build/figures", size=None):
    '''
    plot line, y is a list
    x = [i for i in range(len(y))]
    '''
    if len(x) == 0:
        x = [i for i in range(len(y))]
    if size:
        plt.gcf().set_size_inches(size[0], size[1])
    plt.ion()
    plt.plot(x, y, color='k', label=xlabel)
    # plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc="upper right")
    # plt.title(title)
    # plt.ticklabel_format(useOffset=False)
    plt.savefig(pjoin(save_path, "{}.pdf".format(title)), dpi=300)
    plt.close()


class PlotMultiline():
    linetype = ['silver', 'k', 'g']

    def __init__(self, title="", name="", xlabel="position", ylabel="value", save_path="./build/figures"):
        self.title = title.replace(" ", "_")
        self.save_path = save_path
        self.fmt_index = 0
        self.name = name
        plt.ion()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.gcf().set_size_inches(10, 5)
        # plt.title(title)

    def plot_continuous(self, y, x=[], legend=""):
        if len(x) == 0:
            x = [i for i in range(len(y))]
        plt.plot(x, y, self.linetype[self.fmt_index], label=legend)
        self.fmt_index += 1

    def plot_end(self):
        plt.legend(loc="upper right")
        plt.savefig(
            pjoin(self.save_path, "{}.pdf".format(self.name)), dpi=300)
        plt.close()

    def plt(self, *args, **kargs):
        plt.plot(*args, **kargs)


class PlotPeaks(PlotMultiline):
    def __init__(self, title="", name="", xlabel="参考基因组坐标", ylabel="归一化WPS", save_path="./build/figures"):
        super().__init__(title, name, xlabel, ylabel, save_path)
        plt.gcf().set_size_inches(10, 2.5)

    def plot_continuous(self, y, x=[], legend="", rectangle=False):
        if len(x) == 0:
            x = [i for i in range(len(y))]
        if self.fmt_index == 0:
            plt.plot(x, y, "silver", label=legend)
        else:
            plt.plot(x, y)
            if rectangle:
                pos = (min(x), min(y))
                width = max(x) - min(x)
                height = max(y)-min(y)
                plt.gca().add_patch(Rectangle(pos, width, height, lw=1, ls=':', fill=False))
                plt.gca().text(
                    min(x)+width/2, min(y)-0.01, legend,
                    horizontalalignment='center',
                    verticalalignment='top')
        self.fmt_index += 1

    def plot_end(self):
        plt.legend(loc="best")
        plt.savefig(
            pjoin(self.save_path, "{}.pdf".format(self.name)), dpi=300)
        plt.close()


class PlotPeaks2(PlotPeaks):
    text_pos = []

    def __init__(self, title="", name="", xlabel="参考基因组坐标", ylabel="归一化WPS", save_path="./build/figures", oriData=None):
        super().__init__(title, name, xlabel, ylabel, save_path)
        plt.gcf().set_size_inches(10, 5)
        self.ori = oriData

    def plot_continuous(self, y, x=[], legend="", rectangle=False):

        if len(x) == 0:
            x = [i for i in range(len(y))]
        if self.fmt_index == 0:
            plt.plot(x, y, "silver", label=legend)
        else:
            plt.plot(x, y)
            if rectangle:
                pos = (min(x), min(y))
                width = max(x) - min(x)
                height = max(y)-min(y)
                # plt.gca().add_patch(Rectangle(pos, width, height, lw=1, ls=':', fill=False))
                # plt.gca().text(
                #     min(x)+width/2, min(y)-0.01, legend,
                #     horizontalalignment='center',
                #     verticalalignment='top')
            self.text_pos.append((min(x)+width/2, 1.06))
        self.fmt_index += 1

    def plot_tri(self, index):
        pivot = plt.Polygon(self.ori["triPoints"]
                            [index], closed=True, hatch='//', fill=False, zorder=20)
        plt.gca().text(self.text_pos[index][0], self.text_pos[index][1], "波峰三角形面积",
                       horizontalalignment='center',
                       verticalalignment='top')
        plt.gca().add_patch(pivot)

    def plot_height(self, index):
        data = self.ori["data"]
        plt.gca().text(self.text_pos[index][0], self.text_pos[index][1], "波峰高度",
                       horizontalalignment='center',
                       verticalalignment='top')
        plt.hlines([data[self.ori['indexs'][index]], min(
            data[self.ori['startPos'][index]], data[self.ori['endPos'][index]])],
            xmin=[self.ori['indexs'][index]]*2,
            xmax=[self.ori['indexs'][index]+80]*2,
            lw=1.4, color="k")
        plt.vlines(self.ori['indexs'][index]+40, ymin=min(
            data[self.ori['startPos'][index]], data[self.ori['endPos'][index]]),
            ymax=data[self.ori['indexs'][index]],
            lw=1.4, color="k", ls='--')

    def plot_weight(self, index):
        data = self.ori["data"]
        plt.gca().text(self.text_pos[index][0], self.text_pos[index][1], "波峰宽度",
                       horizontalalignment='center',
                       verticalalignment='top')
        plt.vlines([self.ori['startPos'][index], self.ori['endPos'][index]],
                   ymin=[data[self.ori['indexs'][index]]-0.035]*2,
                   ymax=[data[self.ori['indexs'][index]]+0.035]*2,
                   lw=1.4, color="k")
        plt.hlines(data[self.ori['indexs'][index]],
                   xmin=self.ori['startPos'][index],
                   xmax=self.ori['endPos'][index],
                   lw=1.4, color="k", ls='--')

    def plot_angle(self, index):
        data = self.ori["data"]
        plt.gca().text(self.text_pos[index][0], self.text_pos[index][1], "波峰夹角",
                       horizontalalignment='center',
                       verticalalignment='top')
        a = self.ori['triPoints'][index]
        plt.plot([x for x, _ in a], [x for _, x in a], 'k-', lw=1.4)


class PlotCoverage1(PlotMultiline):
    def __init__(self, title, name, rawDepthList, xlabel="", ylabel="", save_path="./build/figures"):
        super().__init__(title, name, xlabel, ylabel, save_path)
        plt.gcf().set_size_inches(6, 4)
        self.raw = rawDepthList
        self.len = len(rawDepthList)

    def plot_continuous(self, y, x=[], title=""):
        plt.plot(range(self.len), self.raw, 'silver')  # origin coverage line
        plt.plot(x, y, '--k')   # linear model line
        pos_x = [x[0][0], x[-1][0]]
        pos_ymin = [x - 0.15 for x in [y[0][0], y[-1][0]]]
        pos_ymax = [x + 0.15 for x in [y[0][0], y[-1][0]]]
        plt.vlines(pos_x, pos_ymin,
                   pos_ymax, ls='-', colors='k')

    def plot_end(self):
        plt.savefig(
            pjoin(self.save_path, "{}.pdf".format(self.name)), dpi=300)
        plt.close()


class PlotCoverage2(PlotCoverage1):
    '''
    画两直线之间的角度
    '''

    def __init__(self, title, name, rawDepthList, xlabel="", ylabel="", save_path="./build/figures"):
        super().__init__(title, name, rawDepthList, xlabel, ylabel, save_path)

    def plot_continuous(self, y1, y2, x, title, models, k):
        plt.plot(range(self.len), self.raw, 'silver')  # origin coverage line

        A, B, C, D = (x[0], y1[0]), (x[len(y1)-1], y1[-1]
                                     ), (x[-len(y2)+1], y2[0]), (x[-1], y2[-1])

        intersection = line_intersection((A, B), (C, D))

        # plt.plot([x[0] for x in (A, B, C, D)], [x[1]
        #          for x in (A, B, C, D)], 'o')

        int_x = int(intersection[0])
        plt.plot(x[:int_x], y1[:int_x], '--k')
        plt.plot(x[int_x:], y2[-len(x)+int_x:], '--k')
        plt.savefig(
            pjoin(self.save_path, "{}2.pdf".format(self.name)), dpi=300)
        _y = np.linspace(intersection[1][0], intersection[1][0]+2, num=100)
        _b = y1[int_x]-k*x[int_x]
        _x = (_y-_b[0][0])/k[0][0]
        plt.plot(_x, _y, ':k')


class PlotCoverage3(PlotCoverage1):
    index = 0

    def plot_continuous(self, y, x=[], title=""):
        plt.plot(range(self.len), self.raw, 'silver',
                 zorder=0)  # origin coverage line
        if self.index == 0:
            plt.plot(x, y, '--k', label='600bp平均覆盖度最小值')   # linear model line
            self.index += 1
        else:
            plt.plot(x, y, ':k', label='300bp平均覆盖度最小值')   # linear model line

        pos_x = [x[0], x[-1]]
        pos_ymin = [x - 0.10 for x in [y[0], y[-1]]]
        pos_ymax = [x + 0.10 for x in [y[0], y[-1]]]
        plt.vlines(pos_x, pos_ymin, pos_ymax, ls='-', colors='k')
        # if self.index == 0:
        #     plt.gca().text(x[len(x)//2], y[0]+0.2, "600bp区间最小覆盖度", horizontalalignment='center',
        #                    verticalalignment='bottom')
        #     self.index += 1
        # else:
        #     plt.gca().text(x[len(x)//2], y[0]-0.2, "300bp区间最小覆盖度", horizontalalignment='center',
        #                    verticalalignment='top')

    def plot_end(self):
        plt.legend(loc='best')
        return super().plot_end()


def plot_linear_line(modelList, xList, yList, title="linear_line", save_path="./build/figures"):
    figure, axes = plt.subplots(3, 1, figsize=(3 * 3, 3 * 3))
    for i in range(3):
        axes[i].plot(xList[i], yList[i], color='k', label='original')
        axes[i].plot(xList[i], modelList[i].coef_ * xList[i] +
                     modelList[i].intercept_, color='gray', label='linear')
        axes[i].legend()
        if i == 0:
            axes[i].set_title("left 2/3")
        elif i == 1:
            axes[i].set_title("right 2/3")
        else:
            axes[i].set_title("full curve")

    plt.savefig(pjoin(save_path, "{}.png".format(title)))
    plt.close()


def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def plot_feature_line(data, title="", save_path="./build/figures"):
    plt.ion()
    namemap = dict(zip(["depth", "bigDepth", "smallDepth", "height",
                   "width", "angel", "area", "leftK", "rightK",
                        "var_Hight", "var_Width", "var_Angel", "var_Area",
                        "leftKVar", "rightKVar", "m0k", "m1k", "m2k", "m0b", "m1b", "m2b"],
                       ["mean_coverage", "broad_interval_coverage",
                        "narrow_interval_coverage", "height_mean", "width_mean",
                        "degree_mean", "area_mean", "peak_left_coefficient_mean", "peak_right_coefficient_mean",
                        "height_variance", "width_variance", "degree_variance",
                        "area_variance""peak_left_coefficient_variance",
                        "peak_right_coefficient_variance", "left_coefficient", "right_coefficient", "full_coefficient", "left_intercept", "right_intercept", "full_intercept"]))
    class_map = {
        0: '无开放区域',
        1: '开放区域位于中间',
        2: '开放区域位于左侧',
        3: '开放区域位于右侧',
    }

    def get_map(title):
        if title in namemap:
            return namemap[title]
        else:
            return title

    for classType in data['classType'].unique():
        ax = sns.distplot(data[data['classType'] == classType]
                          [title], label=class_map[classType], axlabel=get_map(title))
    fig = ax.get_figure()
    fig.legend(loc="upper left", mode='expand', ncol=4)
    fig.set_size_inches(9, 4)
    fig.savefig(pjoin(save_path, "figures", "{}.pdf".format(get_map(title))),
                bbox_inches='tight', dpi=300)
    print('save to {}'.format(
        pjoin(save_path, "figures", "{}.png".format(get_map(title)))))
    plt.close()
