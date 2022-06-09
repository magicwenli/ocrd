import argparse
from copy import deepcopy
import glob
import math
import random
from os.path import join as pjoin

import numpy as np
import pandas as pd

from utils import draw_utils, feat_utils, log
from utils.Config import Config

logger = log.get_logger(__name__, log.INFO)
conf = Config()
conf.plot = False

parser = argparse.ArgumentParser(description="extract features")
parser.add_argument("-n", "--none-bed", required=True,
                    help="none tss bed file", dest="none_bed")
parser.add_argument("-t", "--tss-bed", required=True,
                    help="tss bed file", dest="tss_bed")
parser.add_argument("-f", "--fraction", required=True,
                    help="coverage fraction: 0.0 ~ 1.0, default: 1.0", dest="fraction", default=1.0, type=float)
parser.add_argument("-o", "--output", help="output dir", dest="output")
args = parser.parse_args()

if not args.output:
    output_dir = './build/'
else:
    output_dir = args.output

fraction = args.fraction

bamfileList = ['/Source/fastq/IA01.bam',
               '/Source/fastq/IC05.bam', ]

# none, mid, left, right
bedDict = {
    0: [args.none_bed],
    1: [args.tss_bed],
    2: [args.tss_bed],
    3: [args.tss_bed]
}

rounds = 0
conting = 1


for classType, bedfileList in bedDict.items():
    # filter for certain type
    # if classType in [0, 2, 3]:
    #     continue

    features_df = pd.DataFrame(columns=(
        'classType', "depth", "bigDepth", "smallDepth", "height", "width", "angel", "area", "leftK", "rightK", "var_Hight", "var_Width", "var_Angel", "var_Area", "leftKVar", "rightKVar", "mTpye", "m0k", "m1k", "m2k", "m0b", "m1b", "m2b", "vector_0", "vector_1", "vector_2", "vector_3", "vector_4", "vector_5", "conting", "start", "end"))
    data_info = {'length': 0}
    c_rounds = 1

    for bed in feat_utils.bedIterator(bedfileList, data_info, conting):
        if classType == 3:  # right
            start = int((int(bed[1]) + int(bed[2])) / 2) + \
                random.randint(-200, 200)
            # end = start + 2000 + random.randint(-1000, 500)
            end = start + 2000
        elif classType == 2:  # left
            end = int((int(bed[1]) + int(bed[2])) / 2) + \
                random.randint(-200, 200)
            # start = end - 2000 + random.randint(-500, 1000)
            start = end - 2000
        else:
            start = int(bed[1])
            end = int(bed[2])

        print("\r*** Type %d, Round %d of %d, Total %d ***" %
              (classType, c_rounds, data_info["length"], rounds), end='')
        rounds += 1
        c_rounds += 1

        # filter for certain round
        # if rounds != 10:
        #     continue

        step = end - start
        length = step + 1
        # lFdepth_Nor isize in [120, 180]
        # sFdepth_Nor isize in [35, 80]
        wpsList_Nor, lFdepth_Nor, sFdepth_Nor = feat_utils.callOneBed(
            bamfileList, str(bed[0]), start, end, win=120, fraction=fraction)
        if conf.plot:
            draw_utils.plot_one(wpsList_Nor, range(
                0, end-start), title="WPS波形（初始）", xlabel='初始WPS波形', ylabel='WPS', size=(10, 2.5))
            draw_utils.plot_one(lFdepth_Nor, range(
                0, end-start), title="测序覆盖度（初始）", xlabel='参考基因组坐标', ylabel='测序深度')
        wpsList_Ori = deepcopy(wpsList_Nor)
        adjustWpsList_Nor = feat_utils.adjustWPS(wpsList_Nor)

        base = feat_utils.getBaseline(adjustWpsList_Nor)
        smoothWpsList_Nor = feat_utils.smoothWPS(adjustWpsList_Nor, base)
        lFdepth_Nor = feat_utils.savgol_filter(lFdepth_Nor, 801, 2)

        if conf.plot:
            p = draw_utils.PlotMultiline(
                "", ylabel="", xlabel="参考基因组坐标", name="sg_wps")
            p.plot_continuous(adjustWpsList_Nor,
                              range(0, end-start), "初始WPS波形")
            p.plot_continuous(smoothWpsList_Nor,
                              range(0, end-start), "S-G滤波后WPS波形")
            p.plot_end()

        peaks, properties = feat_utils.find_peaks(
            smoothWpsList_Nor, height=0.28, distance=25, prominence=0.25, width=[25, 170])
        peakObjectList = feat_utils.getValley(
            smoothWpsList_Nor, adjustWpsList_Nor, peaks, 5)

        if conf.plot:
            p = draw_utils.PlotPeaks(
                "smoothed WPS and peaks", name="wps_peaks_2")
            p.plot_continuous(smoothWpsList_Nor,
                              range(0, end-start), "S-G滤波后WPS波形")
            index = 0
            for peak in peakObjectList:
                index += 1
                p.plot_continuous(smoothWpsList_Nor[peak.startPos:peak.endPos],
                                  range(peak.startPos, peak.endPos), legend="波峰 %d" % index, rectangle=True)
            p.plot_end()

        try:
            feat_utils.getAllPeaksAveHeight(peakObjectList, smoothWpsList_Nor)
        except:
            logger.warn("getALLPeaksAveHeight error")

        # 左侧2/3 右侧2/3 和 全部覆盖度曲线的线性模型
        sublen = int(len(lFdepth_Nor) / 3)
        models = []
        xList = []
        yList = []

        model, new_x, new_y = feat_utils.linearJudgeNDR(
            lFdepth_Nor, 0, 2 * sublen - 1)
        if model.coef_ == 0:
            logger.debug("model.coef_ == 0")
            continue
        models.append(model)
        if conf.plot:
            xList.append(new_x)
            yList.append(new_y)

        model, new_x, new_y = feat_utils.linearJudgeNDR(
            lFdepth_Nor, sublen, 3 * sublen - 1)
        if model.coef_ == 0:
            logger.debug("model.coef_ == 0")
            continue
        models.append(model)
        if conf.plot:
            xList.append(new_x)
            yList.append(new_y)

        model, new_x, new_y = feat_utils.linearJudgeNDR(
            lFdepth_Nor, 0, len(lFdepth_Nor) - 1)
        if model.coef_ == 0:
            logger.debug("model.coef_ == 0")
            continue
        models.append(model)
        if conf.plot:
            xList.append(new_x)
            yList.append(new_y)
        if conf.plot:
            draw_utils.plot_linear_line(models, xList, yList)

        vectors = feat_utils.getMidVector(models)

        meanWidth, meanHeight, meanAngel, meanArea, varWidth, varHeight, varAngel, varArea, peakOri = feat_utils.haveNearContinuouslyPeak(
            smoothWpsList_Nor, adjustWpsList_Nor, peakObjectList)

        if conf.plot:
            p = draw_utils.PlotPeaks2(
                "WPS波形特征", name="wps_peaks_features", oriData=peakOri)
            p.plot_continuous(smoothWpsList_Nor,
                              range(0, end-start), "S-G滤波后WPS波形")
            index = 0
            for peak in peakObjectList:
                p.plot_continuous(smoothWpsList_Nor[peak.startPos:peak.endPos],
                                  range(peak.startPos, peak.endPos), legend="波峰 %d" % index, rectangle=True)
                if index == 1:
                    p.plot_tri(index)
                if index == 2:
                    p.plot_height(index)
                if index == 3:
                    p.plot_weight(index)
                if index == 4:
                    p.plot_angle(index)
                index += 1
            p.plot_end()

        ndrAreaDepth, smallNdrAreaDepth, minIndex = feat_utils.judgeLowDepth(
            lFdepth_Nor, 0, len(lFdepth_Nor) - 1)

        if conf.plot:
            # plot 最小覆盖度
            pc = draw_utils.PlotCoverage3(
                '最小覆盖度区间', 'coverage_features_4', lFdepth_Nor)
            pc.plot_continuous([ndrAreaDepth for _ in range(600)],
                               range(minIndex-300, minIndex+300))
            pc.plot_continuous([smallNdrAreaDepth for _ in range(200)],
                               range(minIndex-100, minIndex+100))
            pc.plot_end()

            # titles = ['左侧2/3区域', '全部区域', "右侧2/3区域"]
            for i, v in enumerate([0, 2, 1]):
                pc = draw_utils.PlotCoverage1(
                    '覆盖度曲线', 'coverage_features_{}'.format(i), lFdepth_Nor)
                _y = models[v].coef_ * xList[v] + models[v].intercept_
                pc.plot_continuous(_y, xList[v])

                pc.plot_end()

            pc = draw_utils.PlotCoverage2(
                '角平分线', 'angle_features_1', lFdepth_Nor)
            _y1 = models[0].coef_ * xList[0] + models[0].intercept_
            _y2 = models[1].coef_ * xList[1] + models[1].intercept_

            # print(models[0].coef_, models[0].intercept_)
            # print(models[1].coef_, models[1].intercept_)

            k1, k2 = models[0].coef_, models[1].coef_
            k_bisector = (
                k1*k2 - 1 - math.sqrt((1-k1*k2)**2+(k1+k2)**2))/(k1+k2)

            pc.plot_continuous(_y1, _y2, xList[2], '角平分线', models, k_bisector)
            pc.plot_end()

        feat = []

        feat.append(classType)
        depth = np.sum(lFdepth_Nor) / length
        feat.append(depth)
        feat.append(ndrAreaDepth)
        feat.append(smallNdrAreaDepth)
        feat.append(meanHeight)
        feat.append(meanWidth)
        feat.append(meanAngel)
        feat.append(meanArea)
        leftK = np.mean([p.leftK for p in peakObjectList])
        feat.append(leftK)
        rightK = np.mean([p.rightK for p in peakObjectList])
        feat.append(rightK)
        feat.append(varHeight)
        feat.append(varWidth)
        feat.append(varAngel)
        feat.append(varArea)
        leftKVar = np.var([p.leftK for p in peakObjectList])
        feat.append(leftKVar)
        rightKVar = np.var([p.rightK for p in peakObjectList])
        feat.append(rightKVar)
        mTpye = feat_utils.modelType(models)
        feat.append(mTpye)
        m0k = models[0].coef_[0, 0]
        m1k = models[1].coef_[0, 0]
        m2k = models[2].coef_[0, 0]
        m0b = models[0].intercept_[0]
        m1b = models[1].intercept_[0]
        m2b = models[2].intercept_[0]
        feat.extend([m0k, m1k, m2k, m0b, m1b, m2b])
        for v in vectors:
            feat.append(v)
        feat.append(int(conting))
        feat.append(start)
        feat.append(end)

        bed_df = pd.DataFrame([feat], columns=(
            'classType', "depth", "bigDepth", "smallDepth", "height", "width", "angel", "area", "leftK", "rightK", "var_Hight", "var_Width", "var_Angel",
            "var_Area", "leftKVar", "rightKVar", "mTpye", "m0k", "m1k", "m2k", "m0b", "m1b", "m2b", "vector_0", "vector_1", "vector_2", "vector_3", "vector_4", "vector_5", "conting", "start", "end"))

        features_df = pd.concat([features_df, bed_df], ignore_index=True)

    features_df.to_csv(pjoin(output_dir, "features_df_{}.csv".format(
        classType)), sep=',', header=True, index=True)

all_filenames = [i for i in glob.glob(pjoin(output_dir, "features_df_*.csv"))]
list_of_files = [pd.read_csv(f, index_col=0) for f in all_filenames]
# combine all files in the list
combined_csv = pd.concat(list_of_files, ignore_index=True)
# export to csv
combined_csv.to_csv(pjoin(output_dir, "features_df_combined.csv"),
                    index=False, encoding='utf-8-sig')
logger.info("csv combined into features_df_combined.csv")
