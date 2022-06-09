import math
import os
import sys
import random

import numpy as np
import pandas as pd
import peakutils
import pysam
from scipy.integrate import simps
from scipy.signal import find_peaks, find_peaks_cwt, savgol_filter
from sklearn import linear_model, preprocessing

from . import log
from .Peak import Peak

logger = log.get_logger(__name__, log.INFO)


def readIterator(filenames, contig, start, end, fraction=1.0):
    '''
    :param filenames: list of bam files
    :param contig: contig name without `chr`. ex: 1
    :param start: start position
    :param end: end position
    :return: iterator of reads
    '''
    for bamfile in filenames:
        logger.debug("Reading %s" % bamfile)
        if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam", ".bai")) or os.path.exists(bamfile+".bai")):
            input_file = pysam.Samfile(bamfile, "rb")

            try:
                input_file.fetch(contig, start-1, end)
            except ValueError:
                contig = 'chr' + str(contig)

            for read in input_file.fetch(contig, start-1, end):
                if random.random() < fraction:
                    yield read
                else:
                    pass
            input_file.close()
        else:
            logger.warning("BAM file %s or BAM index not found" % bamfile)
            sys.exit()


def callOneBed(filenames, contig, bed1, bed2, win, fraction=1.0):
    bed1 = bed1 - win
    bed2 = bed2 + win
    length = bed2 - bed1 + 1
    array = np.zeros(length, dtype=int)
    depth = np.zeros(length, dtype=int)
    depth2 = np.zeros(length, dtype=int)
    if type(contig) is str:
        contig = contig.replace("chr", "")
    for r in readIterator(filenames, contig, bed1, bed2, fraction):
        if (not r.is_reverse) and (not r.is_unmapped) and (not r.mate_is_unmapped) and r.mate_is_reverse:
            logger.debug("{} {} {} {}".format(r.reference_name, r.template_length,
                         r.reference_start, r.reference_start+r.template_length))
            if r.template_length >= 35 and r.template_length <= 80:
                start = r.reference_start - bed1
                end = r.reference_start + r.template_length - bed1
                # depth + 1
                dstart = start
                dend = end
                if dstart < 0:
                    dstart = 0
                if dend > length:
                    dend = length
                d = dstart
                while d < dend:
                    depth2[d] += 1
                    d += 1

            if r.template_length < win or r.template_length > 180:
                continue
            start = r.reference_start - bed1
            end = r.reference_start + r.template_length - bed1
            # depth + 1
            dstart = start
            dend = end
            if dstart < 0:
                dstart = 0
            if dend >= length:
                dend = length
            d = dstart
            while d < dend:
                depth[d] += 1
                d += 1

            # [$start+W/2,$end-W/2] WPS+1
            region1 = start + int(win / 2)
            region2 = end - int(win / 2)
            if region1 < 0:
                region1 = 0
            if region2 > length:
                region2 = length
            array[region1: region2] += 1

            # [$start-w/2,$start-1+w/2] WPS-1
            region1 = start - int(win / 2)
            region2 = start + int(win / 2) + 1
            if region1 < 0:
                region1 = 0
            if region2 > length:
                region2 = length
            array[region1: region2] -= 1

            # [end-w/2+1,$end+w/2] WPS-1
            region1 = end - int(win / 2) + 1
            region2 = end + int(win / 2)
            if region1 < 0:
                region1 = 0
            if region2 > length:
                region2 = length
            array[region1: region2] -= 1

    lenth1 = len(array) - win - 1
    bed1 += win
    bed2 -= win
    array = np.array(array[win: lenth1], dtype=float)
    depth = depth[win: lenth1]
    depth2 = depth2[win: lenth1]
    return array, depth, depth2


def bedIterator(filenames, datainfo={}, chrom="chr1"):
    if type(chrom) is int:
        chrom = "chr" + str(chrom)
    chrom_num = chrom.replace("chr", "")
    for file in filenames:
        if os.path.exists(file):
            data = pd.read_csv(file, sep="\t", header=None)
            data = data[(data[0] == chrom) | (
                data[0] == int(chrom_num)) | (data[0] == chrom_num)]
            datainfo['length'] = data.shape[0]
            for i in range(len(data)):
                yield data.iloc[i]
        else:
            logger.warning("BED file %s not found" % file)
            pass


def adjustWPS(wpsList):
    subarray = wpsList[0: len(wpsList)]
    tmpLength = len(subarray)
    medianA = [0] * tmpLength
    chaA = [0] * tmpLength
    adjustWin = 1000
    start = 0
    end = adjustWin + 1
    while end <= tmpLength:
        # 滑动窗口 找到长1000窗口的中位数、均值和跨度（最大值减最小值）
        # 以窗口中值为索引存放
        tmpArray = subarray[start: end]
        minn, median, maxn = tmpArray.min(), np.median(
            np.array(tmpArray)), tmpArray.max()
        tmploc = int((start + end + 1) / 2)
        medianA[tmploc] = median
        chaA[tmploc] = maxn - minn + 1
        start += 1
        end += 1
    x = 0
    while x < tmpLength:
        loc = x
        if loc < 501:
            loc = 501
        if loc >= tmpLength - 501:
            loc = tmpLength - 501
        # 遍历
        if loc < 0:
            loc = 0
        # wps值 = （wps-所在窗口的均值）/ 窗口跨度
        if chaA[loc] != 0:
            subarray[x] = (subarray[x] - medianA[loc]) / chaA[loc]
        else:
            subarray[x] = 0
        x += 1
    return np.array(subarray)


def getBaseline(adjustWpsList_Nor):
    try:
        base = peakutils.baseline(adjustWpsList_Nor, 8)
    except (ZeroDivisionError, ValueError) as reason:  # 'ZeroDivisionError'除数等于0的报错方式^M
        base = np.zeros(len(adjustWpsList_Nor))
    return base


def smoothWPS(adjustWpsList_Nor, base):
    try:
        adjustWpsList_Nor = np.subtract(adjustWpsList_Nor, base)
        smoothWpsList_Nor = savgol_filter(adjustWpsList_Nor, 51, 1)
        smoothWpsList_Nor = preprocessing.minmax_scale(smoothWpsList_Nor)
    except ValueError:
        logger.warn("smooth error")
    return smoothWpsList_Nor


def trackLeftValley(wpsList, peakIndex, slidWinSize):
    win1 = wpsList[peakIndex - slidWinSize:peakIndex]
    win2 = wpsList[peakIndex - slidWinSize * 2:peakIndex - slidWinSize]
    win3 = wpsList[peakIndex - slidWinSize * 3:peakIndex - slidWinSize * 2]
    i = 3
    length = len(win1) + len(win2) + len(win3)
    while (peakIndex - slidWinSize * i > 0 and slidWinSize * i <= 250 and not (
            np.sum(win1) > np.sum(win2) and np.sum(win2) < np.sum(win3))) or slidWinSize * i < 60:
        win1 = win2
        win2 = win3
        win3 = wpsList[peakIndex - slidWinSize *
                       (i + 1): peakIndex - slidWinSize * i]

        i += 1
    if peakIndex - slidWinSize * i < 0 or slidWinSize * i > 250:
        if peakIndex > 80:
            return peakIndex - 80
        else:
            return 0
    valIndex = peakIndex - slidWinSize * i + int((slidWinSize + 1) / 2)
    curIndex = valIndex + int((slidWinSize + 1) / 2)
    while curIndex < peakIndex - 40 and abs(wpsList[valIndex] - wpsList[curIndex]) < 0.05:
        curIndex += 1

    return curIndex


def trackRightValley(wpsList, peakIndex, slidWinSize):
    win1 = wpsList[peakIndex:peakIndex + slidWinSize]
    win2 = wpsList[peakIndex + slidWinSize:peakIndex + slidWinSize * 2]
    win3 = wpsList[peakIndex + slidWinSize * 2:peakIndex + slidWinSize * 3]
    i = 3
    while ((not (np.sum(win1) > np.sum(win2) and np.sum(win2) < np.sum(win3)) and peakIndex + slidWinSize * i < len(
            wpsList) - 1) and slidWinSize * i <= 250) or slidWinSize * i < 60:
        win1 = win2
        win2 = win3
        win3 = wpsList[peakIndex + slidWinSize *
                       i: peakIndex + slidWinSize * (i + 1)]
        i += 1
    if peakIndex + slidWinSize * i > len(wpsList):
        return len(wpsList) - 1
    valIndex = peakIndex + slidWinSize * i - int((slidWinSize + 1) / 2)
    curIndex = valIndex - int((slidWinSize + 1) / 2)
    while curIndex > peakIndex + 60 and abs(wpsList[valIndex] - wpsList[curIndex]) < 0.05:
        curIndex -= 1

    return curIndex


def getValley(wpsList, rawDataList, peaks, slidWinSize):
    '''
    wpsList: smoothed wps signal
    rawDataList: adjusted wps signal
    peaks: peak index
    slidWinSize: sliding window size
    '''
    leftValleyList = []
    rightValleyList = []
    for peakX in peaks:
        left = trackLeftValley(wpsList, peakX, slidWinSize)
        # start postion of the peak
        right = trackRightValley(wpsList, peakX, slidWinSize)
        # end postion of the peak
        leftValleyList.append(left)
        rightValleyList.append(right)
    peakObjectList = []

    for i in range(len(peaks)):
        leftPeakK = 0
        rightPeakK = 0
        if peaks[i] - leftValleyList[i] != 0 and rawDataList[peaks[i]] != rawDataList[leftValleyList[i]]:
            leftPeakK = (
                rawDataList[peaks[i]] - rawDataList[leftValleyList[i]]) / (peaks[i] - leftValleyList[i])
        if peaks[i] - rightValleyList[i] != 0 and rawDataList[peaks[i]] != rawDataList[rightValleyList[i]]:
            rightPeakK = (
                rawDataList[peaks[i]] - rawDataList[rightValleyList[i]]) / (peaks[i] - rightValleyList[i])
        peakObjectList.append(
            Peak(peaks[i], leftValleyList[i], rightValleyList[i], rightValleyList[i] - leftValleyList[i], leftPeakK, rightPeakK))
    return peakObjectList


def linearJudgeNDR(smoothData, startPos, endPos):
    model = linear_model.LinearRegression()
    x = np.array([i for i in range(startPos, endPos)])[:, np.newaxis]
    y = smoothData[startPos:endPos][:, np.newaxis]
    model.fit(x, y)
    return model, x, y


def getPeaksAveHeight(wpsList, startPos, endPos):
    if endPos > len(wpsList) or startPos == endPos:
        return 0
    x = np.arange(startPos, endPos, 1)
    minNum = np.min(wpsList[startPos:endPos])
    minLine = np.array([minNum for i in range(endPos - startPos)])
    return simps(np.subtract(np.array(wpsList)[startPos:endPos], minLine), x=x)


def getAllPeaksAveHeight(peakObjectList, wpsList):
    for peak in peakObjectList:
        peak.aveHeight = getPeaksAveHeight(wpsList, peak.startPos, peak.endPos)


def getTheta(x1, y1):
    '''
    由(x2,y2)逆时针旋转到(x1,y1)
    '''
    theta = np.arctan2(y1, x1)
    if theta < 0:
        theta = 2 * np.pi + theta
    return np.degrees(theta)


def getMidVector(modelList):
    '''
    获得两个向量的和的弧度
    :param modelList:
    :return:
    '''
    k1 = modelList[0].coef_
    k2 = modelList[1].coef_
    b1 = modelList[0].intercept_
    b2 = modelList[1].intercept_
    intersectionX = (b2 - b1) / (k1 - k2)  # 交点
    intersectionY = (k2 * b1 - k1 * b2) / (k2 - k1)
    #
    point1X = 0
    point1Y = b1
    point2X = 1990
    point2Y = k2 * point2X + b2

    vectorX = point1X - intersectionX + point2X - intersectionX
    vectorY = point1Y - intersectionY + point2Y - intersectionY

    radian = getTheta(vectorX, vectorY)
    radian1 = getTheta(point1X - intersectionX, point1Y - intersectionY)
    radian2 = getTheta(point2X - intersectionX, point2Y - intersectionY)
    if abs(radian1 - radian2) >= 180:
        radianMid = (radian1 + radian2) / 2 + 180
    else:
        radianMid = (radian1 + radian2) / 2
    radianMid %= 360

    return [radianMid[0][0], radian[0][0], intersectionX[0][0], intersectionY[0][0], vectorX[0][0], vectorY[0][0]]


def judgeLowDepth(depth, startPos, endPos):
    '''
    narrow_interval_coverage 和 broad_interval_coverage,寻找平均覆盖度最小的区间，区间大小分别为 200bp(单核小体)和 600bp(三核小体)
    '''
    ndrWin = 300
    smallNdrWin = 100
    depSum = np.zeros(len(depth) + 1)
    depSum[1:] = np.cumsum(depth)  # depth[i - > j] = depSum[j + 1] - depSum[i]
    # depthList = []
    ndrAreaDepth = 1000000
    smallNdrAreaDepth = 1000000

    startPos = max(startPos, ndrWin)
    endPos = min(endPos, len(depth) - ndrWin)
    minIndex = startPos
    for i in range(startPos, endPos):
        if depSum[i + ndrWin] - depSum[i - ndrWin] < ndrAreaDepth:
            ndrAreaDepth = depSum[i + ndrWin] - depSum[i - ndrWin]
            minIndex = i
        smallNdrAreaDepth = min(
            depSum[i + smallNdrWin] - depSum[i - smallNdrWin], smallNdrAreaDepth)
    return ndrAreaDepth / (2 * ndrWin), smallNdrAreaDepth / (2 * smallNdrWin), minIndex


def getTriangleArea(point0, point1, point2):
    '''
    三角形面积
    '''
    line0 = getPointDis(point0, point1)
    line1 = getPointDis(point0, point2)
    line2 = getPointDis(point1, point2)
    p = (line0 + line1 + line2) / 2
    return (p*(p - line0)*(p - line1)*(p - line2)) ** 0.5


def getPointDis(point0, point1):
    return ((point0[0] - point1[0]) ** 2 + (point0[1] - point1[1]) ** 2) ** 0.5


def haveNearContinuouslyPeak(smoothData, rawDataList, peakObjectList):
    '''
    :param smoothData:
    :param rawDataList:
    :param peakObjectList:
    :return: meanWidth, meanHeight, meanAngel, meanArea, varWidth, varHeight, varAngel, varArea
    '''
    peakWidth = []
    peakHeight = []
    peakAngel = []
    peakArea = []
    nearPeakX = []
    peakStart = []
    peakEnd = []
    peakTriPoints = []
    for peak in peakObjectList:
        nearPeakX.append(peak.peakIndex)
        peakStart.append(peak.startPos)
        peakEnd.append(peak.endPos)
        peakWidth.append(peak.width)
        peakHeight.append(
            rawDataList[peak.peakIndex] - min(rawDataList[peak.startPos], rawDataList[peak.endPos]))
        if 1 + peak.leftK * peak.rightK != 0:  # 两直线不垂直
            peakAngel.append(math.atan(
                abs((peak.leftK - peak.rightK) / (1 + peak.leftK * peak.rightK))) * 180 / 3.1415)
            # 角度
        else:
            peakAngel.append(90)
        peakArea.append(getTriangleArea([float(peak.startPos), smoothData[peak.startPos]],
                                        [float(peak.peakIndex),
                                         smoothData[peak.peakIndex]],
                                        [float(peak.endPos), smoothData[peak.endPos]]) / 100)
        peakTriPoints.append([(float(peak.startPos), smoothData[peak.startPos]),
                              (float(peak.peakIndex),
                               smoothData[peak.peakIndex]),
                              [float(peak.endPos), smoothData[peak.endPos]]])
    origin = {
        "indexs": nearPeakX,
        "startPos": peakStart,
        "data": smoothData,
        "endPos": peakEnd,
        "width": peakWidth,
        "height": peakHeight,
        "angel": peakAngel,
        "area": peakArea,
        "triPoints": peakTriPoints
    }
    peakWidth = np.array(peakWidth)
    peakHeight = np.array(peakHeight)
    peakAngel = np.array(peakAngel)
    peakArea = np.array(peakArea)

    meanWidth = np.mean(peakWidth)
    meanHeight = np.mean(peakHeight)
    meanAngel = np.mean(peakAngel)
    meanArea = np.mean(peakArea)

    varWidth = np.var(peakWidth)
    varHeight = np.var(peakHeight)
    varAngel = np.var(peakAngel)
    varArea = np.var(peakArea)
    # print('haveNearContinuouslyPeak done')
    return meanWidth, meanHeight, meanAngel, meanArea, varWidth, varHeight, varAngel, varArea, origin


def modelType(models):
    if models[0].coef_ <= 0 and models[1].coef_ >= 0:
        return 1.0
    elif models[0].coef_ <= 0 and models[1].coef_ <= 0:
        return 2.0
    elif models[0].coef_ >= 0 and models[1].coef_ >= 0:
        return 3.0
    else:
        return 4.0


if __name__ == "__main__":
    pass
    # Test bedIterator
    # for i in bedIterator(['/Source/bed/ENCFF479AQG.bed','/Source/bed/ENCFF479AQG.bed'],1):
    #     print(i.chr,i.start,i.end)

    # Test callOneBed
    # logger.info("Test callOneBed")
    # bamfileList = ['/Source/fastq/IA01.bam']
    # rounds = 0
    # s = e = 0
    # conting = 1
    # peaksList = []
    # peakWdith = []
    # peakDis = []

    # features_df = pd.DataFrame(columns=('index', "depth", "bigDepth", "smallDepth", "height", "width", "angel", "area", "leftK", "rightK", "var_Hight", "var_Width",
    #                            "var_Angel", "var_Area", "leftKVar", "rightKVar", "mTpye", "m0k", "m1k", "m2k", "m0b", "m1b", "m2b", "vectors_", "conting_", "start_", "end_"))

    # for bed in bedIterator(['/Source/OCRDetector/result/TSS.bed'], conting):
    #     logger.info("*** Round %d ***" % rounds)
    #     rounds += 1
    #     step = bed[1] - bed[2]
    #     peaksList = []
    #     length = step + 1
    #     wpsList_Nor, lFdepth_Nor, sFdepth_Nor = callOneBed(
    #         bamfileList, str(bed[0]), bed[1], bed[2], win=120)
    #     # np.savetxt("wps_Nor.txt", wpsList_Nor)
    #     adjustWpsList_Nor = AdjustWPS(wpsList_Nor)
    #     # np.savetxt("wps_Nor_adjust.txt", adjustWpsList_Nor)
    #     base = getBaseline(adjustWpsList_Nor)
    #     smoothWpsList_Nor = smoothWPS(adjustWpsList_Nor, base)
    #     lFdepth_Nor = savgol_filter(lFdepth_Nor, 801, 2)

    #     peakHeight = []

    #     peaks, properties = find_peaks(
    #         smoothWpsList_Nor, height=0.28, distance=25, prominence=0.25, width=[25, 170])
    #     peakObjectList = getValley(
    #         smoothWpsList_Nor, adjustWpsList_Nor, peaks, 5)
    #     peaksList.append([adjustWpsList_Nor, peaks, peakObjectList])
    #     peaksList.append([smoothWpsList_Nor, peaks, peakObjectList])

    #     try:
    #         getAllPeaksAveHeight(peakObjectList, smoothWpsList_Nor)
    #     except:
    #         logger.warn("getALLPeaksAveHeight error")

    #     # 左侧2/3 右侧2/3 和 全部覆盖度曲线的线性模型
    #     sublen = int(len(lFdepth_Nor) / 3)
    #     models = []

    #     model, new_x, new_y = linearJudgeNDR(lFdepth_Nor, 0 , 2 * sublen - 1)
    #     if model.coef_ == 0:
    #         logger.info("model.coef_ == 0")
    #         continue
    #     models.append(model)

    #     model, new_x, new_y = linearJudgeNDR(lFdepth_Nor, sublen , 3 * sublen - 1)
    #     if model.coef_ == 0:
    #         logger.info("model.coef_ == 0")
    #         continue
    #     models.append(model)

    #     model, new_x, new_y = linearJudgeNDR(lFdepth_Nor, 0, len(lFdepth_Nor) - 1)
    #     if model.coef_ == 0:
    #         logger.info("model.coef_ == 0")
    #         continue
    #     models.append(model)

    #     vectors = getMidVector(models)

    #     meanWidth, meanHeight, meanAngel, meanArea, varWidth, varHeight, varAngel, varArea = haveNearContinuouslyPeak(
    #         smoothWpsList_Nor, adjustWpsList_Nor, peakObjectList)
    #     ndrAreaDepth, smallNdrAreaDepth, minIndex = judgeLowDepth(
    #         lFdepth_Nor, 0, len(lFdepth_Nor) - 1)

    #     depth = np.sum(lFdepth_Nor) / length
    #     bigDepth = ndrAreaDepth
    #     smallDepth = smallNdrAreaDepth
    #     height = meanHeight
    #     width = meanWidth
    #     angel = meanAngel
    #     area = meanArea
    #     leftK = np.mean([p.leftK for p in peakObjectList])
    #     rightK = np.mean([p.rightK for p in peakObjectList])
    #     var_Hight = varHeight
    #     var_Width = varWidth
    #     var_Angel = varAngel
    #     var_Area = varArea
    #     leftKVar = np.var([p.leftK for p in peakObjectList])
    #     rightKVar = np.var([p.rightK for p in peakObjectList])
    #     mTpye = modelType(models)
    #     m0k = models[0].coef_
    #     m1k = models[1].coef_
    #     m2k = models[2].coef_
    #     m0b = models[0].intercept_
    #     m1b = models[1].intercept_
    #     m2b = models[2].intercept_
    #     vectors_ = vectors
    #     conting_ = int(conting)
    #     start_ = bed[1]
    #     end_ = bed[2]

    #     bed_df = pd.DataFrame([[rounds, depth, bigDepth, smallDepth, height, width, angel, area, leftK, rightK, var_Hight, var_Width, var_Angel, var_Area, leftKVar, rightKVar, mTpye, m0k, m1k, m2k, m0b, m1b, m2b, vectors_, conting_, start_, end_]], columns=(
    #         'index', "depth", "bigDepth", "smallDepth", "height", "width", "angel", "area", "leftK", "rightK", "var_Hight", "var_Width", "var_Angel", "var_Area", "leftKVar", "rightKVar", "mTpye", "m0k", "m1k", "m2k", "m0b", "m1b", "m2b", "vectors_", "conting_", "start_", "end_"))

    #     features_df = pd.concat([features_df, bed_df], ignore_index=True)

    # features_df.to_csv("features_df.csv", sep=',', header=True, index=True)
