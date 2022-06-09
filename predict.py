import argparse
import copy
import random
from os.path import join as pjoin
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cleanlab.classification import CleanLearning
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, classification_report,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from utils.log import get_logger

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description="prediction")
parser.add_argument("csv", help="combined csv file")
parser.add_argument('-o', '--output', help='output dir',
                    dest='output', default='./build/')
parser.add_argument('-f', '--fraction', help='fraction',
                    dest='fraction', default=1.0, type=float)
args = parser.parse_args()

output_dir = args.output
fraction = args.fraction

logger = get_logger("predict", file_handler=pjoin(
    output_dir, "..", "prediction.log"))


def drawRes(x, y, xNoise, yNoise, title):
    plt.plot(x, y)
    plt.plot(xNoise, yNoise)
    plt.xlabel('Noise Level')
    plt.ylabel('Accuracy')
    plt.ylim(0, 1.0)
    plt.title(title)
    plt.legend(['Normal', 'Confident Learning'])
    # plt.show()
    plt.savefig(pjoin(output_dir, title + '.png'))
    logger.info("save to {}".format(
        pjoin(output_dir, 'figures', title + '.png')))
    plt.close()


def drawResNormol(x, y, title):
    plt.plot(x, y)
    plt.ylim(0, 1.0)
    plt.title(title)
    # plt.show()
    plt.savefig(pjoin(output_dir, 'figures', title + '.png'))
    plt.close()


def transferYPred(y_pred):
    print(y_pred)
    for i in range(len(y_pred)):
        max_value = max(y_pred[i])
        for j in range(len(y_pred[i])):
            if max_value == y_pred[i][j]:
                y_pred[i][j] = 1
            else:
                y_pred[i][j] = 0
    print(y_pred)
    return y_pred


def getResult(testX, testY, clf, title):
    y_pred = clf.predict(testX)
    y_prob = clf.predict_proba(testX)
    accscore = str(accuracy_score(testY, y_pred))[0:6]
    aucscore = str(roc_auc_score(testY, y_prob, multi_class='ovr'))[0:6]
    recallsocre = str(recall_score(testY, y_pred, average='weighted'))[0:6]
    f1score = str(f1_score(testY, y_pred, average='weighted'))[0:6]

    logger.info("\n-------- {} --------\nauc: {}, acc: {}, recall: {}, f1: {}".format(
        title, aucscore, accscore, recallsocre, f1score))

    report = classification_report(testY, y_pred)
    logger.info("\n"+report)
    conf_matrix = confusion_matrix(testY, y_pred)

    return [float(accscore), float(aucscore), float(recallsocre), float(f1score)]


def main():
    data = pd.read_csv(args.csv)
    df_result = pd.DataFrame(columns=(
        "fraction", "isCL", "noiseRatio", "accscore", "aucscore", "recallsocre", "f1score"))

    # depth  bigDepth  smallDepth
    # m0k       m1k       m2k       m0b       m1b       m2b
    X = np.hstack((data.iloc[:, 1:4].values, data.iloc[:, 17:23].values))
    Y = data.iloc[:, 0].values
    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=int)

    # split data into train and test
    (trainX, testX, trainY, testY) = train_test_split(
        X, Y, test_size=0.35, random_state=150)

    # scale data
    std = StandardScaler()
    std = std.fit(trainX)
    # apply scale rule to train and test data
    trainX = std.transform(trainX)
    testX = std.transform(testX)

    # ratio of noise
    ratioList = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    scoreList = []
    noiseScoreist = []

    for ratio in ratioList:
        clf = RandomForestClassifier()
        clfNoise = CleanLearning(clf=RandomForestClassifier())
        newTrainX = trainX
        newTrainY = copy.deepcopy(trainY)
        logger.info("Noise ratio: %f" % ratio)

        for i in range(len(newTrainX)):
            if (random.random() < ratio):
                newTrainY[i] = random.randint(1, 4) - 1

        clf.fit(newTrainX, newTrainY)
        clfNoise.fit(newTrainX, newTrainY)

        # accscore, aucscore, recallsocre, f1score
        scores = getResult(testX, testY, clf, 'Normal - Ratio ' + str(ratio))
        noiseScores = getResult(testX, testY, clfNoise,
                                'CL - Ratio ' + str(ratio))
        scoreList.append(scores)
        noiseScoreist.append(noiseScores)

        df_result.loc[len(df_result)] = [fraction, False, ratio, scores[0],
                                         scores[1], scores[2], scores[3]]
        df_result.loc[len(df_result)] = [fraction, True, ratio, noiseScores[0],
                                         noiseScores[1], noiseScores[2], noiseScores[3]]

    scoreList = np.array(scoreList)
    noiseScoreist = np.array(noiseScoreist)
    titleList = ['accscore', 'aucscore', 'recallsocre', 'f1score']
    # precision (精确度)：正确预测为正的，占全部预测为正的比例。
    # recall（召回率）：正确预测为正的，占全部实际为正的比例。
    # f1-score (f1值)：精确率和召回率的调和平均数。
    # support （各分类样本的数量或测试集样本的总数量）。
    # macro avg (宏平均值)：所有标签结果的平均值。
    # weighted avg（加权平均值）：所有标签结果的加权平均值。

    df_result.to_csv(pjoin(output_dir, "..", 'result.csv'),
                     mode='a', index=False, header=False)
    for i in range(len(scoreList[0])):
        drawRes(ratioList, scoreList[:, i], ratioList, noiseScoreist[:, i], str(
            i+1) + '-' + titleList[i])


if __name__ == '__main__':
    main()
