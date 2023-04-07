import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn.metrics import confusion_matrix
import scipy.sparse as sp
import pandas as pd
import sys
import os

def ROC(y_true, y):
    tn, fp, fn, tp = metrics.confusion_matrix(y_true.toarray().ravel(), y.toarray().ravel()).ravel()
    if not(fp + tn == 0):
        fpr = fp / (fp + tn)
    else:
        fpr = 0

    if not(tp + fn == 0):
        tpr = tp / (tp + fn)
    else:
        tpr = 0
    xy_arr = [fpr, tpr]

    if not(tp + fp == 0):
        precision = tp / (tp + fp)
    else:
        precision = 0

    if not (tp + fn == 0):
        recall = tp / (tp + fn)
    else:
        recall = 0
    xy_arr2 = [recall, precision]

    return fpr, tpr, xy_arr, recall, precision, xy_arr2


def auc(xy):
    au = 0.
    prev_x = 0
    prev_y = 0
    for x,y in xy:
        if x != prev_x:
            au += (x - prev_x) * (y + prev_y) / 2
            prev_x = x
            prev_y = y
    return au


def aupr(xy):
    au = 0.
    prev_x = 0
    prev_y = 1
    for x,y in xy:
        if x != prev_x:
            au += (x - prev_x) * (y + prev_y) / 2
            prev_x = x
            prev_y = y
    return au


def init():
    fprs = [1]
    tprs = [1]
    xy_arrs = [[1, 1]]
    recalls = [0]
    precisions = [1]
    xy_arr2s = [[0, 1]]
    return fprs, tprs, xy_arrs, recalls, precisions, xy_arr2s


def groundTruthtoMatrix(groundTruthFile):
    row = []
    col = []
    data = []
    print("+++ Load groundtruth file from %s"% groundTruthFile)
    with open(groundTruthFile) as f:
        for line in f:
            numbers = line.split()
            row.append(int(numbers[0]))
            col.append(int(numbers[1]))
            data.append(1)
    global matrixSize
    sparse_matrix = sp.csc_matrix((data, (row, col)), shape=(matrixSize, matrixSize))
    print("=== Finally, %s is made to sparse-matrix."% groundTruthFile)
    return sparse_matrix


def datatoMatrix(mode, assessedFile):
    row = []
    col = []
    data = []
    if not(os.path.exists(assessedFile)):
        print("!!! Error: %s is not available" % assessedFile)
        sys.exit()
    else :
        print("+++ Load assessed file from %s" % assessedFile)
    with open(assessedFile) as f:
        for line in f:
            numbers = line.split()
            row.append(int(numbers[0]))
            col.append(int(numbers[1]))
            if(mode == "LSH") :
                data.append(float(numbers[2]))
            else :
                data.append(int(numbers[2]))
    sparse_matrix = sp.csc_matrix((data, (row, col)), shape=(matrixSize, matrixSize))
    return sparse_matrix


def modeWork(mode, y_true, add):
    fprs, tprs, xy_arrs, recalls, precisions, xy_arr2s = init()
    global simList
    if(mode == "LSH") :
        csc = datatoMatrix(moed, add)
        for i in simList:
            new_value = []
            new_row = []
            new_col = []
            for value, row, col in zip(csc.data, csc.indices, csc.indptr[1:]):
                if(value >= i):
                    new_row.append(row)
                    new_col.append(col)
                    new_value.append(1)

            new_csc = sp.csc_matrix((data, (row, col)), shape=(matrixSize, matrixSize))
            fpr, tpr, xy_arr, recall, precision, xy_arr2 = ROC(y_true, new_csc)
            fprs.insert(0, fpr)
            tprs.insert(0, tpr)
            xy_arrs.insert(0, xy_arr)
            recalls.insert(1, recall)
            precisions.insert(1, precision)
            xy_arr2s.insert(1, xy_arr2)
    else :
        for i in simList:
            name = "{:.2f}".format(i)
            address = add + name + ".txt"
            fpr, tpr, xy_arr, recall, precision, xy_arr2 = ROC(y_true, datatoMatrix(mode, address))
            fprs.insert(0, fpr)
            tprs.insert(0, tpr)
            xy_arrs.insert(0, xy_arr)
            recalls.insert(1, recall)
            precisions.insert(1, precision)
            xy_arr2s.insert(1, xy_arr2)

    fprs.insert(0, 0)
    tprs.insert(0, 0)
    xy_arrs.insert(0, [0, 0])
    recalls.append(1)
    precisions.append(0)
    xy_arr2s.append([1, 0])
    dict = {}
    dict[mode] = [fprs, tprs, xy_arrs, recalls, precisions, xy_arr2s]
    return dict


def logSave(logFile, dict):
    with open(logFile, "w") as f:
        print("--- The values of the confusion matrix are being saved...")
        for key, value in dict.items():
            if (key != "LSH"):
                f.write(key + '\n')
                for i in range(len(value[0])):
                    f.write('{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n'
                            .format(value[0][i], value[1][i], value[3][i], value[4][i]))


def logLoad(logFile):
    dirc = {}
    if not(os.path.exists(logFile)):
        print("!!! Error: %s is not available" % logFile)
        sys.exit()
    else:
        print("+++ Load past data file from %s" % logFile)

    with open(logFile) as f:
        for line in f:
            if(line == "Linclust\n" or line == "ALFATClust\n"):
                label = line.strip()
                fprs = []
                tprs = []
                xy_arrs = []
                recalls = []
                precisions = []
                xy_arr2s = []
            else:
                numbers = line.split('\t')
                fprs.append(float(numbers[0]))
                tprs.append(float(numbers[1]))
                xy_arrs.append([float(numbers[0]), float(numbers[1])])
                recalls.append(float(numbers[2]))
                precisions.append(float(numbers[3]))
                xy_arr2s.append([float(numbers[2]), float(numbers[3])])
            dirc[label] = [fprs, tprs, xy_arrs, recalls, precisions, xy_arr2s]

    print("+++ Past data file is loaded!")
    return dirc


def main(ParameterList):
    global simList
    simList = [0.5, 0.6, 0.7, 0.8, 0.9]
    ListLength = len(ParameterList)
    if(ListLength < 7 or ListLength > 8):
        print("!!! Error: Incorrect number (%d) of input parameters" %(ListLength-1))
        sys.exit()

    ParameterNum = int(ParameterList[1])
    global matrixSize
    matrixSize = int(ParameterList[2])

    groudthFile = ParameterList[3]
    if not(os.path.exists(groudthFile)):
        print("!!! Error: GroudthFile is not available")
        sys.exit()
    else :
        y_true = groundTruthtoMatrix(groudthFile)

    labelList = ParameterList[4].split(',')
    if (len(labelList) != ParameterNum):
        print("!!! Error: Incorrect number of Label parameters")
        sys.exit()

    inputFileList = ParameterList[5].split(',')
    if (len(inputFileList) != ParameterNum):
        print("!!! Error: Incorrect number of inputFile parameters")
        sys.exit()

    outputFile = ParameterList[6]

    Matrics = {}

    if (ListLength == 8):
        logFile = ParameterList[7]
        Matrics.update(logLoad(logFile))

    for label, inputFile in zip(labelList, inputFileList):
        Matrics.update(modeWork(label, y_true, inputFile))

    print("--- ROC curve being plotted...")
    colorList = ['#F3D266', '#3979F2', '#999999']
    progress = 0
    for key in Matrics.keys():
        plt.plot(Matrics[key][0], Matrics[key][1], '.-', color = colorList[progress],
                 label = u"AUC in %s=%0.3f" %(key, auc(Matrics[key][2])))
        progress += 1

    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], '--', color='#FC8A61')
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.grid(linestyle='-.')
    plt.grid(True)
    plt.savefig(outputFile+"_ROC.png", dpi=500)
    plt.close()

    print("--- PR curve being plotted...")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot([0, 1], [0, 1], '--', color='#FC8A61')
    progress = 0
    for key in Matrics.keys():
        plt.plot(Matrics[key][3], Matrics[key][4], '.-', color = colorList[progress],
                 label = u"AUPR in %s=%0.3f" %(key, aupr(Matrics[key][5])))
        progress += 1

    plt.legend(loc='lower left')
    plt.grid(linestyle='-.')
    plt.grid(True)
    plt.savefig(outputFile+"_PR.png", dpi=500)

    logSave(outputFile + "_log.txt", Matrics)


main(sys.argv)