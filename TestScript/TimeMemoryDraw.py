import matplotlib.pyplot as plt
import sys
import os

def InfoAdd(label, inputFile):
    global Time, Memory
    if not(os.path.exists(inputFile)):
        print("!!!Error: %s is not available" % inputFile)
        sys.exit()
    else :
        print("+++ Load time information file from %s" % inputFile)

    seqNum = []
    runTime = []
    runMemory = []
    with open(inputFile) as f:
        for line in f:
            numbers = line.split('\t')
            seqNum.append(float(numbers[0]) / 1000)
            runTime.append(float(numbers[1]))
            runMemory.append(float(numbers[2]) / 1024)

    data = list(zip(seqNum, runTime, runMemory))
    data_sorted = sorted(data, key=lambda x: x[0])
    seqNum_sorted, runTime_sorted, runMemory_sorted = zip(*data_sorted)

    Time[label] = [seqNum_sorted, runTime_sorted]
    Memory[label] = [seqNum_sorted, runMemory_sorted]


def TimeDraw(outputFile):
    global Time
    print("--- Runtime versus input set size on linear scales being plotted...")
    plt.figure(figsize=(30, 15))
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.title("Runtime versus input set size on linear scales")
    plt.xlabel("Count of sequences in thousands")
    plt.ylabel("Runtime in seconds")
    plt.yscale('symlog')

    colorList = ['#F3D266', '#3979F2', '#999999']
    pointList = ['o-', 's--', '^:', 'v-.', 'p--']
    progress = 0
    for key in Time.keys():
        plt.plot(Time[key][0], Time[key][1], pointList[progress%5], color=colorList[int(progress/5)], label=key)
        progress += 1
    #plt.legend(frameon='true', fontsize=7, loc='upper right', bbox_to_anchor=(0.93, 1.05))
    plt.legend(frameon='true', fontsize=7)
    plt.savefig(outputFile+"/Runtime.png", dpi=500)


def MemoryDraw(outputFile):
    global Memory
    print("--- Runmemory versus input set size on linear scales being plotted...")
    plt.figure(figsize=(30, 15))
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.title("Runmemory versus input set size on linear scales")
    plt.xlabel("Count of sequences in thousands")
    plt.ylabel("Maximum resident set size in Mbytes")

    colorList = ['#F3D266', '#3979F2', '#999999']
    pointList = ['o-', 's--', '^:', 'v-.', 'p--']
    progress = 0
    for key in Time.keys():
        plt.plot(Memory[key][0], Memory[key][1], pointList[progress%5], color=colorList[int(progress/5)], label=key)
        progress += 1
    #plt.legend(loc = 'center right', frameon=False, fontsize=7)
    plt.legend(frameon=False, fontsize=7)
    plt.savefig(outputFile + "/Runmemory.png", dpi=500)


def main(ParameterList):
    ListLength = len(ParameterList)
    if (ListLength != 5):
        print("!!!Error: Incorrect number (%d) of input parameters" % (ListLength - 1))
        sys.exit()

    ParameterNum = int(ParameterList[1])

    labelList = ParameterList[2].split(',')
    if (len(labelList) != ParameterNum):
        print("!!!Error: Incorrect number of Label parameters")
        sys.exit()

    inputFileList = ParameterList[3].split(',')
    if (len(inputFileList) != ParameterNum):
        print("!!!Error: Incorrect number of inputFile parameters")
        sys.exit()

    OutputFile = ParameterList[4]

    global Time, Memory
    Time = {}
    Memory = {}
    for label, inputFile in zip(labelList, inputFileList):
        InfoAdd(label, inputFile)

    TimeDraw(OutputFile)
    plt.close()
    MemoryDraw(OutputFile)


main(sys.argv)