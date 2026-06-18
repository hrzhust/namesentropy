##This is the script to calculate entropy values from raw txt files
import math
from collections import Counter
from typing import List

path = r'your\path\to\raw\txt'
path = r'C:\D\Names25\618\usraw.txt'
allnames = ''

def freqCount(tokens):
    freqdict = {}
    i = 0
    lw = len(tokens)
    for word in tokens:
        i += 1
        if word not in freqdict.keys():
            freqdict[word] = 1
        else:
            freqdict[word] += 1
    freqlist = sorted(freqdict.items(), key=lambda e: e[1], reverse=True)
    fv = []
    for x,y in freqlist:
        # print(x,y)
        fv.append(y)
    return fv

def freqCount_nofunc(tokens):
    freqdict = {}
    i = 0
    lw = len(tokens)
    for word in tokens:
        i += 1
        if word not in freqdict.keys():
            freqdict[word] = 1
        else:
            freqdict[word] += 1
    freqlist = sorted(freqdict.items(), key=lambda e: e[1], reverse=True)
    fv = []
    for x,y in freqlist:
        # print(x,y)
        fv.append(y)
    return fv

def calentropy(freqlist):
    sumfreq=sum(freqlist)
    sumH=0
    for fl in freqlist:
        if fl>0:
            p=fl/sumfreq
            sumH += -p*math.log(p,2)
    return sumH
def calentropy_rel(freqlist):
    sumfreq=sum(freqlist)
    sumH=0
    for fl in freqlist:
        if fl>0:
            p=fl/sumfreq
            sumH += -p*math.log(p,2)
    V = 0
    for value in freqlist:
        if value != 0:
            V += 1
    return sumH/math.log(V,2)


def zhang_entropy(freq: List[int]) -> float:
    if not freq:
        return 0.0
    T = sum(freq)
    freq_counts = Counter(freq)
    H_z = 0.0
    for f, n_f in freq_counts.items():
        p_hat = f / T
        inner_sum = 0.0
        prod = 1.0
        for v in range(1, T - f + 1):
            if v > 1:
                j = v - 2
                prod *= (1 + (1 - f) / (T - 1 - j))
            inner_sum += prod / v
        H_z += n_f * p_hat * inner_sum
    return H_z

with open(path,'r',encoding='utf-8',errors='ignore') as of:
    flines = of.readlines()
    for line in flines:
        # print()
        allnames += line[:-1]+' '

    # print(allnames)

words = allnames.lower().split()

fl = freqCount_nofunc(words)
print("Entropy:")
print(calentropy(fl))
print("Relative Entropy:")
print(calentropy_rel(fl))
print("Zhang\'s Estimator:")
print(zhang_entropy(fl))
