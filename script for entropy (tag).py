##This is the script to calculate entropy values from tagged txt files

import math
from collections import Counter
from typing import List


def freqCount(tokens):
    freqdict = {}
    for word in tokens:
        freqdict[word] = freqdict.get(word, 0) + 1
    freqlist = sorted(freqdict.items(), key=lambda e: e[1], reverse=True)
    fv = []
    for x, y in freqlist:
        fv.append(y)
    return fv

def TTR(freqlist):
    sumfreq = sum(freqlist)
    types = len(freqlist)
    return types / sumfreq

def calentropy(freqlist):
    sumfreq = sum(freqlist)
    sumH = 0
    for fl in freqlist:
        if fl > 0:
            p = fl / sumfreq
            sumH += -p * math.log(p, 2)
    return sumH

def calentropy_rel(freqlist):
    sumfreq = sum(freqlist)
    sumH = 0
    for fl in freqlist:
        if fl > 0:
            p = fl / sumfreq
            sumH += -p * math.log(p, 2)
    V = len([v for v in freqlist if v > 0])
    return sumH / math.log(V, 2)

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

txtPath = r'your\path\to\tagged\txt'

txtPath = r'C:\D\Names25\618\cntg.txt'
txtPath = r'C:\D\Names25\618\ustg.txt'

resd = dict()
rowpatterns = []

with open(txtPath,'r') as ot:
    lines = ot.readlines()
    for line in lines[:-1]:
        rowpattern = ''

        words = line.split('\t')
        for w in words:
            if w.startswith('#'):
                rowpattern += w[1:] + '+'
                if w not in resd.keys():
                    resd[w] = 1
                else:
                    resd[w] += 1
        rowpattern = rowpattern.rstrip('+')
        rowpatterns.append(rowpattern)

print("***************************")
fv = freqCount(rowpatterns)

print("Entropy:", calentropy(fv))
print("Relative Entropy:", calentropy_rel(fv))
print("Zhang Entropy:", zhang_entropy(fv))



