#!/usr/bin/python

import re
import matplotlib.pyplot as plt
import numpy as np

with open("genome.txt", "r") as gfile:
    dnaString = gfile.readline()
gfile.close()

with open("read.txt", "r") as infile:
    dnaArray = []
    for line in infile:
        dnaArray.append(line.strip('\n'))
infile.close()

countArray = [0] * len(dnaString)
for substr in dnaArray:
    regex = re.compile(substr)
    for match in regex.finditer(dnaString):
        countArray[match.start()] += 1

N = len(dnaString)
x = range(N)
width = 1/1.5
plt.bar(x, countArray, width, color="blue")

plt.title("Read count Histogram")
plt.xlabel("Position")
plt.ylabel("Frequency")

plt.show()
