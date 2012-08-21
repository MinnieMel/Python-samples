#-------------------------------------------------------------------------------
# Name:        kNN classifier
# Purpose:
#
# Author:      Michael
#
# Created:     02/12/2011
# Copyright:   (c) Michael 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
import math
from operator import itemgetter

AML_INDEX = 0
ALL_INDEX = 1

k = 3
training_offset = 27
testing_offset = 21

# takes two columns
def euclidean_distance(a, b):
    if len(a) != len(b):
        print 'Dimension length of a and b do not match!'
    sum = 0
    for x in range(len(a)):
        sum = sum + math.pow(a[x]-b[x], 2)
    return math.sqrt(abs(sum))


def kNN(sample, training_data, training_offset):
    distances = {}
    for x in range(len(training_data)):
        distances[x] = 0
    for idx, column in enumerate(training_data):
        distances[idx] = euclidean_distance(sample, column)
    sorted_distances = sorted(distances.items(), key=itemgetter(1))
    #print sorted_distances
    aml_score = 0
    all_score = 0
    for x in range(k):
        #print sorted_distances[x]
        if sorted_distances[x][0] < training_offset:
            aml_score = aml_score + 1
        else:
            all_score = all_score + 1
    if aml_score > all_score:
        return AML_INDEX
    else:
        return ALL_INDEX

def parse(training, testing):
    training_data = readDataFromFile(training)
    testing_data = readDataFromFile(testing)
    print '=== Classifying testing data using training data ==='
    print 'Testing Errors: ' + repr(classifyData(training_data, testing_data, training_offset, testing_offset))
    print '=== Classifying training data using training data ==='
    print 'Training Errors: ' + repr(classifyData(training_data, training_data, training_offset, training_offset))

def getSampleLabel(index, offset):
    if index < offset:
        return AML_INDEX
    else:
        return ALL_INDEX

def classifyData(training, testing, training_offset, testing_offset):
    error_count = 0
    for idx, sample in enumerate(testing):
        prediction = kNN(sample, training, training_offset)
        if (prediction != getSampleLabel(idx, testing_offset)):
            print 'Error Column: ' + repr(idx + 1)
            error_count = error_count + 1
    return error_count

def readDataFromFile(filepath):
    f = open(filepath, 'r')
    data = []
    for line in f:
        vectors = line.split()
        if len(data) != len(vectors):
            for x in range(len(vectors)):
                data.append([])
        for idx, vector in enumerate(vectors):
            data[idx].append(int(vector))
    f.close()
    return data

def main():
    if (len(sys.argv) == 3):
        parse(sys.argv[1], sys.argv[2])
    else:
        training = raw_input("Specify input filepath for training data:")
        testing = raw_input("Specify input filepath for testing data:")
        parse(training, testing)
    pass

if __name__ == '__main__':
    main()