#-------------------------------------------------------------------------------
# Name:        Viterbi algorithm
# Purpose:
#
# Author:      Michael Vu
#
# Created:     19/11/2011
# Copyright:   (c) Mike 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
import random
import math

p = 0.98
q = 0.999

alphabet = ('A', 'C', 'G', 'T')

states = (0, 1, 2, 3, 4, 5, 6, 7)

trans = ((0.180*p, 0.274*p, 0.426*p, 0.120*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4),
        (0.171*p, 0.368*p, 0.274*p, 0.188*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4),
        (0.161*p, 0.339*p, 0.375*p, 0.125*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4),
        (0.079*p, 0.355*p, 0.384*p, 0.182*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4),
        ((1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.300*q, 0.205*q, 0.285*q, 0.210*q),
        ((1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.322*q, 0.298*q, 0.078*q, 0.302*q),
        ((1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.248*q, 0.246*q, 0.298*q, 0.208*q),
        ((1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.177*q, 0.239*q, 0.292*q, 0.292*q))

def readSeqsFromFile(filepath):
    f = open(filepath, 'r')
    sqs = [];
    for line in f:
        if line[:1] == '>':
            sqs.append('')
        elif(len(sqs) == 0):
            continue
        else:
            sqs[-1] += line[:len(line)-1]
    f.close()
    #print sqs
    return sqs

def getStates(x):
    states = []
    for i, y in enumerate(alphabet):
        if (y == x):
            states.append(i)
            states.append(i+4)
            return states

def viterbi(seq,trans):
    V = [{}]
    path = {}
    islands = []
    island = None
    L = len(seq)

    states = getStates(seq[0])
    for y in states:
        V[0][y] = math.log(0.5)
        path[y] = [y]
        #path.append([states[idx]])

    for t in range(1, L):
        V.append({})
        newpath = {}

        next_states = getStates(seq[t])
        for next in next_states:
            probs = {}
            prev_states = getStates(seq[t-1])
            for prev in prev_states:
                prob = V[t-1][prev] + math.log(trans[prev][next])
                probs[prev] = prob
            top = max(probs, key = probs.get)
            V[t][next] = probs[top]
            newpath[next] = path[top] + [next]

        path = newpath

    (prob, state) = max([(V[L - 1][y], y) for y in path.keys()])

    print 'Optimal alignment score (log): ' + repr(prob)
    #check = math.log(0.5)
    #for t in range(1, len(path[state])):
    #    check = check + math.log(trans[path[state][t-1]][path[state][t]])
    #print 'Post-score (log): ' + repr(check)

    print 'Optimal path:'
    path_repr = ''
    for idx, x in enumerate(path[state]):
        if (island == None) & (x < 4):
            island = [idx]
        elif (island != None) & (x >= 4):
            island.append(idx - 1)
            islands.append(island)
            island = None
        path_repr = path_repr + repr(x)
    if (island != None):
        island.append(len(path[state]) - 1)
        islands.append(island)
    print path_repr

    print 'Positions of CpG islands inferred by Viterbi algorithm:'
    for island in islands:
        print repr(island[0]) + ' ' + repr(island[1])


def parse(filepath):
    for idx, sq in enumerate(readSeqsFromFile(filepath)):
        print '>seq ' + repr(idx)
        viterbi(sq, trans)


def main():
    if (len(sys.argv) > 1):
        parse(sys.argv[1])
    else:
        parse(raw_input("Specify input filepath:"))
    pass

if __name__ == '__main__':
    main()