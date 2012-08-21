#-------------------------------------------------------------------------------
# Name:        Posterior decoding algorithm
# Purpose:
#
# Author:      Mike
#
# Created:     21/11/2011
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

# to prevent arithmetic underflow
def logplus(p, q):
    if p == float('-inf'):
        return q
    if q == float('-inf'):
        return p
    x = max(p, q)
    if (x == p):
        diff = q - p
    else:
        diff = p - q
    if diff < -50:
        return x + 0
    else:
        return x + math.log(1+math.exp(diff))

def forward(seq, trans):
    F = [{}]
    L = len(seq)

    states = getStates(seq[0])
    for y in states:
        F[0][y] = math.log(0.5)

    for t in range(1, L):
        F.append({})

        next_states = getStates(seq[t])
        for next in next_states:
            prev_states = getStates(seq[t-1])
            prob = float('-inf')
            for prev in prev_states:
                prob = logplus(prob, F[t-1][prev] + math.log(trans[prev][next]))
            F[t][next] = prob
        #print repr(t) + ' : ' + repr(F[t])

    return F

def backward(seq, trans):
    B = []
    L = len(seq)
    for x in range(L+1):
        B.append({})

    init_states = getStates(seq[L-1])
    for k in init_states:
        B[L][k] = math.log(1)
    prob = 0
    for t in range(L-1, 0, -1):
        next_states = getStates(seq[t-1])
        for next in next_states:
            prev_states = getStates(seq[t])
            prob = float('-inf')
            for prev in prev_states:
                prob = logplus(prob, B[t+1][prev] + math.log(trans[prev][next]))
            B[t][next] = prob
        #print repr(t) + ' : ' + repr(B[t])
    B.pop(0)
    #print B
    return B

def otherState(state):
    if state < 4:
        return state + 4
    else:
        return state - 4

def posterior(seq, trans):
    print 'Posterior prob of state variables:'
    f = forward(seq, trans)
    b = backward(seq, trans)
    path = ''
    islands = []
    island = None
    for idx, x in enumerate(seq):
        states = getStates(x)
        probs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        for state in states:
            probs[state] = f[idx][state] + b[idx][state]
        (score, state) = min([(x, state) for state, x in enumerate(probs)])
        score = math.exp(score - probs[otherState(state)])
        score = score/(1+score)
        otherScore = 1 - score
        probs[state] = score
        probs[otherState(state)] = otherScore
        print repr(idx) + ': ' + repr(probs)
        path = path + repr(otherState(state))
        if (island == None) & (otherState(state) < 4):
            island = [idx]
        elif (island != None) & (otherState(state) >= 4):
            island.append(idx - 1)
            islands.append(island)
            island = None
    if (island != None):
        island.append(len(seq) - 1)
        islands.append(island)
        island = None
    print 'Posterior path:'
    print path
    print 'CpG islands inferred by posterior decoding algorithm:'
    for island in islands:
        print repr(island[0]) + ' ' + repr(island[1])

def parse(filepath):
    for idx, sq in enumerate(readSeqsFromFile(filepath)):
        print '>seq ' + repr(idx)
        posterior(sq, trans)

def main():
    if (len(sys.argv) > 1):
        parse(sys.argv[1])
    else:
        parse(raw_input("Specify input filepath:"))
    pass

if __name__ == '__main__':
    main()
