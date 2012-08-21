#-------------------------------------------------------------------------------
# Name:       CS284A/184A Problem Set 2
# Purpose:  Sequence Alignment
#
# Author:      Michael Vu
#
# Created:     22/10/2011
# Copyright:   (c) Mike 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import sys


# matrix[row][column] == matrix[b][a]
# b = row, a = column
def alignSequences(a, b):
    # Build matrix of scores. 2 points for equal characters, -1 for different
    matrix = [2 if i == j else -1 for i in b for j in a]
    for x in range(0, len(matrix), len(a)):
        row = [matrix.pop() for i in range(0, len(a))]
        # each row must be reversed because pop() is LIFO
        row.reverse()
        matrix.insert(0, row)
    # calculate optimal alignment scores for 1st row
    prev = 0
    for column, x in enumerate(matrix[0]):
        matrix[0][column] = x + prev
        prev = matrix[0][column]
    # calculate optimal alignment scores for 1st column
    prev = 0
    for row, x in enumerate(matrix):
        matrix[row][0] = x[0] + prev
        prev = matrix[row][0]

    # diagonal parent
    d = lambda row, column: matrix[row-1][column-1]
    # horizontal parent
    h = lambda row, column: matrix[row][column-1]
    # vertical parent
    v = lambda row, column: matrix[row-1][column]


    # calculate optimal alignment scores row-by-row
    prev = 0
    for row, x in enumerate(matrix[1:]):
        for column, value in enumerate(matrix[row+1][1:]):
            matrix[row+1][column+1] = max(d(row+1, column+1), h(row+1, column+1), v(row+1, column+1)) + value

    # matrix completed
    #for row in matrix:
    #    print repr(row)
    score = matrix[len(b)-1][len(a)-1]

    # helper method to find a path back to start from the last element
    # expects an empty list for path
    # recursive solution
    def trace_path(path, row, column):
        if row < 0 or column < 0:
            return path
        dv = d(row, column) if row != 0 and column != 0 else -sys.maxint - 1
        hv = h(row, column) if column != 0 else -sys.maxint - 1
        vv = v(row, column) if row != 0 else -sys.maxint - 1

        if max(dv,hv,vv) <= matrix[row][column]:
            path.append((row,column))
            return trace_path(path, row-1, column-1)
        #diagonal?
        elif dv >= hv and dv >= vv:
            path.append((row,column))
            return trace_path(path, row-1, column-1)
        #horizontal?
        elif hv >= dv and hv >= vv:
            return trace_path(path, row, column-1)
        #vertical?
        elif vv >= dv and vv >= hv:
            return trace_path(path, row-1, column)

    # compute alignment
    alignment = trace_path([], len(b)-1, len(a)-1)
    alignment.reverse()
    #print alignment

    # this block is used to format the strings with appropriate amounts of '-'
##    total_difference = 0
##    for x in alignment:
##        difference = x[1]-x[0]-total_difference
##        if difference > 0:
##            #difference = difference - total_difference
##            b = b[:x[0]] + '-'*abs(difference) + b[x[0]:]
##        elif difference < 0:
##            #difference = difference + total_difference
##            a = a[:x[1]] + '-'*abs(difference) + a[x[1]:]
##        else:
##            continue
##        print 'alignment: ' + repr(x)
##        print 'difference: ' + repr(difference)
##        total_difference = total_difference + difference
##        print 'resulting total_difference: ' + repr(total_difference)



    # to account for string extension when you add '-'
    prev = 0
    a_offset = 0
    b_offset = 0
    for x in alignment:
        difference = x[1]-x[0]-prev
        #print 'difference:' + repr(difference)
        prev = difference + prev
        if difference > 0:
            b = b[:x[0]+b_offset] + '-'*abs(difference) + b[x[0]+b_offset:]
            b_offset = b_offset + abs(difference)
        elif difference < 0:
            a = a[:x[1]+a_offset] + '-'*abs(difference) + a[x[1]+a_offset:]
            a_offset = a_offset + abs(difference)
        else:
            continue

    # top it off
    difference = len(a) - len(b)
    if difference < 0:
        a = a[:len(a)] + '-'*-difference
    if difference > 0:
        b = b[:len(b)] + '-'*difference

    length = len(a)
    matches = 0
    for x in xrange(length):
        if a[x] == b[x]:
            matches = matches + 1

    print 'Optimal alignment score: ' + repr(matches*2-(length-matches))
    print 'Length: ' + repr(length)
    print 'Nucleotide matches: ' + repr(matches)
    print 'Similarity measure: ' + repr(matches/length)

    print a
    print b

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
    return sqs

def parse(filepath):
    sqs = readSeqsFromFile(filepath)
    p = sqs[0]
    for i, sq in enumerate(sqs[1:]):
        print '=== #1 and #' + repr(i+2) + ' ==='
        alignSequences(p, sq)

def main():
    # A recursive function is used to trace a path through the matrix.
    # Since the matrix can get very large, more recursive calls are needed to
    # the point where it may exceed the system's recursion limit.
    # This is a hacky workaround. Maybe in the future an iterative algorithm
    # should be used
    sys.setrecursionlimit(99999)
    if (len(sys.argv) > 1):
        parse(sys.argv[1])
    else:
        parse(raw_input("Specify input filepath:"))
    pass

if __name__ == '__main__':
    main()

##class node:
##    #diagonal points in direction of (0,0) to (m,n)
##    #top left, bottom right, left, right, top, bottom
##    tl, br, l, r, t, b, score = None, None, None, None, None, None, 0
##
##    def __init__(self, data):
##        self.score = data


### left, top left, top
##def score(l, tl, t, points):
##    return min(l, d, t) + points

##def alignSequences(a, b):
##    # start building a table
##    head = node(0)
##    tail = head
##    for row, i in enumerate(b):
##        first_in_lastrow = None
##        for column, j in enumerate(a):
##            buffer = tail
##            l, tl, t = None, None, None
##            if (row == 0 and column == 0):
##                l = -column-1
##                l = -row-1
##                tl = 0
##                head.br = buffer
##            elif (row == 0):
##                t = -column-1
##                l = tail
##                tail.r = buffer
##            elif (column == 0):
##                l = -row-1
##                t = first_in_lastrow
##                first_in_lastrow.b = buffer
##            else:
##                tl = buffer.tl
##                l = buffer.l
##                t = buffer.t
##            if (i == j):
##                buffer.score = score(l, tl, t, 2)
##            else:
##                buffer.score = score(l, tl, t, -1)
##            buffer.l = l
##            buffer.tl = tl
##            buffer.t = t
##            tail = buffer
##            if (first_in_lastrow == None):
##                first_in_lastrow = tail