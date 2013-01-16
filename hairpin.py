#!/usr/bin/env python


#
# Screening hairpin
#
# The primer folds itself.
# The screeing continues until the strand is completely reversed.
#
#  e.g) 5' GTTTTGCTCGCCAACGGTTTGG 3'
#
#       G                           -- top
#      (
#       TTTTGCTCGCCAACGGTTTGG       -- bottom
#
#       G
#      T
#       TTTGCTCGCCAACGGTTTGG
#
#       TG
#      (
#       TTTGCTCGCCAACGGTTTGG
#
#      ....
#
#       TTTGGCAACCGCTCGTTTTG
#      G
#       G
#
#       GTTTGGCAACCGCTCGTTTTG
#      (
#       G
#
# When the length of top is n there are n possible overlaps to be considered,
# but we do not need to calculate all n overlaps.
#  e.g) when length of top is three,
#    1)   TTG
#         TTGCTCGCCAACGGTTTGG
#
#    2)  TTG
#         TTGCTCGCCAACGGTTTGG
#
#    3) TTG
#         TTGCTCGCCAACGGTTTGG
#
# We do not need to consider 3), since 3) is equivalent to
# the overlap of top with length 2.
# (The 'TT' itself overlaps each other.)
#
#       TG
#       TTTGCTCGCCAACGGTTTGG
#
# This case would be considered when the top length is 2.
# Therefore, when we consider the overlaps, we only need to just two cases.
# 1) Top fully overlaps with bottom.
# 2) Top without first nucleotide overlaps with bottom.
#


def isComplement(x, y):
    cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return cs[x] == y

mycomparision = 0


def complement_positions(top, bottom):
    positions = []
    global mycomparision
    for i, (t, b) in enumerate(zip(top, bottom)):
        mycomparision += 1
        if isComplement(t, b):
            positions.append(i)

    return positions

def get_score(positions):
    # Find the longest consecutive matches
    score = 0
    maxscore = 0
    if len(positions) == 1:
        score = 1
    elif len(positions) > 1:
        score = 1
        lastposition = positions[0]
        for position in positions[1:]:
            if position == lastposition + 1:
                score += 1
            else:
                if score > maxscore:
                    maxscore = score
                score = 1
            lastposition = position
    if score > maxscore:
        maxscore = score
    return maxscore

class Hairpin(object):
    def __init__(self, sequence):
        self.sequence = sequence

    def score(self):

        def detect_hairpins(toplength):
            top = self.sequence[:toplength]
            bottom = self.sequence[toplength:max(toplength*2,
                                                 len(self.sequence))]
            top = top[::-1]
            hairpins = []
            for startIndex in [0, 1]:
                positions = complement_positions(top[startIndex:], bottom)
                if positions:
                    hairpins.append((top, bottom, startIndex, positions))
            return hairpins

        def find_maxscore(hairpins):
            # maxscore = (score, top, bottom, startIndex, positions)
            maxscore = (0, None, None, None, None)
            for top, bottom, index, positions in hairpins:
                #print_hairpin(top, bottom, index, positions)
                score = get_score(positions)
                if score > maxscore[0]:
                    maxscore = (score, top, bottom, index, positions)
            return maxscore

        def print_hairpin(top, bottom, startIndex, positions):
            if positions:
                matches = [s for s in ' '*min(len(top), len(bottom))]
                for i in positions:
                    matches[i] = '|'
                indentation = startIndex

                print '-'*(len(top + bottom) + startIndex)
                print ' '+''.join(top) + " <5'"
                print ' '*indentation, ''.join(matches)
                print ' '*indentation, ''.join(bottom) + " >3'"
                print '-'*(len(top + bottom) + startIndex)

        # maxscore = (score, top, bottom, startIndex, positions)
        maxscore = (0, None, None, None, None)
        for toplength in range(1, len(self.sequence)):
            hairpins = detect_hairpins(toplength)
            score = find_maxscore(hairpins)
            if score[0] > maxscore[0]:
                maxscore = score
        print 'Max scored hairpin is', maxscore
        print_hairpin(maxscore[1], maxscore[2], maxscore[3], maxscore[4])
        return maxscore[0]


primers = ['GTTTTGCTCGCCAACGGTTTGG',
           'GTTTTGCCAACACCCGGTTTGG',
           'GCCTTGCTCGCACCCGGTAAGG',
           'GTTTTGCTCGCCGACGGTTCGG',
           'GTTTTAATCGCCCCCGGTGGGG',
           'GTTTTGCTCGAGACCGATCTCG',
           'GTTTCGAGAGCACCCGGTCTCG',
           'TTCGTGGCTCGCACCCGGTCACGA',
           'GTTTCGAGCGCACCCGGGCTCG',
           'GTTTCTATAGCACCCGGTATAG',
           'GTTCCACTCGCACCCGGTCTGG',
           'GTTTTGCTCGCACCCGGTCTGG',
           'GAGGAAGTGGACACGGGTTAG',
           'TTTACTCGCAGCCGGTCTG',
           'TTTGCTCGCAGCCGGTCTG',
           'TTTGATCGCAGACGGTCTG',
           'TTTGATCGCAGACGCTGCG',
           'TTTGCTCGCAGCCGGTCTG',
           'ATTAGGCAGAGGTGAAAAAG']

s = 'GAGGAAGTGGACACGGGTTAG'
Hairpin(s).score()
#print 'num of comparisions = %d' % mycomparision
#print '----------------------'

#for n, primer in enumerate(primers):
#    print 'For the primer:', primer
#    Hairpin(primer).score()
