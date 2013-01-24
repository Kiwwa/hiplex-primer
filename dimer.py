#!/usr/bin/env python

#
# Screening dimer
#
# Slide bottom strand +1/-1 while finding matches against Top strand.
# By both forward screening and backward screening,
# we can find all possible matches.
#
#  Forward screening
#                       5' CCCAGTTTTAATATTTG 3'
#                          |     || ||     |
#                       3' GTTTATAATTTTGACCC 5'
#
#                       5' CCCAGTTTTAATATTTG 3'
#                           | | | |||| |   |
#                        3' GTTTATAATTTTGACCC 5'  ->
#
#                       5' CCCAGTTTTAATATTTG 3'
#                            ||  | | | |   |
#                         3' GTTTATAATTTTGACCC 5' -->
#
#                       ....
#
#                       5' CCCAGTTTTAATATTTG 3'
#
#                                       3' GTTTATAATTTTGACCC 5'
#
#
#  Backward screening
#                       5' CCCAGTTTTAATATTTG 3'
#                               ||  ||
#                 <-   3' GTTTATAATTTTGACCC 5'
#
#                       5' CCCAGTTTTAATATTTG 3'
#                               |   | |
#                <--  3' GTTTATAATTTTGACCC 5'
#
#                       ....
#
#                       5' CCCAGTTTTAATATTTG 3'
#
#       3' GTTTATAATTTTGACCC 5'
#
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


class Dimer(object):
    def __init__(self, top, bottom=None):
        self._top = top
        self._bottom = bottom if bottom is not None else top[::-1]
        self._maxscore = None  # (score, direction, startIndex,  positions)

    def score(self):

        def detect_dimers():
            matches = []
            # forward screening
            for i in range(0, len(self._top)):
                top = self._top[i:]
                positions = complement_positions(top, self._bottom)
                if positions:
                    matches.append(('forward', i, positions))
            # backward screening
            # We do not need to start from index 0,
            # since it is already convered by forward screening.
            for i in range(1, len(self._bottom)):
                bottom = self._bottom[i:]
                positions = complement_positions(self._top, bottom)
                if positions:
                    matches.append(('backward', i, positions))
            return matches

        def find_maxscore(matches):
            # maxscore = (score, direction, startIndex, positions)
            maxscore = (0, None, None, None)
            for direction, index, positions in matches:
                #print_dimer(direction, index, positions)
                score = get_score(positions)
                if score > maxscore[0]:
                    maxscore = (score, direction, index, positions)
            return maxscore

        def print_dimer(direction, index, positions):
            if positions:
                matchSigns = [s for s in ' '*min(len(self._top),
                                                 len(self._bottom))]
                for i in positions:
                    matchSigns[i] = '|'
                matchSigns = ''.join(matchSigns)
                indentation = index
                if direction == 'backward':
                    top = ' '*indentation + "5'> " + self._top + " >3'"
                    bottom = "3'< " + self._bottom + " <5'"
                    matchSigns = ' '*(indentation+4) + matchSigns
                if direction == 'forward':
                    top = "5'> " + self._top + " >3'"
                    bottom = ' '*indentation + "3'< " + self._bottom + " <5'"
                    matchSigns = ' '*(indentation+4) + matchSigns
                print '-'*2*(max(len(self._top), len(self._bottom)) + 10)
                print ''.join(top)
                print ''.join(matchSigns)
                print ''.join(bottom)
                print '-'*2*(max(len(self._top), len(self._bottom)) + 10)

        matches = detect_dimers()
        maxscore = find_maxscore(matches)
        print 'Max scored dimer is', maxscore
        print_dimer(maxscore[1], maxscore[2], maxscore[3])
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

s = 'GGATTGATAATGTAATAGG'
t = 'CATTATGGGTGGTATGTTGG'
Dimer(s, t[::-1]).score()
#print 'num of comparisions = %d' % mycomparision
#print '----------------------'

#for primer in primers:
#    print 'For the primer:', primer
#    Dimer(primer).score()
