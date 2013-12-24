# hairpin code

# was in find_best_primer
if options.maxhairpinsize is not None:
    filtered = filter_primer(options,
                             'hairpin',
                             primer_suffix)
primer_score = score_primer(options.melt, primer_suffix)
if filtered:
    logging.info(("Score: %4d, %s5'> %s <3', Filtered (%s)"
                  % (primer_score,
                     ' ' * suffix_start,
                     primer_suffix,
                     'hairpin')))
    continue
else:
    None # code just continue the find_best_primer function
    

# was in get_best_window
filtered, reason = check_filter(options, scoredBlocks)
if filtered:
    logging.info('Window: %d is filtered due to %s.'
                 % (window, reason))
    continue
None

# was only referenced once throughout Primer_design.py, in the above code
def check_filter(options, scoredBlocks):
    filtered = False
    reason = None
    for block in scoredBlocks:
        forward = block.primer_forward
        reverse = block.primer_reverse
        # If forward or reverse is None, that means
        # the primer was filtered at the stage where the primer is selected
        # among the primers within the allowed primer variation.
        if forward is None or \
           reverse is None:
            return True, 'hairpin'
        # primer-dimer and difference in Tm of both forward and reverse
        # have to be considered with two primers together.
        # if filter_primer(options, 'dimer', forward.bases, reverse.bases):
        #    return True, 'dimer'
        # if filter_primer(options, 'tmdiff', forward.bases, reverse.bases):
        #    return True, 'tmdiff'
    return filtered, reason
    

# whole class, lol (deals with detection and reporting of hairpins)
class Hairpin(object):

    def __init__(self, sequence):
        self.sequence = sequence
        self._score = None
        self._top = None
        self._bottom = None
        self._index = None
        self._positions = None

    def print_hairpin(self):

        def hairpin_string(top, bottom, startIndex, positions):
            if positions:
                matches = [s for s in ' ' * min(len(top), len(bottom))]
                for i in positions:
                    matches[i] = '|'
                indentation = startIndex
                top = ''.join(top) + " <5'"
                matches = ' ' * indentation + ''.join(matches)
                bottom = ' ' * indentation + ''.join(bottom) + " >3'"
                hairpin = '%s\n%s\n%s' % (top, matches, bottom)
                return hairpin
                # print '-'*(len(top + bottom) + startIndex)
                # print ' '+''.join(top) + " <5'"
                # print ' '*indentation, ''.join(matches)
                # print ' '*indentation, ''.join(bottom) + " >3'"
                # print '-'*(len(top + bottom) + startIndex)

        if self._score is None:
            self.score()
        return hairpin_string(self._top, self._bottom,
                              self._index, self._positions)

    def score(self):

        def detect_hairpins(toplength):
            top = self.sequence[:toplength]
            bottom = self.sequence[toplength:max(toplength * 2,
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
                score = get_score(positions)
                if score > maxscore[0]:
                    maxscore = (score, top, bottom, index, positions)
            return maxscore

        if self._score is None:
            # maxscore = (score, top, bottom, startIndex, positions)
            maxscore = (0, None, None, None, None)
            for toplength in range(1, len(self.sequence)):
                hairpins = detect_hairpins(toplength)
                score = find_maxscore(hairpins)
                if score[0] > maxscore[0]:
                    maxscore = score
            self._score = maxscore[0]
            self._top = maxscore[1]
            self._bottom = maxscore[2]
            self._index = maxscore[3]
            self._positions = maxscore[4]
        return self._score


# pretty sure this wasn't being used at all; code for dimer-detection?        
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
                matchSigns = [s for s in ' ' * min(len(self._top),
                                                   len(self._bottom))]
                for i in positions:
                    matchSigns[i] = '|'
                matchSigns = ''.join(matchSigns)
                indentation = index
                if direction == 'backward':
                    top = ' ' * indentation + "5'> " + self._top + " >3'"
                    bottom = "3'< " + self._bottom + " <5'"
                    matchSigns = ' ' * (indentation + 4) + matchSigns
                if direction == 'forward':
                    top = "5'> " + self._top + " >3'"
                    bottom = ' ' * indentation + "3'< " + self._bottom + " <5'"
                    matchSigns = ' ' * (indentation + 4) + matchSigns
                print '-' * 2 * (max(len(self._top), len(self._bottom)) + 10)
                print ''.join(top)
                print ''.join(matchSigns)
                print ''.join(bottom)
                print '-' * 2 * (max(len(self._top), len(self._bottom)) + 10)

        matches = detect_dimers()
        maxscore = find_maxscore(matches)
        # print 'Max scored dimer is', maxscore
        #print_dimer(maxscore[1], maxscore[2], maxscore[3])
        return maxscore[0]