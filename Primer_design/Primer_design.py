#!/usr/bin/env python

'''
Primer design program for amplicon sequencing.

Authors: Bernie Pope, Danny Park, Tu Nguyen-Dumont, Sori Kang, Luke Shillabeer
'''

from Bio import SeqIO
import sys
import math
import os
import csv
import itertools
import argparse
import logging

parser = argparse.ArgumentParser(description='A primer design tool for '
                                             'amplicon sequencing')
parser.add_argument('--refdir',
                    metavar='REF',
                    type=str,
                    required=True,
                    help='directory containing reference files for each '
                         'chromosome in fasta format')
parser.add_argument('--genes',
                    metavar='GENES',
                    type=str,
                    required=True,
                    help='gene list file in CSV format')
parser.add_argument('--blocksize',
                    metavar='B',
                    type=int,
                    required=True,
                    help='block size')
parser.add_argument('--maxprimersize',
                    metavar='M',
                    type=int,
                    required=True,
                    help='maximum primer size')
parser.add_argument('--primervar',
                    metavar='V',
                    type=int,
                    required=True,
                    help='primer size variance')
parser.add_argument('--splicebuffer',
                    metavar='S',
                    type=int,
                    required=True,
                    help='splice site buffer size')
parser.add_argument('--melt',
                    metavar='M',
                    type=int,
                    required=True,
                    help='ideal melting temperature Tm')
parser.add_argument('--log',
                    metavar='L',
                    type=str,
                    help='optional log file')
parser.add_argument('--idtfile',
                    metavar='FILE',
                    type=str,
                    required=True,
                    help='CSV output file for IDT')
parser.add_argument('--roverfile',
                    metavar='RFILE',
                    type=str,
                    help='TSV output for ROVER tool')
parser.add_argument('--blocksizevar',
                    metavar='B',
                    type=int,
                    help='block size variance')
parser.add_argument('--scale',
                    metavar='S',
                    type=str,
                    default='25 nmole',
                    help='ordering scale for the oligos')
parser.add_argument('--purification',
                    default='Standard Desalting',
                    choices=('Standard Desalting', 'HPLC'),
                    help='purification')
parser.add_argument('--senseheelseq',
                    metavar='SH',
                    type=str,
                    help='Sense strand primer heel sequence')
parser.add_argument('--antisenseheelseq',
                    metavar='AH',
                    type=str,
                    help='Sense strand primer heel sequence')

# Global Vars
banner_width = 80   # width of banner message separator
block_sizes = []    # XXX Why is this outside a def or class? (was above                         print_blocksize_distribution)
complements = {'N': 'N', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def main():
    options = parser.parse_args()
    
    with open(options.idtfile, 'w') as itdfileopen, \
         open(options.roverfile, 'w') as roverfileopen:
            
        if options.log:
            logging.basicConfig(filename=options.log,
                                level=logging.DEBUG,
                                filemode='w',
                                format='%(message)s')

        logging.info('command line: {}'.format(' '.join(sys.argv)))

        gene_file = GeneFile(options.genes)
        for (gene_name,
             exon_id,
             chromosome,
             exon_start,
             exon_end) \
        in gene_file.process():
            # score all the sliding windows over this exon
            window_scores = score_exon_windows(options,
                                               chromosome,
                                               gene_name,
                                               exon_id,
                                               exon_start,
                                               exon_end)
            # find the best scoring window for this exon
            # best_window = get_best_window(options, window_scores)
            best_blocks = get_optimal_primer_combination(options, 
                                                         window_scores)
            # print out the primers for the best window
            if best_blocks is not None:
                print_best_primers(options,
                                   gene_name,
                                   exon_id,
                                   chromosome,
                                   exon_start,
                                   exon_end,
                                   best_blocks,
                                   itdfileopen,
                                   roverfileopen)
        gene_file.close()


class GeneFile(object):

    def __init__(self, filename):
        self.file = open(filename, 'U')
        # mapping from gene name to most recently used block number
        # block numbers start at 1
        self.gene_blocks = {}

    def process(self):
        for line in self.file:
            fields = line.split()
            chromosome, start, end, name = fields[0:4]
            start = int(start)
            end = int(end)
            gene_name_parts = name.split(',')
            if len(gene_name_parts) >= 2:
                gene_name = gene_name_parts[0]
            else:
                gene_name = name
            if gene_name in self.gene_blocks:
                current_block_number = self.gene_blocks[gene_name]
                self.gene_blocks[gene_name] = \
                    block_number = \
                    current_block_number + 1
            else:
                self.gene_blocks[gene_name] = block_number = 1
            yield gene_name, block_number, chromosome, start, end

    def close(self):
        self.file.close()


class ScoredBlock(object):

    def __init__(self, block_num, start, end, primer_forward, primer_reverse):
        # sequential number of the block in the window
        self.block_num = block_num
        self.start = start  # first coordinate in the block
        self.end = end     # last coordinate in the block
        self.primer_forward = primer_forward
        self.primer_reverse = primer_reverse


class ScoredPrimer(object):

    def __init__(self, score, direction, start, end, bases):
        self.score = score
        self.direction = direction
        self.start = start
        self.end = end
        self.bases = bases

    def __lt__(self, other):
        return self.end < other.end

# Compute the best primers for an exon, using a sliding window.
#
# |.................Region...................|
#
#   MP   SL   SB      Exon      SB   SL   MP
# |----|----|----|------------|----|----|----|
#
#      ^========^========^=========^
#        block1   block2   block3
#
#      |..........Window...........|
#
# MP = max primer size
# SL = slack
# SB = splice buffer
#
# Window is num_blocks * block_size.
#
# Region = exon + (2 * slack) + (2 * max_primer)
#
# A region is the total chunk of DNA that we read from the fasta file.
#
# We slide the window along slack+1 steps. For each window position
# we find the best primers for each block in that position.


def score_exon_windows(options,
                       chromosome,
                       gene_name,
                       exon_id,
                       exon_start,
                       exon_end):
    exon_size = (exon_end - exon_start) + 1

    # add the splice buffer on to the start and end of the exon coordinates
    exon_buffer_start = exon_start - options.splicebuffer
    exon_buffer_end = exon_end + options.splicebuffer
    exon_buffer_size = (exon_buffer_end - exon_buffer_start) + 1
    num_blocks = (int(math.ceil
                      (float(exon_buffer_size) / float(options.blocksize))))
    window_size = num_blocks * options.blocksize
    slack = window_size - exon_buffer_size

    # the region must include the exon, plus the splice buffer,
    # plus slack on either side, plus maximum primers on either side
    region_start = exon_buffer_start - slack - options.maxprimersize
    region_end = exon_buffer_end + slack + options.maxprimersize

    logging.info('*' * banner_width)
    logging.info('chrom:\t\t%s' % chromosome)
    logging.info('exon:\t\t%s' % exon_id)
    logging.info('exon start:\t%d' % exon_start)
    logging.info('exon buff start:%d' % exon_buffer_start)
    logging.info('exon end:\t%d' % exon_end)
    logging.info('exon buff end:\t%d' % exon_buffer_end)
    logging.info('exon size:\t%d' % exon_size)
    logging.info('exon buff size:\t%d' % exon_buffer_size)
    logging.info('block_size:\t%d' % options.blocksize)
    logging.info('num_blocks:\t%d' % num_blocks)
    logging.info('window_size:\t%d' % window_size)
    logging.info('slack:\t\t%d' % slack)
    logging.info('region_start:\t%d' % region_start)
    logging.info('region_end:\t%d' % region_end)

    region = get_region(options, chromosome, region_start, region_end)

    # get the first and last 5 bases in the region so we can print it
    # out for diagnostic purposes.
    if len(region.bases) > 10:
        region_start_bases = region.bases[0:5]
        region_end_bases = region.bases[-5:]
        logging.info(('region starts with: %s, region ends with: %s'
                      % (region_start_bases, region_end_bases)))

    window_scores = {}
    window_count = 1

    # slide the window over the exon, starting "slack" number of positions
    # before the start of the exon (with splice buffer added). The last
    # position of the window is when it lines up with the start of the exon.
    for window_start \
            in range(exon_buffer_start - slack, exon_buffer_start + 1):
        logging.info('=' * banner_width)
        logging.info(('Scoring window starting at %d (%d/%d)'
                      % (window_start,
                         window_count,
                         slack + 1)))
        window_scores[window_start] = []

        # score each block in the window, trying the different primer
        # sizes at each end.
        for block_num in range(0, num_blocks):
            block_start = window_start + (block_num * options.blocksize)
            block_end = block_start + options.blocksize - 1

            logging.info('-' * banner_width)
            logging.info('block:\t\t%d' % block_num)
            logging.info('block start:\t%d' % block_start)
            logging.info('block end:\t%d' % block_end)

            # find the best primer for either the start or the end of a block.
            def find_best_primer(direction, primer_end):

                logging.info('~' * (banner_width / 2))
                logging.info('Primer dir:\t%s' % direction)
                if direction == 'forward':
                    primer_start = (primer_end - options.maxprimersize) + 1
                    primer = get_primer(region, primer_start, primer_end)
                elif direction == 'reverse':
                    primer_start = (primer_end + options.maxprimersize) - 1
                    # the reverse primer comes from the opposite strand
                    primer = get_primer(region, primer_end, primer_start)
                    primer = reverse_complement(primer)

                logging.info('Primer start:\t%d' % primer_start)
                logging.info('Primer end:\t%d' % primer_end)

                primer_len = len(primer)

                # check that the primer length is what we believe it to be
                assert primer_len == (abs(primer_end - primer_start) + 1)

                # score all suffixes of the primer (we try the first N
                # suffixes, where N = primer_variance)
                best_primer = None
                for suffix_start in range(0, options.primervar):
                    primer_suffix = primer[suffix_start:]
                    primer_score = score_primer(options.melt, primer_suffix)
                    filtered = False
                    logging.info(("Score: %4d, %s5'> %s <3'"
                                  % (primer_score,
                                     ' ' * suffix_start,
                                     primer_suffix)))
                    # calculate the start position of this primer candidate
                    if direction == 'forward':
                        candidate_start = primer_start + suffix_start
                    elif direction == 'reverse':
                        candidate_start = primer_start - suffix_start
                    candidate = ScoredPrimer(primer_score,
                                             direction,
                                             candidate_start,
                                             primer_end,
                                             primer_suffix)
                    if best_primer == None:
                        # this is the first one we've seen
                        best_primer = candidate
                    elif primer_score == best_primer.score:
                        # we have a tie for best primer, choose the shortest
                        # one
                        if len(candidate.bases) < len(best_primer.bases):
                            best_primer = candidate
                    elif primer_score < best_primer.score:
                        # this is the new best score
                        best_primer = candidate
                if best_primer:
                    logging.info('Best primer is: %s' % best_primer.bases)
                else:
                    logging.info('There is no best primer we can choose.')
                return best_primer

            # find the best forward and reverse primers for this block
            # using the coordinate of the last position in the primer to
            # indicate the primer position
            primer_forward = find_best_primer('forward', block_start - 1)
            primer_reverse = find_best_primer('reverse', block_end + 1)

            # Add primer even though it is None
            block_score = ScoredBlock(block_num,
                                      block_start,
                                      block_end,
                                      primer_forward,
                                      primer_reverse)
            window_scores[window_start].append(block_score)
        window_count += 1
    return window_scores


def get_optimal_primer_combination(options, window_scores):

    def get_best_scored_primer(block, direction, var_start, var_end):
        # Among the primers within variance boundary,
        # find the best scored primers.
        # The primers could be multiple if their score is the same.
        primers = block[direction][var_start:var_end]
        best_primers = []
        best_score = None
        for n, primer in enumerate(primers):
            if primer is not None:
                if best_score is None:
                    best_primers = [(n + var_start, primer)]
                    best_score = primer.score
                elif best_score > primer.score:
                    best_primers = [(n + var_start, primer)]
                    best_score = primer.score
                elif best_score == primer.score:
                    best_primers.append((n + var_start, primer))
        return best_primers

    def get_block_primers(window_scores):
        # Make the list of primers that belong to
        # the same block and the same direction.
        # The list of primer should be ordered on the basis of window.
        # <window_scores>
        # Each window has a list of ScoredBlock.
        # Each ScoredBlock has two ScoredPrimer for forward and reverse.
        # e.g, window_scores = {1234:[ScoredBlock, ScoredBlock],
        #                       1235:[ScoredBlock, ScoredBlock]}
        # => <blocks>
        # Each block has two types of primers. forward and reverse.
        # e.g, blocks = {0: {'forward':[ScoredPrimer, ..],
        #                      'reverse':[ScoredPrimer, ..]},
        #                1: {'forward':[ScoredPrimer, ..],
        #                      'reverse':[ScoredPrimer, ..]}
        num_blocks = 0
        blocks = {}
        # Get window range and order
        windows = [window for window in window_scores]
        windows.sort()
        # Get block numbers
        for window, scoredBlocks in window_scores.items():
            num_blocks = len(scoredBlocks)
            break

        for num in range(0, num_blocks):
            blocks[num] = {}
            blocks[num]['forward'] = []
            blocks[num]['reverse'] = []

        for window in windows:
            scoredBlocks = window_scores[window]
            for block in scoredBlocks:
                forward = block.primer_forward
                reverse = block.primer_reverse
                blocks[block.block_num]['forward'].append(forward)
                blocks[block.block_num]['reverse'].append(reverse)

        return windows, blocks

    def get_optimal_distance_pair(forward_primers, reverse_primers,
                                  block_size):
        # Among the all the possible combination of the best primers,
        # find the primer pairs so that the block size could be
        # very close to the option parameter; blocksize
        optimal_distance = None
        optimal_forward = None
        optimal_reverse = None
        for f_pos, f_primer in forward_primers:
            for r_pos, r_primer in reverse_primers:
                distance = abs(r_primer.end - f_primer.end + 1 - block_size)
                if optimal_distance is None:
                    optimal_distance = distance
                    optimal_forward = (f_pos, f_primer)
                    optimal_reverse = (r_pos, r_primer)
                elif optimal_distance > distance:
                    optimal_distance = distance
                    optimal_forward = (f_pos, f_primer)
                    optimal_reverse = (r_pos, r_primer)
        return optimal_forward, optimal_reverse

    def get_subset_combination(blocks, num_blocks, start_index,
                               block_size, maxvar):
        # Find the best combination of primers within
        # (start_index ~ start_index + maxvar)
        # If there is no primers for any block and any direction,
        # which means there is no candidate primer to be combined,
        # return None.
        var_start = start_index
        var_end = var_start + maxvar + 1
        total_score = 0
        subset = []
        positions = []
        # no_candidates = False
        for block_num in range(0, num_blocks):
            best_forward_primers = get_best_scored_primer(blocks[block_num],
                                                          'forward',
                                                          var_start,
                                                          var_end)
            best_reverse_primers = get_best_scored_primer(blocks[block_num],
                                                          'reverse',
                                                          var_start,
                                                          var_end)
            if not best_forward_primers or not best_reverse_primers:
                # Could not find any primer within block variance
                # no_candidate = True
                break
            best_forward, best_reverse = get_optimal_distance_pair(
                                                        best_forward_primers,
                                                        best_reverse_primers,
                                                        block_size)
            best_forward_index, forward_primer = best_forward
            best_reverse_index, reverse_primer = best_reverse
            positions += [best_forward_index, best_reverse_index]
            block_start = forward_primer.end + 1
            block_end = reverse_primer.end - 1
            optimal_block = ScoredBlock(block_num,
                                        block_start,
                                        block_end,
                                        forward_primer,
                                        reverse_primer)
            subset.append(optimal_block)
            total_score += (forward_primer.score + reverse_primer.score)
            
            # Need to adjust the var_start and var_end, since
            # for the following block, the var_end should not be greater than
            # the index of the previous block's reverse primer.
            var_end = best_reverse_index + 1
            var_start = max(0, best_reverse_index - maxvar)
        if len(subset) == num_blocks:
            return total_score, subset, positions
        else:
            return None, None, None

    maxvar = options.blocksizevar if options.blocksizevar is not None else 0
    windows, blocks = get_block_primers(window_scores)
    num_windows = len(window_scores)
    num_blocks = len(blocks)
    optimal_score = None
    optimal_blocks = []
    optimal_subset_start = 0
    # Scan the primers within the maximum block size variance, and
    # select the best scored primers of each direction.
    for start_index in range(0, num_windows):
        score, subset, positions = get_subset_combination(
            blocks, num_blocks, start_index,
            options.blocksize, maxvar)
        if subset:
            window_selected = [windows[i] for i in positions]
            logging.info('Subset start: %d, score: %d, window selected: %s'
                         % (windows[start_index],
                            score,
                            ','.join([str(w) for w in window_selected])))
            if optimal_score is None:
                optimal_score = score
                optimal_blocks = subset
                optimal_subset_start = start_index
            elif optimal_score > score:
                optimal_score = score
                optimal_blocks = subset
                optimal_subset_start = start_index
            elif optimal_score == score:
                logging.info('%s subset start (%d) and subset start (%d)'
                             % ('Tie for the optimal combination between',
                                windows[start_index],
                                windows[optimal_subset_start]))
        else:
            logging.info('Subset start: %d, Could not find primer.'
                         % windows[start_index])
    if optimal_score is None:
        logging.info('Could not find primer pairs because there may not be '
                     'any candidate primers.')
    else:
        logging.info('Best combination: subset start (%d), score (%d)'
                     % (windows[optimal_subset_start], optimal_score))
    return optimal_blocks


# given all the scores for each position of 
# the sliding window, find the best one

def get_best_window(options, window_scores):
    best_window = None
    best_score = None
    for window, scoredBlocks in window_scores.items():
        total_score = get_total_window_score(scoredBlocks)
        logging.info('Window: %d, total score: %d' % (window, total_score))
        if best_window == None:
            # this is the first one we've seen
            best_window = window
            best_score = total_score
        elif total_score == best_score:
            # XXX we have a tie for the best, what to do?
            logging.info('Tie for best window between %d and %d'
                         % (best_window, window))
        elif total_score < best_score:
            best_window = window
            best_score = total_score
    if best_window is None:
        logging.info('Could not find the best window.')
    else:
        logging.info('Best window: %d, score: %d' % (best_window, best_score))
        check_potential_issues(window_scores[best_window])
    return best_window


def get_total_window_score(scoredBlocks):
    score = 0
    for block in scoredBlocks:
        score += block.primer_forward.score + block.primer_reverse.score
    return score


def check_potential_issues(scoredBlocks):
    problems = None
    for block in scoredBlocks:
        forward_hairpin = Hairpin(block.primer_forward.bases)
        reverse_hairpin = Hairpin(block.primer_reverse.bases)
        if forward_hairpin.score() >= 5 or reverse_hairpin.score() >= 5:
            logging.info('+' * banner_width)
            logging.info('Potential issues.')
            logging.info('block num: %s, start: %s, end: %s' %
                         (block.block_num, block.start, block.end))
        if forward_hairpin.score() >= 5:
            logging.info('Forward primer %s has hairpin with length of %s.'
                         % (block.primer_forward.bases,
                            forward_hairpin.score()))
            logging.info(forward_hairpin.print_hairpin())
        if reverse_hairpin.score() >= 5:
            logging.info('Reverse primer %s has hairpin with length of %s.'
                         % (block.primer_reverse.bases,
                            reverse_hairpin.score()))
            logging.info(reverse_hairpin.print_hairpin())


def print_blocksize_distribution():
    if not block_sizes:
        return
    block_sizes.sort()
    min = block_sizes[0]
    max = block_sizes[-1]
    max_count = 0
    distribution = {}
    for size in block_sizes:
        if size not in distribution:
            distribution[size] = 1
        else:
            distribution[size] += 1
        if distribution[size] > max_count:
            max_count = distribution[size]
    for y in range(max_count, 0, -1):
        ys = []
        for key, value in distribution.items():
            if value >= y:
                ys.append('--')
            else:
                ys.append('')
        print '\t'.join(ys)
    print '\t'.join([str(size) for size in distribution])


def print_best_primers(options, gene_name, exon_id, chromosome,
                       exon_start, exon_end, scored_blocks, 
                       idtfile, roverfile):
    csv_file = idtfile
    rover_file = roverfile
    
    primer_name_prefix = gene_name + '_' + str(exon_id) + '_'
    print('-' * banner_width)
    print('gene: %s, exon: %s, %s:%d-%d' % (gene_name, exon_id, chromosome,
                                            exon_start, exon_end))
    print('-' * banner_width)
    print #formatting the output
    for block in scored_blocks:

        block_size = block.end - block.start + 1
        forward = block.primer_forward
        reverse = block.primer_reverse

        # required for hairpin scores (currently requires upper bases)
        fwd_upper_bases = forward.bases.upper()
        rev_upper_bases = reverse.bases.upper()

        if options.senseheelseq is not None:
            forward.bases = options.senseheelseq + str(forward.bases)
        if options.antisenseheelseq is not None:
            reverse.bases = options.antisenseheelseq + str(reverse.bases)

        block_size = block.end - block.start + 1
        # block_sizes.append(block_size)
        print('block %d, %d-%d, block size: %d'
              % (block.block_num, block.start, block.end, block_size))
        print('forward: %d-%d, %s'
              % (forward.start, forward.end, forward.bases))
        print('reverse: %d-%d, %s'
              % (reverse.start, reverse.end, reverse.bases))
        primer_name_forward = (primer_name_prefix + 'F'
                                                  + str(block.block_num + 1))
        primer_name_reverse = (primer_name_prefix + 'R'
                                                  + str(block.block_num + 1))
        csv_file.write(','.join([primer_name_forward,
                                 str(forward.bases),
                                 options.scale,
                                 options.purification])
                       + '\n')
        csv_file.write(','.join([primer_name_reverse,
                                 str(reverse.bases),
                                 options.scale,
                                 options.purification])
                       + '\n')
        print #formatting the output

        # generate the ROVER compatible input file (tab delimited format)
        # (compliant to BED file format)
        rover_file.write((chromosome + '\t'
                                     + str(block.start)
                                     + '\t'
                                     + str(block.end)
                                     + '\n'))


class Region(object):
    # a (possibly large) chunk of DNA from the reference
    
    def __init__(self, chromosome, start, end, file, bases):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.file = file
        self.bases = bases




def get_region(options, chromosome, start, end):
    # SeqIO uses 0 based indexes, but our input data uses 1 based indexes,
    # so we need to convert here
    
    ref_filename = os.path.join(options.refdir, chromosome + '.fa')
    logging.info('Reading region on %s, start: %d, end: %d from file %s'
                 % (chromosome, start, end, ref_filename))
    with open(ref_filename) as ref_file:
        # XXX should really cache this
        bases_list = list(SeqIO.parse(ref_file, 'fasta'))
        if len(bases_list) == 1:
            bases = bases_list[0]
            if start > 0 and end <= len(bases):
                # here is where we do index conversion from 1 based to 0 based
                segment = bases[start - 1:end].seq
                normalised_segment = segment.upper()
                validate_sequence(normalised_segment)
                return Region(chromosome,
                              start,
                              end,
                              ref_filename,
                              normalised_segment)
            else:
                exit('chromosome %s does not span %d %d'
                     % (chromosome, start, end))
        exit('got wrong number of sequences in bases_list: %d'
             % len(bases_list))


def get_primer(region, start, end):
    new_start = start - region.start
    new_end = end - region.start
    return region.bases[new_start:new_end + 1]


def validate_sequence(sequence):
    # check our sequence only has A,T,G,C
    
    for base in sequence:
        if base not in 'ATGC':
            exit('bad base found: %s' % base)
    return sequence


def reverse_complement(sequence):
    return ''.join(map(complement, sequence[::-1]))


def complement(base):
    if base in complements:
        return complements[base]
    else:
        exit('found bad base %s' % base)


def score_primer(ideal_melting_temp, sequence):
    t_m = melting_temp(sequence)
    return (ideal_melting_temp - t_m) ** 2

base_temp = {'A': 2, 'T': 2, 'G': 4, 'C': 4}


def melting_temp(sequence):
    t_m = 0
    for base in sequence:
        t_m += base_temp[base]
    return t_m


def isComplement(x, y):
    cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return cs[x] == y


def complement_positions(top, bottom):
    positions = []
    for i, (t, b) in enumerate(zip(top, bottom)):
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


def filter_primer(options, filter, primer_a, primer_b=None):

    def gc_count(primer):
        count = 0
        for n in primer:
            if n == 'G' or n == 'C':
                count += 1
        return count

    filtered = False
    if filter == 'hairpin':
        hairpin = Hairpin(primer_a).score()
        if hairpin > options.maxhairpinsize:
            filtered = True
    if filter == 'dimer':
        dimer = Dimer(primer_a, primer_b).score()
        if dimer > 4:
            filtered = True
    if filter == 'endwithgc':
        if not (primer_a.endswith('G') or
                primer_a.endswith('C')):
            filtered = True
        elif gc_count(primer_a[-5:]) > 3:
            filtered = True
    if filter == 'tmdiff':
        tm_a = melting_temp(primer_a)
        tm_b = melting_temp(primer_b)
        if abs(tm_a - tm_b) > 5:
            filtered = True
    return filtered


if __name__ == '__main__':
    main()
