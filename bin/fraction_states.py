#!/usr/bin/env python
#
#  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
#
#  This file is part of Pergola.
#
#  Pergola is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Pergola is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.

################################################################
### Jose A Espinosa. CSN/CB-CRG Group. June 2018             ###
################################################################
### Script creates takes bed with the segmentation anotated  ###
### by chromHMM and calculates the fraction of time spend in ###
### each of the states.                                      ###
################################################################

from argparse import ArgumentParser
from sys import stderr
from os import path
import pybedtools

parser = ArgumentParser(description='File input and options to process data')
parser.add_argument('-b', '--bed_file', help='Bed file containing segmentation annotated with chromHMM', required=True)
parser.add_argument('-p', '--phases_bed_file', help='Bed file containing phases file', required=True)
parser.add_argument('-n', '--n_states', help='Number of states of the HMM', required=False, type=int, default=4)
parser.add_argument('-c', '--chrom_sizes_file', help='File containing length of chromosomes, longest trajectory', required=False)
parser.add_argument('-t', '--tag', help='tag for output file naming', required=False)


args = parser.parse_args()

print >> stderr, "Bed file with chromHMM segmentation: ", args.bed_file
print >> stderr, "Bed file with phases: ", args.phases_bed_file
print >> stderr, "Number of states: ", args.n_states


bed_file_seg = args.bed_file
bed_file_phases = args.phases_bed_file
n_states = args.n_states
light_phases_file = "light_phases.bed"

if args.chrom_sizes_file:
    chrom_sizes = args.chrom_sizes_file
    print >> stderr, "Chrom_sizes file: ", args.chrom_sizes_file
    # get complement of phases, i.e. light phases bed
    pybedtools.BedTool(bed_file_phases).complement(g=chrom_sizes).saveas(light_phases_file)

if args.tag:
    tag = args.tag
    print >> stderr, "Tag output files: ", args.tag
else:
    tag = "file"

bed_seg = pybedtools.BedTool(bed_file_seg)

def featuretype_filter(feature, featuretype):
    if feature[3] == featuretype:
        return True
    return False

def subset_featuretypes(featuretype):
    result = bed_seg.filter(featuretype_filter, featuretype).saveas()
    return pybedtools.BedTool(result.fn)


def coverage_in_features(features_fn, bed_file_ph):
    """
    Calculates coverage of meals on phase bed file
    """
    return pybedtools.BedTool(bed_file_ph).coverage(
                             features_fn,
                             stream=True)

def add_id(f, id="NA", state="NA", phase="NA"):
    """
    adds mice id to the end of the file
    """
    if id_mice % 2 == 1:
        group = "Ctrl"
    else:
        group = "HF"

    f.strand = id
    f.score = group
    f.chrom = state
    f.append(phase)

    return f

for state in range(1, n_states+1):
    bed_state = subset_featuretypes (str(state)).saveas().fn

    id_mice = int(bed_file_seg.split("_")[1])

    # if path.isfile(light_phases_file):
    #     coverage_in_features(bed_state, bed_file_phases).each(add_id, id_mice, state, "dark").saveas(bed_file_seg + ".state_" + str(state) + "_dark" + ".cov")
    #     coverage_in_features(bed_state, light_phases_file).each(add_id, id_mice, state, "light").saveas(bed_file_seg + ".state_" + str(state) + "_light" + ".cov")
    # else:
    #     coverage_in_features(bed_state, bed_file_phases).each(add_id, id_mice, state, "NA").saveas(bed_file_seg + ".state_" + str(state) + ".cov")

    coverage_in_features(bed_state, bed_file_phases).each(add_id, id_mice, state, tag).saveas(bed_file_seg + ".state_" + str(state) + "_" + tag + ".cov")

