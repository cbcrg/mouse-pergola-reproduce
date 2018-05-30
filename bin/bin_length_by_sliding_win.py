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
### Jose A Espinosa. CSN/CB-CRG Group. May 2018              ###
################################################################
### Script creates several bed files corresponding to a bin  ###
### of meals duration for a single period of time (intersect)###
### e.g. week                                                ###
################################################################

from argparse import ArgumentParser
from os import path
from sys import stderr
import pybedtools
from pybedtools.featurefuncs import greater_than, less_than

parser = ArgumentParser(description='File input and options to process data')
parser.add_argument('-b', '--bed_file', help='Bed file containing meals of a single animal/track', required=True)
parser.add_argument('-ct', '--coordinates_time_to_bed', help='Coordinates to generate a bed of a given time period',
                    required=True, nargs='+')
parser.add_argument('-bins', '--bins', help='Lengths used to bin the meals by meal duration ', required=True, \
                    nargs='+', type=str)

args = parser.parse_args()

print >> stderr, "Bed file feeding behavior: %s", args.bed_file
print >> stderr, "Coordinates for time period intersection: %s" % ' '.join(str(c) for c in args.coordinates_time_to_bed)
print >> stderr, "Values to perform the binning: %s" % ' '.join(str(b) for b in args.bins)

bed_file_meals = args.bed_file
bed_meals = pybedtools.BedTool(bed_file_meals)

coord_period = args.coordinates_time_to_bed
bed_period = [("chr1", coord_period[0], coord_period[1], " ", 0, 0)]

bed_meals_by_period = bed_meals.intersect(pybedtools.BedTool(bed_period))

bins = [ int(i) for i in args.bins ]

for i, b in enumerate(bins):
    if i == 0:
        bed_meals_by_period.filter(greater_than, 0).filter(less_than, b+1).saveas("0_" + str(b) + "_" + path.basename(bed_file_meals)) # < 31  <=30, include 30
    elif i == len(bins) - 1:
        bed_meals_by_period.filter(greater_than, bins[i-1]).filter(less_than, b+1).saveas(str(bins[i-1]) + "_" + str(b) + "_" + path.basename(bed_file_meals)) # > 30, 30 not included, <121, <=120 include 120
        bed_meals_by_period.filter(greater_than, b).saveas(str(b) + "_" + path.basename(bed_file_meals)) # >120
    else:
        bed_meals_by_period.filter(greater_than, bins[i - 1]).filter(less_than, b+1).saveas(str(bins[i-1]) + "_" + str(b) + "_" + path.basename(bed_file_meals)) # > 30 # < 121 <= 120
