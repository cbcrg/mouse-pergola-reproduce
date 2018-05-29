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
### Jose A Espinosa. CSN/CB-CRG Group. Jan 2018              ###
################################################################
### Script creates BedTools objects from a file containing   ###
### mice feeding behavior and uses tools from pybedtools to  ###
### intersect them with day phases (light/dark).             ###
### Generates a bed file for each track with the result of   ###
### the above described operations.                          ###
################################################################

from argparse import ArgumentParser, FileType
from os import path, getcwd, makedirs, chdir
from shutil import rmtree
from sys import stderr, stdout
import subprocess
import pybedtools
from pergola import mapping
from pergola import intervals
from pergola import tracks

_stats_available = ['mean', 'count', 'sum', 'max', 'min', 'median']
_behaviors_available = ['feeding', 'drinking']

parser = ArgumentParser(description='File input and options to process data')
parser.add_argument('-f', '--file_mice_behavior', help='Files containing feeding/drinking behaviors of mice',
                    required=True, nargs='+')
parser.add_argument('-m', '--mapping_file', help='Mapping file to read mice behavioral files', required=True)
parser.add_argument('-p', '--phases_exp_file', help='Experiment phases file', required=True)
parser.add_argument('-pl', '--phases_exp_file_long', help='Experiment phases file whole trajectory', required=False)
parser.add_argument('-mp', '--mapping_phases', help='Mapping file of experiment phases file', required=True)
parser.add_argument('-s', '--statistic', help='Choose one of the possible statistical available on Bedtools map option',
                    required=True, choices=_stats_available)
parser.add_argument('-b', '--behavioral_type',
                    help='Choose whether to work with drinking or feeding mice behavioral data',
                    required=True, choices=_behaviors_available)

args = parser.parse_args()

print >> stderr, "Statistic to be calculated: %s" % args.statistic
print >> stderr, "Working with mice %s behavioral data" % args.behavioral_type

# Statistic to calculate
statistic = args.statistic

## Dictionary to set colors of each type of food
# food_sc    orange
# food_fat    black
# water    blue
# saccharin    red

### Feeding data
if args.behavioral_type == "feeding":
    data_type_1 = "food_sc"
    data_type_2 = "food_fat"
    data_type_joined = ["food_sc", "food_fat"]
    data_type_col = {data_type_1: 'blue', data_type_2: 'red'}

### Drinking data
elif args.behavioral_type == 'drinking':
    data_type_1 = "water"
    data_type_2 = "saccharin"
    data_type_col = {data_type_1: 'blue', data_type_2: 'red'}
else:
    print >> stderr, "Behavioral data type not available in script, please try again with \"drinking\" or \"feeding\""

mapping_data = mapping.MappingInfo(args.mapping_file)

behavior_mice = dict()
end_time = -10000
data_read_all_batches = None

for f in args.file_mice_behavior:
    int_data = intervals.IntData(f, map_dict=mapping_data.correspondence)

    data_read = int_data.read(relative_coord=True)

    chr_file_n = "chrom"
    mapping.write_chr_sizes(data_read, file_n=chr_file_n)
    chr_file = chr_file_n + ".sizes"

    if end_time < int_data.max - int_data.min:
        end_time = int_data.max - int_data.min

    if not data_read_all_batches:
        data_read_all_batches = data_read
    else:
        data_read_all_batches = tracks.merge_tracks(data_read_all_batches, data_read)

stdout.write(str(end_time))

list_ctrl = [1, 3, 5, 7, 9, 11, 13, 15, 17]
list_hf = [2, 4, 8, 10, 12, 14, 16, 18]

bed_dict = dict()

bed_dict['ctrl'] = {}
bed_dict['hf'] = {}

bed_dict['ctrl']['food'] = data_read_all_batches.convert(mode="bed", data_types=data_type_joined, data_types_actions="all",
                                                              color_restrictions=data_type_col, tracks=list_ctrl)
bed_dict['hf']['food'] = data_read_all_batches.convert(mode="bed", data_types=data_type_joined, data_types_actions="all",
                                                            color_restrictions=data_type_col, tracks=list_hf)

# Generating sequence of days without light and dark phases #del
mapping.write_period_seq(end=end_time, delta=86400, tag="day", name_file="days_seq")
# mapping.write_cytoband(start=0, end=end_time, start_phase="dark", track_line=False, lab_bed=False)
mapping.write_cytoband(start=26953, end=end_time, start_phase="light", track_line=False, lab_bed=False)
days_bed_f = "days_seq.bed"
days_bed = pybedtools.BedTool(days_bed_f)

## weeks
mapping.write_period_seq(end=end_time, delta=604800, tag="week", name_file="weeks_seq", lab_bed=True, track_line=True)
weeks_bed_f = "weeks_seq.bed"

## Reading experimental phases from csv file
mapping_data_phases = mapping.MappingInfo(args.mapping_phases)

int_exp_phases = intervals.IntData(args.phases_exp_file, map_dict=mapping_data_phases.correspondence)
data_read_exp_phases = int_exp_phases.read(relative_coord=True)
d_exp_phases_bed2file = data_read_exp_phases.convert(mode="bed", data_types_actions="all")
d_exp_phases_bed2file[d_exp_phases_bed2file.keys()[0]].save_track(bed_label="True", name_file="exp_phases")
d_exp_phases_bed = data_read_exp_phases.convert(mode="bed", data_types_actions='one_per_channel')

int_exp_phases_long = intervals.IntData(args.phases_exp_file_long, map_dict=mapping_data_phases.correspondence)
data_read_exp_phases_long = int_exp_phases_long.read(relative_coord=True)
# generates a bed for habituation and one for development
d_exp_phases_bed_long = data_read_exp_phases_long.convert(mode="bed", data_types_actions='one_per_channel')




def length_bed (b):
    """
    calculates length of bed items and dumps it on the score field
    """
    b.score = b.end - b.start
    return b

def rate_bed (b):
    """
    calculates the eating rate by divinding by score by the length of bed items and dumps it on the score field
    """
    b.score = str(float(b.score) / (b.end - b.start))
    return b

for exp_group, dict_exp_gr in bed_dict.iteritems():

    for data_type, dict_bed in dict_exp_gr.iteritems():
        for tr, bed in dict_bed.iteritems():
            bed_BedTools = bed.create_pybedtools()

            ### bouts per week
            weeks_bed = pybedtools.BedTool(weeks_bed_f)
            weeks_bed.intersect(bed_BedTools, c=True).moveto("tr_"+ tr[0] + '.' + exp_group + '.' + data_type + '.' + "n_bouts.tbl")

            ### mean value
            statistic="mean"
            weeks_bed.map(bed_BedTools, c=5, o=statistic, null=0).moveto("tr_"+ tr[0] + '.' + exp_group + '.' + data_type + '.' + 'mean.tbl')

            ### eating rate
            weeks_bed.map(bed_BedTools.each(rate_bed), c=5, o=statistic, null=0).moveto("tr_"+ tr[0] + '.' + exp_group + '.' + data_type + '.' + 'rate.tbl')

            ### meal duration
            weeks_bed.map(bed_BedTools.each(length_bed), c=5, o=statistic, null=0).moveto("tr_"+ tr[0] + '.'+ exp_group+ '.' + data_type + '.' + 'duration.tbl')

            ### intermeal interval
            weeks_bed.map(weeks_bed.intersect(bed_BedTools.complement(g=chr_file)).each(length_bed).sort(), c=5, o=statistic, null=0)\
                .moveto("tr_"+ tr[0] + '.' + exp_group + '.' + data_type + '.' + 'intermeal_time.tbl')

            for key, bed_phase in d_exp_phases_bed.iteritems():
                exp_phase = key[1]

                if not path.isfile(exp_phase + ".bed"):
                    bed_phase.create_pybedtools().saveas(exp_phase + ".bed")

                # bouts per experimental phase
                exp_phase_events_bed = bed_BedTools.intersect(pybedtools.BedTool(exp_phase + ".bed"))

                # dark phases per experimental phase
                pybedtools.BedTool(exp_phase + ".bed").intersect(pybedtools.BedTool("phases_dark.bed")).sort().saveas(exp_phase + "_dark.bed")

                ###################
                # Generate mean value of the whole record after intersecting with phase
                if exp_phase_events_bed.count() == 0:
                    # When there is any interval we set the mean to zero
                    list_no_intervals = [("chr1", 0, 1, "no_hits", 0, 0)]
                    pybedtools.BedTool(list_no_intervals).saveas(
                        'tr_' + exp_group + '.' + '.'.join(tr) + ".day." + exp_phase + ".bed")
                else:
                    days_bed.map(exp_phase_events_bed, c=5, o=statistic, null=0).intersect(
                        pybedtools.BedTool(exp_phase + ".bed")).saveas(
                        'tr_' + exp_group + '.' + '.'.join(tr) + ".day." + exp_phase + ".bed")


for key, bed_phase in d_exp_phases_bed_long.iteritems():
    exp_phase = key[1]

    if not path.isfile(exp_phase + "_long.bed"):
        bed_phase.create_pybedtools().saveas(exp_phase + "_long.bed")
        # dark phases per experimental phase
    pybedtools.BedTool(exp_phase + "_long.bed").intersect(pybedtools.BedTool("phases_dark.bed")).sort().saveas(exp_phase + "_dark_long.bed")

        # data_read_exp_phases