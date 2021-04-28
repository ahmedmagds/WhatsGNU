#!/usr/bin/env python3
# Copyright (C) 2019 Ahmed M. Moustafa
#########################################################################################
# LICENSE
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# DATE CREATED: December 01, 2020

# AUTHOR: Ahmed M Moustafa

# CONTACT1: moustafaam@chop.edu
# CONTACT2: ahmedmagdy2009@hotmail.com

# AFFILIATION: Pediatric Infectious Disease Division, Children’s Hospital of Philadelphia,
# Abramson Pediatric Research Center, University of Pennsylvania, Philadelphia,
# Pennsylvania, 19104, USA

import os
import sys
import time
import argparse
from collections import Counter
from collections import OrderedDict

PARSER = argparse.ArgumentParser(
    prog="Complete_Reference_selector.py",
    description="find complete genomes from NCBI",
)
PARSER.add_argument(
    "NCBI_csv",
    type=str,
    help="GenBank csv for a species",
)
PARSER.add_argument(
    "genomes", type=str, help="A folder or text file of genomes to analyze"
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
#####Mash files to be processed######
TOPGEN_REPORTS = ARGS.genomes
TOPGEN_LIST = []
try:
    for file in os.listdir(TOPGEN_REPORTS):
        if file.endswith(".txt"):
            TOPGEN_LIST.append(TOPGEN_REPORTS + file)
    print(
        "You provided folder of {} files".format(len(TOPGEN_LIST))
    )
    if len(TOPGEN_LIST) == 0:
        PARSER.exit(
            status=0,
            message="The directory did not have any file.txt\n",
        )
except:
    TOPGEN_LIST.append(TOPGEN_REPORTS)
    print("You provided one file")
############################
GCA_complete = {}
LIST_FILE_OBJECT = open(ARGS.NCBI_csv, 'r')
for line in LIST_FILE_OBJECT:
    if line.startswith("#"):
        continue
    line_list = line.split(",")
    GCA_complete[line_list[5].split('"')[1].split('.')[0]] = line_list[6]
#print(GCA_complete)

#############################
cwd = os.getcwd()
TIMESTR = time.strftime("%m%d%Y_%H%M")
Name_Time = 'Completeness_report_{}.txt'.format(TIMESTR)
file_report = os.path.join(cwd, Name_Time)
output_file = open(file_report, "w")
output_file.write("Strain\tCompleteness\n")
#############################
CNT = Counter()
SCORES = Counter()
for TOPGEN_REPORT in TOPGEN_LIST:
    TOPGEN_REPORT_OBJECT = open(TOPGEN_REPORT, "r")
    FILE_NAME = TOPGEN_REPORT.rsplit(OS_SEPARATOR, 1)[-1]
    for line in TOPGEN_REPORT_OBJECT:
        line = line.rstrip()
        GCA_acc = line.split('.')[0]
        output_file.write(GCA_acc+'\t'+GCA_complete[GCA_acc]+'\n')
