#!/usr/bin/env python3
# PROGRAM: WhatsGNU is a Python3 program that ranks protein sequences in a genome
# faa file generated from annotation programs based on the number of observed
# exact protein matches in a public or private database.

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

# DATE CREATED: February 28, 2022

# AUTHOR: Ahmed M Moustafa

# CONTACT1: moustafaam@email.chop.edu
# CONTACT2: ahmedmagdy2009@hotmail.com

# AFFILIATION: Pediatric Infectious Disease Division, Children’s Hospital of Philadelphia,
# Abramson Pediatric Research Center, University of Pennsylvania, Philadelphia,
# Pennsylvania, 19104, USA

# CITATION1: Ahmed M Moustafa and Paul J Planet
# WhatsGNU: a tool for identifying proteomic novelty
# Genome Biology(2020)21:58, doi: https://doi.org/10.1186/s13059-020-01965-w.

import sys
import os
import argparse
import tempfile
import subprocess
from collections import defaultdict


PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_alleles_strain_mapper.py", description="WhatsGNU_alleles_strain_mapper script for WhatsGNU v1.4.",)
#PARSER.add_argument("traits", type=str, help="traits file")
PARSER.add_argument("ids_hits", type=str, help="path to ids_hits WhatsGNU output file")
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
####################################
ids_hits_obj = open(ARGS.ids_hits,'r')
op_file = ARGS.ids_hits.split('.')[0] + '_mapping.csv'
op_obj = open(op_file,'w')
ids_hits_obj.readline()
alleles = []
Strains_Dict = defaultdict(list)
for line in ids_hits_obj:
    line_lst = line.rstrip().split('\t')
    allele_name = line_lst[0]
    alleles.append(allele_name)
    strain_names = [i.split('|')[0] for i in line_lst[1:]]
    for j in strain_names:
        Strains_Dict[j].append(allele_name)
#print(Strains_Dict)
final_dict = {}
op_obj.write('Isolate,{}\n'.format(','.join(alleles)))
for strain in Strains_Dict:
    strain_alleles = Strains_Dict[strain]
    tmp_lst = []
    for allele in alleles:
        if allele in strain_alleles:
            hit = '1'
        else:
            hit = '0'
        tmp_lst.append(hit)
    final_dict[strain] = tmp_lst
#print(final_dict)
summary_dict = defaultdict(list)
for x in final_dict:
    op_obj.write('{},{}\n'.format(x,','.join(final_dict[x])))
    summary_dict[''.join(final_dict[x])].append(x)
op_file2 = ARGS.ids_hits.split('.')[0] + '_summary.txt'
op_obj2 = open(op_file2,'w')
op_obj2.write('alleles_order: {}\n'.format(','.join(alleles)))
for comb in summary_dict:
    #print(comb,len(summary_dict[comb]))
    op_obj2.write('{}\t{}\n'.format(comb,len(summary_dict[comb])))
