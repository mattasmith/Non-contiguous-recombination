#! /usr/local/bin/python
"""Script for designing a set of non-contiguous recombination libraries for site-directed, structure-guided homologous recombination.

    ******************************************************************
    Copyright (C) 2011  Matt Smith, California Institute of Technology

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    *******************************************************************


SCHEMA and non-contiguous recombination were developed in the laboratory of Frances H. Arnold at the California Institute of Technology.

References:

Smith, M.A. et al., Chimeragenesis of distantly-related proteins by non-contiguous recombination, Protein Science 22(2):231-238 (2013).
Voigt, C.A. et al., Protein building blocks preserved by recombination, Nature Structural Biology 9(7):553-558 (2002).
Karypis, G. et al., Multilevel Hypergraph Partitioning: Applications in VLSI Domain, 34th Design and Automation Conference, 526-529, (1997).
Karypis, G. et al., Multilevel k-way Hypergraph Partitioning, 36th Design Automation Conference, 343-348, 1999.
Edgar, R.C., MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97.
"""

from subprocess import Popen
from sys import exit
from sys import path
path.append('./tools')
import os
from make_alignment_and_contacts2 import make_alignment_and_contacts2
from run_shmetis import run_shmetis
from PDB_tools import split_into_chains


## information about non-contiguous recombination
print('\n********************* Non-contiguous recombination *********************\n')
print('Written by Matt Smith, 2011.\n')
print('SCHEMA and non-contiguous recombination were developed in the \nlaboratory of Frances H. Arnold at the California Institute of Technology.\n')

## check the files are there
if os.path.isfile('alignment.fasta') == False:
    exit('Error: cannot find \'alignment.fasta\' alignment file')
if os.path.isfile('init.txt') == False:
    exit('Error: cannot find \'init.txt\' setup file.')
if os.path.isdir('./tools') == False:
    exit('Error: the non-contiguous recombination tools are missing')
if len([f for f in os.listdir('./tools/muscle') if 'muscle' in f])==0:
    exit('Error: cannot find MUSCLE in \'tools/muscle\'')
if len([f for f in os.listdir('./tools/muscle') if 'muscle' in f])>1:
    exit('Error: please provide just one MUSCLE executable in \'tools/muscle\'')
if len([f for f in os.listdir('./tools') if 'hmetis' in f])==0:
    exit('Error: cannot find hmetis package in \'tools\'')
if len([f for f in os.listdir('./tools') if 'hmetis-1.5' in f])>1:
    exit('Error: please provide just one hmetis package in \'tools\'')

## load in the initial file
data = [s for s in open('init.txt').read().split('\n') if (len(s)>0 and s[0]!='#')]
for i,datum in enumerate(data):
    if 'Number of blocks' in datum.split(' = ')[0]:
        numberofblocks_str = datum.split(' = ')[1]
        if '-' in numberofblocks_str:
            numberofblocks_min = int(numberofblocks_str.split('-')[0])
            numberofblocks_max = int(numberofblocks_str.split('-')[1])
        else:
            numberofblocks_min = int(numberofblocks_str)
            numberofblocks_max = int(numberofblocks_str)
    if 'Find all PDB structures' in datum.split(' = ')[0]:
        searchPDB = int(datum.split(' = ')[1])
# end for i, datum


## find the muscle version
muscle_file = [f for f in os.listdir('./tools/muscle') if 'muscle' in f]
muscle_version = muscle_file[0]

## find the hmetis version
hmetis_file = [f for f in os.listdir('./tools') if 'hmetis-1.5' in f]
hmetis_version = hmetis_file[0]


## download all available structures or check user pdb files
if searchPDB == 1:
    Popen('python ./tools/search_download_save.py',shell=True).wait()
else:
    if os.path.isdir('./structures') == False:
        exit('Error: you need to provide at least one pdb in a folder called \'structures\'')
    elif len([f for f in os.listdir('./structures') if os.path.splitext(f)[-1].lower() == '.pdb'])==0:
        exit('Error: there are no pdbs in \'structures\'')
    else:
        print('Structures provided by user:')
        structurefilelist = os.listdir('./structures')
        for filename in structurefilelist:
            if os.path.splitext(filename)[-1].lower() == '.pdb':
                print filename
                split_into_chains(filename,'./structures/')


## create the contact maps - one for each parent (if the parent has a structure)
num_contact_maps = make_alignment_and_contacts2(muscle_version)

## formulate and solve with graph partitioning
print ('\nDesigning libraries...')
run_success = run_shmetis(num_contact_maps, numberofblocks_min, numberofblocks_max, hmetis_version)

## done!
