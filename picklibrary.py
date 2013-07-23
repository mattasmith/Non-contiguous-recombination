#! /usr/local/bin/python
"""Script for building the sequences of a non-contiguous library, given a library design.

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
"""


import pickle
from sys import path
import sys
import os
path.append('./tools')
import sequence_tools
import PDB_tools
import chimera_tools


a = sys.argv[1]

## check the files are there
if os.path.isdir('./output') == False:
    exit('Error:  cannot find \''+a+'.output\' output file')
if os.path.isfile('./output/'+a+'.output') == False:
    exit('Error: cannot find \''+a+'.output\' output file')
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

nblocks = int(a.split('library')[1].split('_')[0])

## find the muscle version
muscle_file = [f for f in os.listdir('./tools/muscle') if 'muscle' in f]
muscle_version = muscle_file[0]

## load the alignment
alignment, contacts = pickle.load(open('./alignment_and_contacts.pkl'))
numberofparents = len(alignment[0])

## load the assignments file
assignments_file = './output/'+a+'.output'
## read in the library file - this is residue number and block letter
assignments_line = [l for l in open(assignments_file).read().split('\n') if len(l)>0 and l[0]!='#']
assignment = [ord(l.split('\t')[2]) - ord('A') for l in assignments_line if l.split('\t')[2] !='-']
nodes_outputfile = [int(l.split('\t')[1])-1 for l in assignments_line if l.split('\t')[2] !='-'] # -1 because counting 0,1,2...

## Create a folder to store the data for this specific library
library_path = './picked_libraries/'+a
if os.path.isdir('./picked_libraries') == False:
    os.mkdir('./picked_libraries')
if os.path.isdir(library_path) == False:
    os.mkdir(library_path)

## load in the number of contact maps used
num_contact_maps = pickle.load(open('num_contact_maps.pkl'))

## give warning if residues have no SCHEMA contacts!!
for res,align in enumerate(alignment):
    ## if residue is not conserved, but not in nodes
    if (len(set(align)) != 1) and (res not in nodes_outputfile):
        print ('Warning: residue ' + str(res+1) + ' is not conserved but has no SCHEMA contacts. By default it is been assigned to block A.')
        nodes_outputfile.append(res)
        nodes_outputfile.sort()
        assignment.insert(nodes_outputfile.index(res),0) # assign unspecified residues to block A by default
# end for res, align

print('\nCurrent block sizes:')
for i in range(nblocks):
    print(chr(ord('A')+i)+'\t'+str(sum([1 for b in assignment if b==i])))


## calculate the SCHEMA E for each contact
SCHEMA_E = chimera_tools.calculate_contact_SCHEMA_E_2(alignment,contacts)


#unique_contacts = list(set(contacts)) # unique list of all the contacts

(library_chimerablocks, library_seq) = chimera_tools.build_library(assignment, alignment, nodes_outputfile)
library_seq_concat = [ ''.join(seq) for seq in library_seq ]


## calculate the E and m value for each chimera
m_chi = {}
E_chi = {}
for i, chimera in enumerate(library_seq):
    (m_chi[i],par) = chimera_tools.calculate_m(zip(*alignment),chimera)
    E_chi[i] = chimera_tools.calculate_SCHEMA_E(contacts, alignment, chimera)/num_contact_maps

## check that the average E is the same as the average E_chi, same for m!!
average_E_chi = sum(E_chi.values())/float(len(E_chi))
average_m_chi = sum(m_chi.values())/float(len(m_chi))
average_E = chimera_tools.score_library(assignment, SCHEMA_E, nodes_outputfile, numberofparents)/num_contact_maps
print '<E> = ', average_E
print '<m> = ', average_m_chi

## output E, m, library_chimerablocks, library sequences to text file for Matlab to find the optimal subset for synthesis
f = open(library_path+'/chimeras.output','w')
for i, chimera in enumerate(library_chimerablocks):
    out = ''.join([str(block+1) for block in chimera])+' '+str(E_chi[i])+' '+str(m_chi[i])+' '+library_seq_concat[i] 
    f.write(out+'\n')
f.close()


