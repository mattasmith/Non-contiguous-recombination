Non-contiguous-recombination
============================

This is a software package for protein engineers. It uses protein structure and sequence information to aid researchers in designing chimeric proteins and recombination libraries. 

These tools identify elements of structure that can be shuffled between two or more proteins. 

SCHEMA and non-contiguous recombination were developed in the laboratory of Frances H. Arnold at the California Institute of Technology.


Contents
============================

1. References
2. Installation
3. Setup
4. Execution
5. Output
6. Things to be careful about



1. References
============================

Smith, M. et al., Chimeragenesis of distantly-related proteins by non-contiguous recombination, Protein Science 22(2):231-238 (2013).
Voigt, C. et al., Protein building blocks preserved by recombination, Nature Structural Biology 9(7):553-558 (2002).

This script uses George Karypis' hmetis graph partitioning algorithm:

Karypis, G. et al. Multilevel Hypergraph Partitioning: Applications in VLSI Domain, 34th Design and Automation Conference, 526-529, (1997)
Karypis, G. et al. Multilevel k-way Hypergraph Partitioning, 36th Design Automation Conference, 343-348, 1999.

And Robert Edgar's MUSCLE alignment program:

Edgar, R.C. MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97.


2. Installation
============================

Non-contiguous recombination is written in python 2.6 (www.python.org). 
It uses two external programs: MUSCLE and hmetis.

MUSCLE can be downloaded from: http://www.drive5.com/muscle/downloads.htm
It should be unpacked and the executable placed in tools/muscle

hmetis 1.5 can be downloaded from: http://glaros.dtc.umn.edu/gkhome/metis/hmetis/download
It should be unpacked and the hmetis folder placed in tools


3. Setup
============================

The init.txt file specifies the parameters to be used by the non-contiguous recombination script.

Number of blocks - the number of blocks in the designed libraries. It can either be a number (e.g. 8) or a range of numbers (e.g. 2-6) for designing a range of libraries with different block sizes.

Find all PDB structures - the script should search the PDB database for structures ( = 1 ) or the user will provide structures ( = 0 ). Structures should be in a folder called 'structures'. Separate PDB files should be used for separate chains.


4. Execution
============================

Run: 
python ncr.py

For a list of the chimeric proteins in a given library, run:
python picklibrary.py libraryX
where libraryX is the name of the given library (e.g. library12_2)


5. Output
============================

ncr.py will create a set of library designs in a directory 'output'. Residues (numbered by  the alignment) are assigned to a block by the letters A,B,C,D, etc. Conserved residues (and residues with no SCHEMA contacts) are represented by '-'.

The script will also create a .csv file that contains the average SCHEMA E and average number of mutations for each library design, along with the estimated fraction of the library that is folded and the number of mutations in each block.


6. Things to be careful about
============================

* A good sequence alignment is essential. Structural information can help alignments if the sequence identity is not that high. 

* Residues that do not have any SCHEMA contacts will not be assigned to a block, even if the residues are not conserved. There is no structural information for such residues, so the block choice is left up to the user.

* Blocks may not be contiguous in 3D space. Indeed, certain sets of residues in surface loops can be disconnected from the rest of the contacts and these residues will be assigned to the smallest block to balance the mutations in each block.

