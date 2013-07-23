from sys import path
import os
path.append('./tools')
import sequence_tools
import PDB_tools


## read the parents from the fasta alignment file
par_names,parents,parent_alignment = sequence_tools.read_fasta_alignment('./alignment.fasta')



## do PDB search on each parent, and merge results
print 'searching the PDB for related structures'
pdb_ids = []
for par in parents:
    pdb_ids.extend(PDB_tools.pdb_search(par))

pdb_ids = sorted(set(pdb_ids))



## download each PDB, parse, and save each chain
structure_path = './structures/'
if os.path.isdir(structure_path) == False:
    os.mkdir(structure_path)
for id in pdb_ids:
    PDB_tools.download_parse_save(id,structure_path)

