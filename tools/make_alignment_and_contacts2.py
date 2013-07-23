import pickle
from os import listdir
from sys import path
path.append('./tools')
import sequence_tools
import PDB_tools
import contact_tools

def make_alignment_and_contacts2(muscle_version):
    ## read the parents and their alignment from the fasta alignment file
    par_names,parents,parent_alignment = sequence_tools.read_fasta_alignment('./alignment.fasta')
    
    ## cycle through every pdb, keep if it has significant sequence identity
    structure_path = './structures/'
    pdbs = [f for f in listdir(structure_path) if f.split('.')[-1]=='pdb']
    structure_list = []
    structure_list_noinsert = []
    par_identity = []
    for par in range(len(parents)):
        par_identity.append([])
        structure_list.append([])
    
    for pdb in pdbs:
        structure = PDB_tools.PDB_structure(structure_path+pdb)
        str_align = sequence_tools.muscle_add_sequence(parent_alignment,structure.sequence,muscle_version)
        for par in range(len(parents)):
            iden = len([pos for pos in str_align if pos[par]==pos[-1]])/float(len(str_align))
            if iden > 0.9: ## the structure has to have at least 90% identity with one of the parents
                structure.get_contacts(['H'])  ## atom types to ignore from search
                structure_list[par].append(structure)
                par_identity[par].append(iden)
    
    ### Treat each parent separately, find the structure for each parent and create separate contact matrices for each one.
    ## Beware of different numberings of the different contact matrices
    ##  Sum the contact matrices for the graph partition, then divide <E> by the number of contact matrices used.
    ## Contacts are effectively being weighted by the number of parents that they appear in.
    
    alignment = []
    numbered_alignment = []
    alignment_contacts = []
    contacts = []
    
    ## for each parent
    ## add the structures from structure_list to the parent_alignment
    ## make a numbered_alignment of parents and structures
    ## remove parents and gaps in all parents from alignment
    for par in range(len(parents)):
        alignment.append([])
        numbered_alignment.append([])
        alignment_contacts.append([])
        contacts.append([])
        ## add the PDB sequences to the parent_alignment
        alignment[par] = parent_alignment
        for struct in structure_list[par]:
            alignment[par] = sequence_tools.muscle_add_sequence(alignment[par],struct.sequence,muscle_version)
        ## print the alignment
        #sequence_tools.print_alignment(alignment[par]) 
        ## make an alignment of residue numbers (so we can map alignment postions back to structure postions)
        numbered_alignment[par] = sequence_tools.number_alignment(alignment[par])
        # remove the first len(parent) sequences, and remove postions that have gaps in all parents    - this will make numbered_alignment same length as parent_alignment!
        numbered_alignment[par] = [pos[len(parents):] for pos in numbered_alignment[par] if not all([p=='-' for p in pos[:len(parents)]])]
        if len(numbered_alignment[par][0]) == 0:
            print('Parent '+par_names[par]+' does not have any structures.  No contact map exists for this parent.')
        else:
            print('Parent '+par_names[par]+' has at least one structure.  A contact map will be used for this parent. Building contact map...')
        
        ## count contacts
        alignment_contacts[par] = contact_tools.find_alignment_contacts(structure_list[par],numbered_alignment[par]) # takes a long time
        contacts[par] = sorted([k for k in alignment_contacts[par].iterkeys() if alignment_contacts[par][k]/float(len(structure_list[par])) >= 0.5 ]) # only keep the contacts that are observed in at least 50% of the structures 
        
    ## convert a summed contacts file, total_contacts
    total_contacts = []
    no_struct_count = 0
    for i in range(len(parents)):
        total_contacts = total_contacts + contacts[i]
        if len(contacts[i]) == 0:
            no_struct_count = no_struct_count + 1
    
    ## calculate the number of contact maps used. Divide <E> by this number to normalize it to SCHEMA E.
    num_contact_maps = len(parents) - no_struct_count
    # print(str(no_struct_count)+' parent structures are missing.  Divide <E> by '+str(len(parents)-no_struct_count))
    
    
    
    ##################################### save to a pickle ########################################
    ## take parent alignment - numbered alignment is the same length as the parent alignment
    ## take total contacts - want all contacts 
    pickle.dump((parent_alignment,total_contacts,),open('alignment_and_contacts.pkl','w'))
    
    ## save the contacts[] file
    ## for contact venn diagram and calculating <E> from each contact map
    pickle.dump(contacts,open('contacts_each_parent.pkl','w'))
    
    ## return the num_contact_maps
    ## for normalizing <E> later
    pickle.dump(num_contact_maps,open('num_contact_maps.pkl','w'))
    return num_contact_maps

