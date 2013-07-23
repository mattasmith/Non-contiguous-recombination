

def calculate_contact_SCHEMA_E_2(alignment,contacts):
    '''All SCHEMA contacts as dictionary with E value'''
    # update 1: multiple contacts
    # update 2: backbone continuity in loops, add contacts to loops to stop them being cut
    SCHEMA_E = {}
    for contact in contacts:
        WT_contacts = tuple((alignment[contact[0]][par],alignment[contact[1]][par]) for par in range(len(alignment[0]))) # what is seen in the three parents
        E = 0
        for p1 in alignment[contact[0]]:
            for p2 in alignment[contact[1]]:
                if (p1,p2) not in WT_contacts:
                    E+=1
        if E>0:
            if SCHEMA_E.has_key(contact):
                SCHEMA_E[contact] += E
            else:
                SCHEMA_E[contact] = E
    
    # backbone continuity in loops
    # insertions must be contained within a loop
    # find insertion indices
    insertion_index = {}
    for i,align in enumerate(alignment):
        if '-' in align:
            insertion_index[i] = align
    # find continuous parts of insertions
    cont_insertion_index = [c for c in insertion_index.keys() if insertion_index.has_key(c+1)]
    # add these to the SCHEMA_E, with a high energy value
    for i in cont_insertion_index:
        if SCHEMA_E.has_key((i,i+1)):
            SCHEMA_E[(i,i+1)] = SCHEMA_E[(i,i+1)]+100  # large number - will not break in partitioning
        else:
            SCHEMA_E[(i,i+1)] = 100
        
    return SCHEMA_E


def score_library(library, SCHEMA_E, nodes, numberofparents):
    '''given a library, return average SCHEMA E'''
    score = 0
    for contact in SCHEMA_E.iterkeys():
        pos1 = nodes.index(contact[0])
        pos2 = nodes.index(contact[1])    
        if library[pos1]!=library[pos2]:
            score += SCHEMA_E[contact]
    ave_E = score/float(numberofparents*numberofparents)
    return ave_E


def calculate_SCHEMA_E(contacts, alignment, sequence): # check matches chimera_tools function - yes
    '''given a chimera, return SCHEMA E'''
    # update 1: backbone continuity in loops, add contacts to loops to stop them being cut
    score = 0
    numberofparents = len(alignment[0])
    parents = zip(*alignment)
    ## do not use SCHEMA_E.keys() when counting contacts - mult contacts are counted only once!
    for contact in contacts: 
        pos1 = contact[0]
        pos2 = contact[1]
        # if the residue pair is not found in any of the parents
        if len([p for p in parents if not ((p[contact[0]]==sequence[contact[0]]) & (p[contact[1]]==sequence[contact[1]]))]) == numberofparents:
            score += 1
    # backbone continuity in loops
    # insertions must be contained within a loop
    # find insertion indices
    insertion_index = {}
    for i,align in enumerate(alignment):
        if '-' in align:
            insertion_index[i] = align
    # find continuous parts of insertions
    cont_insertion_index = [c for c in insertion_index.keys() if insertion_index.has_key(c+1)]
    # add these to the SCHEMA_E, with a high energy value
    for i in cont_insertion_index:
         # if the residue pair is not found in any of the parents
        if len([p for p in parents if not ((p[i]==sequence[i]) & (p[i]==sequence[i]))]) == numberofparents:
            score += 100
    E = score
    return E


def calculate_SCHEMA_E_broken_contacts(parents,contacts,sequence):
    # first need to calculate the contacts that are seen in the parents
    parent_contacts = {}
    for contact in contacts:
        res_id = set()
        for parent in parents:
            res_id.add(tuple([parent[i] for i in contact]))
        parent_contacts[contact] = tuple(res_id)
    # scan through the sequence and calculate the number of broken contacts
    SCHEMA_E = 0
    broken_contacts = []
    for contact in contacts:
        if not tuple([sequence[i] for i in contact]) in parent_contacts[contact]:
            SCHEMA_E+=1
            broken_contacts.append(tuple([sequence[i]+str(i) for i in contact]))
    return SCHEMA_E,broken_contacts


def build_library(library, alignment, nodes):
    '''given a library, return the sequences of the members'''
    import itertools
    ## build the library using alignment, assignment and nodes!
    numberofparents = len(alignment[0])
    numberofblocks = max(library)+1
    library_seq = []
    library_chimerablocks = []
    ## construct all the chimeras.  
    for chimerablocks in itertools.product( range(numberofparents), repeat=numberofblocks ):
        library_chimerablocks.append(chimerablocks)
        ## set chimera equal to one of the parents
        ## so conserved/not contacting residues are assigned
        chimera = [ alignment[i][0] for i in range(len(alignment)) ]
        ## adjust chimera for residues that differ
        for i in range(len(nodes)):
            chimera[nodes[i]] = alignment[nodes[i]][chimerablocks[int(library[i])]]
        library_seq.append(chimera)

    return library_chimerablocks, library_seq
 

def calculate_m(parents,sequence):
    dist = [len([1 for p in zip(par,sequence) if p[0]!=p[1]]) for par in parents]
    M = min(dist)
    par = dist.index(M)
    return M,par


def calculate_average_m(library, alignment, nodes):
    '''given a library, return average m'''
    library_chimerablocks, library_seq = build_library(library, alignment, nodes)
    ## find <m> for the library
    ## take each chimera and find distance from closest parent 
    average_m = 0.0
    for chi in library_seq:
        m = calculate_m(zip(*alignment),chi)
        average_m = average_m + m[0]
 
    average_m = average_m / len(library_seq)
    return average_m


def calculate_block_alignment(library, alignment, nodes):
    '''create a block alignment list for each residue of block it belongs to, and the parent alignment'''
    # make all residues equal to block -1
    # this will assign conserved residues to their own block
    full_assignment = [ -1 for res in range(len(alignment)) ]
    # assign residues not conserved to their blocks
    for i,res in enumerate(nodes):
        full_assignment[res] = library[i]
    # add full alignment to front of alignment
    block_alignment = zip(*alignment)
    block_alignment.insert(0, tuple(full_assignment))
    block_alignment = zip(*block_alignment)
    return block_alignment


def write_solution(outputfile, alignment, assignment, nodes):
    '''save the graph partition output to a text file'''
    f = open('output/'+outputfile,'w')
    f.write('# Non-contiguous recombination solution, Matt Smith 2011\n')
    f.write('# Residue\tBlock')
    for a,res in enumerate(alignment):
        if a in nodes:
            f.write('\n\t'+str(a+1)+'\t'+chr(ord('A')+assignment[nodes.index(a)]))
        if a not in nodes:
            f.write('\n\t'+str(a+1)+'\t'+'-')
    f.close()
