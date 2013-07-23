from sys import path
import subprocess
import os
path.append('./tools')
import chimera_tools
import pickle
from random import random

def run_shmetis(num_contact_maps, numberofblocks_min, numberofblocks_max, hmetis_version):
    ## create folder for solutions
    output_path = './output/'
    if os.path.isdir(output_path) == False:
        os.mkdir(output_path)
    
    ## load in an alignment and structural contacts
    alignment,contacts = pickle.load(open('alignment_and_contacts.pkl'))
    numberofparents = len(alignment[(1)])
        
    ## calculate the SCHEMA E for each contact
    SCHEMA_E = chimera_tools.calculate_contact_SCHEMA_E_2(alignment,contacts)
    
    ## get the alignment indices of all residues that have at least one SCHEMA contact
    ## these are residues that have a SCHEMA contact
    nodes = sorted(set([c[0] for c in SCHEMA_E.keys()]+[c[1] for c in SCHEMA_E.keys()]))
    
    ## construct and write the hmetis HGraphFile
    # first line is the number of edges, nodes, and 1 (which means weighted edges)
    graph_str = '%d %d %s\n' % (len(SCHEMA_E),len(nodes),'1')
    for residue in nodes:
        for contact in [c for c in SCHEMA_E.iterkeys() if residue in c]:
            pos1 = nodes.index(residue)+1
            # each res that is contacting residue (note: indices start at 1)
            pos2 = nodes.index([c for c in contact if c!=residue][0])+1
            if pos2 > pos1:
                # format is: weight1 node1_1 node1_2; weight2 node2_1 node2_2; weight3 node3_1 node3_2; ...
                graph_str += '%i %i %i\n' % (SCHEMA_E[contact], pos1, pos2) ##
    
    ## Create a HGraphFile from the contact matrix
    open('SCHEMA.hgr','w').write(graph_str)
    HGraphFile = 'SCHEMA.hgr'

    
    ## try a range of different number of blocks
    for nblocks in range(numberofblocks_min,numberofblocks_max+1):
        
        E = []
        block_count = []
        hmetis_curve_m = []
        hmetis_curve_E = []
        hmetis_curve_FF = []
        outputfilelist = []
        
        Numberofblocks = str(nblocks)
                            
        for UB in range(30):
            UBfactor = str(UB+1)
            
            for repeat in range(1):
                
                E = []
                output = subprocess.Popen('./tools/' +hmetis_version+ '/shmetis '+HGraphFile+' '+Numberofblocks
                                      +' '+UBfactor,shell=True,stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()
                
                ## import the assignments from the output file
                assignment = [int(l) for l in open('SCHEMA.hgr.part.'+Numberofblocks).read().split('\n') if len(l)>0]
                E = (chimera_tools.score_library(assignment, SCHEMA_E, nodes, numberofparents))/num_contact_maps
                blocks = sorted(set(assignment))
                block_count = ([assignment.count(b) for b in blocks])
                av_m = chimera_tools.calculate_average_m(assignment,alignment,nodes)  ## calculate <m>
                out = tuple( [Numberofblocks, UBfactor, E, av_m ] + block_count )
                print ('  library%s_%s: E = %.3f, m = %.3f, blocks = ['+'%i '*len(block_count)+']'  ) % out
                hmetis_curve_E.append(E)
                hmetis_curve_m.append(av_m)
                ## save the assignments in an output file
                outputfile ='library'+Numberofblocks+'_'+UBfactor+'.output' # name of the output file
                outputfilelist.append(outputfile);
                #open('output/'+outputfile,'w').write('\n'.join([str(s) for s in assignment]))
                chimera_tools.write_solution(outputfile, alignment, assignment, nodes)
                
            # end for repeat
        # end for UB
        f = open('library'+Numberofblocks+'_result_list.csv','w')
        for i in range(len(hmetis_curve_E)):
            f.write(outputfilelist[i]+', '+str(hmetis_curve_E[i])+', '+str(hmetis_curve_m[i])+',\n')
        f.close()
        
        os.remove('SCHEMA.hgr.part.'+Numberofblocks)
    # end for nblocks
    os.remove('SCHEMA.hgr')
    
    ## give warning if residues have no SCHEMA contacts!!
    for res,align in enumerate(alignment):
        ## if residue is not conserved, but not in nodes
        if (len(set(align)) != 1) and (res not in nodes):
            print ('Warning: residue ' + str(res+1) + ' is not conserved but has no SCHEMA contacts, so it has not been assigned to a block.')
    # end for res, align
    return 1

