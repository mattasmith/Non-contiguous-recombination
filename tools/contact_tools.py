

def find_alignment_contacts(structure_list,numbered_alignment):
    alignment_contacts = {}
    for i,pos1 in enumerate(numbered_alignment):
        for j,pos2 in enumerate(numbered_alignment):
            if i<j:
                for struc in range(len(pos1)):
                    if (pos1[struc],pos2[struc]) in structure_list[struc].contacts:
                        if alignment_contacts.has_key((i,j)):
                            alignment_contacts[(i,j)]+=1
                        else:
                            alignment_contacts[(i,j)] = 1
    return alignment_contacts


