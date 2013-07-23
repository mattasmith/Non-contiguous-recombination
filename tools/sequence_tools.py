from subprocess import Popen
from random import choice
import os


def rand_tag():
    '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
    alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
    tag = ''.join([choice(alpha) for i in range(15)])
    return tag


# functions to run  muscle and read output
def read_MSA_file(filename):
    "clustalW, muscle, and procons output fasta files"
    out_fasta = open(filename,'r').read()
    sequences = [''.join(seq.split('\n')[1:]) for seq in out_fasta.split('>')[1:]]
    alignment = zip(*sequences)
    return alignment


def muscle_add_sequence(alignment,sequence,muscle_version):
    tag = rand_tag() # make a rand tag to make unique files
    new_seq = 'muscle_new_seq_'+tag+'.fasta'
    open(new_seq,'w').write('>new_seq\n'+sequence+'\n\n')
    existing_align = 'muscle_existing_align_'+tag+'.fasta'
    sequences = [''.join(s) for s in zip(*alignment)]
    open(existing_align,'w').write(''.join(['>seq'+str(i+1)+'\n'+sequences[i]+'\n\n' for i in range(len(sequences))]))
    muscle_filename = 'muscle_align_'+tag+'.out'
    muscle_command = ['./tools/muscle/' +muscle_version+ ' ',
                      '-profile ',
                      '-in1 '+existing_align,
                      '-in2 '+new_seq,
                      '-out '+muscle_filename] # output file 
    cmd = ' '.join(muscle_command)
    fnull = open(os.devnull, 'w')
    Popen(cmd,shell=True,stdout = fnull, stderr = fnull).wait() # runs muscle silently
    fnull.close()
    alignment = read_MSA_file(muscle_filename)
    os.remove(muscle_filename)
    os.remove(new_seq)
    os.remove(existing_align)
    return alignment


## fucntions to read, write, and print alignments
def read_alignment(filename):
    file = open(filename).read()
    data = [line for line in file.split('\n') if len(line) > 0 and line[0]!='#']
    if '>seq_names' in file:
        seq_names = data[data.index('>seq_names')+1:data.index('>alignment')]
    else:
        seq_names = []
    ali_data = data[data.index('>alignment')+1:]
    alignment =  [pos.split()[1:] for pos in ali_data]
    return alignment,seq_names


def write_alignment(alignment,seq_names=[],filename='new_alignment.aln'):
    num_digits = len(str(len(alignment)))
    alignment_file = open(filename,'w')
    if len(seq_names)>0:
        alignment_file.write('>seq_names\n')
        for name in seq_names:
            alignment_file.write(name+'\n')
    alignment_file.write('>alignment\n')
    for i,pos in enumerate(alignment):
        line = str(i).ljust(num_digits)+'  '+'  '.join(pos)
        alignment_file.write(line+'\n')
    alignment_file.close()


def print_alignment(alignment,seq_names=[]):
    if seq_names == []:
        seq_names=len(alignment[0])*['']
    name_length =  max([len(name) for name in seq_names])
    screen_width = 200
    num_lines = len(alignment)/screen_width    
    for i in range(num_lines+1):
        align_seg = alignment[(screen_width*i):(screen_width*(i+1))]
        conservation = []
        for pos in align_seg:
            if all([s==pos[0] for s in pos]):
                conservation.append('*')
            else:
                conservation.append(' ')
        seqs = zip(*align_seg)
        for i,seq in enumerate(seq_names):
            print seq.ljust(name_length)+':'+''.join(seqs[i])
        print ''.ljust(name_length)+':'+''.join(conservation)+'\n\n'


def number_alignment(alignment):
    '''This numbers an alignment so the alignment position can be mapped back to the original sequence (and its features) that were aligned'''
    numbered_sequences = []
    for sequence in zip(*alignment):
        num_seq = []
        k = 0
        for pos in sequence:
            if pos=='-':
                num_seq.append('-')
            else:
                num_seq.append(k)
                k+=1
        numbered_sequences.append(num_seq)
    numbered_alignment = zip(*numbered_sequences)
    return numbered_alignment


def read_fasta(filename):
    data = [s for s in open(filename).read().split('>') if len(s)>0]
    seq_names = [d.split('\n')[0] for d in data]
    sequences = [''.join([r for r in ''.join(d.split('\n')[1:]) if r.isalpha()==True]) for d in data]
    return seq_names,sequences


def read_fasta_alignment(filename):
    '''Read a fasta file of an alignment of parents in PROMALS3D'''
    data = [s for s in open(filename).read().split('>') if len(s)>0]
    par_names = [d.split('\n')[0] for d in data]
    sequences = [''.join(d.split('\n')[1:]) for d in data]
    parents = [''.join([r for r in ''.join(d.split('\n')[1:]) if r.isalpha()==True]) for d in data]
    parent_alignment = zip(*sequences)  # format used by muscle align
    return par_names,parents,parent_alignment

