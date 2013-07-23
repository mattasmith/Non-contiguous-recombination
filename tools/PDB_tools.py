import urllib
import urllib2
import os

threetoone={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
            'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
onetothree = dict([(threetoone[x],x) for x in threetoone.keys()])


def quick_sequence(pdbfile):
    return ''.join([threetoone[l[17:20]] for l in open(pdbfile).read().split('\n') if l[13:15]=='CA'])


def pdb_search(sequence):
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = """
    <orgPdbQuery>
        <queryType>org.pdb.query.simple.SequenceQuery</queryType>
        <description>Sequence Search (Structure:Chain = 4HHB:A, Expectation Value = 10.0, Search Tool = blast)</description>
        <structureId></structureId>
        <chainId></chainId>
        <sequence>%s</sequence>
        <eCutOff>1e-10</eCutOff>
        <searchTool>blast</searchTool>
    </orgPdbQuery>""" % sequence

    req = urllib2.Request(url, data=queryText)
    f = urllib2.urlopen(req)
    result = f.read()
    PDBs = tuple(sorted(set([p for p in result.replace('\r','').split('\n') if len(p)>0])))
    return PDBs


def download_parse_save(pdb,path=''):
    pdb_file = urllib.urlopen('http://www.rcsb.org/pdb/files/'+pdb+'.pdb').read().split('\n')
    atoms = [l.strip() for l in pdb_file if l[:4]=='ATOM']
    nonhydrogen = [l for l in atoms if l[-1]!='H']
    nonaltconf = [l for l in nonhydrogen if l[16:17]=='A' or l[16:17]==' ']
    files = []
    for chain in set([l[21:22] for l in nonaltconf]):
        print 'writing '+path+pdb+'.'+chain+'.pdb'
        open(path+pdb+'.'+chain+'.pdb','w').write('\n'.join([l for l in nonaltconf if l[21:22]==chain]))
        files.append(path+pdb+'.'+chain+'.pdb')
    return files

def parse_save(pdb,path=''):
    pdb_file = open(path+pdb+'.pdb','r').read().split('\n')
    atoms = [l.strip() for l in pdb_file if l[:4]=='ATOM']
    nonhydrogen = [l for l in atoms if l[-1]!='H']
    nonaltconf = [l for l in nonhydrogen if l[16:17]=='A' or l[16:17]==' ']
    files = []
    for chain in set([l[21:22] for l in nonaltconf]):
        print 'writing '+path+pdb+'.'+chain+'.pdb'
        open(path+pdb+'.'+chain+'.pdb','w').write('\n'.join([l for l in nonaltconf if l[21:22]==chain]))
        files.append(path+pdb+'.'+chain+'.pdb')
    return files

def split_into_chains(filename,path=''):
    pdb_fullname = os.path.splitext(filename)[0]
    pdb_firstname = os.path.splitext(pdb_fullname)[0]
    ## if the name of the pdb does not have any '.'s in it - need to split into chains
    if pdb_fullname == pdb_firstname:
        ## split into chains
        files = parse_save(pdb_fullname,'./structures/')
        os.remove('./structures/'+filename)
        return files
    else:
        pass

## HOW TO USE THE pdb structure class
#structure = PDB_structure()
#structure.read('1jpz.pdb')
#structure.filter('include',[['ATOM'],'all','all',['','A'],'all',['A'],'all']) #[record_types, atom_number, atom_type, alt_loc, residue_type, chain, residue_number]
#structure.filter('exclude',['none','none',['OXT'],'none','none','none','none']) #[record_types, atom_number, atom_type, alt_loc, residue_type, chain, residue_number]
class PDB_structure:
    '''represents a PDB as a list of lines'''
    def __init__(self,filename=''):
        if len(filename)>0: self.read(filename)

    def initialize(self):
        '''run this method after updating the structure in any way'''
        self.parse_sequence() # adds sequence and numbered sequence
        self.split_structure() # adds split residues
        self.get_backbone()
        
    def parse_sequence(self):
        seq = []
        for line in self.structure:
            if line[0:4]=='ATOM' and line[12:16].strip()=='CA':
                seq.append(threetoone[line[17:20]]+line[22:26].strip())
        self.sequence = ''.join([p[0] for p in seq])
        self.num_seq = seq

    def split_structure(self):
        residue_numbers = [str(i) for i in list(set([int(line[22:26].strip()) for line in self.structure if line[0:4]=='ATOM']))]
        protein = []
        for res in residue_numbers:
            residue = []
            for line in self.structure:
                residue_number = line[22:26].strip()
                if residue_number==res:
                    residue.append(line)
            protein.append(residue)
        self.residues = protein

    def get_backbone(self):
        self.backbone = [tuple([float(c) for c in [l[30:38],l[39:46],l[47:54]] ]) for l in self.structure if (l[13:15]=='CA' and l[0:4]=='ATOM')]

    def get_contacts(self,exclude = []):
        coordinates = []
        for residue in self.residues:
            coor = [[float(c) for c in [l[30:38],l[39:46],l[47:54]] ] for l in residue if (l[12:16].strip() not in exclude and l[0:4]=='ATOM')]
            coordinates.append(coor)
        contacts = []
        for i,res1 in enumerate(coordinates):
            for j,res2 in enumerate(coordinates):
                if i<j:
                    coor_pairs = [(c1,c2) for c1 in res1 for c2 in res2]
                    for cp in coor_pairs:
                        coor1,coor2 = cp
                        distance = sum([(coor1[k]-coor2[k])**2 for k in range(3)])**0.5
                        if distance < 4.5:
                            contacts.append((i,j))
                            break
        self.contacts = tuple(contacts)


    # methods for reading/writing files and PDB strings
    def read(self,filename):
        pdb_file = open(filename,'r').read().split('\n')
        self.name = filename
        self.structure = pdb_file
        self.initialize()
        
    def read_str(self,pdb_str):
        pdb_file = pdb_str.split('\n')
        self.name = 'from string'
        self.structure = pdb_file
        self.initialize()
        
    def write(self,filename):
        out_str = '\n'.join(self.structure)
        open(filename,'w').write(out_str)


    ## methods for updating PDBs:
    def remove_residues(self,residue_list):
        for res in residue_list:
            self.residues[res]=[]

        ## update everything
        new_pdb = []
        for res in self.residues:
            for line in res:
                new_pdb.append(line)
        self.structure = new_pdb
        self.initialize()


    def renumber(self):
        pdb = self.structure
        #first, renumber the atoms
        renum_atoms = []
        k = 1
        for line in pdb:
            renum_atoms.append(line[:6]+str(k).rjust(5)+line[11:])
            k = k+1
        #seconds, renumber the residues
        renum_res = []
        k = 0
        res_type_prev = ''
        for line in renum_atoms:
            if line[12:16].strip()=='N' or res_type_prev != line[17:20]: #if there is a nitrogen, or the residue name changes
                k = k+1
            renum_res.append(line[:22]+str(k).rjust(4)+line[26:])
            res_type_prev = line[17:20]
        self.structure = renum_res
        self.initialize()


    def filter(self,filter_type,records_list): # [record_types, atom_number, atom_type, alt_loc, residue_type, chain, residue_number]
        '''a very gernal filtering method for PDBs'''
        
        pdb = [line for line in self.structure if len(line)>30] #prefilter

        # record type - cols 0:6
        if records_list[0]=='all':
            record_types = list(set([line[0:6].strip() for line in pdb]))
        elif records_list[0]=='none':
            record_types = []
        else:
            record_types = records_list[0]

        # atom number - cols 6:11
        if records_list[1]=='all':
            atom_numbers = list(set([line[6:11].strip() for line in pdb]))
        elif records_list[1]=='none':
            atom_numbers = []
        else:
            atom_numbers = records_list[1]

        # atom type - cols 12:16
        if records_list[2]=='all':
            atom_types = list(set([line[12:16].strip() for line in pdb]))
        elif records_list[2]=='none':
            atom_types = []
        else:
            atom_types = records_list[2]

        # alt_loc - cols 16
        if records_list[3]=='all':
            alt_locs = list(set([line[16].strip() for line in pdb]))
        elif records_list[3]=='none':
            alt_locs = []
        else:
            alt_locs = records_list[3]

        # residue type - cols 17:20
        if records_list[4]=='all':
            residue_types = list(set([line[17:20].strip() for line in pdb]))
        elif records_list[4]=='none':
            residue_types = []
        else:
            residue_types = records_list[4]

        # chain - cols 21
        if records_list[5]=='all':
            chains = list(set([line[21].strip() for line in pdb]))
        elif records_list[5]=='none':
            chains = []
        else:
            chains = records_list[5]

        # residue number - cols 22:26
        if records_list[6]=='all':
            residue_numbers = list(set([line[22:26].strip() for line in self.structure]))
        elif records_list[6]=='none':
            residue_numbers = []
        else:
            residue_numbers = records_list[6]

        if filter_type=='include':
            new_pdb = []
            for line in pdb:
                record_type = line[0:6].strip()
                atom_number = line[6:11].strip()
                atom_type = line[12:16].strip()
                alt_loc = line[16].strip()
                residue_type = line[17:20].strip()
                chain = line[21].strip()
                residue_number = line[22:26].strip()
                if record_type in record_types and atom_number in atom_numbers and atom_type in atom_types and alt_loc in alt_locs and residue_type in residue_types and chain in chains and residue_number in residue_numbers:
                    new_pdb.append(line)

        if filter_type=='exclude':
            new_pdb = []
            for line in pdb:
                record_type = line[0:6].strip()
                atom_number = line[6:11].strip()
                atom_type = line[12:16].strip()
                alt_loc = line[16].strip()
                residue_type = line[17:20].strip()
                chain = line[21].strip()
                residue_number = line[22:26].strip()
                if record_type not in record_types and atom_number not in atom_numbers and atom_type not in atom_types and alt_loc not in alt_locs and residue_type not in residue_types and chain not in chains and residue_number not in residue_numbers:
                    new_pdb.append(line)

        self.structure = new_pdb
        self.initialize()



