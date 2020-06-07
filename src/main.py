import Bio
import pickle
from Bio.PDB import PDBList
import sys
import os
package_dir = os.path.abspath('../packages')
sys.path.append(package_dir)

import aadict as aa
import make_http_requests as mhr
import textwrap
from collections import defaultdict

SEQ_DIR = "../data/uniprot/"
DATA_DIR = "../data/"
pdblist = PDBList()


pdb_id = '2ga6'
annotations = defaultdict(list)
URL = 'https://www.uniprot.org/uniprot/'
pdblist.retrieve_pdb_file(pdb_id,file_format='pdb',pdir='../data/PDB', overwrite=False)#download pdb file
with open(DATA_DIR+'PDB/pdb'+pdb_id+'.ent') as filehandle:       
        records = filehandle.readlines()
        uniprot_id = ""
        for line in records:
            if "UNP" in line:
                unp_index = line.split()
                index = unp_index.index('UNP')+1
                uniprot_id = unp_index[index]
                print(uniprot_id)
                if uniprot_id != None:
                    with open(DATA_DIR+"mapping.txt","a") as f:
                        f.write(pdb_id+' '+uniprot_id+'\n')
                    break
                with open(DATA_DIR+"no_uniprot_id.txt","a") as fh:
                    fh.write(pdb_id+'\n')
        response = mhr.make_request(URL+uniprot_id+".txt")
        response = response.content.decode('utf-8').strip()
     
        annotation_dict = defaultdict(list)
        metals = response.split('FT   METAL           ')

        annotations = metals[1:]
        annotations = [int(x.split('\nFT')[0]) for x in annotations]
        for annotation in annotations:
            try:
                    ann = response.split('''FT   METAL           '''+str(annotation)+'''
FT                   /note="''')[1].split(' ')[0].strip()
                    annotation_dict[ann].append(annotation)
                   
            except:
                a =1
        sequence = response.strip()
        
        sequence = sequence[(sequence.rindex(';')+1) : (-2)].strip().replace(' ','')
        
        for key, value in annotation_dict.items():
            with open(SEQ_DIR+key+".fasta", "a") as file_handle:
                sequences = aa.insert_annotations(sequence, value,'#')
                file_handle.write('>'+pdb_id +' : '+uniprot_id+"\n")
                file_handle.write(sequences)
                print(sequences)
      
