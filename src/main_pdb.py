import Bio
import pickle
from Bio.PDB import PDBList
import sys
import os
package_dir = os.path.abspath('../packages')
sys.path.append(package_dir)
import json
import pprint as pp
import aadict as aa
import make_http_requests as mhr
import textwrap
import re
from collections import defaultdict

SEQ_DIR = "../data/pdb/"
DATA_DIR = "../data/"


pdb_ids = []
with open(SEQ_DIR+"pdb_ids.txt","r") as file_handle:
    pdb_ids = file_handle.read().split('\n')

pdb_ids = list(set(pdb_ids))
for pdb_id in pdb_ids:
    EB_URL = 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/'

    resp = mhr.make_request(EB_URL+pdb_id)

    if resp:
            uniprot_id = resp.text.split('{"')[3].replace('"','').replace(":",'')
            print(uniprot_id)
            annotation_dict = defaultdict(list)
            if uniprot_id:
                    URL = 'https://www.uniprot.org/uniprot/'

                    annotations = defaultdict(list)
                    uniprot_id = uniprot_id.strip()
                    print("processing "+uniprot_id)

                    response = mhr.make_request(URL+uniprot_id+".txt")

                    if response:
                        response = response.content.decode('utf-8').strip()
                        metals = response.split('FT   METAL           ')

                        annotations = metals[1:]
                        annotations = [int(x.split('\nFT')[0]) for x in annotations]
                        

                        if len(annotations) < 1:
                            print("no metal site")
                            with open("no_metal_site.txt", "a") as f:
                                    f.write(uniprot_id+"\n")
                    for annotation in annotations:
                                try:
                                        ann = response.split(
                                            '''FT   METAL           '''+str(annotation)+'''
FT                   /note="'''
                                        )[1].split(' ')[0].strip()
                                        ann = re.compile('[^a-zA-Z]').sub('', ann).replace('FT', "").replace('"','')
                                        if annotation not in annotation_dict[ann]:
                                            annotation_dict[ann].append(annotation)
                                    
                                except:
                                    with open(SEQ_DIR+"error.txt","a") as f:
                                        f.write(uniprot_id+"\n")
                    
                    if response:             
                        sequence = response.strip()

                        sequence = sequence[(sequence.rindex(';')+1) : (-2)].strip().replace(' ','')

                        sequences = ''
                        for key, value in annotation_dict.items():
                            print(key)
                            with open(SEQ_DIR+"sequences/"+key+".fasta", "a") as file_handle:
                                sequence = sequence.replace("\n",'')
                                sequences = aa.insert_annotations(sequence, value,'#')
                                sequences = textwrap.fill(sequences,60)
                                if len(sequences) > 0:
                                    file_handle.write('>'+uniprot_id +"\n")
                                    file_handle.write(sequences+"\n")

        
            