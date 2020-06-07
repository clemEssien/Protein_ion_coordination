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
import re
from collections import defaultdict

SEQ_DIR = "../data/uniprot/"
DATA_DIR = "../data/"



annotations = defaultdict(list)
URL = 'https://www.uniprot.org/uniprot/'

uniprot_id = 'P45452'
response = mhr.make_request(URL+uniprot_id+".txt")
response = response.content.decode('utf-8').strip()
# print(response)

#Zinc, Calcium, Mangnesium, Calcium, Sodium, Iron, Manganese, 
annotation_dict = defaultdict(list)
metals = response.split('FT   METAL           ')


annotations = metals[1:]
annotations = [int(x.split('\nFT')[0]) for x in annotations]
print(annotations)
if len(annotations) < 1:
    print("no metal site")
    with open("no_metal_site.txt", "a") as f:
        f.write(uniprot_id+"\n")

for annotation in annotations:
            try:
                    ann = response.split('''FT   METAL           '''+str(annotation)+'''
FT                   /note="''')[1].split(' ')[0].strip()
                    ann = re.compile('[^a-zA-Z]').sub('', ann).replace('FT', "")
                    if annotation not in annotation_dict[ann]:
                        annotation_dict[ann].append(annotation)
                   
            except:
                a =1
sequence = response.strip()

sequence = sequence[(sequence.rindex(';')+1) : (-2)].strip().replace(' ','')


for key, value in annotation_dict.items():
    print(key)
    with open(SEQ_DIR+key+".fasta", "a") as file_handle:
        sequence = sequence.replace("\n",'')
        sequences = aa.insert_annotations(sequence, value,'#')
        sequences = textwrap.fill(sequences,60)
        file_handle.write('>'+uniprot_id +"\n")
        file_handle.write(sequences+"\n")

