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
records = []
#catch all exceptions
try:
    with open(SEQ_DIR+"uniprot_data1.txt") as filehandle:
        records = filehandle.read().split("//")
except expression as identifier:
    print("error reading file")

if len(records) > 0:
    for record in records:
        
        if record: 
            if 'FT   METAL           ' in record:
                annotation_text = record.split('FT   METAL           ')
                annotation_dict = defaultdict(list)
                
                for i in range(1,len(annotation_text)):
                    annotation = int(annotation_text[i].split('FT   METAL           ')[0].split("\n")[0])
                    metal = annotation_text[i].split(str(annotation)+'''
FT                   /note="''')[1].split('"')[0].split(';')[0]
                    annotation_dict[metal].append(annotation)
                print((annotation_dict))

                sequence = annotation_text[i].strip()
                sequence = sequence[(sequence.rindex(';')+1) : (-2)].strip().replace(' ','')

                sequences = ''
                for key, value in annotation_dict.items():
                    # print(key)
                    # with open(SEQ_DIR+"sequences/"+key+".fasta", "a") as file_handle:
                    sequence = sequence.replace("\n",'')
                    sequences = aa.insert_annotations(sequence, value,'#')
                    sequences = textwrap.fill(sequences,60)
                    if len(sequences) > 0:
                        #print('>'+uniprot_id +"\n")
                        print(sequences+"\n")

