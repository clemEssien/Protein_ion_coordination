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
        records = filehandle.read().split("//\n")
except expression as identifier:
    print("error reading file")

if len(records) > 0:
    for record in records:
        
        if record: 
            if 'FT   METAL           ' in record:
                annotation_text = record.split('FT   METAL           ')
                annotation_dict = defaultdict(list)
                uniprot_id = re.split('\nAC[\s]+',record)[1].split('\n')[0]
                print('processing..... '+uniprot_id)
                for i in range(1,len(annotation_text)):
                    annotation = int(re.split('FT[\s]METAL[\s]',annotation_text[i])[0].split("\n")[0])
                    metal = re.split('[\d]+\nFT[\s]+/note="',annotation_text[i])
                    if metal:
                        metal = re.split('[^a-zA-Z]+',metal[1])[0]
                        annotation_dict[metal].append(annotation)
                    
                if annotation_dict:
                    sequence = annotation_text[i].strip()
                    sequence = re.split(';\n',sequence)[1].split('\n//')[0].replace(' ','')
                    sequences = ''
                    # for key, value in annotation_dict.items():
                    #     sequence = sequence.replace("\n",'')
                    #     sequences = aa.insert_annotations(sequence, value,'#')
                    #     sequences = textwrap.fill(sequences,60)
                    #     if len(sequences) > 0:
                    #         print('>'+uniprot_id +"\n")
                    #         print(sequences+"\n")
                    for key, value in annotation_dict.items():
                        with open(SEQ_DIR+"new_current_sequences/"+key+".fasta", "a") as file_handle:
                            sequence = sequence.replace("\n",'')
                            sequences = aa.insert_annotations(sequence, value,'#')
                            sequences = textwrap.fill(sequences,60)
                            if len(sequences) > 0:
                                file_handle.write('>'+uniprot_id +"\n")
                                file_handle.write(sequences+"\n")
