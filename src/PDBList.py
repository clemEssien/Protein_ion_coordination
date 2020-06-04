import Bio
import pickle
from Bio.PDB import PDBList
import sys
import os
package_dir = os.path.abspath('../packages')
sys.path.append(package_dir)

import aadict as aa
import textwrap
from collections import defaultdict

SEQ_DIR = "../data/sequences/"
pdblist = PDBList()
# entry_list = pdblist.get_all_entries()
## pdblist.download_entire_pdb(file_format='pdb')

# with open('../data/pdb_entries.txt', 'w') as filehandle:
#     filehandle.write("\n".join(entry_list))
# for pdb_id in entry_list:   
#     pdb_id = pdb_id.strip()
#     pdblist.retrieve_pdb_file(pdb_id,file_format='pdb',pdir='../data/PDB', overwrite=False)
#     with open('../data/PDB/pdb/'+pdb_id+'.ent') as filehandle:
#         for line in filehandle:




count = 0
res_list = []
amino_acid_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
metals_list = ['CS', 'K.', 'LI', 'NA', 'RB', 'BA', 'BE', 'CA', 'MG', 'SR', 'SC', 'Y.',
       'HF', 'LA', 'TI', 'ZR', 'CE', 'TA', 'V.', 'CR', 'MO', 'PR', 'W.', 'MN',
       'RE', 'U.', 'FE', 'OS', 'RU', 'CO', 'IR', 'RH', 'SM', 'EU', 'NI', 'PD',
       'PT', 'AG', 'AU', 'CU', 'GD', 'CD', 'HG', 'TB', 'ZN', 'AL', 'DY', 'GA',
       'IN', 'TL', 'HO', 'PB', 'SN', 'AS', 'BI', 'ER', 'SB', 'LU', 'YB']

pdb_id = '830c copy'
annotations = defaultdict(list)

# pdblist.retrieve_pdb_file(pdb_id,file_format='pdb',pdir='../data/PDB', overwrite=False)#download pdb file
with open('../data/PDB/pdb'+pdb_id+'.ent') as filehandle:
         metals = []
         res = {}
         site = []
         records = filehandle.readlines()
         chain_info = ""
         sequences = ""
         
         for i in range(1,len(records)):
            index = 0
            if ('REMARK 620 COORDINATION ANGLES FOR:  M RES CSSEQI METAL') in (records[i]) :
                metal = records[i+1].split()[-1]
                metals.append(records[i+1].split()[-1])
                
                j=i+1

                try:
                    while True: 
                        if records[j].strip() == 'REMARK 620':
                            break
                        else:
                            if records[j].strip() == 'REMARK 800':
                                break
                            #print(records[j])
                            metal = records[j].split() 
                            if len(metal) == 6:
                                if metal[5] in metals_list : 
                                    #print(metal[5]+'_'+metal[4]+'_'+metal[3]+'------')
                                    index +=1
                                    if metal[5]+'_'+metal[3] not in site:
                                        site.append(metal[5]+'_'+metal[3])
                                        #print(metal[5]+'_'+metal[4]+'_'+metal[3]+'------')


                            else: 
                                if (records[j].strip() != 'REMARK 620 N RES CSSEQI ATOM ' and  
                                    'REMARK 620 N' not in records[j].strip()) :
                                        residues = records[j].split()
                                        features = records[i+1].split()
                                        # print(features[2]+'_'+features[3]+'_'+features[4]+'.........')
                                        
                                        if aa.aa_residue(residues[3]) != None:
                                            residues = aa.aa_residue(residues[3])+'_'+residues[4]+'_'+residues[5]
                                            # print(residues)
                                            annotations[features[2]+'_'+features[3]+'_'+features[4]].append(residues)
                                        # with open("one.txt", "a") as f:
                                        
                                        #     f.write(records[j] +'-'+residues) 
                        j+=1
                except IndexError:
                    print("nothing")

            if ('SEQRES') in records[i]:
                sequences += records[i]

aa_acid = aa.sort_aa_dict(annotations)
seq_rec = sequences.split('\n')   

for keys, values in aa_acid.items():
    chain = keys.split('_')
    print(pdb_id+'_'+chain[1]+':')
    seq = ""
    for line in seq_rec:
        seq_info = line.split()
        if len(seq_info) > 0:
            if seq_info[2] == chain[1].strip():
                seq = seq+(aa.aa_residue_list(' '.join(seq_info[4:]).strip()))
    if len(seq) > 0:
        seq = textwrap.fill(seq,60)
        position = [int(x.split('_')[2]) for x in values]
        position.sort()
        if not os.path.exists(SEQ_DIR):
            os.mkdir(SEQ_DIR)
        with open("../data/sequences/"+chain[0]+'.fasta',"a") as filehandle:
            seq = aa.insert_annotations(seq, position,'#')
            filehandle.write(pdb_id+'_'+chain[0]+'_'+chain[1]+':'+'\n'+seq+'\n')
        
        print((values)) 
        
