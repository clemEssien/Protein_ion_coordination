import urllib.parse
import urllib.request
import sys
import os
package_dir = os.path.abspath('../../packages')
sys.path.append(package_dir)
import make_http_requests as mhr

pdb_ids = []
DATA_DIR = "../../data/"
with open(DATA_DIR +"pdb_entries.txt") as filehandle:
       p_ids = filehandle.readlines()
       for pdb_id in p_ids:
              pdb_ids.append(pdb_id.strip())

for i in range(58,len(pdb_id)+58,58):
       # print(p_ids[(i-1000):(i)])
       arr = p_ids[(i-58):(i)]
       url = 'https://www.uniprot.org/uploadlists/'
       # pdb_ids = ['5k21','830c']
       # print(' '.join(pdb_ids))
       params = {
       'from': 'PDB_ID',
       'to': 'ACC',
       'format': 'tab',
       'query': ' '.join(arr)
       }

       data = mhr.make_request(url, params)
       print(data.content.decode('utf-8'))


       result = data.content.decode('utf-8')

       with open(DATA_DIR +"mapping.txt", "a") as file_handle:
              file_handle.write(result)



'''
MTGKTKPAII GGVVIAALAA AGLGVWLFTD GRGGRSTTEP VTMTLDVKND 
        60         70         80         90        100
QVAKHDFGKP GMDVGDMDIF SDILSVDGKQ VGYDGGACFF TNVTPDNPMT 
       110        120        130        140        150
YCELTIHLDA GEIFARSLTP HTLAPFTMAI TGGTGEYANS KGELTVSGVA 
       160 
TPDEKYELKL TK     
'''

'''
MTGKTKPAII GGVVIAALAA AGLGVWLFTD GRGGRSTTEP VTMTLDVKND
QVAKHDFGKP GMDVGDMDIF SDILSVDGKQ VGYDGGACFF TNVTPDNPMT
YCELTIHLDA GEIFARSLTP HTLAPFTMAI TGGTGEYANS KGELTVSGVA
TPDEKYELKL TK
'''
#https://www.uniprot.org/uniprot/P32179.txt
