# OS: Ubuntu, 18.04.1 LTS
# Python: Python 2.7.15
# Mongodb: v3.2.21 
# Siteng Cai
# import feedparser
from datetime import datetime as dt
import sys
if sys.version_info[0] >= 3:
    from urllib.request import urlretrieve
else:
    # Not Python 3 - today, it is most likely to be Python 2
    # But note that this might need an update when Python 4
    # might be around one day
    from urllib import urlretrieve
    
import gzip
import shutil
# import pymongo
# from pymongo import MongoClient
import sys
import os.path
import argparse
# from crontab import CronTab
# import configparser
import re
import json
import itertools
import os
import getpass
import time
#GLOBAL VAR
USER_NAME = getpass.getuser()
PARENT_DIR = os.path.dirname(os.getcwd())


def getUniprot():
    creattime = '-'.join(time.ctime(os.path.getctime('../data/uniprot/uniprot.txt')).split())
    print("zip old file created on "+creattime)
    with gzip.open('../data/uniprot/uniprot'+str(creattime)+'.txt.gz', 'wb') as f_out:
        with open('../data/uniprot/uniprot.txt', 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)
    
    print("download new file")
    uniprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
    urlretrieve(uniprot_url, 'uniprot.txt.gz')
    print("Unzip file...")
    with gzip.open('uniprot.txt.gz', 'rb') as f_in:
        with open('../data/uniprot/uniprot.txt', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("File name:uniprot.txt")
getUniprot()