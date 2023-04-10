#!/usr/bin/env python3
import os
import sys
import argparse
import urllib.request
import zipfile
from tqdm import tqdm
PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_db_download.py",
    description="This script will download databases for WhatsGNU for WhatsGNU v1.4",)
PARSER.add_argument("db_name", type=str, help="database name[Sau, TB, Pa, Kp, Staphopia, S.enterica, all] for WhatsGNU v1.4")
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
################################
class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)
#################################
Script_Path = os.path.realpath(__file__)
DB_Folder_Path = os.path.join(Script_Path.rsplit(OS_SEPARATOR,2)[0], "db")
if os.path.exists(DB_Folder_Path):
    print("Found database Folder in {}".format(DB_Folder_Path))
else:
    os.mkdir(DB_Folder_Path)
DB_Dict = {}
DB_file = os.path.join(Script_Path.rsplit(OS_SEPARATOR,2)[0], "databases_available.csv")
DB_file_obj = open(DB_file, 'r')
for line in DB_file_obj:
    line_lst = line.rstrip().split(',')
    DB_Dict[line_lst[0]] = line_lst[1:]

if ARGS.db_name == 'all':
    for db in DB_Dict:
        DB_list = DB_Dict[db]
        DB_FILE_Name = DB_list[0]
        COMP_DB_file = os.path.join(DB_Folder_Path,'{}.zip'.format(db))#compressed DB
        if os.path.exists(DB_FILE_Name):
            print("Found database {}".format(DB_FILE_Name))
        else:
            url = DB_list[1]
            download_url(url, COMP_DB_file)
            #urllib.request.urlretrieve(url, COMP_DB_file)
            os.mkdir(DB_FILE_Name)
            DB_DIRECTORY = os.path.join(DB_Folder_Path,DB_FILE_Name)#compressed DB
            with zipfile.ZipFile(COMP_DB_file, 'r') as zip_ref:
                zip_ref.extractall(DB_DIRECTORY)
            os.remove(COMP_DB_file)
            print("Downloaded and decompressed database {}".format(DB_FILE_Name))
else:
    DB_list = DB_Dict[ARGS.db_name]
    DB_FILE_Name = DB_list[0]
    COMP_DB_file = os.path.join(DB_Folder_Path,'{}.zip'.format(ARGS.db_name))#compressed DB
    if os.path.exists(DB_FILE_Name):
        print("Found database {}".format(DB_FILE_Name))
    else:
        url = DB_list[1]
        download_url(url, COMP_DB_file)
        #urllib.request.urlretrieve(url, COMP_DB_file)
        os.mkdir(DB_FILE_Name)
        DB_DIRECTORY = os.path.join(DB_Folder_Path,DB_FILE_Name)#compressed DB
        with zipfile.ZipFile(COMP_DB_file, 'r') as zip_ref:
            zip_ref.extractall(DB_DIRECTORY)
        os.remove(COMP_DB_file)
        print("Downloaded and decompressed database {}".format(DB_FILE_Name))
