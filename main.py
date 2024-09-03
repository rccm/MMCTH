
import warnings
warnings.filterwarnings("ignore", message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library.")

import sys
import os
import time
import logging
import sqlite3
import numpy as np
from mpi4py import MPI
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.utils import fetch_methods
from src.data_processors.process_main import MainProcessor
from src.data_processors.process_wMOD06_woTC import WMOD06woTC
from src.data_processors.process_woMOD06_wTC import WoMOD06wTC

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def create_connection(db_file):
    """ Create a database connection to the SQLite database specified by db_file """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except sqlite3.Error as e:
        print(e)
        return None

# Fetch the group files
def determine_process_method(mod06_file, tccloud_file, inputfile_list):
    conditions_to_methods = {
        (True, False): (WoMOD06wTC, inputfile_list),
        (False, True): (WMOD06woTC, inputfile_list),
        (False, False): (MainProcessor, inputfile_list),
        (True, True):  (None, None) 
    }
    default_method = (MainProcessor, inputfile_list)
    condition = (os.path.basename(mod06_file) == 'X', os.path.basename(tccloud_file) == 'X')
    process_class, args = conditions_to_methods.get(condition, default_method)

    # Print out which process method is being used
    print(f"Using process method: {process_class.__name__}")
    
    return process_class, args

def fetch_and_process_group(inputfile_list):
    # Fetch the input files 
    mod21_file = inputfile_list[0]
    mod06_file = inputfile_list[1]
    mod03_file = inputfile_list[2]
    tccloud_file = inputfile_list[3]
    agp_file = inputfile_list[4]
    era5_single_file = inputfile_list[5]
    era5_multi_file = inputfile_list[6]
    print(f'There are {len(inputfile_list)} input files')
    print('------------Input File List--------------')
    print(f'MOD021KM File >> {mod21_file}') 
    print(f'MOD06 File >>  {mod06_file}') 
    print(f'MOD03 File >> {mod03_file}') 
    print(f'TC_Cloud File >> {tccloud_file}') 
    print(f'AGP File >> {agp_file}') 
    print(f'ERA5 Single-Level File >> {era5_single_file}') 
    print(f'ERA5 Pressure-level File >> {era5_multi_file}') 
    print(f'----------------------------------------') 
    process_class, args = determine_process_method(mod06_file, tccloud_file, inputfile_list)
    if process_class is None:
        # Skip this set of input files
        print(f"Skipping this set of input files: {inputfile_list}")
        return 
    processor = process_class(args, logger)
    start_time = time.time()
    processor.run_process()
    end_time = time.time()
    duration = end_time - start_time
    print(f"Process {process_class.__name__} took {duration:.2f} seconds to run.")

def main():
    date = '2002-12-03'
    database_path = "/data/keeling/a/gzhao1/f/Database/inputfiles_2002.sqlite"
    conn = create_connection(database_path)
    if conn is None:
        logger.error("Error! cannot create the database connection.")
        return
    if rank == 0:
        # Fetch all group_ids to distribute
        
        groups = fetch_methods.fetch_files_by_date(conn, date)
        groups_ids = groups.keys()
    
        group_ids = list(groups_ids)
        print(f'<<<<There are {len(group_ids)} groups of input files>>>')
       # Distribute the jobs evenly through the cores
        chunks = [[] for _ in range(size)]
        for i, group_id in enumerate(group_ids):
            chunks[i % size].append(group_id)
    else:
        chunks = None
    group_ids_chunk = comm.scatter(chunks, root=0)
    # Each process fetches and processes its chunk of group_ids
    for group_id in group_ids_chunk[3:]:
        print(f'input file list for group_id: {group_id}')
        fetch_and_process_group(groups[group_id])
        break
    conn.close()

if __name__ == "__main__":
    main()

# retrieve all the variables 


