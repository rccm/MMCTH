import warnings
warnings.filterwarnings("ignore", message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library.")

import sys
import os
import time
import logging
import sqlite3
import numpy as np
from dask import delayed, compute
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.utils import fetch_methods
from src.data_processors.process_main import MainProcessor
from src.data_processors.process_wMOD06_woTC import WMOD06woTC
from src.data_processors.process_woMOD06_wTC import WoMOD06wTC

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

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

@delayed
def fetch_and_process_group(inputfile_list):
    # Fetch the input files 
    mod21_file = inputfile_list[0]
    mod06_file = inputfile_list[1]
    mod03_file = inputfile_list[2]
    tccloud_file = inputfile_list[3]
    agp_file = inputfile_list[4]
    era5_single_file = inputfile_list[5]
    era5_multi_file = inputfile_list[6]

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
    return f"Completed processing group with MOD021KM File {mod21_file}"

def main():
    # Initialize Dask client
    client = Client()

    cluster = SLURMCluster(
        cores=16,                   # Number of cores per job
        memory="64GB",              # Memory per job
        processes=1,                # Number of processes per job
        walltime="01:00:00",        # Walltime per job
    )

    # Scale the cluster to use 4 jobs (nodes)
    cluster.scale(jobs=4)

    # Connect to the cluster
    client = Client(cluster)

    # Print the number of workers and cores
    num_workers = len(client.scheduler_info()['workers'])
    total_cores = sum(worker_info['nthreads'] for worker_info in client.scheduler_info()['workers'].values())
    print(f"Number of workers: {num_workers}")
    print(f"Total number of cores: {total_cores}")
    date = '2002-12-03'
    database_path = "/data/keeling/a/gzhao1/f/mmcth/lists/inputfiles_2002.sqlite"
    conn = create_connection(database_path)
    if conn is None:
        logger.error("Error! cannot create the database connection.")
        return

    # Fetch all group_ids to distribute
    groups = fetch_methods.fetch_files_by_date(conn, date)
    group_ids = list(groups.keys())
    
    print(f'<<<<There are {len(group_ids)} groups of input files>>>')

    # Create a list of delayed tasks
    tasks = [fetch_and_process_group(groups[group_id]) for group_id in group_ids]

    # Execute all tasks in parallel using Dask
    results = compute(*tasks)

    # Print the results
    for result in results:
        print(result)
    
    conn.close()
    client.close()

if __name__ == "__main__":
    main()
