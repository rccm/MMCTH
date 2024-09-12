import argparse
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
from src.data_processors.process_wMOD06_woTC import WMOD06WoTC
from src.data_processors.process_woMOD06_wTC import WoMOD06WTC

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
        logger.error(e)
        return None

# Fetch the group files
def determine_process_method(mod06_file, tccloud_file, inputfile_list):
    conditions_to_methods = {
        (True, False): (WoMOD06WTC, inputfile_list),
        (False, True): (WMOD06WoTC, inputfile_list),
        (False, False): (MainProcessor, inputfile_list),
        (True, True):  (None, None) 
    }
    default_method = (MainProcessor, inputfile_list)
    condition = (os.path.basename(mod06_file) == 'X', os.path.basename(tccloud_file) == 'X')
    process_class, args = conditions_to_methods.get(condition, default_method)

    # Print out which process method is being used
    logger.debug(f"Using process method: {process_class.__name__}")
    
    return process_class, args

def fetch_and_process_group(inputfile_list):
    # Fetch the input files 
    mod21_file = inputfile_list[0]
    mod06_file = inputfile_list[1]
    mod03_file = inputfile_list[2]
    tccloud_file = inputfile_list[3]
    agp_file = inputfile_list[4]
    era5_single_file_current = inputfile_list[5]
    era5_multi_file_current = inputfile_list[6]
    
    logger.debug(f'There are {len(inputfile_list)} input files')
    logger.debug(f'MOD021KM File >> {mod21_file}') 
    logger.debug(f'MOD06 File >>  {mod06_file}') 
    logger.debug(f'MOD03 File >> {mod03_file}') 
    logger.debug(f'TC_Cloud File >> {tccloud_file}') 
    logger.debug(f'AGP File >> {agp_file}') 
    logger.debug(f'ERA5 Single-Level File >> {era5_single_file_current}') 
    logger.debug(f'ERA5 Pressure-level File >> {era5_multi_file_current}') 

    process_class, args = determine_process_method(mod06_file, tccloud_file, inputfile_list)
    if process_class is None:
        logger.debug(f"Skipping this set of input files: {inputfile_list}")
        return 
    processor = process_class(args, logger)
    start_time = time.time()
    processor.run_process()
    end_time = time.time()
    duration = end_time - start_time
    logger.debug(f"Process {process_class.__name__} took {duration:.2f} seconds to run.")

def main():
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description='Process satellite data for a given year, month, or date.')
    parser.add_argument('-y', '--year', type=int, required=True, help='Year to process (e.g., 2016)')
    parser.add_argument('-m', '--month', type=int, help='Month to process (optional, e.g., 3 for March)')
    parser.add_argument('-d', '--date', type=str, help='Day to process in MM-DD format (optional)')
    parser.add_argument('-o', '--orbit', type=str, help='Day to process in MM-DD format (optional)')
    
    args = parser.parse_args()
    year = args.year
    month = args.month
    date = args.date
    orbit = args.orbit

    # Construct database path
    db_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_{year}.sqlite"
    conn = create_connection(db_path)
    if conn is None:
        logger.error("Cannot connect to the database.")
        return
    
    # Fetch group data based on input parameters
    if date:
        # If --date is provided
        month, day = map(int, date.split('-'))
        full_date = f"{year}-{month:02d}-{day:02d}"
        logger.info(f'Fetching files for {full_date}')
        groups = fetch_methods.fetch_files_by_date(conn, full_date)
    elif month:
        # If only --month is provided
        logger.info(f'Fetching files for {year}-{month:02d}')
        groups = fetch_methods.fetch_files_by_month(conn, year, month)
    elif orbit:
        logger.info(f'Fetching files for {year}-{month:02d}')
        groups = fetch_methods.fetch_files_by_orbit(conn, orbit)
    else:
        # Default: process entire year
        logger.info(f'Fetching files for the entire year {year}')
        groups = fetch_methods.fetch_files_by_year(conn, year)

    group_ids = list(groups.keys())
    logger.info(f'Found {len(group_ids)} groups.')

    if rank == 0:
    # Split group_ids into roughly equal chunks
        group_ids_chunk = np.array_split(group_ids, size)
    else:
        group_ids_chunk = None

    # Scatter job chunks to processes
    group_ids_chunk = comm.scatter(group_ids_chunk, root=0)

    # Each process handles its assigned group IDs
    for group_id in group_ids_chunk:
        logger.debug(f'Processing group ID: {group_id}')
        fetch_and_process_group(groups[group_id])
        break 
    conn.close()

if __name__ == "__main__":
    main()
