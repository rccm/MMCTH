import os
import fnmatch
import sqlite3
import sys
from datetime import datetime
from find_groupfiles import find_nite_files
sys.path.append(os.path.abspath('../'))
from data_readers.mod_read import MODL1Granule
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def process_files(file_list):
    night_files = []
    for file_path in file_list:
        result = find_nite_files(file_path)
        if result:
            night_files.append(result)
    return night_files

if __name__ == '__main__':
    # Read the file list
    if rank == 0:
        with open('/data/gdi/f/gzhao1/Database/MOD021KM_completelist.txt', 'r') as f:
            file_list = f.read().splitlines()
        # Split the file list into chunks for each process
        chunk_size = len(file_list) // size
        chunks = [file_list[i * chunk_size:(i + 1) * chunk_size] for i in range(size)]
        # Distribute the remainder if the file list isn't evenly divisible
        if len(file_list) % size != 0:
            remainder = file_list[size * chunk_size:]
            chunks[-1].extend(remainder)
    else:
        chunks = None
    # Scatter the chunks to all processes
    chunk = comm.scatter(chunks, root=0)

    # Each process processes its chunk
    processed_chunk = process_files(chunk)

    # Gather all the processed chunks back to the root process
    gathered_data = comm.gather(processed_chunk, root=0)

    if rank == 0:
        # Flatten the list of lists
        all_night_files = [file for sublist in gathered_data for file in sublist]
        # Save the results to a file or database if needed
        with open('/data/gdi/f/gzhao1/Database/processed_night_files.txt', 'w') as out_file:
            for night_file in all_night_files:
                out_file.write(f"{night_file}\n")
