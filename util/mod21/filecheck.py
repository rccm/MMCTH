from mpi4py import MPI
from pyhdf.SD import SD, HDF4Error
import numpy as np
import os

# Function to read file paths and filter out already checked ones
def read_and_filter_files(file_path, log_file):
    checked_files = set()
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            checked_files = set(f.read().splitlines())
    
    with open(file_path, 'r') as f:
        all_files = [line.strip() for line in f if line.strip() not in checked_files]
    
    return all_files

# Function to check readability of a single HDF4 file
def check_hdf4_readability(file_path):
    try:
        SD(file_path)  # Try opening the file
        return None  # File is readable
    except (HDF4Error, Exception):
        print(file_path)
        return file_path  # Return if unreadable

def log_checked_file(file_path, log_file):
    with open(log_file, 'a') as f:
        f.write(f"{file_path}\n")


def log_unreadable_file(file_path, rank):
    unreadable_log_file = f'unreadable_mod21_{rank}.txt'
    with open(unreadable_log_file, 'a') as f:
        f.write(f"{file_path}\n")

# Main process
file_to_check = '/data/keeling/a/gzhao1/f/Database/mod021km_2022.list'
log_file = './checked_files.log'

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# On rank 0, read and filter file paths
if rank == 0:
    file_paths = read_and_filter_files(file_to_check, log_file)
    file_chunks = np.array_split(file_paths, size)  # Evenly distribute file paths
else:
    file_chunks = None

# Scatter file chunks to processes
file_chunk = comm.scatter(file_chunks, root=0)

# Check files and log them immediately
for file_path in file_chunk:
    result = check_hdf4_readability(file_path)
    log_checked_file(file_path, log_file)  # Log the file as soon as it's checked
    if result:
        log_unreadable_file(result, rank)  # Log unreadable file immediately

# No need to gather unreadable files at rank 0
