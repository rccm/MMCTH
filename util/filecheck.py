from mpi4py import MPI
from pyhdf.SD import SD, HDF4Error
import numpy as np

# Function to read file paths from a text file
def read_file_paths(file_path):
    with open(file_path, 'r') as file:
        file_paths = [line.strip() for line in file]
    return file_paths

# Function to check readability of a single HDF4 file
def check_hdf4_readability(file_path):
    try:
        # Attempt to open the file
        hdf = SD(file_path)
        return None  # If file is readable, return None (ignore it)
    except HDF4Error:
        print(file_path)
        return file_path  # Return the file path if it's unreadable
    except Exception:
        print(file_path)
        return file_path  # Handle unexpected errors similarly by marking file as unreadable



file_to_check = '/data/keeling/a/gzhao1/f/Database/mod06_full.list'
# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# On rank 0, read file paths and distribute
if rank == 0:
    file_paths = read_file_paths(file_to_check)
    file_paths = np.array_split(file_paths, size)  # Evenly distribute file paths across processes
else:
    file_paths = None

# Scatter file paths to each process
file_chunk = comm.scatter(file_paths, root=0)

# Each process checks the readability of its assigned files and only keeps unreadable ones
local_unreadable_files = [check_hdf4_readability(f) for f in file_chunk if check_hdf4_readability(f) is not None]

# Gather unreadable files from all processes
unreadable_files = comm.gather(local_unreadable_files, root=0)

# On rank 0, write unreadable files to a text file
if rank == 0:
    # Flatten the list of unreadable files
    unreadable_files = [file for sublist in unreadable_files for file in sublist if file]
    with open('unreadable_mod06.txt', 'w') as f:
        for file in unreadable_files:
            f.write(f"{file}\n")
    print(f"Unreadable files have been written to unreadable_mod06.txt")
