from pyhdf.SD import SD, HDF4Error

# Function to read file paths from a text file
def read_file_paths(file_path):
    with open(file_path, 'r') as file:
        file_paths = [line.strip() for line in file]
    return file_paths

# Function to check readability of HDF4 files and dump unreadable files to a text file
def check_hdf4_readability(file_paths):
    unreadable_files = []
    
    for file_path in file_paths:
        try:
            # Attempt to open the file
            hdf = SD(file_path)
            print(f"File {file_path} is readable.")
        except HDF4Error as e:
            print(f"File {file_path} is not readable. Error: {e}")
            unreadable_files.append(file_path)
        except Exception as e:
            print(f"An unexpected error occurred with file {file_path}. Error: {e}")
            unreadable_files.append(file_path)
    
    # Write unreadable files to a text file
    with open('unreadable_files_maiac_2.txt', 'w') as f:
        for file in unreadable_files:
            f.write(f"{file}\n")
    print(f"Unreadable files have been written to unreadable_files.txt")

# Read file paths from file_paths.txt
file_paths = read_file_paths('../lists/MOD21KM_2019.list')
file_paths = read_file_paths('./unreadable_files_maiac.txt')

# Run the check
check_hdf4_readability(file_paths)
