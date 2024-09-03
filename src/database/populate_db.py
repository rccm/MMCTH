import os
import fnmatch
import sqlite3
import sys
from datetime import datetime
from find_groupfiles import find_group_files
sys.path.append(os.path.abspath('../'))
from data_readers.mod_read import MODL1Granule
from mpi4py import MPI
def read_file_paths(file_path):
    with open(file_path, 'r') as file:
        file_paths = [line.strip() for line in file]
    return file_paths

def create_connection(db_file):
    """ Create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
        return conn
    except Exception as e:
        print(e)
        return None

def create_tables(conn):
    """ Create tables in the SQLite database """
    cursor = conn.cursor()
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS groups (
        group_id INTEGER PRIMARY KEY AUTOINCREMENT,
        group_name TEXT UNIQUE NOT NULL,
        group_datetime TEXT
    )
    ''')
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS files (
        file_id INTEGER PRIMARY KEY AUTOINCREMENT,
        group_id INTEGER,
        file_path TEXT,
        orbit_number INTEGER,
        date TEXT,
        modis_id TEXT,
        path_number INTEGER,
        FOREIGN KEY (group_id) REFERENCES groups(group_id)
        UNIQUE(group_id, file_path)           
    )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_groups_group_name ON groups (group_name)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_files_group_id ON files (group_id)')
    conn.commit()

def insert_group_files_batch(conn, groups):
    cursor = conn.cursor()
    for group_name, group_data in groups.items():
        cursor.execute('SELECT group_id FROM groups WHERE group_name = ?', (group_name,))
        existing_group = cursor.fetchone()
        if existing_group:
            group_id = existing_group[0]
            current_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            cursor.execute("UPDATE groups SET group_datetime = ? WHERE group_id = ?", (current_datetime, group_id))
        else:
            cursor.execute('INSERT INTO groups (group_name) VALUES (?)', (group_name,))
            group_id = cursor.lastrowid
            current_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            cursor.execute("UPDATE groups SET group_datetime = ? WHERE group_id = ?", (current_datetime, group_id))

        files_to_insert = [(group_id, fp, group_data['orbit_number'], group_data['date'], group_data['modis_id'], group_data['path_number']) for fp in group_data['file_paths']]
        cursor.executemany('''
            INSERT INTO files (group_id, file_path, orbit_number, date, modis_id, path_number) 
            VALUES (?, ?, ?, ?, ?, ?)
        ''', files_to_insert)
    conn.commit()

def generate_groups():
    groups = {}
    mod_l1b_fullpaths = read_file_paths('../../lists/MOD21KM_2019.list')
    for mod_l1b_fullpath in mod_l1b_fullpaths:
        mod_06, mod_03, tc_cloud, agp_file,orbit, mod_id, path_number, dd = find_group_files(mod_l1b_fullpath)
        if mod_06 and mod_03 and tc_cloud:
            group_name = f"group_{mod_id}"
            groups[group_name] = {
                "file_paths": [mod_06, mod_03, tc_cloud,agp_file],
                "orbit_number": orbit,
                "date": dd,
                "modis_id": mod_id,
                "path_number": path_number
            }
    return groups

def process_files(file_paths):
    groups = {}
    for mod_l1b_fullpath in file_paths:
        group_data = find_group_files(mod_l1b_fullpath)
        if group_data is None:
            continue
        group_name = f"group_{group_data['mod_id']}"
        groups[group_name] = {
            "file_paths": [mod_l1b_fullpath,group_data['mod_06'], group_data['mod_03'], group_data['tc_cloud'], group_data['agp'], \
                           group_data['era5_single'],group_data['era5_multi'],group_data['era5_single_extra'],group_data['era5_multi_extra']],
            "orbit_number": group_data['orbit'],
            "date": group_data['date'],
            "modis_id": group_data['mod_id'],
            "path_number": group_data['path_number']
        }
    return groups

def main_single_core():
    database_path = "/data/keeling/a/gzhao1/f/mmcth/lists/inputfiles_2002.sqlite"
    conn = create_connection(database_path)
    if conn is not None:
        create_tables(conn)
        groups = generate_groups()
        batch_size = 1000  # Adjust based on performance and memory constraints
        group_items = list(groups.items())
        for i in range(0, len(group_items), batch_size):
            batch = group_items[i:i + batch_size]
            for group_name, group_data in batch:
                files_data = [{"file_path": fp} for fp in group_data['file_paths']]
                insert_group_files_batch(conn, group_name, group_data, files_data)
        
        conn.close()
    else:
        print("Error! cannot create the database connection.")
def vacuum_database(conn):
    """Apply the VACUUM command to the SQLite database to optimize and reduce file size."""
    cursor = conn.cursor()
    cursor.execute('VACUUM')
    conn.commit()

def main_parrallel():
    year = '2016'

    list_filepath = f'/data/keeling/a/gzhao1/f/Database/MOD21KM_{year}.list'  

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        mod_l1b_fullpaths = read_file_paths(list_filepath)
        # mod_l1b_fullpaths = read_file_paths('valid_file_list.txt')
        print(len(mod_l1b_fullpaths))
        chunk_size = (len(mod_l1b_fullpaths) + size - 1) // size  # ensures all chunks are distributed
        chunks = [mod_l1b_fullpaths[i:i + chunk_size] for i in range(0, len(mod_l1b_fullpaths), chunk_size)]
        if len(chunks) < size:  # ensure that we have exactly `size` chunks
            for _ in range(size - len(chunks)):
                chunks.append([])
    else:
        chunks = None

    file_chunk = comm.scatter(chunks, root=0)
    local_groups = process_files(file_chunk)
    all_groups = comm.gather(local_groups, root=0)

    if rank == 0:
        groups = {}
        for group in all_groups:
            groups.update(group)
        database_path = f"/data/gdi/f/gzhao1/Database/inputfiles_{year}.sqlite"
        conn = create_connection(database_path)
        if conn is not None:
            create_tables(conn)
            insert_group_files_batch(conn, groups)
            vacuum_database(conn)  
            conn.close()
        else:
            print("Error! cannot create the database connection.")

if __name__ == "__main__":
    main_parrallel()