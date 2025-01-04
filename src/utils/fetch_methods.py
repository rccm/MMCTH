import sqlite3
import argparse

def create_connection(db_file):
    """ Create a database connection to the SQLite database specified by db_file """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except sqlite3.Error as e:
        print(e)
        return None

def fetch_files_by_group(conn, group_id):
    cursor = conn.cursor()
    cursor.execute('''
        SELECT file_path FROM files WHERE group_id = ?
    ''', (group_id,))
    files = cursor.fetchall()
    return [file[0] for file in files]

def fetch_files_by_orbit(conn, orbit_number,groud_id_only = False):
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE orbit_number = ?
    ''', (orbit_number,))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def fetch_files_by_path(conn, path_number,groud_id_only = False):
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE path_number = ?
    ''', (path_number,))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def fetch_files_by_date(conn, date):
    date_time = date + ' 00:00:00'
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE date = ?
    ''', (date_time,))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def fetch_files_by_month(conn, year, month):
    cursor = conn.cursor()
    start_date = f"{year}-{month:02d}-01"
    end_date = f"{year}-{month:02d}-31"  # Assumes all months have 31 days, SQLite handles overflow dates correctly.
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE date BETWEEN ? AND ?
    ''', (start_date, end_date))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def fetch_files_by_modisid(conn, modis_id):
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE modis_id LIKE ?
    ''', ('%' + modis_id + '%',))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def fetch_files_by_year(conn, year):
    cursor = conn.cursor()
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"
    cursor.execute('''
        SELECT group_id, file_path FROM files WHERE date BETWEEN ? AND ?
    ''', (start_date, end_date))
    rows = cursor.fetchall()
    files_by_group = {}
    for group_id, file_path in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
    return files_by_group

def process_groups(groups):
    for group_id, file_paths in groups.items():
        print(f"Processing group {group_id} {file_paths}")

def fetch_all_groups(conn):
    """ Fetch all distinct group IDs in the database """
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT group_id FROM files
    ''')
    groups = cursor.fetchall()
    return [group[0] for group in groups]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process satellite data for a given year, month, or date.')
    parser.add_argument('-y', '--year', type=int, required=True, help='Year to process (e.g., 2016)')
    parser.add_argument('-m', '--month', type=int, help='Month to process (optional, e.g., 3 for March)')
    parser.add_argument('-d', '--date', type=str, help='Day to process in MM-DD format (optional)')
    parser.add_argument('-o', '--orbit', type=str, help='Orbit to process (optional)')
    parser.add_argument('-i', '--modisid', type=str, help='MODIS ID (optonal e.g. A2001.1234')
    parser.add_argument('-p', '--path', type=int, help='Path number (optional)')
    
    args = parser.parse_args()
    year = args.year
    month = args.month
    date = args.date
    orbit = args.orbit
    modisid = args.modisid

    # Construct database path
    database_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_{year}.sqlite"
    
    conn = create_connection(database_path)
    if conn is None:
        print("Error! cannot create the database connection.")
    if date:
        # If --date is provided
        month, day = map(int, date.split('-'))
        full_date = f"{year}-{month:02d}-{day:02d}"
        groups = fetch_files_by_date(conn, full_date)
    elif month:
        # If only --month is provided
        groups = fetch_files_by_month(conn, year, month)
    elif orbit:
        groups = fetch_files_by_orbit(conn, orbit)
    elif modisid:
        groups = fetch_files_by_modisid(conn, modisid) 
    else:
        # Default: process entire year
        groups = fetch_files_by_year(conn, year)
 
    process_groups(groups)