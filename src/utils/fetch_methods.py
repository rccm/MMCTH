import sqlite3

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
        print(f"Processing group {group_id}")

def fetch_all_groups(conn):
    """ Fetch all distinct group IDs in the database """
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT group_id FROM files
    ''')
    groups = cursor.fetchall()
    return [group[0] for group in groups]


if __name__ == "__main__":
    
    database_path = "/data/keeling/a/gzhao1/f/mmcth/lists/inputfiles.sqlite"
    conn = create_connection(database_path)
    if conn is None:
        print("Error! cannot create the database connection.")
    
    # Example: Fetching and processing files by group
    group_id = 10
    date = '2002-12-03'
    groups = fetch_files_by_date(conn, date)
    # process_groups(groups)
    groups = fetch_files_by_month(conn, 2002, 12)
    process_groups(groups)