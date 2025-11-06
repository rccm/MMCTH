import sqlite3
import argparse
import subprocess
import re
from datetime import datetime
import calendar

def _range_bounds(start_date_str, end_date_str):
    """Return inclusive timestamp bounds 'YYYY-MM-DD 00:00:00' to 'YYYY-MM-DD 23:59:59'."""
    # Validates YYYY-MM-DD and builds safe bounds
    datetime.strptime(start_date_str, "%Y-%m-%d")
    datetime.strptime(end_date_str, "%Y-%m-%d")
    start_ts = f"{start_date_str} 00:00:00"
    end_ts   = f"{end_date_str} 23:59:59"
    return start_ts, end_ts

def _add_months(yyyy_mm_dd, n):
    """Add n months to a YYYY-MM-DD date, clamping day to end-of-month if needed."""
    dt = datetime.strptime(yyyy_mm_dd, "%Y-%m-%d")
    y = dt.year + (dt.month - 1 + n) // 12
    m = (dt.month - 1 + n) % 12 + 1
    last_day = calendar.monthrange(y, m)[1]
    d = min(dt.day, last_day)
    return f"{y:04d}-{m:02d}-{d:02d}"

def fetch_files_by_date_range(conn, start_date, end_date):
    """
    Inclusive date range (YYYY-MM-DD to YYYY-MM-DD).
    Matches all times on end day by extending to 23:59:59.
    """
    cursor = conn.cursor()
    start_ts, end_ts = _range_bounds(start_date, end_date)
    cursor.execute('''
        SELECT group_id, file_path, orbit_number
        FROM files
        WHERE date >= ? AND date <= ?
    ''', (start_ts, end_ts))
    rows = cursor.fetchall()
    files_by_group, orbit_numbers = {}, set()
    for group_id, file_path, orbit_number in rows:
        files_by_group.setdefault(group_id, []).append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, sorted(orbit_numbers)

def fetch_files_next_n_months(conn, start_date, n_months):
    """
    From start_date (YYYY-MM-DD) through the **last day** of the month that is
    (n_months-1) months later. For example, start=2001-01-01, n=6 → through 2001-06-30.
    """
    if n_months <= 0:
        return {}, []
    # Compute inclusive end date = last day of target month
    end_anchor = _add_months(start_date, n_months - 1)  # same day in target month
    y, m = map(int, end_anchor.split("-")[:2])
    last_day = calendar.monthrange(y, m)[1]
    end_date = f"{y:04d}-{m:02d}-{last_day:02d}"
    return fetch_files_by_date_range(conn, start_date, end_date)

def get_year_from_orbit(orbit_number):
    """
    Use MtkOrbitToTimeRange to get the date for an orbit number and extract the year
    Returns (year, full_date_string) or (None, None) on error
    """
    try:
        # Run MtkOrbitToTimeRange command
        result = subprocess.run(
            ['MtkOrbitToTimeRange', str(orbit_number)],
            capture_output=True, 
            text=True, 
            check=True
        )
        # Parse the output to get the date
        lines = result.stdout.strip().split('\n')
        if len(lines) >= 1:
            timestamp = lines[0].strip()
            # Use regex to extract date in YYYY-MM-DD format
            match = re.search(r'(\d{4}-\d{2}-\d{2})', timestamp)
            if match:
                date_str = match.group(1)
                year = int(date_str.split('-')[0])
                return year, date_str
        return None, None
        
    except FileNotFoundError:
        print("Error: MtkOrbitToTimeRange command not found. Make sure it's in your PATH.")
        return None, None
    except subprocess.CalledProcessError as e:
        print(f"Error running MtkOrbitToTimeRange for orbit {orbit_number}: {e}")
        print(f"Stderr: {e.stderr}")
        return None, None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None, None

def get_orbit_from_date(date_string):
    """
    Use MtkTimeToOrbit to get orbit number from a date string
    Returns orbit_number or None on error
    """
    try:
        # Convert date to format expected by MtkTimeToOrbit (YYYY-MM-DDTHH:MM:SSZ)
        # If only date is provided, use noon as default time
        if 'T' not in date_string:
            date_string = date_string + 'T12:00:00Z'
        
        # Run MtkTimeToOrbit command
        result = subprocess.run(
            ['MtkTimeToOrbit', date_string],
            capture_output=True, 
            text=True, 
            check=True
        )
        
        # Parse the output to get orbit number
        orbit_match = re.search(r'(\d+)', result.stdout.strip())
        if orbit_match:
            return int(orbit_match.group(1))
        return None
        
    except FileNotFoundError:
        print("Error: MtkTimeToOrbit command not found. Make sure it's in your PATH.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error running MtkTimeToOrbit for date {date_string}: {e}")
        print(f"Stderr: {e.stderr}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

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
        SELECT file_path, orbit_number FROM files WHERE group_id = ?
    ''', (group_id,))
    files = cursor.fetchall()
    return [file[0] for file in files], files[0][1] if files else None

def fetch_files_by_orbit(conn, orbit_number, groud_id_only=False):
    cursor = conn.cursor()
    try:
        cursor.execute('''
            SELECT group_id, file_path FROM files WHERE orbit_number = ?
        ''', (orbit_number,))
        rows = cursor.fetchall()
        files_by_group = {}
        for group_id, file_path in rows:
            if group_id not in files_by_group:
                files_by_group[group_id] = []
            files_by_group[group_id].append(file_path)
        return files_by_group, orbit_number
    except sqlite3.Error as e:
        print(f"SQLite error in fetch_files_by_orbit: {e}")
        return {}, None

def fetch_files_by_path(conn, path_number, groud_id_only=False):
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path, orbit_number FROM files WHERE path_number = ?
    ''', (path_number,))
    rows = cursor.fetchall()
    files_by_group = {}
    orbit_numbers = set()
    for group_id, file_path, orbit_number in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, list(orbit_numbers)

def fetch_files_by_date(conn, date):
    date_time = date + ' 00:00:00'
    cursor = conn.cursor()
    cursor.execute('''
        SELECT group_id, file_path, orbit_number FROM files WHERE date = ?
    ''', (date_time,))
    rows = cursor.fetchall()
    files_by_group = {}
    orbit_numbers = set()
    for group_id, file_path, orbit_number in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, list(orbit_numbers)

def fetch_files_by_month(conn, year, month):
    cursor = conn.cursor()
    start_date = f"{year}-{month:02d}-01"
    end_date = f"{year}-{month:02d}-31"
    cursor.execute('''
        SELECT group_id, file_path, orbit_number FROM files WHERE date BETWEEN ? AND ?
    ''', (start_date, end_date))
    rows = cursor.fetchall()
    files_by_group = {}
    orbit_numbers = set()
    for group_id, file_path, orbit_number in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, list(orbit_numbers)

def fetch_files_by_modisid(conn, modis_id):
    cursor = conn.cursor()
    
    cursor.execute('''
        SELECT group_id, file_path, orbit_number FROM files WHERE modis_id LIKE ?
    ''', ('%' + modis_id + '%',))
    rows = cursor.fetchall()
    files_by_group = {}
    orbit_numbers = set()
    for group_id, file_path, orbit_number in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, list(orbit_numbers)

def fetch_files_by_year(conn, year):
    cursor = conn.cursor()
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"
    cursor.execute('''
        SELECT group_id, file_path, orbit_number FROM files WHERE date BETWEEN ? AND ?
    ''', (start_date, end_date))
    rows = cursor.fetchall()
    files_by_group = {}
    orbit_numbers = set()
    for group_id, file_path, orbit_number in rows:
        if group_id not in files_by_group:
            files_by_group[group_id] = []
        files_by_group[group_id].append(file_path)
        orbit_numbers.add(orbit_number)
    return files_by_group, list(orbit_numbers)

def process_groups(groups, orbit_info):
    """Enhanced to show orbit information"""
    if isinstance(orbit_info, list) and len(orbit_info) > 0:
        print(f"Found {len(orbit_info)} unique orbit(s): {orbit_info}")
    elif orbit_info is not None:
        print(f"Orbit: {orbit_info}")
    
    for group_id, file_paths in groups.items():
        print(f"Processing group {group_id}:")
        for file_path in file_paths:
            print(f"  {file_path}")

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
    parser.add_argument('-y', '--year', type=int, help='Year to process (e.g., 2016)')
    parser.add_argument('-m', '--month', type=int, help='Month to process (optional, e.g., 3 for March)')
    parser.add_argument('-d', '--date', type=str, help='Day to process in MM-DD format (optional)')
    parser.add_argument('-o', '--orbit', type=str, help='Orbit to process (optional)')
    parser.add_argument('-i', '--modisid', type=str, help='MODIS ID (optional e.g. A2001.1234')
    parser.add_argument('-p', '--path', type=int, help='Path number (optional)')
    parser.add_argument('-s', '--start-date', type=str, help='Start date YYYY-MM-DD (range mode)')
    parser.add_argument('-e', '--end-date', type=str, help='End date YYYY-MM-DD (inclusive, range mode)')
    parser.add_argument('-n', '--months', type=int, help='If set with --start-date, span this many months (inclusive)')
    
    args = parser.parse_args()
    year = args.year
    month = args.month
    date = args.date
    orbit = args.orbit
    modisid = args.modisid
    path = args.path

    # Construct database path
    database_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_all.sqlite"
    
    conn = create_connection(database_path)
    if conn is None:
        print("Error! cannot create the database connection.")
        exit(1)
    
    # Get orbit information for different query types
    orbit_info = None
    if args.start_date and args.end_date:
        groups, orbit_info = fetch_files_by_date_range(conn, args.start_date, args.end_date)
    # Next: start date + N months (inclusive to last day of the Nth month)
    elif args.start_date and args.months:
        groups, orbit_info = fetch_files_next_n_months(conn, args.start_date, args.months)
    elif date:
        # If --date is provided
        month_num, day = map(int, date.split('-'))
        full_date = f"{year}-{month_num:02d}-{day:02d}"
        groups, orbit_info = fetch_files_by_date(conn, full_date)
        
        # Also get orbit from MtkTimeToOrbit for comparison
        calculated_orbit = get_orbit_from_date(full_date)
        if calculated_orbit:
            print(f"Calculated orbit from date {full_date}: {calculated_orbit}")
            
    elif month:
        # If only --month is provided
        groups, orbit_info = fetch_files_by_month(conn, year, month)
    elif orbit:
        # Auto-detect year from orbit if needed
        orbit_year, orbit_date = get_year_from_orbit(orbit)
        if orbit_year and orbit_year != year:
            print(f"Warning: Orbit {orbit} corresponds to year {orbit_year}, but you specified year {year}")
            print(f"Using database for year {orbit_year} instead")
            database_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_{orbit_year}.sqlite"
            conn = create_connection(database_path)
        groups, orbit_info = fetch_files_by_orbit(conn, orbit)
    elif modisid:
        groups, orbit_info = fetch_files_by_modisid(conn, modisid)
    elif path:
        groups, orbit_info = fetch_files_by_path(conn, path)
    else:
        # Default: process entire year
        groups, orbit_info = fetch_files_by_year(conn, year)
 
    process_groups(groups, orbit_info)