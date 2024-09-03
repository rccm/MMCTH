import os
import fnmatch
import sqlite3
import sys
from datetime import datetime


def create_connection(db_file):
    """ Create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
        return conn
    except Exception as e:
        print(e)
        return None

def query_groups_by_year(conn, year):
    cursor = conn.cursor()
    query = '''
    SELECT g.group_name, f.file_path, f.orbit_number, f.date, f.modis_id, f.path_number
    FROM groups g
    JOIN files f ON g.group_id = f.group_id
    WHERE f.date LIKE ?
    '''
    cursor.execute(query, (f'{year}%',))
    rows = cursor.fetchall()
    
    # Process the results to group by group_name
    results = {}
    for row in rows:
        group_name = row[0]
        file_info = {
            'file_path': row[1],
            'orbit_number': row[2],
            'date': row[3],
            'modis_id': row[4],
            'path_number': row[5]
        }
        if group_name not in results:
            results[group_name] = []
        results[group_name].append(file_info)
    
    return results

def print_groups_files(groups):
    for group_name, files in groups.items():
        print(f"Group: {group_name}")
        for file_info in files:
            print(f"  File Path: {file_info['file_path']}")
            print(f"  Orbit Number: {file_info['orbit_number']}")
            print(f"  Date: {file_info['date']}")
            print(f"  MODIS ID: {file_info['modis_id']}")
            print(f"  Path Number: {file_info['path_number']}")
        print()

def main():
    year = 2002
    database_path = "/data/keeling/a/gzhao1/f/mmcth/lists/inputfiles.sqlite"
    conn = create_connection(database_path)
    if conn is not None:
        groups = query_groups_by_year(conn, year)
        print_groups_files(groups)
        conn.close()
    else:
        print("Error! cannot create the database connection.")

if __name__ == "__main__":
    main()