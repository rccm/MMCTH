#!/usr/bin/env python3

import sqlite3

def verify_merge():
    merged_db = "merged_database.sqlite"
    conn = sqlite3.connect(merged_db)
    cursor = conn.cursor()
    
    # Count groups and files
    cursor.execute("SELECT COUNT(*) FROM groups")
    group_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM files")
    file_count = cursor.fetchone()[0]
    
    print(f"Groups in merged database: {group_count}")
    print(f"Files in merged database: {file_count}")
    
    # Check for any orphaned files (files without groups)
    cursor.execute("""
        SELECT COUNT(*) FROM files 
        WHERE group_id NOT IN (SELECT group_id FROM groups)
    """)
    orphaned_files = cursor.fetchone()[0]
    print(f"Orphaned files: {orphaned_files}")
    
    # Check sample data
    cursor.execute("SELECT * FROM groups LIMIT 5")
    sample_groups = cursor.fetchall()
    print(f"\nSample groups: {sample_groups}")
    
    cursor.execute("SELECT * FROM files LIMIT 5")
    sample_files = cursor.fetchall()
    print(f"Sample files: {sample_files}")
    
    conn.close()

if __name__ == "__main__":
    verify_merge()