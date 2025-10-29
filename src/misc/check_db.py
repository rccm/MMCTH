import sqlite3
import os

def inspect_database(db_path):
    if not os.path.exists(db_path):
        print(f"Database file does not exist: {db_path}")
        return
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # List all tables
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()
    print(f"Tables: {tables}")
    
    # If files table exists, show its structure and some data
    if any('files' in table for table in tables):
        cursor.execute("PRAGMA table_info(files)")
        columns = cursor.fetchall()
        print("Files table structure:")
        for col in columns:
            print(f"  {col}")
        
        cursor.execute("SELECT * FROM files LIMIT 5")
        sample_data = cursor.fetchall()
        print("Sample data:")
        for row in sample_data:
            print(f"  {row}")
    
    conn.close()

# Test with your specific year
year = 2003  # change to the year you're testing
database_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_{year}.sqlite"
inspect_database(database_path)