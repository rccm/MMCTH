#!/bin/bash

merged_db="merged_database.sqlite"
rm -f "$merged_db"

# Create the merged database with schema
first_db="/data/keeling/a/gzhao1/f/Database/inputfiles_2000.sqlite"
sqlite3 "$first_db" ".schema groups" | sqlite3 "$merged_db"
sqlite3 "$first_db" ".schema files" | sqlite3 "$merged_db"

# Get list of all database files sorted by year
dbs=$(ls /data/keeling/a/gzhao1/f/Database/inputfiles_2*.sqlite | sort)

for db in $dbs; do
    year=$(echo "$db" | grep -oE '[0-9]{4}')
    echo "Processing $db (year: $year)"
    
    sqlite3 "$merged_db" << EOF
ATTACH DATABASE '$db' AS source;

-- Insert groups, ignoring duplicates
INSERT OR IGNORE INTO main.groups (group_name, group_datetime)
SELECT group_name, group_datetime FROM source.groups;

-- Insert files with proper group_id mapping, ignoring duplicates
INSERT OR IGNORE INTO main.files (group_id, file_path, orbit_number, date, modis_id, path_number)
SELECT 
    main_groups.group_id,
    source_files.file_path, 
    source_files.orbit_number, 
    source_files.date, 
    source_files.modis_id, 
    source_files.path_number
FROM source.files AS source_files
JOIN source.groups AS source_groups ON source_files.group_id = source_groups.group_id
JOIN main.groups AS main_groups ON source_groups.group_name = main_groups.group_name;

DETACH DATABASE source;
EOF

    echo "Completed $db"
done

# Create indexes after all data is inserted for better performance
sqlite3 "$merged_db" << EOF
CREATE INDEX IF NOT EXISTS idx_groups_group_name ON groups (group_name);
CREATE INDEX IF NOT EXISTS idx_files_group_id ON files (group_id);
EOF

echo "Merging complete: $merged_db"