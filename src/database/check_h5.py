import os, sqlite3
from pathlib import Path

for year in range(2000,2023):

    db = f"/data/keeling/a/gzhao1/f/Database/inputfiles_{str(year)}.sqlite"
    con = sqlite3.connect(db)
    cur = con.cursor()

    rows = cur.execute("""
        SELECT file_id, file_path
        FROM files
        WHERE file_path LIKE '%MISR_AM1_TC_CLOUD%' AND file_path LIKE '%.h5'
    """).fetchall()

    fixed = 0
    missing = []
    for fid, path in rows:
        cand = Path(path).with_suffix(".hdf")
        if cand.exists():
            cur.execute("UPDATE files SET file_path = ? WHERE file_id = ?", (str(cand), fid))
            fixed += 1
        else:
            missing.append((fid, path, str(cand)))

    con.commit()
    con.close()

    print(f"updated: {fixed}, missing: {len(missing)}")
    for m in missing[:10]:
        print("No .hdf sibling:", m)
