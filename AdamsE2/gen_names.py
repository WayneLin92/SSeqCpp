"""doc"""
import argparse
import os
import subprocess
import sqlite3

gen_names = {
0: "h_0",
1: "h_1",
2: "h_2",
3: "h_3",
7: "h_4",
18: "h_5",
69: "h_6",
324: "h_7",
}

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('--db', default=R"C:\Users\lwnpk\Documents\Projects\AdamsSSDB\AdamsSS.db", help='open the script in vscode')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    conn = sqlite3.connect(args.db)
    with conn:
        c = conn.cursor()
        sql = f"Update AdamsE2_generators SET name=?2 where id=?1"
        c.executemany(sql, gen_names.items())
        c.close()
    conn.close()