"""doc"""
import argparse
import os
import subprocess

import json
import sqlite3

MARKDOWN = """
# Diagram status
"""

MARKDOWN = """
# {spectra}
- *t_max*: {t_max}
- *certain/uncertain differential: {num_diff1}/{num_diff2}
- *certain/uncertain homotopy relations: {num_htpy1}/{num_htpy2}
- *first uncertain differential: {deg_diff} d_{r}{x}=?
- *first uncertain relations: {deg_rel} {rel}=0
"""

PATH_SS_JSON = "C:/Users/lwnpk/OneDrive/Projects/algtop_cpp/ss/ss.json"


def load_db(path):
    conn_ring = sqlite3.connect(path)
    c_ring = conn_ring.cursor()

    c_ring.close()
    conn_ring.close()


if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edit", action="store_true", help="open the script in vscode")
    parser.add_argument("diagram", help="Name of the diagram")
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    with open(PATH_SS_JSON) as file_ss_json:
        ss_json = json.load(file_ss_json)
        print(ss_json)
