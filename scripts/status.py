"""doc"""
import argparse
import os
import subprocess

import json
import sqlite3
from rich import print
from rich.text import Text

TITLE = """
[bold]Diagram status[bold]
"""

colors = ["red", "blue", "green", "orange", "purple"]

TEMPLATE = """
### {spectra}
- t_max: [{c}]{t_max}[/{c}]
- certain/uncertain differentials: [{c}]{num_diff1}/{num_diff2} [/{c}]
- certain/uncertain homotopy relations: [{c}]{num_htpy1}/{num_htpy2} [/{c}]
- first uncertain differential: [{c}]{deg_diff} d_{r}({x})=? [/{c}]
- first uncertain relations: [{c}]{deg_rel} {rel}=0 [/{c}]
"""

PATH_SS_JSON = "C:/Users/lwnpk/OneDrive/Projects/algtop_cpp/ss/ss.json"
DIR_DB_ROOT = "C:/Users/lwnpk/Documents/Projects/algtop_cpp_build/bin/Release"


def load_db(path, name):
    conn = sqlite3.connect(os.path.join(DIR_DB_ROOT, path))
    c = conn.cursor()

    result = {}

    res = c.execute(f"select max(t) from {name}_AdamsE2_ss")
    result["t_max"] = res.fetchone()[0]

    res = c.execute(
        f"select count(*) from {name}_AdamsE2_ss where level>9800 and diff is not NULL"
    )
    result["num_diff1"] = res.fetchone()[0]

    res = c.execute(
        f"select count(*) from {name}_AdamsE2_ss where level>9800 and diff is NULL"
    )
    result["num_diff2"] = res.fetchone()[0]

    res = c.execute(
        f'select count(*) from {name}_pi_relations where not instr(name, "O")'
    )
    result["num_htpy1"] = res.fetchone()[0]

    res = c.execute(f'select count(*) from {name}_pi_relations where instr(name, "O")')
    result["num_htpy2"] = res.fetchone()[0]

    res = c.execute(
        f"select t-s, s, 10000-level, base from {name}_AdamsE2_ss where level>9800 and diff is NULL ORDER BY t-s,s"
    )
    if res_fetchone := res.fetchone():
        result["deg_diff"] = (res_fetchone[0], res_fetchone[1])
        result["r"] = res_fetchone[2]
        result["x"] = res_fetchone[3]
    else:
        result["deg_diff"] = None
        result["r"] = None
        result["x"] = None

    res = c.execute(
        f'select t-s, s, name from {name}_pi_relations where instr(name, "O") ORDER BY t-s,s'
    )
    if res_fetchone := res.fetchone():
        result["deg_rel"] = (res_fetchone[0], res_fetchone[1])
        result["rel"] = res_fetchone[2]
    else:
        result["deg_rel"] = None
        result["rel"] = None

    c.close()
    conn.close()
    return result


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

    print(TITLE)
    spectra = ss_json["diagrams"][args.diagram]["spectra"]
    dir_db = ss_json["diagrams"][args.diagram]["dir"]
    for i, name in enumerate(spectra):
        print(spectra[name])
        result = load_db(os.path.join(dir_db, spectra[name]), name)
        out = TEMPLATE.format(
            spectra=name,
            t_max=result["t_max"],
            num_diff1=result["num_diff1"],
            num_diff2=result["num_diff2"],
            num_htpy1=result["num_htpy1"],
            num_htpy2=result["num_htpy2"],
            deg_diff=result["deg_diff"],
            r=result["r"],
            x=result["x"],
            deg_rel=result["deg_rel"],
            rel=result["rel"],
            c=colors[i],
        )

        print(out)
