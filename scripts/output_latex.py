import sqlite3
import os
from collections import defaultdict


def get_complex_name(path):
    path = os.path.basename(path)
    names = ["S0", "C2", "Ceta", "Cnu", "Csigma", "RP10", "RPinf", "X2"]
    for name in names:
        if path.startswith(name):
            return name
    raise ValueError(f"{path=} is not recognized")


def use_name(name) -> bool:
    return len(name) <= 100 and "+" not in name  # TODO: improve


def new_row(
    indices, table_head, row, nRows, nTables, hLine: bool, greyRow: bool
) -> str:
    result = ""
    if indices["row"] == 1:
        result += table_head
    elif hLine:
        result += "\hline\n"
    if greyRow:
        result += "\\rowcolor{Gainsboro!60}\n"
    result += row
    if indices["table"] <= nTables:
        nRows -= 4
    if indices["row"] == nRows:
        result += "\\end{tabular}\\hspace{1pt}\\hfill\n"
        if indices["table"] % nTables == 0:
            result += "\n"
        indices["row"] = 0
        indices["table"] += 1
    indices["row"] += 1
    return result


def empty_rows(index, nRow, nColumn) -> str:
    result = ""
    if index % nRow != 1:
        while index % nRow != 1:
            index += 1
            result += " ~ ".join(("~" for _ in range(nColumn))) + "\\\\\n"
        result += "\\end{tabular}\\hspace{1pt}\\hfill\n"
    return result


############################### Indecomposables ##########################################
def str_gen_names(c, table, letter):
    nRows = 52
    nColumn = 5
    nTables = 4
    table_head = R"""\begin{tabular}{lrr}\hline
name & stem & s \\\hline
"""

    sql = f"SELECT name, s, t FROM {table} order by s, t"
    gens_group_by_st = defaultdict(int)
    result = ""
    indices = {"row": 1, "table": 1}
    prev_s = -1
    for name, s, t in c.execute(sql):
        name1 = name
        if not name or name.startswith("x_{") or not use_name(name):
            gens_group_by_st[(s, t)] += 1
            if gens_group_by_st[(s, t)] == 1:
                name1 = f"{letter}_{{{t-s}, {s}}}"
            else:
                name1 = f"{letter}_{{{t-s}, {s}, {gens_group_by_st[(s, t)]}}}"
        row = f"${name1}$ & {t-s} & {s}\\\\\n"
        result += new_row(
            indices, table_head, row, nRows, nTables, hLine=s != prev_s, greyRow=False
        )
        prev_s = s
    result += empty_rows(indices["row"], nRows, nColumn)
    return result


def str_gen_alias(c, table, letter):
    nRows = 48
    nColumn = 5
    nTables = 1
    table_head = R"""\begin{tabular}{lrrl}\hline
name & stem & s & alias \\\hline
"""
    sql = f"SELECT name, s, t FROM {table} order by s, t"
    gens_group_by_st = defaultdict(int)
    result = ""
    indices = {"row": 1, "table": 1}
    for name, s, t in c.execute(sql):
        name1 = name
        if name and not name.startswith("x_{") and not use_name(name):
            gens_group_by_st[(s, t)] += 1
            if gens_group_by_st[(s, t)] == 1:
                name1 = f"{letter}_{{{t-s}, {s}}}"
            else:
                name1 = f"{letter}_{{{t-s}, {s}, {gens_group_by_st[(s, t)]}}}"
            row = f"${name1}$ & {t-s} & {s} & ${name}$ \\\\\n"
            result += new_row(
                indices, table_head, row, nRows, nTables, hLine=False, greyRow=False
            )
    result += empty_rows(indices["row"], nRows, nColumn)
    return result


def str_gen_names_from_db(cw, path_ring):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    result = R"\subsection{Indecomposables}" + "\n"
    result += str_gen_names(c_ring, cw + "_AdamsE2_generators", "x")
    result += R"\subsection{Other names of some of the indecomposables}" + "\n"
    result += str_gen_alias(c_ring, cw + "_AdamsE2_generators", "x")

    c_ring.close()
    conn_ring.close()

    return result


############################### Basis ##########################################
def load_gen_names(c, table, letter):
    sql = f"SELECT name, s, t FROM {table} order by id"
    gens_group_by_st = defaultdict(int)
    result = []
    for name, s, t in c.execute(sql):
        name1 = name
        if not name or name.startswith("x_{") or not use_name(name):
            gens_group_by_st[(s, t)] += 1
            if gens_group_by_st[(s, t)] == 1:
                name1 = f"{letter}_{{{t-s}, {s}}}"
            else:
                name1 = f"{letter}_{{{t-s}, {s}, {gens_group_by_st[(s, t)]}}}"
        result.append(name1)
    return result


def str2array(str_array: str):
    return tuple(int(i) for i in str_array.split(","))


def latex_mon(str_mon, gen_names):
    result = ""
    mon = str2array(str_mon)
    for i in range(0, len(mon), 2):
        base = gen_names[mon[i]]
        if "^" in base:
            base = f"{{{base}}}"
        if mon[i + 1] == 1:
            result += f"{base}"
        elif mon[i + 1] < 10:
            result += f"{base}^{mon[i + 1]}"
        else:
            result += f"{base}^{{{mon[i + 1]}}}"
    if result == "":
        result = "1"
    return result


def str_basis(c, table, gen_names):
    nRows = 52
    nColumns = 4
    nTables = 3
    table_head = R"""\begin{tabular}{lrrr}\hline
mon & stem & s & i\\\hline
"""
    sql = f"SELECT mon, s, t FROM {table} where t-s>0 and s<=(t-s)/2+1 order by t-s, s, id"
    result = ""
    indices = {"row": 1, "table": 1, "no": 1}
    prev_s, prev_t = -1, -1
    color = 0
    for str_mon, s, t in c.execute(sql):
        str_stem = str(t - s)
        str_i = None
        if (s, t) != (prev_s, prev_t):
            color = 1 - color
            str_stem = f"\\hypertarget{{basis:{s}:{t}}}{{{t-s}}}"
            str_i = f"\\hyperlink{{rel:{s}:{t}}}{{{1}}}"
            indices["no"] = 1
        str_i = str_i or indices["no"]
        row = f"${latex_mon(str_mon, gen_names)}$ & {str_stem} & {s} & {str_i}\\\\\n"
        result += new_row(
            indices,
            table_head,
            row,
            nRows,
            nTables,
            hLine=t - s != prev_t - prev_s,
            greyRow=color,
        )
        indices["no"] += 1
        prev_s, prev_t = s, t

    result += empty_rows(indices["row"], nRows, nColumns)
    return result


def str_basis_from_db(cw, path_ring):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    gen_names = load_gen_names(c_ring, cw + "_AdamsE2_generators", "x")
    result = str_basis(c_ring, cw + "_AdamsE2_basis", gen_names)

    c_ring.close()
    conn_ring.close()

    return result


############################### Relations ##########################################
def load_basis(c, table, gen_names):
    count = defaultdict(int)
    result = {}
    sql = f"SELECT mon, s, t FROM {table} where t-s>0 order by t-s, s, id"
    for str_mon, s, t in c.execute(sql):
        count[(s, t)] += 1
        result[latex_mon(str_mon, gen_names)] = count[(s, t)]
    return result


def latex_poly(str_mon_s: list[str], gen_names, basis, abbreviate_min=2):
    if len(str_mon_s) < abbreviate_min:
        result = "+".join(latex_mon(mon, gen_names) for mon in str_mon_s)
    else:
        result = "+".join(
            f"m_{{{basis[latex_mon(mon, gen_names)]}}}" for mon in str_mon_s
        )
    if result == "":
        result = "0"
    return result


def str_relations(c, table, gen_names, basis, remove_single=True, abbreviate_min=2):
    nRows = 52
    nColumns = 5
    nTables = 10
    table_head = R"""\begin{tabular}{r@{\hspace{3.3pt}}c@{\hspace{3.3pt}}lrr}\hline
LM & = & basis & stem & s \\\hline
"""
    sql = f"SELECT rel, s, t FROM {table} order by t-s, s"
    result = ""
    indices = {"row": 1, "table": 1}
    prev_s, prev_t = -1, -1
    color = 0
    for str_rel, s, t in c.execute(sql):
        str_mons = str_rel.split(";")
        if remove_single and len(str_mons) == 1:
            continue
        str_stem = str(t - s)
        str_s = str(s)
        if (s, t) != (prev_s, prev_t):
            color = 1 - color
            str_stem = f"\\hypertarget{{rel:{s}:{t}}}{{{t-s}}}"
            str_s = f"\\hyperlink{{basis:{s}:{t}}}{{{s}}}"
        row = f"${latex_mon(str_mons[0], gen_names)}$ & = & ${latex_poly(str_mons[1:], gen_names, basis, abbreviate_min)}$ & {str_stem} & {str_s} \\\\\n"
        result += new_row(
            indices,
            table_head,
            row,
            nRows,
            nTables,
            hLine=t - s != prev_t - prev_s,
            greyRow=color,
        )
        prev_s, prev_t = s, t

    result += empty_rows(indices["row"], nRows, nColumns)
    return result


def str_relations_monomials(c, table, gen_names):
    nRows = 52
    nColumns = 3
    nTables = 3
    table_head = R"""\begin{tabular}{lrr}\hline
LM & stem & s \\\hline
"""
    sql = f"SELECT rel, s, t FROM {table} order by t-s, s"
    result = ""
    indices = {"row": 1, "table": 1}
    prev_s, prev_t = -1, -1
    color = 0
    for str_rel, s, t in c.execute(sql):
        str_mons = str_rel.split(";")
        if len(str_mons) != 1:
            continue
        row = f"${latex_mon(str_mons[0], gen_names)}$ & {t-s} & {s} \\\\\n"
        if (s, t) != (prev_s, prev_t):
            color = 1 - color
        result += new_row(
            indices,
            table_head,
            row,
            nRows,
            nTables,
            hLine=t - s != prev_t - prev_s,
            greyRow=color,
        )
        prev_s, prev_t = s, t

    result += empty_rows(indices["row"], nRows, nColumns)
    return result


def str_relations_from_db(cw, path_ring, remove_single=True, abbreviate_min=2):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    gen_names = load_gen_names(c_ring, cw + "_AdamsE2_generators", "x")
    basis = load_basis(c_ring, cw + "_AdamsE2_basis", gen_names)
    result = "\\subsection{Relations with multiple summands}\n"
    result += str_relations(
        c_ring,
        cw + "_AdamsE2_relations",
        gen_names,
        basis,
        remove_single,
        abbreviate_min,
    )
    result += "\\subsection{Monomials that equal to zero}\n"
    result += str_relations_monomials(c_ring, cw + "_AdamsE2_relations", gen_names)

    c_ring.close()
    conn_ring.close()

    return result
