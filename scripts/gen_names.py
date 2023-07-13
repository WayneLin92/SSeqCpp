"""doc"""
import json
import sqlite3
import os
from collections import defaultdict

PATH_JSON = "gen_names.json"


def get_complex_name(path):
    path = os.path.basename(path)
    names = ["S0", "C2", "Ceta", "Cnu", "Csigma"]
    for name in names:
        if path.startswith(name):
            return name
    raise ValueError(f"{path=} is not recognized")


def tex_outside_delimiter(text: str, symbol: str):
    """Return if symbol appears in text and is outside any pair of delimiters including ()[]{}"""
    left, right = "([{", ")]}"
    left_minus_right = 0
    for c in text:
        if c in left:
            left_minus_right += 1
        elif c in right:
            left_minus_right -= 1
        elif c == symbol and left_minus_right == 0:
            return True
    return False


def tex_pow(base, exp: int) -> str:
    """Return base^exp in latex."""
    if type(base) != str:
        base = str(base)
    if exp == 1:
        return base
    else:
        if tex_outside_delimiter(base, "^"):
            base = "(" + base + ")"
        return f"{base}^{exp}" if len(str(exp)) == 1 else f"{base}^{{{exp}}}"


def get_mon_name(m, gen_names):
    it = map(int, m.split(","))
    return "".join(tex_pow(gen_names[i], e) for i, e in zip(it, it))


def get_poly_name(p, gen_names):
    coeff = "+".join(get_mon_name(m, gen_names) for m in p.split(";"))
    if tex_outside_delimiter(coeff, "+"):
        coeff = "(" + coeff + ")"
    return "(" + coeff + R"i_1)"


def set_AdamsE2_name(cw: str, path: str):
    with open(PATH_JSON) as fp:
        gen_names = json.load(fp)
    gen_names = {int(k): v for k, v in gen_names[cw].items()}

    conn = sqlite3.connect(path)
    c = conn.cursor()

    sql = f"select id, s, t from {cw}_AdamsE2_generators"
    names = defaultdict(int)
    for id, s, t in c.execute(sql):
        if id not in gen_names:
            if names[(t - s, s)] == 0:
                gen_names[id] = f"x_{{{t-s},{s}}}"
            else:
                gen_names[id] = f"x_{{{t-s},{s},{names[(t-s, s)]}}}"
        names[(t - s, s)] += 1

    sql = f"Update {cw}_AdamsE2_generators SET name=?2 where id=?1"
    c.executemany(sql, gen_names.items())

    c.close()
    conn.commit()
    conn.close()


def set_E2_name_two_cell(path):
    gen_names_C = {}
    complex = get_complex_name(path)
    sql = f"select id, to_S0 from {complex}_AdamsE2_generators"
    for id, toS0 in c.execute(sql):
        if id == 0:
            gen_names_C[id] = R"i_0"
        else:
            gen_names_C[id] = get_poly_name(toS0, gen_names)

    sql = f"Update {complex}_AdamsE2_generators SET name=?2 where id=?1"
    c.executemany(sql, gen_names_C.items())


def set_pi_name(path_S0):
    with open("pi_gen_names.json") as fp:
        pi_gen_names = json.load(fp)
    pi_gen_names = {
        (int((s := k.split(","))[0]), int(s[1])): v for k, v in pi_gen_names.items()
    }

    conn_S0 = sqlite3.connect(
        os.path.join(
            R"C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release", path_S0
        )
    )
    c_S0 = conn_S0.cursor()

    sql = f"select id, s, t from S0_pi_generators"
    name_i = defaultdict(int)
    id_to_name = {}
    for id, s, t in c_S0.execute(sql):
        if (t - s, s) in pi_gen_names:
            id_to_name[id] = pi_gen_names[(t - s, s)]
        else:
            if name_i[t - s] > 0:
                id_to_name[id] = f"\\rho_{{{t-s}, {name_i[t-s]}}}"
            else:
                id_to_name[id] = f"\\rho_{{{t-s}}}"
        name_i[t - s] += 1

    sql = f"Update S0_pi_generators SET name=?2 where id=?1"
    c_S0.executemany(sql, id_to_name.items())

    c_S0.close()
    conn_S0.commit()
    conn_S0.close()


def reset_pi_name(path_S0):
    conn_S0 = sqlite3.connect(
        os.path.join(
            R"C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release", path_S0
        )
    )
    c_S0 = conn_S0.cursor()

    sql = f"Update S0_pi_generators SET name=NULL"
    c_S0.execute(sql)

    c_S0.close()
    conn_S0.commit()
    conn_S0.close()
