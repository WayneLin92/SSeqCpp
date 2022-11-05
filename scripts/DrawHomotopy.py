"""Generate html file for a spectral sequence"""
import sqlite3
import math
from collections import defaultdict
import argparse
import os
import subprocess

path_tpl = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\index_tpl.html"
path_html = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\index.html"
path_js = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\data.js"

class Config:
    bullets_tilt_angle = -17 / 180 * math.pi
    bullets_radius = 0.08

cosA = math.cos(Config.bullets_tilt_angle)
sinA = math.sin(Config.bullets_tilt_angle)

def str2mon(str_mon: str):
    list_mon = str_mon.split(",")
    it = iter(list_mon)
    return tuple((int(g), int(e)) for g, e in zip(it, it))

def str2array(str_array: str):
    return tuple(int(i) for i in str_array.split(","))

def get_complex_name(db):
    db = os.path.basename(db)
    if (db[1] == '0'):
        return "S0"
    elif (db[1] == '2'):
        return "C2"
    elif (db[1] == 'e'):
        return "Ceta"
    elif (db[1] == 'n'):
        return "Cnu"
    elif (db[1] == 's'):
        return "Csigma"
    else:
        raise ValueError(f"{db=} is not recognized")

def smoothen(radius):
    radius_ub = defaultdict(lambda: 0.1)
    b_changed = True
    while b_changed:
        for x, y in radius:
            r = radius[(x, y)] * 1.085
            r1 = radius[(x, y)] * 1.13
            radius_ub[(x + 1, y)] = min(radius_ub[(x + 1, y)], r)
            radius_ub[(x - 1, y)] = min(radius_ub[(x - 1, y)], r)
            radius_ub[(x, y + 1)] = min(radius_ub[(x, y + 1)], r)
            radius_ub[(x, y - 1)] = min(radius_ub[(x, y - 1)], r)
            radius_ub[(x + 1, y + 1)] = min(radius_ub[(x + 1, y + 1)], r1)
            radius_ub[(x + 1, y - 1)] = min(radius_ub[(x + 1, y - 1)], r1)
            radius_ub[(x - 1, y + 1)] = min(radius_ub[(x - 1, y + 1)], r1)
            radius_ub[(x - 1, y - 1)] = min(radius_ub[(x - 1, y - 1)], r1)
        
        b_changed = False
        for x, y in radius:
            if radius[(x, y)] > radius_ub[(x, y)] * 1.005:
                radius[(x, y)] = radius_ub[(x, y)]
                b_changed = True

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-i', default=R'C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\benchmark\S0_AdamsSS_t100.db', help='the database of the spectral sequence')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    table = get_complex_name(args.i)
    conn = sqlite3.connect(args.i)
    c = conn.cursor()
    
    path_plot = args.i.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()
    
    path_S0 = R'C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\benchmark\S0_AdamsSS_t100.db'
    conn_S0 = sqlite3.connect(path_S0)
    c_S0 = conn_S0.cursor()

    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()
    tpl_title = "Adams Spectral Sequence"
    
    # Load `basis`, `bullets` and `bullets_ind`
    basis = []
    bullets = defaultdict(list)
    bullets_ind = defaultdict(set)
    bullets_from_S0 = defaultdict(set)
    index = 0
    sql = f"SELECT mon, s, t FROM {table}_pi_basis ORDER BY id"
    for str_mon, s, t in c.execute(sql):
        mon = list(map(int, str_mon.split(",") if len(str_mon) > 0 else []))
        basis.append(mon)
        bullets[(t - s, s)].append(index)
        if table.startswith("S0"):
            if len(mon) == 2 and mon[1] == 1:
                bullets_ind[(t - s, s)].add(index)
        else:
            if len(mon) == 1:
                bullets_ind[(t - s, s)].add(index)
            if mon[-1] == 0:
                bullets_from_S0[(t - s, s)].add(index)
        index += 1

    # Load `ids`
    ids = defaultdict(int)
    sql = f"SELECT s, t, min(id) FROM {table}_ss GROUP BY t, s;"
    for s, t, id in c.execute(sql):
        ids[(s, t)] = id

    # Load `AdamsE2_ss`
    basis_ss = []
    primitive_candidates = {}
    sql = f"SELECT base, level, diff, s, t, id FROM {table}_ss ORDER BY id"
    for base, level, diff, s, t, id in c.execute(sql):
        if level > 5000:
            r = 10000 - level
            s1, t1 = s + r, t + r - 1
        else:
            s1, t1 = s - level, t - level + 1
        basis_ss.append((base, level, diff, ids[(s, t)], ids[(s1, t1)]))
        primitive_candidates[(base, level, ids[(s, t)])] = id
        
    # Load `AdamsE2_ss_primitives`
    primitive_bullets = set()  # set of (base, level, ids[(s,t)])
    primitive_diffs = set()  # set of (src, r)
    if c.execute(f"SELECT name FROM sqlite_master WHERE name='{table}_ss_primitives'").fetchone() is not None:
        sql = f"SELECT base, level, diff, s, t FROM {table}_ss_primitives WHERE t-s<=110 ORDER BY t-s,s"
        for base, level, diff, s, t in c.execute(sql):
            if diff is None:
                primitive_bullets.add((base, level, ids[(s, t)]))
            else:
                k = (base, level, ids[(s, t)])
                if k in primitive_candidates:
                    primitive_diffs.add(primitive_candidates[k])
                else:
                    print(f"Migration input ({base=}, {level=}, {diff=}) changed!")
    else:
        print("ignoring AdamsE2_ss_primitives because it is not found")

    # Load `stable_levels`
    stable_levels = {}
    sql = f"SELECT s, t, l FROM {table}_ss_stable_levels GROUP BY t, s;"
    for s, t, l in c_plot.execute(sql):
        stable_levels[(t - s, s)] = l

    # Determine the maximum radius
    radius = defaultdict(lambda: Config.bullets_radius)
    for x, y in bullets:
        len_bullets_xy = len(bullets[(x, y)])
        length_world = (len_bullets_xy - 1) * Config.bullets_radius * 3
        if (length_world > 1 - Config.bullets_radius * 6):
            length_world = 1 - Config.bullets_radius * 6
        sep_world = length_world / (len_bullets_xy - 1) if len_bullets_xy > 1 else Config.bullets_radius * 3
        r = sep_world / 3

        if r < radius[(x, y)]:
            radius[(x, y)] = r
    
    # Make the radius smooth
    smoothen(radius)

    # Plot the bullets. Compute `index2xyr` for lines
    tpl_bullets = ""
    index2xyrp = {}
    for x, y in bullets:
        r = radius[(x, y)]
        len_bullets_xy = len(bullets[(x, y)])
        sep_world = r * 3
        length_world = (len_bullets_xy - 1) * sep_world

        bottom_right_x = x + length_world / 2 * cosA
        bottom_right_y = y + length_world / 2 * sinA
        for i, index in enumerate(bullets[(x, y)]):
            cx = bottom_right_x - i * sep_world * cosA
            cy = bottom_right_y - i * sep_world * sinA
            base, level, diff, offset, offsetD = basis_ss[index]
            list_base = list(base.split(","))
            if table.startswith("S0"):
                str_fill = ' fill="blue"' if (offset + int(base.split(",")[0])) in bullets_ind[(x, y)] else ''
            else:
                index_first_bullet = offset + int(base.split(",")[0])
                if index_first_bullet in bullets_ind[(x, y)]:
                    str_fill = ' fill="Maroon"'
                elif index_first_bullet in bullets_from_S0[(x, y)]:
                    str_fill = ''
                else:
                    str_fill = ' fill="purple"'
            if level >= stable_levels[(x, y)]:
                if diff:
                    page = 10000 - level
                else:
                    page = 200
            elif level >= 5000:
                page = 200
            elif diff:
                page = level
            else:
                page = 200
            str_class = '"b prim_b"' if (base, level, offset) in primitive_bullets else 'b'
            tpl_bullets += f'<circle id=b{index} class={str_class} cx={cx:.6g} cy={cy:.6g} r={r:.6g}{str_fill} data-b={base} data-l={level} data-d="{diff}" data-i={offset} data-j={offsetD} data-g={int(len(str_fill) > 0)} data-page={page}>'
            tpl_bullets += f'<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
            index2xyrp[index] = (cx, cy, r, page)

    # Store gen_names
    sql = f"SELECT id, name, s FROM S0_AdamsE2_generators order by id"
    gens = defaultdict(int)
    content_js = "gen_names = [\n"
    offset_names = 0
    for id, name, s in c_S0.execute(sql):
        if name:
            name = name.replace("\\", "\\\\")
            content_js += f' "{name}",\n'
        else:
            content_js += f' "a_{{{s},{gens[s]}}}",\n'
        gens[s] += 1
        offset_names += 1
    if not table.startswith("S0"):
        sql = f"SELECT id, name, s FROM {table}_generators order by id"
        gens = defaultdict(int)
        for id, name, s in c.execute(sql):
            if name:
                name = name.replace("\\", "\\\\")
                content_js += f' "{name}",\n'
            else:
                content_js += f' "x_{{{s},{gens[s]}}}",\n'
    content_js += "];\n\n"

    # Store basis
    content_js += "const basis = [\n"
    if table.startswith("S0"):
        for mon in basis:
            content_js += f" {mon},\n"
    else:
        for mon in basis:
            mon[-1] += offset_names
            mon.append(1)
            content_js += f" {mon},\n"
    content_js += f"];\n\n"

    # Store multiplicative structure in path_js
    arr_factors = {1, 3, 7, 15, 23, 29, 33, 42}
    lines = [[], [], []] # h0,h1,h2
    b2g = {1: 0, 3: 1, 7: 2} # basis_id to gen_id
    sql = f"SELECT id1, id2, prod FROM {table}_ss_products"
    content_js += "const basis_prod = {\n"
    for id1, id2, prod in c_plot.execute(sql):
        if len(prod) > 0 and id1 in arr_factors:
            content_js += f' "{id1},{id2}": [{prod}],\n'

            if id1 in b2g:
                for id3 in map(int, prod.split(",")):
                    lines[b2g[id1]].append((id2, id3))
    content_js += "};\n\n"

    diff_lines = [[], [], [], [], []] # d2,d3,d4,d5,d6
    sql = f"SELECT src, r, tgt FROM {table}_ss_diffs"
    for src, r, tgt in c_plot.execute(sql):
        index = min(r - 2, 4)
        for id_tgt in map(int, tgt.split(",")):
            diff_lines[index].append((src, id_tgt, ))

    null_diff_lines = [[], [], [], [], []] # d2,d3,d4,d5,d6
    sql = f"SELECT src, r, tgt FROM {table}_ss_nd"
    for src, r, tgt in c_plot.execute(sql):
        index = min(r - 2, 4)
        if len(tgt) > 0:
            null_diff_lines[index].append((src, int(tgt.split(",")[-1])))

    with open(path_js, "w") as file:
        file.write(content_js)
    
    # Plot the struct lines
    tpl_lines = ["", "", ""]
    for i in range(3):
        for i1, i2 in lines[i]:
            tpl_lines[i] += f'<line class=strt_l x1={index2xyrp[i1][0]:.6g} y1={index2xyrp[i1][1]:.6g} x2={index2xyrp[i2][0]:.6g} y2={index2xyrp[i2][1]:.6g} stroke-width={min(index2xyrp[i1][2], index2xyrp[i2][2]) / 4:.6g} data-page={min(index2xyrp[i1][3], index2xyrp[i2][3])} />'
            
    # Plot the diff lines
    tpl_diff_lines = ["", "", "", "", ""]
    for i in range(5):
        for i1, i2 in diff_lines[i]:
            str_class = '"diff_l prim_l"' if i1 in primitive_diffs else 'diff_l'
            tpl_diff_lines[i] += f'<line class={str_class} x1={index2xyrp[i1][0]:.6g} y1={index2xyrp[i1][1]:.6g} x2={index2xyrp[i2][0]:.6g} y2={index2xyrp[i2][1]:.6g} stroke-width={min(index2xyrp[i1][2], index2xyrp[i2][2]) / 4:.6g} data-page={min(index2xyrp[i1][3], index2xyrp[i2][3])} />\n'

    # Plot the null diff lines
    for i in range(5):
        for i1, i2 in null_diff_lines[i]:
            if round(index2xyrp[i1][0]) <= 127:
                tpl_diff_lines[i] += f'<line class=dashed_l x1={index2xyrp[i1][0]:.6g} y1={index2xyrp[i1][1]:.6g} x2={index2xyrp[i2][0]:.6g} y2={index2xyrp[i2][1]:.6g} stroke-width={min(index2xyrp[i1][2], index2xyrp[i2][2]) / 4:.6g} stroke-dasharray="0.2,0.2" data-r={round(index2xyrp[i2][1]) - round(index2xyrp[i1][1])} />\n'

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("bullets:b8999a4e", tpl_bullets)

    content = content.replace("h0_lines:839e707e", tpl_lines[0])
    content = content.replace("h1_lines:2b707f82", tpl_lines[1])
    content = content.replace("h2_lines:e766a4f6", tpl_lines[2])

    content = content.replace("d2_lines:d9dec788", tpl_diff_lines[0])
    content = content.replace("d3_lines:1a07bf89", tpl_diff_lines[1])
    content = content.replace("d4_lines:7d935e97", tpl_diff_lines[2])
    content = content.replace("d5_lines:f6d7f992", tpl_diff_lines[3])
    content = content.replace("d6_lines:05c49c4a", tpl_diff_lines[4])
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)

    c.close()
    conn.close()
    c_S0.close()
    conn_S0.close()
    c_plot.close()
    conn_plot.close()
