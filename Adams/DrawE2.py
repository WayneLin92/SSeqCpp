"""Generate html file for a spectral sequence"""
import sqlite3
import math
from collections import defaultdict
import argparse
import os
import subprocess
from turtle import st

path_tpl = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\index_tpl.html"
path_html = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\index.html"
path_js = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\scripts\data.js"

class Config:
    bullets_tilt_angle = -17 / 180 * math.pi
    bullets_radius = 0.08

cosA = math.cos(Config.bullets_tilt_angle)
sinA = math.sin(Config.bullets_tilt_angle)

def str2array(str_array: str):
    return tuple(int(i) for i in str_array.split(","))

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-i', default=R'C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\C2_AdamsSS.db', help='the database of the spectral sequence of a complex')
    parser.add_argument('--table', default='C2_AdamsE2', help='the database of the spectral sequence of a complex')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    conn = sqlite3.connect(args.i)
    c = conn.cursor()

    path_plot = args.i.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()

    path_S0 = R'C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\S0_AdamsSS_t249.db'
    conn_S0 = sqlite3.connect(path_S0)
    c_S0 = conn_S0.cursor()

    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()
    tpl_title = "Adams Spectral Sequence"
        
    # Load `basis`, `bullets'
    basis = []
    bullets = defaultdict(list)
    bullets_color = defaultdict(list)
    index = 0
    sql = f"SELECT mon, s, t FROM {args.table}_basis ORDER BY id"
    for str_mon, s, t in c.execute(sql):
        mon = list(map(int, str_mon.split(",")))
        basis.append(mon)
        bullets[(t - s, s)].append(index)
        if len(mon) == 1:
            bullets_color[(t - s, s)].append("darkorange")
        elif len(mon) == 3 and mon[-1] == 0 and mon[-2] == 1:
            bullets_color[(t - s, s)].append("blue")
        elif mon[-1] == 0:
            bullets_color[(t - s, s)].append("black")
        else:
            bullets_color[(t - s, s)].append("Maroon")
        index += 1

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
            str_fill = f' fill="{bullets_color[(x, y)][i]}"' if bullets_color[(x, y)][i] != "black" else ''
            page = 200
            level = 9800
            str_class = 'b'
            tpl_bullets += f'<circle id=b{index} class={str_class} cx={cx:.6g} cy={cy:.6g} r={r:.6g}{str_fill} data-b={index} data-l={level} data-d="{None}" data-i={0} data-j={0} data-g={int(len(str_fill) > 0)} data-page={page}>'
            tpl_bullets += f'<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
            index2xyrp[index] = (cx, cy, r, page)

    # Store gen_names
    sql = f"SELECT id, name, s FROM AdamsE2_generators order by id"
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
    sql = f"SELECT id, name, s FROM {args.table}_generators order by id"
    gens = defaultdict(int)
    for id, name, s in c.execute(sql):
        if name:
            name = name.replace("\\", "\\\\")
            content_js += f' "{name}",\n'
        else:
            content_js += f' "x_{{{s},{gens[s]}}}",\n'
        gens[s] += 1
    
    content_js += "];\n\n"

    # Store basis
    content_js += "const basis = [\n"
    for mon in basis:
        mon[-1] += offset_names
        mon.append(1)
        content_js += f" {mon},\n"
    content_js += f"];\n\n"

    # Store multiplicative structure in path_js
    arr_factors = {1, 2, 5, 13, 21, 28, 33, 42}
    lines = [[], [], []] # h0,h1,h2
    b2g = {1: 0, 2: 1, 5: 2} # basis_id to gen_id
    # b2g = {1: 0, 2: 1, 558: 2}
    sql = f"SELECT id1, id2, prod FROM {args.table}_basis_products"
    content_js += "const basis_prod = {\n"
    for id1, id2, prod in c_plot.execute(sql):
        if len(prod) > 0 and id1 in arr_factors:
            content_js += f' "{id1},{id2}": [{prod}],\n'

            if id1 in b2g:
                for id3 in map(int, prod.split(",")):
                    lines[b2g[id1]].append((id2, id3))
    content_js += "};\n\n"

    with open(path_js, "w") as file:
        file.write(content_js)
    
    # Plot the struct lines
    tpl_lines = ["", "", ""]
    for i in range(3):
        for i1, i2 in lines[i]:
            tpl_lines[i] += f'<line class=strt_l x1={index2xyrp[i1][0]:.6g} y1={index2xyrp[i1][1]:.6g} x2={index2xyrp[i2][0]:.6g} y2={index2xyrp[i2][1]:.6g} stroke-width={min(index2xyrp[i1][2], index2xyrp[i2][2]) / 4:.6g} data-page={min(index2xyrp[i1][3], index2xyrp[i2][3])} />'
            
    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("bullets:b8999a4e", tpl_bullets)

    content = content.replace("h0_lines:839e707e", tpl_lines[0])
    content = content.replace("h1_lines:2b707f82", tpl_lines[1])
    content = content.replace("h2_lines:e766a4f6", tpl_lines[2])

    content = content.replace("d2_lines:d9dec788", "")
    content = content.replace("d3_lines:1a07bf89", "")
    content = content.replace("d4_lines:7d935e97", "")
    content = content.replace("d5_lines:f6d7f992", "")
    content = content.replace("d6_lines:05c49c4a", "")
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)

    c.close()
    conn.close()
    c_S0.close()
    conn_S0.close()