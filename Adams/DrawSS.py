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

def str2mon(str_mon: str):
    list_mon = str_mon.split(",")
    it = iter(list_mon)
    return tuple((int(g), int(e)) for g, e in zip(it, it))

def str2array(str_array: str):
    return tuple(int(i) for i in str_array.split(","))

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-i', default=R'C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\S0_AdamsSS_t245.db', help='the database of the spectral sequence')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions

    conn = sqlite3.connect(args.i)
    with conn:
        c = conn.cursor()
        with open(path_tpl, encoding="utf-8") as fp:
            content_tpl = fp.read()
        tpl_title = "Adams Spectral Sequence"
        
        # Load `basis`, `bullets` and `bullets_ind`
        basis = []
        bullets = defaultdict(list)
        bullets_ind = defaultdict(set)
        index = 0
        sql = f"SELECT mon, s, t FROM AdamsE2_basis ORDER BY id"
        for str_mon, s, t in c.execute(sql):
            mon = str2mon(str_mon)
            basis.append(str_mon)
            bullets[(t - s, s)].append(index)
            if len(mon) == 1 and mon[0][1] == 1:
                bullets_ind[(t - s, s)].add(index)
            index += 1

        # Load `ids`
        ids = defaultdict(int)
        sql = "SELECT s, t, min(id) FROM AdamsE2_ss GROUP BY t, s;"
        for s, t, id in c.execute(sql):
            ids[(s, t)] = id

        # Load `AdamsE2_ss`
        basis_ss = []
        primitive_candidates = {}
        sql = f"SELECT base, level, diff, s, t, id FROM AdamsE2_ss ORDER BY id"
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
        if c.execute("SELECT name FROM sqlite_master WHERE name='AdamsE2_ss_primitives'").fetchone() is not None:
            sql = f"SELECT base, level, diff, s, t FROM AdamsE2_ss_primitives WHERE t-s<=110 ORDER BY t-s,s"
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
        sql = "SELECT s, t, l FROM AdamsE2_ss_stable_levels GROUP BY t, s;"
        for s, t, l in c.execute(sql):
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
                base, level, diff, offset, offsetD = basis_ss[index]
                str_fill = ' fill="blue"' if "," not in base and (offset + int(base)) in bullets_ind[(x, y)] else ''
                if level >= stable_levels[(x, y)]:
                    if diff:
                        page = 10000 - level
                    else:
                        page = 200
                elif level >= 5000:
                    page = 200
                else:
                    page = level
                str_class = '"b prim_b"' if (base, level, offset) in primitive_bullets else 'b'
                tpl_bullets += f'<circle id=b{index} class={str_class} cx={cx:.6g} cy={cy:.6g} r={r:.6g}{str_fill} data-b={base} data-l={level} data-d="{diff}" data-i={offset} data-j={offsetD} data-g={int(len(str_fill) > 0)} data-page={page}>'
                tpl_bullets += f'<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
                index2xyrp[index] = (cx, cy, r, page)

        # Store gen_names
        sql = f"SELECT id, name, s FROM AdamsE2_generators order by id"
        gens = defaultdict(int)
        content_js = "gen_names = [\n"
        for id, name, s in c.execute(sql):
            if name:
                name = name.replace("\\", "\\\\")
                content_js += f' "{name}",\n'
            else:
                content_js += f' "x_{{{gens[s]},{s}}}",\n'
            gens[s] += 1
        content_js += "];\n\n"

        # Store basis
        content_js += "const basis = [\n"
        for str_mon in basis:
            content_js += f" [{str_mon}],\n"
        content_js += f"];\n\n"

        # Store multiplicative structure in path_js
        arr_factors = {1, 2, 5, 13, 21, 28, 33, 42}
        lines = [[], [], []] # h0,h1,h2
        b2g = {1: 0, 2: 1, 5: 2} # basis_id to gen_id
        # b2g = {1: 0, 2: 1, 558: 2}
        sql = f"SELECT id1, id2, prod FROM AdamsE2_ss_products"
        content_js += "const basis_prod = {\n"
        for id1, id2, prod in c.execute(sql):
            if len(prod) > 0 and ({id1, id2} & arr_factors):
                content_js += f' "{id1},{id2}": [{prod}],\n'

                if id1 in b2g:
                    for id3 in map(int, prod.split(",")):
                        lines[b2g[id1]].append((id2, id3))
                if id2 in b2g:
                    for id3 in map(int, prod.split(",")):
                        lines[b2g[id2]].append((id1, id3))
        content_js += "};\n\n"

        diff_lines = [[], [], [], [], []] # d2,d3,d4,d5,d6
        sql = f"SELECT src, r, tgt FROM AdamsE2_ss_diffs"
        for src, r, tgt in c.execute(sql):
            index = min(r - 2, 4)
            for id_tgt in map(int, tgt.split(",")):
                diff_lines[index].append((src, id_tgt, ))

        null_diff_lines = [[], [], [], [], []] # d2,d3,d4,d5,d6
        sql = f"SELECT src, r, tgt FROM AdamsE2_ss_nd"
        for src, r, tgt in c.execute(sql):
            index = min(r - 2, 4)
            for id_tgt in map(int, tgt.split(",")):
                null_diff_lines[index].append((src, id_tgt))
                break

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
