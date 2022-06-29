"""doc"""
import sqlite3
import math
from collections import defaultdict

path_db = (
    R"C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\AdamsE2Export.db"
)
path_tpl = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\index_tpl.html"
path_html = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\index.html"
path_js = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss\AdamsSS_tmp\scripts\data.js"

class Config:
    bullets_tilt_angle = -20 / 180 * math.pi
    bullets_radius = 0.1

cosA = math.cos(Config.bullets_tilt_angle)
sinA = math.sin(Config.bullets_tilt_angle)

def str2mon(str_mon: str):
    list_mon = str_mon.split(",")
    it = iter(list_mon)
    return tuple((int(g), int(e)) for g, e in zip(it, it))

conn = sqlite3.connect(path_db)
with conn:
    c = conn.cursor()
    with open(path_tpl) as fp:
        content_tpl = fp.read()
    tpl_title = "Adams Spectral Sequence"
    bullets = defaultdict(list)
    bullets_ind = defaultdict(set)

    content_js = "const basis = [\n"
    sql = f"SELECT mon, s, t FROM AdamsE2_basis ORDER BY id"
    mon2index = {}
    index = 0
    for str_mon, s, t in c.execute(sql):
        mon = str2mon(str_mon)
        mon2index[mon] = index
        bullets[(t - s, s)].append(index)
        content_js += f" [{str_mon}],\n"
        if len(mon) == 1 and mon[0][1] == 1:
            bullets_ind[(t - s, s)].add(index)
        dict_mon = dict(mon)
        for i in range(3):
            dict_mon1 = dict_mon.copy()
            if i in dict_mon1:
                if dict_mon1[i] == 1:
                    del dict_mon1[i]
                else:
                    dict_mon1[i] -= 1
            mon1 = tuple(sorted(dict_mon1.items()))
        index += 1
    content_js += f"];\n\n"

    tpl_bullets = ""
    radius = defaultdict(lambda: 0.1)

    # Determine the maximum radius
    for x, y in bullets:
        len_bullets_xy = len(bullets[(x, y)])
        length_world = (len_bullets_xy - 1) * Config.bullets_radius * 3
        if (length_world > 1 - Config.bullets_radius * 3):
            length_world = 1 - Config.bullets_radius * 3
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
            radius_ub[(x + 1, y)] = min(radius_ub[(x + 1, y)], r)
            radius_ub[(x - 1, y)] = min(radius_ub[(x - 1, y)], r)
            radius_ub[(x, y + 1)] = min(radius_ub[(x, y + 1)], r)
            radius_ub[(x, y - 1)] = min(radius_ub[(x, y - 1)], r)
        
        b_changed = False
        for x, y in radius:
            if radius[(x, y)] > radius_ub[(x, y)] * 1.005:
                radius[(x, y)] = radius_ub[(x, y)]
                b_changed = True

    # Plot the bullets
    index2xyr = {}
    for x, y in bullets:
        r = radius[(x, y)]
        len_bullets_xy = len(bullets[(x, y)])
        sep_world = r * 3
        length_world = (len_bullets_xy - 1) * sep_world

        top_left_x = x - length_world / 2 * cosA
        top_left_y = y - length_world / 2 * sinA
        for i, index in enumerate(bullets[(x, y)]):
            cx = top_left_x + i * sep_world * cosA
            cy = top_left_y + i * sep_world * sinA
            if index in bullets_ind[(x, y)]:
                tpl_bullets += f'<circle id="b{index}" class="b" cx="{cx:.6g}" cy="{cy:.6g}" r="{r:.6g}" fill="blue"><title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
            else:
                tpl_bullets += f'<circle id="b{index}" class="b" cx="{cx:.6g}" cy="{cy:.6g}" r="{r:.6g}"><title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
            index2xyr[index] = (cx, cy, r)
    

    # Multiplicative structure
    lines = [[], [], []] # h0,h1,h2
    b2g = {1: 0, 2: 1, 5: 2}
    sql = f"SELECT id1, id2, prod FROM AdamsE2_basis_products"
    content_js += "const basis_prod = {\n"
    for id1, id2, prod in c.execute(sql):
        if len(prod) > 0:
            content_js += f' "{id1},{id2}": [{prod}],\n'

            if id1 in b2g:
                for id3 in map(int, prod.split(",")):
                    lines[b2g[id1]].append((id2, id3))
            if id2 in b2g:
                for id3 in map(int, prod.split(",")):
                    lines[b2g[id2]].append((id1, id3))
    content_js += "};\n\n"

    sql = f"SELECT id, name, s FROM AdamsE2_generators order by id"
    gens = defaultdict(int)
    content_js += "gen_names = [\n"
    for id, name, s in c.execute(sql):
        if name:
            content_js += f' "{name}",\n'
        else:
            content_js += f' "x_{{{s},{gens[s]}}}",\n'
        gens[s] += 1
    content_js += "];\n\n"

    with open(path_js, "w") as file:
        file.write(content_js)
    
    # Plot the lines
    tpl_lines = ["", "", ""]
    for i in range(3):
        for i1, i2 in lines[i]:
            tpl_lines[i] += f'<line x1="{index2xyr[i1][0]:.6g}" y1="{index2xyr[i1][1]:.6g}" x2="{index2xyr[i2][0]:.6g}" y2="{index2xyr[i2][1]:.6g}" stroke-width="{min(index2xyr[i1][2], index2xyr[i2][2]) / 5:.6g}"></line>\n'
    
    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("bullets:b8999a4e", tpl_bullets)

    content = content.replace("h0_lines:839e707e", tpl_lines[0])
    content = content.replace("h1_lines:2b707f82", tpl_lines[1])
    content = content.replace("h2_lines:e766a4f6", tpl_lines[2])
    with open(path_html, "w") as file:
        file.write(content)

    c.close()
conn.close()
