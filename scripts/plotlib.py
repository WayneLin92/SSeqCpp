import sqlite3
import math
from collections import defaultdict
import os

class Config:
    bullets_tilt_angle = -17 / 180 * math.pi
    bullets_radius = 0.08
cosA = math.cos(Config.bullets_tilt_angle)
sinA = math.sin(Config.bullets_tilt_angle)

path_tpl = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\index_tpl.html"
path_html = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\index.html"
path_js = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\data.js"

########################### Read #################################
def str2array(str_array: str):
    return [int(i) for i in str_array.split(",")] if len(str_array) > 0 else []

def get_complex_name(path):
    path = os.path.basename(path)
    names = ["S0", "C2", "Ceta", "Cnu", "Csigma"]
    for name in names:
        if path.startswith(name):
            return name
    raise ValueError(f"{path=} is not recognized")

def get_complex_t(path_mod):
    """ Return the dimension of the top cell """
    path_mod = os.path.basename(path_mod)
    names_to_t = {"C2" : 1, "Ceta" : 2, "Cnu" : 4, "Csigma" : 8}
    for name in names_to_t:
        if path_mod.startswith(name):
            return names_to_t[name]
    raise ValueError(f"{path_mod=} is not recognized")

def load_pi_basis(c, table, is_ring: bool):
    result = {"basis": [], "bullets": defaultdict(list), "bullets_ind": set(), "bullets_color": defaultdict(list)}
    index = 0
    sql = f"SELECT mon, s, t FROM {table} ORDER BY id"
    for str_mon, s, t in c.execute(sql):
        if is_ring:
            mon = list(map(int, str_mon.split(","))) if len(str_mon) > 0 else []
            mon1 = [[mon[i] // 2, mon[i + 1]] for i in range(0, len(mon) - 2, 2)]
            mon1.sort()
            mon = sum(mon1, [])
        else:
            mon = list(map(int, str_mon.split(",")))
            mon1 = [[mon[i] // 2, mon[i + 1]] for i in range(0, len(mon) - 3, 2)]
            mon1.sort()
            mon = sum(mon1, []) + [mon[-1], 1]

        result["basis"].append(mon)
        result["bullets"][(t - s, s)].append(index)
        if len(mon) == 2 and mon[1] == 1:
            color = "blue"
            result["bullets_ind"].add(index)
        else:
            color = "black"
        result["bullets_color"][(t - s, s)].append(color)
        index += 1
    return result

def load_pi_gen_names(c, table, letter):
    sql = f"SELECT id, name, t - s FROM {table} order by id"
    gens_stem = defaultdict(int)
    result = []
    for id, name, stem in c.execute(sql):
        if name:
            result.append(name)
        # elif gens_stem[stem] == 0:
        #     if stem < 10:
        #         result.append(f"{letter}_{stem}")
        #     else:
        #         result.append(f"{letter}_{{{stem}}}")
        # else:
        #     result.append(f"{letter}_{{{stem},{gens_stem[stem]}}}")
        else:
            if id < 10:
                result.append(f"{letter}_{id}")
            else:
                result.append(f"{letter}_{{{id}}}")
        gens_stem[stem] += 1
    return result

def load_pi_multiplications(c, table):
    sql = f"SELECT id1, id2, prod, O FROM {table}"
    result = []
    for id1, id2, prod, O in c.execute(sql):
        result.append((id1, id2, str2array(prod), O))
    return result

def load_pi_basis_map(c, table, target):
    sql = f"SELECT id, to_{target}, O_{target} FROM {table} WHERE to_{target} IS NOT NULL"
    result = []
    for id, target, O in c.execute(sql):
        result.append((id, str2array(target), O))
    return result

def load_ring(path_ring):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    path_plot = path_ring.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()

    complex_ring = get_complex_name(path_ring)

    is_ring = True
    result = load_pi_basis(c_ring, complex_ring + "_pi_basis", is_ring)
    result["gen_names"] = load_pi_gen_names(c_ring, complex_ring + "_pi_generators", "\\mu")
    result["products"] = load_pi_multiplications(c_plot, complex_ring + "_pi_basis_products")

    c_ring.close()
    conn_ring.close()
    c_plot.close()
    conn_plot.close()

    return result

def load_mod(path_ring, path_mod):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    conn_mod = sqlite3.connect(path_mod)
    c_mod = conn_mod.cursor()

    path_plot_ring = path_ring.replace(".db", "_plot.db")
    conn_plot_ring = sqlite3.connect(path_plot_ring)
    c_plot_ring = conn_plot_ring.cursor()

    path_plot_mod = path_mod.replace(".db", "_plot.db")
    conn_plot_mod = sqlite3.connect(path_plot_mod)
    c_plot_mod = conn_plot_mod.cursor()

    complex_ring = get_complex_name(path_ring)
    complex_mod = get_complex_name(path_mod)
    
    is_ring = True
    result_ring = load_pi_basis(c_ring, complex_ring + "_pi_basis", is_ring)
    result_ring["gen_names"] = load_pi_gen_names(c_ring, complex_ring + "_pi_generators", "\\mu")
    result_ring["products"] = load_pi_multiplications(c_plot_ring, complex_ring + "_pi_basis_products")
    result_ring["basis_map"] = load_pi_basis_map(c_plot_ring, complex_ring + "_pi_basis_maps", complex_mod)

    is_ring = False
    result_mod = load_pi_basis(c_mod, complex_mod + "_pi_basis", is_ring)
    result_mod["gen_names"] = load_pi_gen_names(c_mod, complex_mod + "_pi_generators", "\\iota")
    result_mod["products"] = load_pi_multiplications(c_plot_mod, complex_mod + "_pi_basis_products")
    result_mod["basis_map"] = load_pi_basis_map(c_plot_mod, complex_mod + "_pi_basis_maps", complex_ring)
    result_mod["t"] = get_complex_t(path_mod)

    c_ring.close()
    conn_ring.close()
    c_mod.close()
    conn_mod.close()
    c_plot_ring.close()
    conn_plot_ring.close()
    c_plot_mod.close()
    conn_plot_mod.close()

    return result_ring, result_mod

########################### Process #################################
def smoothen(radius):
    """Make the distribution of radius smooth """
    radius_ub = defaultdict(lambda: 0.1)  # upper bound of radius
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

def get_radius(data):
    """Set the radius of bullets in each coordinate (x, y) """
    radius = defaultdict(lambda: Config.bullets_radius)
    for x, y in data["bullets"]:
        len_bullets_xy = len(data["bullets"][(x, y)])
        length_world = (len_bullets_xy - 1) * Config.bullets_radius * 3
        if (length_world > 1 - Config.bullets_radius * 6):
            length_world = 1 - Config.bullets_radius * 6
        sep_world = length_world / (len_bullets_xy - 1) if len_bullets_xy > 1 else Config.bullets_radius * 3
        r = sep_world / 3

        if r < radius[(x, y)]:
            radius[(x, y)] = r
    smoothen(radius)
    return radius

########################### Write #################################
def export_bullets(data1d, radius1d):
    """ Input `data["bullets"]`, `data["bullets_color"]` and `data["basis"]` """
    tpl_bullets = ""
    index2xyrp = {}
    offset_basis = 0
    for i_data, (data, radius) in enumerate(zip(data1d, radius1d)):
        offset_x = i_data
        for x, y in data["bullets"]:
            r = radius[(x, y)]
            len_bullets_xy = len(data["bullets"][(x, y)])
            sep_world = r * 3
            length_world = (len_bullets_xy - 1) * sep_world

            bottom_right_x = x + length_world / 2 * cosA
            bottom_right_y = y + length_world / 2 * sinA
            for i, index in enumerate(data["bullets"][(x, y)]):
                cx = bottom_right_x - i * sep_world * cosA + offset_x
                cy = bottom_right_y - i * sep_world * sinA
                str_fill = f'fill="{data["bullets_color"][(x, y)][i]}" ' if data["bullets_color"][(x, y)][i] != "black" else ''
                page = 200
                level = 9800
                str_class = 'b'
                if len(data1d) > 0:
                    str_class += f' b{i_data}'
                tpl_bullets += f'<circle id=b{offset_basis + index} class="{str_class}" cx={cx:.6g} cy={cy:.6g} r={r:.6g} {str_fill}'\
                    f'data-b={index} data-l={level} data-d="{None}" data-i={offset_basis} data-j={0} data-g={int(len(str_fill) > 0)} data-page={page}>'\
                        f'<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n'
                index2xyrp[offset_basis + index] = (cx, cy, r, page)
        offset_basis += len(data["basis"])
    return tpl_bullets, index2xyrp

def export_basis_prod(data, index2xyrp, offset_basis=0):
    lines = [[], [], []]
    dashed_lines = [[], [], []]
    b2g = {1: 0, 3: 1, 7: 2}  # basis_id to gen_id
    for id_gen, id1, prod, O in data["products"]:
        if id_gen not in b2g:
            continue
        gen = b2g[id_gen]
        for id2 in prod:
            lines[gen].append((offset_basis + id1, offset_basis + id2))
        if O != -1:
            x1 = index2xyrp[offset_basis + id1][0] - (1 if offset_basis > 0 else 0)
            x2, y2 = round(x1) + 2 ** gen - 1, O
            if (x2, y2) in data["bullets"]:
                id2 = data["bullets"][(x2, y2)][0]
                dashed_lines[gen].append((offset_basis + id1, offset_basis + id2))
            else:
                dashed_lines[gen].append((offset_basis + id1, (x2 + (1 if offset_basis > 0 else 0), y2)))
    return lines, dashed_lines

def export_basis_map(data1, data2, index2xyrp, offset1_basis, offset2_basis, dx):
    lines = []
    dashed_lines = []
    for id1, target, O in data1["basis_map"]:
        for id2 in target:
            lines.append((offset1_basis + id1, offset2_basis + id2))
        if O != -1:
            x1 = index2xyrp[offset1_basis + id1][0] - (1 if offset1_basis > 0 else 0)
            x2, y2 = round(x1) + dx, O
            if (x2, y2) in data2["bullets"]:
                id2 = data2["bullets"][(x2, y2)][0]
                dashed_lines.append((offset1_basis + id1, offset2_basis + id2))
            else:
                dashed_lines.append((offset1_basis + id1, (x2 + (1 if offset2_basis > 0 else 0), y2)))
    return lines, dashed_lines

def element_line(x1, y1, x2, y2, *, w, p: int, class_: str, color: str = None, dashed=False):
    attr_more = ""
    if color: 
        attr_more += f' stroke="{color}"'
    if dashed:
        attr_more += ' stroke-dasharray="0.2,0.2"'


    range_x = f'data-x1={round(min(x1, x2))} data-x2={round(max(x1, x2))}'
    if round(y2) - round(y1) <= 1:
        return f'<line class="{class_}" x1={x1:.6g} y1={y1:.6g} x2={x2:.6g} y2={y2:.6g} {range_x} stroke-width={w:.6g} data-page={p}{attr_more} />'
    else:
        if round(x1) == round(x2):
            return f'<path class="{class_}" d="M {x1:.6g} {y1:.6g} C {x1 + 0.3:.6g} {y1 * 0.7 + y2 * 0.3:.6g}, {x2 + 0.3:.6g} {y1 * 0.3 + y2 * 0.7:.6g}, {x2:.6g} {y2:.6g}" '\
                f'{range_x} stroke-width={w:.6g}{attr_more} />\n'
        else:
            x_diff = abs(round(x2) - round(x1))
            sign = 1 if x2 > x1 else -1
            return f'<path class="{class_}" d="M {x1:.6g} {y1:.6g} C {x1 + 0.5 * sign:.6g} {y1 + 0.5 / x_diff:.6g}, {x2 - x_diff * sign / 4 / (y2 - y1):.6g} {y2 - 0.5:.6g}, {x2:.6g} {y2:.6g}" '\
                f'{range_x} stroke-width={w:.6g}{attr_more} />\n'

def export_h_lines(index2xyrp, lines, dashed_lines, class_: str = None):
    """ Plot the h0, h1, h2 lines """
    tpl_lines = ["", "", ""]
    extend_colors = ["red", "blue", "green"]
    if class_ is None:
        class_ = ""
    else:
        class_ = " " + class_
    for i in range(3):
        for i1, i2 in lines[i]:
            x1, y1, r1, p1 = index2xyrp[i1]
            x2, y2, r2, p2 = index2xyrp[i2]
            if(round(x2) < round(x1)): #########
                print(x1, y1, x2, y2, i1, i2, class_)
                raise ValueError("incorrect line")
            color = None if round(y2) - round(y1) == 1 else extend_colors[i]
            tpl_lines[i] += element_line(x1, y1, x2, y2, w=min(r1, r2) / 4, p=min(p1, p2), class_=f"strt_l{class_}", dashed=False, color=color)
    for i in range(3):
        for i1, O in dashed_lines[i]:
            x1, y1, r1, p1 = index2xyrp[i1]
            if type(O) is tuple:
                x2, y2, r2, p2 = *O, r1, p1
            else:
                x2, y2, r2, p2 = index2xyrp[O]
            tpl_lines[i] += element_line(x1, y1, x2, y2, w=min(r1, r2) / 4, p=min(p1, p2), class_=f"strt_l dashed_l{class_}", dashed=True, color=extend_colors[i])
    return tpl_lines

def export_lines(index2xyrp, lines, dashed_lines, class_: str = None):
    """ Plot lines """
    tpl_lines = ""
    if class_ is None:
        class_ = ""
    else:
        class_ = " " + class_
    for i1, i2 in lines:
        x1, y1, r1, p1 = index2xyrp[i1]
        x2, y2, r2, p2 = index2xyrp[i2]
        tpl_lines += element_line(x1, y1, x2, y2, w=min(r1, r2) / 4, p=min(p1, p2), class_=f"strt_l{class_}", dashed=False)
    for i1, O in dashed_lines:
        x1, y1, r1, p1 = index2xyrp[i1]
        if type(O) is tuple:
            x2, y2, r2, p2 = *O, r1, p1
        else:
            x2, y2, r2, p2 = index2xyrp[O]
        tpl_lines += element_line(x1, y1, x2, y2, w=min(r1, r2) / 4, p=min(p1, p2), class_=f"strt_l dashed_l{class_}", dashed=True)
    return tpl_lines

def export_gen_names_to_js(data1d):
    content_js = "const gen_names = [\n"
    for data in data1d:
        for name in data["gen_names"]:
            name1 = name.replace('\\', '\\\\')
            content_js += f' "{name1}",\n'
        content_js += '\n'
    content_js += "];\n\n"
    return content_js

def export_basis_to_js(data1d):
    content_js = "const basis = [\n"
    offset = 0
    for data in data1d:
        for mon in data["basis"]:
            if offset > 0:
                mon[-2] += offset
            content_js += f" {mon},\n"
        content_js += '\n'
        offset += len(data["gen_names"])
    content_js += f"];\n\n"
    return content_js

def export_basis_prod_to_js(data1d):
    content_js = "const basis_prod = {\n"
    offset = 0
    for data in data1d:
        for id1, id2, prod, _ in data["products"]:
            if len(prod) > 0:
                prod = [i + offset for i in prod]
                content_js += f' "{id1},{id2 + offset}": {prod},\n'
        content_js += '\n'
        offset += len(data["basis"])
    content_js += "};\n\n"
    return content_js

def export_pi_ring(data):
    tpl_title = "Adams Spectral Sequence"
    
    # js
    content_js = ""
    content_js += export_gen_names_to_js([data])
    content_js += export_basis_to_js([data])
    content_js += export_basis_prod_to_js([data])

    # html
    radius = get_radius(data)
    tpl_bullets, index2xyrp = export_bullets([data], [radius])
    lines, dashed_lines = export_basis_prod(data, index2xyrp, 0)
    tpl_lines = export_h_lines(index2xyrp, lines, dashed_lines)
    
    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("<!-- bullets:b8999a4e -->", tpl_bullets)

    content = content.replace("<!-- h0_lines:839e707e -->", tpl_lines[0])
    content = content.replace("<!-- h1_lines:2b707f82 -->", tpl_lines[1])
    content = content.replace("<!-- h2_lines:e766a4f6 -->", tpl_lines[2])

    with open(path_js, "w") as file:
        file.write(content_js)
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)

def export_pi_mod(data_ring, data_mod):
    tpl_title = "Adams Spectral Sequence"
    
    # js
    content_js = 'const MODE="DualSS";\n'
    content_js += f'const sep_max_width={data_mod["t"] + 1};\n'
    content_js += 'var sep_right=10;\n'
    content_js += 'var sep_width=1;\n'
    content_js += export_gen_names_to_js([data_ring, data_mod])
    content_js += export_basis_to_js([data_ring, data_mod])
    content_js += export_basis_prod_to_js([data_ring, data_mod])

    # html
    radius = get_radius(data_ring)
    radius1 = get_radius(data_mod)
    tpl_bullets, index2xyrp = export_bullets([data_ring, data_mod], [radius, radius1])
    lines, dashed_lines = export_basis_prod(data_ring, index2xyrp, 0)
    lines1, dashed_lines1 = export_basis_prod(data_mod, index2xyrp, len(data_ring["basis"]))
    lines_bc, dashed_lines_bc = export_basis_map(data_ring, data_mod, index2xyrp, 0, len(data_ring["basis"]), 0)
    lines_tc, dashed_lines_tc = export_basis_map(data_mod, data_ring, index2xyrp, len(data_ring["basis"]), 0, -data_mod["t"])
    tpl_lines = export_h_lines(index2xyrp, lines, dashed_lines, "l0")
    tpl_lines1 = export_h_lines(index2xyrp, lines1, dashed_lines1, "l1")
    tpl_lines_bc = export_lines(index2xyrp, lines_bc, dashed_lines_bc, "lbc")
    tpl_lines_tc = export_lines(index2xyrp, lines_tc, dashed_lines_tc, "ltc")
    for i in range(3):
        tpl_lines[i] += tpl_lines1[i]
    
    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("<!-- bullets:b8999a4e -->", tpl_bullets)

    content = content.replace("<!-- h0_lines:839e707e -->", tpl_lines[0])
    content = content.replace("<!-- h1_lines:2b707f82 -->", tpl_lines[1])
    content = content.replace("<!-- h2_lines:e766a4f6 -->", tpl_lines[2])

    content = content.replace("<!-- g_maroon_lines:7aa92c76 -->", tpl_lines_bc)
    content = content.replace("<!-- g_purple_lines:114fd993 -->", tpl_lines_tc)

    with open(path_js, "w") as file:
        file.write(content_js)
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)