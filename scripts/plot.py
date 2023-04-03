import sqlite3
import math
from collections import defaultdict
import os
import copy


class Config:
    bullets_tilt_angle = -17 / 180 * math.pi
    bullets_radius = 0.08


cosA = math.cos(Config.bullets_tilt_angle)
sinA = math.sin(Config.bullets_tilt_angle)

path_tpl = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\index_tpl.html"
path_html_tmp = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\index.html"
path_js_tmp = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\ss-fb42729d\AdamsSS_tmp\data.js"

########################### Read #################################
def has_table(c, table_name):
    if (
        c.execute(
            f"SELECT name FROM sqlite_master WHERE name='{table_name}_ss_primitives'"
        ).fetchone()
        is not None
    ):
        return True
    else:
        return False


def has_column(c, table_name, column_name):
    if (
        c.execute(
            f"SELECT * FROM pragma_table_info('{table_name}') WHERE name='{column_name}'"
        ).fetchone()
        is not None
    ):
        return True
    else:
        return False


def str2array(str_array: str):
    return [int(i) for i in str_array.split(",")] if len(str_array) > 0 else []


def get_complex_name(path):
    path = os.path.basename(path)
    names = ["S0", "C2", "Ceta", "Cnu", "Csigma", "RP10", "RPinf", "X2"]
    for name in names:
        if path.startswith(name):
            return name
    raise ValueError(f"{path=} is not recognized")


def get_complex_t(path_mod):
    """Return the dimension of the top cell"""
    path_mod = os.path.basename(path_mod)
    names_to_t = {"C2": 1, "Ceta": 2, "Cnu": 4, "Csigma": 8}
    for name in names_to_t:
        if path_mod.startswith(name):
            return names_to_t[name]
    raise ValueError(f"{path_mod=} is not recognized")


def load_basis(c, complex, is_ring: bool):
    result = {
        "basis": [],
        "bullets": defaultdict(list),
        "bullets_ind": set(),
        "bullets_from_S0": set(),
    }
    index = 0
    sql = f"SELECT mon, s, t FROM {complex}_AdamsE2_basis ORDER BY id"
    for str_mon, s, t in c.execute(sql):
        if is_ring:
            mon = list(map(int, str_mon.split(","))) if len(str_mon) > 0 else []
            mon1 = [[mon[i], mon[i + 1]] for i in range(0, len(mon), 2)]
            mon1.sort()
            mon = sum(mon1, [])
        else:
            mon = list(map(int, str_mon.split(",")))
            if mon[-1] == 0:
                result["bullets_from_S0"].add(index)
            mon1 = [[mon[i], mon[i + 1]] for i in range(0, len(mon) - 1, 2)]
            mon1.sort()
            mon = sum(mon1, []) + [mon[-1], 1]

        result["basis"].append(mon)
        result["bullets"][(t - s, s)].append(index)
        if len(mon) == 2 and mon[1] == 1:
            result["bullets_ind"].add(index)
        index += 1
    return result


def load_basis_from_res(c, complex):
    result = {
        "basis": [],
        "bullets": defaultdict(list),
        "bullets_color": defaultdict(list),
        "bullets_ind": set(),
    }
    sql = f"SELECT id, s, t FROM {complex}_Adams_res_generators ORDER BY id"
    for id, s, t in c.execute(sql):
        result["basis"].append([id, 1])
        result["bullets"][(t - s, s)].append(id)
        result["bullets_color"][(t - s, s)].append("black")
    return result


def load_pi_basis(c, complex, is_ring: bool):
    result = {
        "basis": [],
        "bullets": defaultdict(list),
        "bullets_color": defaultdict(list),
        "bullets_ind": set(),
    }

    id_ind_def = {}
    sql = f"SELECT id, def FROM {complex}_pi_generators_def"
    for id_, def_ in c.execute(sql):
        id_ind_def[id_] = def_

    index = 0
    sql = f"SELECT mon, s, t FROM {complex}_pi_basis ORDER BY id"
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
            result["bullets_ind"].add(index)
            def_type = id_ind_def.get(mon[0], 0)
            if def_type == 1:
                color = "#0000ff"
            elif def_type == 2:
                color = "#0080a0"
            else:
                color = "#00c0ff"
        else:
            color = "black"
        result["bullets_color"][(t - s, s)].append(color)
        index += 1
    return result


def load_basis_ss(c, complex):
    """Return a list of (base, level, diff, id_start_src, id_start_tgt)"""
    ids = defaultdict(int)
    sql = f"SELECT s, t, min(id) FROM {complex}_AdamsE2_ss GROUP BY t, s;"
    for s, t, id in c.execute(sql):
        ids[(s, t)] = id

    sql = f"SELECT base, level, diff, s, t, id FROM {complex}_AdamsE2_ss ORDER BY id"
    result = []
    for base, level, diff, s, t, id in c.execute(sql):
        if level > 5000:
            r = 10000 - level
            s1, t1 = s + r, t + r - 1
        else:
            s1, t1 = s - level, t - level + 1
        result.append((base, level, diff, ids[(s, t)], ids[(s1, t1)]))
    return result


def load_ss_stable_levels(c, complex):
    result = {}
    sql = f"SELECT s, t, l FROM {complex}_AdamsE2_ss_stable_levels GROUP BY t, s;"
    for s, t, l in c.execute(sql):
        result[(t - s, s)] = l
    return result


def load_diffs(c, complex):
    diff_lines = [[], [], [], [], []]
    sql = f"SELECT src, r, tgt FROM {complex}_AdamsE2_ss_diffs;"
    for src, r, tgt in c.execute(sql):
        index = min(r - 2, 4)
        for id_tgt in map(int, tgt.split(",")):
            diff_lines[index].append((src, id_tgt))

    null_diff_lines = [[], [], [], [], []]
    sql = f"SELECT src, r, tgt FROM {complex}_AdamsE2_ss_nd"
    for src, r, tgt in c.execute(sql):
        index = min(r - 2, 4)
        if len(tgt) > 0:
            null_diff_lines[index].append((src, int(tgt.split(",")[-1])))
    return diff_lines, null_diff_lines


def load_gen_names(c, table, letter):
    sql = f"SELECT id, name FROM {table} order by id"
    result = []
    for id, name in c.execute(sql):
        if name:
            result.append(name)
        else:
            if id < 10:
                result.append(f"{letter}_{id}")
            else:
                result.append(f"{letter}_{{{id}}}")
    return result


def load_multiplications(c, table):
    sql = f"SELECT id1, id2, prod FROM {table}"
    result = []
    result_factors = set()
    for id1, id2, prod in c.execute(sql):
        result_factors.add(id1)
        result.append((id1, id2, str2array(prod)))
    return result, result_factors


def load_multiplications_from_res(c, table):
    sql = f"SELECT id, id_ind, prod_h_glo FROM {table}"
    prod = defaultdict(list)
    factors = set()
    for id1, id2, prod_h_glo in c.execute(sql):
        factors.add(id2)
        for id3 in str2array(prod_h_glo):
            prod[(id2, id3)].append(id1)
    prod = [k + (v,) for k, v in prod.items()]
    return prod, factors


def load_pi_multiplications(c, table):
    sql = f"SELECT id1, id2, prod, O FROM {table}"
    result = []
    result_factors = set()
    for id1, id2, prod, O in c.execute(sql):
        result_factors.add(id1)
        result.append((id1, id2, str2array(prod), O))
    return result, result_factors


def load_pi_basis_map(c, table, target):
    sql = (
        f"SELECT id, to_{target}, O_{target} FROM {table} WHERE to_{target} IS NOT NULL"
    )
    result = []
    for id, target, O in c.execute(sql):
        result.append((id, str2array(target), O))
    return result


def load_ss_ring(path_ring):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    path_plot = path_ring.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()

    complex = get_complex_name(path_ring)

    data = load_basis(c_ring, complex, is_ring=True)
    data["gen_names"] = load_gen_names(c_ring, complex + "_AdamsE2_generators", "x")
    data["basis_ss"] = load_basis_ss(c_ring, complex)
    data["stable_levels"] = load_ss_stable_levels(c_plot, complex)
    data["products"], data["products_factors"] = load_multiplications(
        c_plot, complex + "_AdamsE2_ss_products"
    )
    data["diff_lines"], data["null_diff_lines"] = load_diffs(c_plot, complex)

    c_ring.close()
    conn_ring.close()
    c_plot.close()
    conn_plot.close()

    data["name"] = complex
    return data


def load_ss_mod(path_ring, path_mod):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()
    conn_mod = sqlite3.connect(path_mod)
    c_mod = conn_mod.cursor()

    path_plot = path_mod.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()

    complex_ring = get_complex_name(path_ring)
    complex_mod = get_complex_name(path_mod)

    data_ring = {"basis": []}
    data_ring["gen_names"] = load_gen_names(
        c_ring, complex_ring + "_AdamsE2_generators", "x"
    )

    is_ring = False
    data_mod = load_basis(c_mod, complex_mod, is_ring)
    data_mod["gen_names"] = load_gen_names(
        c_mod, complex_mod + "_AdamsE2_generators", "i"
    )
    data_mod["basis_ss"] = load_basis_ss(c_mod, complex_mod)
    data_mod["stable_levels"] = load_ss_stable_levels(c_plot, complex_mod)
    data_mod["products"], data_mod["products_factors"] = load_multiplications(
        c_plot, complex_mod + "_AdamsE2_ss_products"
    )
    data_mod["diff_lines"], data_mod["null_diff_lines"] = load_diffs(
        c_plot, complex_mod
    )

    c_ring.close()
    conn_ring.close()
    c_mod.close()
    conn_mod.close()
    c_plot.close()
    conn_plot.close()

    data_ring["name"] = complex_ring
    data_mod["name"] = complex_mod
    return data_ring, data_mod


def load_ss_from_res(path):
    conn = sqlite3.connect(path)
    c = conn.cursor()

    path_prod = path.replace(".db", "_prod.db")
    conn_prod = sqlite3.connect(path_prod)
    c_prod = conn_prod.cursor()

    complex = get_complex_name(path)

    data = load_basis_from_res(c, complex)
    data["gen_names"] = []
    for i in range(len(data["basis"])):
        if i < 10:
            data["gen_names"].append(f"x_{i}")
        else:
            data["gen_names"].append(f"x_{{{i}}}")
    data["products"], data["products_factors"] = load_multiplications_from_res(
        c_prod, complex + "_Adams_res_products"
    )

    c.close()
    conn.close()
    c_prod.close()
    conn_prod.close()

    data["name"] = complex
    return data


def load_pi_ring(path_ring):
    conn_ring = sqlite3.connect(path_ring)
    c_ring = conn_ring.cursor()

    path_plot = path_ring.replace(".db", "_plot.db")
    conn_plot = sqlite3.connect(path_plot)
    c_plot = conn_plot.cursor()

    complex_ring = get_complex_name(path_ring)

    is_ring = True
    data = load_pi_basis(c_ring, complex_ring, is_ring)
    data["gen_names"] = load_gen_names(c_ring, complex_ring + "_pi_generators", "\\mu")
    data["products"], data["products_factors"] = load_pi_multiplications(
        c_plot, complex_ring + "_pi_basis_products"
    )

    c_ring.close()
    conn_ring.close()
    c_plot.close()
    conn_plot.close()

    data["name"] = complex_ring
    return data


def load_pi_mod(path_ring, path_mod):
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
    data_ring = load_pi_basis(c_ring, complex_ring, is_ring)
    data_ring["gen_names"] = load_gen_names(
        c_ring, complex_ring + "_pi_generators", "\\rho"
    )
    data_ring["products"], data_ring["products_factors"] = load_pi_multiplications(
        c_plot_ring, complex_ring + "_pi_basis_products"
    )
    data_ring["basis_map"] = load_pi_basis_map(
        c_plot_ring, complex_ring + "_pi_basis_maps", complex_mod
    )

    is_ring = False  # TODO: remove unnecessary ones
    data_mod = load_pi_basis(c_mod, complex_mod, is_ring)
    data_mod["gen_names"] = load_gen_names(
        c_mod, complex_mod + "_pi_generators", "\\iota"
    )
    data_mod["products"], data_mod["products_factors"] = load_pi_multiplications(
        c_plot_mod, complex_mod + "_pi_basis_products"
    )
    data_mod["basis_map"] = load_pi_basis_map(
        c_plot_mod, complex_mod + "_pi_basis_maps", complex_ring
    )
    data_mod["t"] = get_complex_t(path_mod)

    c_ring.close()
    conn_ring.close()
    c_mod.close()
    conn_mod.close()
    c_plot_ring.close()
    conn_plot_ring.close()
    c_plot_mod.close()
    conn_plot_mod.close()

    data_ring["name"] = complex_ring
    data_mod["name"] = complex_mod
    return data_ring, data_mod


########################### Process #################################
def smoothen(radius):
    """Make the distribution of radius smooth"""
    radius_ub = defaultdict(lambda: 0.1)  # upper bound of radius
    b_changed = True
    while b_changed:
        for x, y in radius:
            r = radius[(x, y)] + 0.01
            r1 = radius[(x, y)] + 0.01 * 1.4
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
    """Set the radius of bullets in each coordinate (x, y)"""
    radius = defaultdict(lambda: Config.bullets_radius)
    for x, y in data["bullets"]:
        len_bullets_xy = len(data["bullets"][(x, y)])
        length_world = (len_bullets_xy - 1) * Config.bullets_radius * 3
        if length_world > 1 - Config.bullets_radius * 6:
            length_world = 1 - Config.bullets_radius * 6
        sep_world = (
            length_world / (len_bullets_xy - 1)
            if len_bullets_xy > 1
            else Config.bullets_radius * 3
        )
        r = sep_world / 3

        if r < radius[(x, y)]:
            radius[(x, y)] = r
    smoothen(radius)
    return radius


########################### Write #################################
def export_ss_bullets(data, radius, is_ring: bool):
    """Input `data["bullets"]`, `data["basis_ss"]`, `data["stable_levels"]` and `data["bullets_from_S0"]`"""
    tpl_bullets = ""
    index2xyrp = {}

    for x, y in data["bullets"]:
        r = radius[(x, y)]
        len_bullets_xy = len(data["bullets"][(x, y)])
        sep_world = r * 3
        length_world = (len_bullets_xy - 1) * sep_world

        bottom_right_x = x + length_world / 2 * cosA
        bottom_right_y = y + length_world / 2 * sinA
        for i, index in enumerate(data["bullets"][(x, y)]):
            cx = bottom_right_x - i * sep_world * cosA
            cy = bottom_right_y - i * sep_world * sinA
            base, level, diff, offset_x, offset_dx = data["basis_ss"][index]

            # color
            if is_ring:
                str_fill = (
                    ' fill="blue"'
                    if (offset_x + int(base.split(",")[0])) in data["bullets_ind"]
                    else ""
                )
            else:
                index_first_bullet = offset_x + int(base.split(",")[0])
                if index_first_bullet in data["bullets_ind"]:
                    str_fill = ' fill="Maroon"'
                elif index_first_bullet in data["bullets_from_S0"]:
                    str_fill = ""
                else:
                    str_fill = ' fill="#aa00aa"'

            # page
            if level >= data["stable_levels"][(x, y)]:
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

            str_class = "b"
            tpl_bullets += (
                f'<circle id=b{index} class="{str_class}" cx={cx:.6g} cy={cy:.6g} r={r:.6g} {str_fill}'
                f'data-b={base} data-l={level} data-d="{diff}" data-i={offset_x} data-j={offset_dx} data-g={int(len(str_fill) > 0)} data-page={page}>'
                f"<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n"
            )
            index2xyrp[index] = (cx, cy, r, page)
    return tpl_bullets, index2xyrp


def export_bullets(data1d, radius1d, offsets_x):
    """Input `data["bullets"]`, `data["bullets_color"]` and `data["basis"]`"""
    tpl_bullets = ""
    index2xyrp = {}
    offset_basis = 0
    for i_data, (data, radius) in enumerate(zip(data1d, radius1d)):
        offset_x = offsets_x[i_data]
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
                str_fill = (
                    f'fill="{data["bullets_color"][(x, y)][i]}" '
                    if data["bullets_color"][(x, y)][i] != "black"
                    else ""
                )
                page = 200
                level = 9800
                str_class = "b"
                if len(data1d) > 0:
                    str_class += f" b{i_data}"
                tpl_bullets += (
                    f'<circle id=b{offset_basis + index} class="{str_class}" cx={cx:.6g} cy={cy:.6g} r={r:.6g} {str_fill}'
                    f'data-b={index} data-l={level} data-d="{None}" data-i={offset_basis} data-j={0} data-g={int(len(str_fill) > 0)} data-page={page}>'
                    f"<title>({x:.6g}, {y:.6g}) id: {index}</title></circle>\n"
                )
                index2xyrp[offset_basis + index] = (cx, cy, r, page)
        offset_basis += len(data["basis"])
    return tpl_bullets, index2xyrp


def export_basis_ss_prod(data):
    lines = [[], [], []]
    b2g = {1: 0, 3: 1, 7: 2}  # basis_id to gen_id
    for id_gen, id1, prod in data["products"]:
        if id_gen not in b2g:
            continue
        gen = b2g[id_gen]
        for id2 in prod:
            lines[gen].append((id1, id2))
    return lines


def export_pi_basis_prod(
    data, index2xyrp, offset_basis=0, scale=1
):  # TODO: add para dx
    lines = [[], [], [], []]
    dashed_lines = [[], [], [], []]
    b2g = {1: 0, 3: 1, 7: 2, 15: 3}  # basis_id to gen_id
    for id_gen, id1, prod, O in data["products"]:
        if id_gen not in b2g:
            continue
        gen = b2g[id_gen]
        for id2 in prod:
            lines[gen].append((offset_basis + id1, offset_basis + id2))
        if O != -1:
            x1 = index2xyrp[offset_basis + id1][0] - (1 if offset_basis > 0 else 0)
            x2, y2 = round(x1) + (2**gen - 1) * scale, O
            if (x2, y2) in data["bullets"]:
                id2 = data["bullets"][(x2, y2)][0]
                dashed_lines[gen].append((offset_basis + id1, offset_basis + id2))
            else:
                dashed_lines[gen].append(
                    (offset_basis + id1, (x2 + (1 if offset_basis > 0 else 0), y2))
                )
    return lines, dashed_lines


def export_pi_basis_map(data1, data2, index2xyrp, offset1_basis, offset2_basis, dx):
    lines = []
    dashed_lines = []
    for id1, target, O in data1["basis_map"]:
        for id2 in target:
            lines.append((offset1_basis + id1, offset2_basis + id2))
        if O != -1:
            x1 = index2xyrp[offset1_basis + id1][0]
            x2, y2 = round(x1) + dx, O
            if (x2, y2) in data2["bullets"]:
                id2 = data2["bullets"][(x2, y2)][0]
                dashed_lines.append((offset1_basis + id1, offset2_basis + id2))
            else:
                dashed_lines.append((offset1_basis + id1, (x2, y2)))
    return lines, dashed_lines


def element_line(
    x1,
    y1,
    x2,
    y2,
    *,
    width: float,
    page: int,
    class_: str,
    color: str = None,
    dashed=False,
    straight=False,
    r=False,
):
    if dashed and (x1 >= 126.5 or x2 >= 126.5):
        return ""
    attr_more = ""
    if color:
        attr_more += f' stroke="{color}"'
    if dashed:
        attr_more += ' stroke-dasharray="0.2,0.2"'
    if r:
        r = round(y2) - round(y1)
        attr_more += f' data-r="{r}"'

    range_x = f"data-x1={round(min(x1, x2))} data-x2={round(max(x1, x2))}"
    if straight or round(y2) - round(y1) <= 1:
        return f'<line class="{class_}" x1={x1:.6g} y1={y1:.6g} x2={x2:.6g} y2={y2:.6g} {range_x} stroke-width={width:.6g} data-page={page}{attr_more} />'
    else:
        if round(x1) == round(x2):
            return (
                f'<path class="{class_}" d="M {x1:.6g} {y1:.6g} C {x1 + 0.3:.6g} {y1 * 0.7 + y2 * 0.3:.6g}, {x2 + 0.3:.6g} {y1 * 0.3 + y2 * 0.7:.6g}, {x2:.6g} {y2:.6g}" '
                f"{range_x} stroke-width={width:.6g}{attr_more} />\n"
            )
        else:
            x_diff = abs(round(x2) - round(x1))
            sign = 1 if x2 > x1 else -1
            return (
                f'<path class="{class_}" d="M {x1:.6g} {y1:.6g} C {x1 + 0.5 * sign:.6g} {y1 + 0.5 / x_diff:.6g}, {x2 - x_diff * sign / 4 / (y2 - y1):.6g} {y2 - 0.5:.6g}, {x2:.6g} {y2:.6g}" '
                f"{range_x} stroke-width={width:.6g}{attr_more} />\n"
            )


def export_h_lines(index2xyrp, lines, class_: str = None):
    """Plot the h0, h1, h2 lines"""
    tpl_lines = ["", "", ""]
    if class_ is None:
        class_ = ""
    else:
        class_ = " " + class_
    for i in range(3):
        for i1, i2 in lines[i]:
            x1, y1, r1, p1 = index2xyrp[i1]
            x2, y2, r2, p2 = index2xyrp[i2]
            tpl_lines[i] += element_line(
                x1,
                y1,
                x2,
                y2,
                width=min(r1, r2) / 4,
                page=min(p1, p2),
                class_=f"strt_l{class_}",
            )
    return tpl_lines


def export_pi_h_lines(index2xyrp, lines, dashed_lines, class_: str = None):
    """Plot the h0, h1, h2, h3 lines"""
    tpl_lines = ["", "", "", ""]
    extend_colors = ["red", "blue", "green", "DeepSkyBlue"]
    if class_ is None:
        class_ = ""
    else:
        class_ = " " + class_
    for i in range(len(tpl_lines)):
        for i1, i2 in lines[i]:
            x1, y1, r1, p1 = index2xyrp[i1]
            x2, y2, r2, p2 = index2xyrp[i2]
            color = None if round(y2) - round(y1) == 1 else extend_colors[i]
            tpl_lines[i] += element_line(
                x1,
                y1,
                x2,
                y2,
                width=min(r1, r2) / 4,
                page=min(p1, p2),
                class_=f"strt_l{class_}",
                dashed=False,
                color=color,
            )
    for i in range(len(tpl_lines)):
        for i1, O in dashed_lines[i]:
            x1, y1, r1, p1 = index2xyrp[i1]
            if type(O) is tuple:
                x2, y2, r2, p2 = *O, r1, p1
            else:
                x2, y2, r2, p2 = index2xyrp[O]
            tpl_lines[i] += element_line(
                x1,
                y1,
                x2,
                y2,
                width=min(r1, r2) / 4,
                page=min(p1, p2),
                class_=f"strt_l dashed_l{class_}",
                dashed=True,
                color=extend_colors[i],
            )
    return tpl_lines


def export_pi_lines(index2xyrp, lines, dashed_lines, class_: str = None):
    """Plot lines"""
    tpl_lines = ""
    if class_ is None:
        class_ = ""
    else:
        class_ = " " + class_
    for i1, i2 in lines:
        x1, y1, r1, p1 = index2xyrp[i1]
        x2, y2, r2, p2 = index2xyrp[i2]
        tpl_lines += element_line(
            x1,
            y1,
            x2,
            y2,
            width=min(r1, r2) / 4,
            page=min(p1, p2),
            class_=f"strt_l{class_}",
            dashed=False,
        )
    for i1, O in dashed_lines:
        x1, y1, r1, p1 = index2xyrp[i1]
        if type(O) is tuple:
            x2, y2, r2, p2 = *O, r1, p1
        else:
            x2, y2, r2, p2 = index2xyrp[O]
        tpl_lines += element_line(
            x1,
            y1,
            x2,
            y2,
            width=min(r1, r2) / 4,
            page=min(p1, p2),
            class_=f"strt_l dashed_l{class_}",
            dashed=True,
        )
    return tpl_lines


def export_diffs(data, index2xyrp):
    tpl_diff_lines = ["", "", "", "", ""]
    for index, lines in enumerate(data["diff_lines"]):
        for i1, i2 in lines:
            x1, y1, r1, p1 = index2xyrp[i1]
            x2, y2, r2, p2 = index2xyrp[i2]
            tpl_diff_lines[index] += element_line(
                x1,
                y1,
                x2,
                y2,
                width=min(r1, r2) / 4,
                page=min(p1, p2),
                class_="diff_l",
                dashed=False,
                straight=True,
                r=True,
            )
    for index, lines in enumerate(data["null_diff_lines"]):
        for i1, i2 in lines:
            x1, y1, r1, p1 = index2xyrp[i1]
            x2, y2, r2, p2 = index2xyrp[i2]
            tpl_diff_lines[index] += element_line(
                x1,
                y1,
                x2,
                y2,
                width=min(r1, r2) / 4,
                page=min(p1, p2),
                class_="dashed_l",
                dashed=True,
                straight=True,
                r=True,
            )

    return tpl_diff_lines


def export_gen_names_to_js(data1d):
    content_js = "const gen_names = [\n"
    for data in data1d:
        for name in data["gen_names"]:
            name1 = name.replace("\\", "\\\\")
            content_js += f' "{name1}",\n'
        content_js += "\n"
    content_js += "];\n\n"
    return content_js


def export_basis_to_js(data1d):
    content_js = "const basis = [\n"
    offset = 0
    for data in data1d:
        for mon in data["basis"]:
            if data["name"] != "S0":
                mon[-2] += offset
            content_js += f" {mon},\n"
        content_js += "\n"
        if data["name"] == "S0":
            offset = len(data["gen_names"])
    content_js += f"];\n\n"
    return content_js


def export_ss_prod_to_js(data):
    content_js = f"const arr_factors = {sorted(data['products_factors'])};\n"
    content_js += "const basis_prod = {\n"
    for id1, id2, prod in data["products"]:
        if len(prod) > 0:
            prod = [i for i in prod]
            content_js += f' "{id1},{id2}": {prod},\n'
    content_js += "\n"
    content_js += "};\n\n"
    return content_js


def export_pi_prod_to_js(data1d):
    content_js = f"const arr_factors = {sorted(data1d[0]['products_factors'])};\n"
    content_js += "const basis_prod = {\n"
    offset = 0
    for data in data1d:
        for id1, id2, prod, O in data["products"]:
            if len(prod) > 0 or O >= 0:
                prod = [i + offset for i in prod]
                content_js += (
                    f' "{id1},{id2 + offset}": [{prod}, {O if O > 0 else 300}],\n'
                )
        content_js += "\n"
        offset += len(data["basis"])
    content_js += "};\n\n"
    return content_js


def export_ss_ring(data_ring, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Adams Spectral Sequence of {data_ring['name']}"

    # js
    content_js = 'const MODE="SS";\n'
    content_js += export_gen_names_to_js([data_ring])
    content_js += export_basis_to_js([data_ring])
    content_js += export_ss_prod_to_js(data_ring)

    # html
    radius = get_radius(data_ring)
    tpl_bullets, index2xyrp = export_ss_bullets(data_ring, radius, is_ring=True)
    lines = export_basis_ss_prod(data_ring)
    tpl_lines = export_h_lines(index2xyrp, lines)
    tpl_diff_lines = export_diffs(data_ring, index2xyrp)

    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("<!-- bullets:b8999a4e -->", tpl_bullets)

    content = content.replace("<!-- h0_lines:839e707e -->", tpl_lines[0])
    content = content.replace("<!-- h1_lines:2b707f82 -->", tpl_lines[1])
    content = content.replace("<!-- h2_lines:e766a4f6 -->", tpl_lines[2])

    content = content.replace("<!-- d2_lines:d9dec788 -->", tpl_diff_lines[0])
    content = content.replace("<!-- d3_lines:1a07bf89 -->", tpl_diff_lines[1])
    content = content.replace("<!-- d4_lines:7d935e97 -->", tpl_diff_lines[2])
    content = content.replace("<!-- d5_lines:f6d7f992 -->", tpl_diff_lines[3])
    content = content.replace("<!-- d6_lines:05c49c4a -->", tpl_diff_lines[4])

    with open(path_js, "w") as file:
        file.write(content_js)
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)


def export_ss_mod(data_ring, data_mod, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Adams Spectral Sequence of {data_mod['name']}"

    # js
    content_js = 'const MODE="SS";\n'
    content_js += export_gen_names_to_js([data_ring, data_mod])
    content_js += export_basis_to_js([data_ring, data_mod])
    content_js += export_ss_prod_to_js(data_mod)

    # html
    radius = get_radius(data_mod)
    tpl_bullets, index2xyrp = export_ss_bullets(data_mod, radius, is_ring=False)
    lines = export_basis_ss_prod(data_mod)
    tpl_lines = export_h_lines(index2xyrp, lines)
    tpl_diff_lines = export_diffs(data_mod, index2xyrp)

    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("<!-- bullets:b8999a4e -->", tpl_bullets)

    content = content.replace("<!-- h0_lines:839e707e -->", tpl_lines[0])
    content = content.replace("<!-- h1_lines:2b707f82 -->", tpl_lines[1])
    content = content.replace("<!-- h2_lines:e766a4f6 -->", tpl_lines[2])

    content = content.replace("<!-- d2_lines:d9dec788 -->", tpl_diff_lines[0])
    content = content.replace("<!-- d3_lines:1a07bf89 -->", tpl_diff_lines[1])
    content = content.replace("<!-- d4_lines:7d935e97 -->", tpl_diff_lines[2])
    content = content.replace("<!-- d5_lines:f6d7f992 -->", tpl_diff_lines[3])
    content = content.replace("<!-- d6_lines:05c49c4a -->", tpl_diff_lines[4])

    with open(path_js, "w") as file:
        file.write(content_js)
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)


def export_ss_from_res(data, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Adams Spectral Sequence of {data['name']}"

    # js
    content_js = 'const MODE="FromRes";\n'
    content_js += export_gen_names_to_js([data])
    content_js += export_basis_to_js([data])
    content_js += export_ss_prod_to_js(data)

    # html
    radius = get_radius(data)
    tpl_bullets, index2xyrp = export_bullets([data], [radius], [0])
    lines = export_basis_ss_prod(data)
    tpl_lines = export_h_lines(index2xyrp, lines)

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


def export_pi_ring(data, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Homotopy of {data['name']}"

    # js
    content_js = 'const MODE="Pi";\n'
    content_js += export_gen_names_to_js([data])
    content_js += export_basis_to_js([data])
    content_js += export_pi_prod_to_js([data])

    # html
    radius = get_radius(data)
    tpl_bullets, index2xyrp = export_bullets([data], [radius], [0])
    lines, dashed_lines = export_pi_basis_prod(data, index2xyrp, 0)
    tpl_lines = export_pi_h_lines(index2xyrp, lines, dashed_lines)

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


def export_pi_mod(data_ring, data_mod, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Homotopy of {data_mod['name']}"

    # js
    content_js = 'const MODE="DualPi";\n'
    content_js += f'const SEP_MAX_WIDTH={data_mod["t"] + 1};\n'
    content_js += "var SEP_RIGHT=1;\n"
    content_js += "var SEP_WIDTH=1;\n"
    content_js += export_gen_names_to_js([data_ring, data_mod])
    content_js += export_basis_to_js([data_ring, data_mod])
    content_js += export_pi_prod_to_js([data_ring, data_mod])

    # html
    radius = get_radius(data_ring)
    radius1 = get_radius(data_mod)
    tpl_bullets, index2xyrp = export_bullets(
        [data_ring, data_mod], [radius, radius1], [0, 1]
    )
    lines, dashed_lines = export_pi_basis_prod(data_ring, index2xyrp, 0)
    lines1, dashed_lines1 = export_pi_basis_prod(
        data_mod, index2xyrp, len(data_ring["basis"])
    )
    lines_bc, dashed_lines_bc = export_pi_basis_map(
        data_ring, data_mod, index2xyrp, 0, len(data_ring["basis"]), 0
    )
    lines_tc, dashed_lines_tc = export_pi_basis_map(
        data_mod, data_ring, index2xyrp, len(data_ring["basis"]), 0, -data_mod["t"] - 1
    )
    tpl_lines = export_pi_h_lines(index2xyrp, lines, dashed_lines, "l0")
    tpl_lines1 = export_pi_h_lines(index2xyrp, lines1, dashed_lines1, "l1")
    tpl_lines_bc = export_pi_lines(index2xyrp, lines_bc, dashed_lines_bc, "lbc")
    tpl_lines_tc = export_pi_lines(index2xyrp, lines_tc, dashed_lines_tc, "ltc")
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


def export_pi_exact(data_ring, data_mod, path_html=None, path_js=None):
    path_html = path_html or path_html_tmp
    path_js = path_js or path_js_tmp

    tpl_title = f"Exact sequence for {data_mod['name']}"
    T = data_mod["t"]

    # js
    content_js = 'const MODE="Exact";\n'
    content_js += "const arr_factors=[];\n"
    content_js += f"const T_TOP_CELL = {T};\n"
    content_js += export_gen_names_to_js([data_ring, data_mod])
    content_js += export_basis_to_js([data_ring, data_ring, data_mod])

    # html
    data_ring1 = copy.deepcopy(data_ring)
    data_ring2 = copy.deepcopy(data_ring)
    data_mod0 = copy.deepcopy(data_mod)
    data_ring1["bullets"] = {
        (3 * k[0] + 1, k[1]): v for k, v in data_ring["bullets"].items()
    }
    data_ring1["bullets_color"] = {
        (3 * k[0] + 1, k[1]): v for k, v in data_ring["bullets_color"].items()
    }
    data_ring2["bullets"] = {
        (3 * (k[0] + T - 1) + 2, k[1]): v for k, v in data_ring["bullets"].items()
    }
    data_ring2["bullets_color"] = {
        (3 * (k[0] + T - 1) + 2, k[1]): v for k, v in data_ring["bullets_color"].items()
    }
    data_mod0["bullets"] = {(3 * k[0], k[1]): v for k, v in data_mod["bullets"].items()}
    data_mod0["bullets_color"] = {
        (3 * k[0], k[1]): v for k, v in data_mod["bullets_color"].items()
    }

    radius1 = get_radius(data_ring1)
    radius2 = get_radius(data_ring2)
    radius0 = get_radius(data_mod0)
    tpl_bullets, index2xyrp = export_bullets(
        [data_ring1, data_ring2, data_mod0], [radius1, radius2, radius0], [0, 0, 0]
    )

    lines, dashed_lines = export_pi_basis_prod(data_ring1, index2xyrp, 0, 3)
    for i in range(4):
        if 2**i == T:
            lines[i] = [[i1 + len(data_ring1["basis"]), i2] for i1, i2 in lines[i]]
            dashed_lines[i] = [
                [i1 + len(data_ring1["basis"]), i2] for i1, i2 in dashed_lines[i]
            ]
        else:
            lines[i] = []
            dashed_lines[i] = []

    lines_bc, dashed_lines_bc = export_pi_basis_map(
        data_ring1,
        data_mod0,
        index2xyrp,
        offset1_basis=0,
        offset2_basis=len(data_ring1["basis"]) + len(data_ring2["basis"]),
        dx=-1,
    )
    lines_tc, dashed_lines_tc = export_pi_basis_map(
        data_mod0,
        data_ring2,
        index2xyrp,
        offset1_basis=len(data_ring1["basis"]) + len(data_ring2["basis"]),
        offset2_basis=len(data_ring1["basis"]),
        dx=-1,
    )
    tpl_lines = export_pi_h_lines(index2xyrp, lines, dashed_lines)
    tpl_lines_bc = export_pi_lines(index2xyrp, lines_bc, dashed_lines_bc, "lbc")
    tpl_lines_tc = export_pi_lines(index2xyrp, lines_tc, dashed_lines_tc, "ltc")

    with open(path_tpl, encoding="utf-8") as fp:
        content_tpl = fp.read()

    content = content_tpl.replace("title:13dd3d15", tpl_title)
    content = content.replace("<!-- bullets:b8999a4e -->", tpl_bullets)

    content = content.replace("<!-- h0_lines:839e707e -->", tpl_lines[0])
    content = content.replace("<!-- h1_lines:2b707f82 -->", tpl_lines[1])
    content = content.replace("<!-- h2_lines:e766a4f6 -->", tpl_lines[2])
    content = content.replace("<!-- h3_lines:0x2b309dc4U -->", tpl_lines[3])

    content = content.replace("<!-- g_maroon_lines:7aa92c76 -->", tpl_lines_bc)
    content = content.replace("<!-- g_purple_lines:114fd993 -->", tpl_lines_tc)

    with open(path_js, "w") as file:
        file.write(content_js)
    with open(path_html, "w", encoding="utf-8") as file:
        file.write(content)
