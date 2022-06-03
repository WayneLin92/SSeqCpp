"""doc"""
import sqlite3
import textwrap
import struct


path_db = R"C:\Users\lwnpk\Documents\Projects\algtop_cpp_build\bin\Release\AdamsE2Export.db"
js_file_path = R"C:\Users\lwnpk\OneDrive\Projects\HTML\WayneLin92.github.io\programs\ssplot\scripts\data.js"
kLevelMax = 10000
conn = sqlite3.connect(path_db)
with conn:
    c = conn.cursor()
    with open(js_file_path, "w") as file:
        file.write("var data_bullets = [\n")
        index = 0
        arrows = []
        sql = f"SELECT mon, s, t FROM AdamsE2_basis ORDER BY id"

        id_mons = {}
        index = 0
        for b_mon, s, t in c.execute(sql):
            ell = len(b_mon)
            assert len(b_mon) % 8 == 0
            mon = []
            for i in range(0, ell, 8):
                g = struct.unpack('<i', b_mon[i:i+4])[0]
                e = struct.unpack('<i', b_mon[i+4:i+8])[0]
                mon.append((g, e))
            mon = tuple(p for p in mon)

            id_mons[mon] = index
            index += 1

            if len(mon) == 1 and mon[0][1] == 1:
                file.write(f"[{t - s}, {s}, \"\", {1}],\n")  # indecomposables
            else:
                file.write(f"[{t - s}, {s}, \"\", {0}],\n")
            
        file.write(textwrap.dedent("""\
            ];
            var data_lines = [
            """))
        for mon in id_mons:
            dict_mon = dict(mon)
            if 0 in dict_mon:
                dict_mon_div_h0 = dict_mon.copy()
                
                if dict_mon_div_h0[0] == 1:
                    del dict_mon_div_h0[0]
                else:
                    dict_mon_div_h0[0] -= 1

                mon_div_h0 = tuple(sorted(dict_mon_div_h0.items()))
                file.write(f"[{id_mons[mon_div_h0]}, {id_mons[mon]}],\n")
            
            if 1 in dict_mon:
                dict_mon_div_h1 = dict_mon.copy()
                
                if dict_mon_div_h1[1] == 1:
                    del dict_mon_div_h1[1]
                else:
                    dict_mon_div_h1[1] -= 1

                mon_div_h1 = tuple(sorted(dict_mon_div_h1.items()))
                file.write(f"[{id_mons[mon_div_h1]}, {id_mons[mon]}],\n")
            
            if 2 in dict_mon:
                dict_mon_div_h2 = dict_mon.copy()
                
                if dict_mon_div_h2[2] == 1:
                    del dict_mon_div_h2[2]
                else:
                    dict_mon_div_h2[2] -= 1

                mon_div_h2 = tuple(sorted(dict_mon_div_h2.items()))
                file.write(f"[{id_mons[mon_div_h2]}, {id_mons[mon]}],\n")
                

        file.write(textwrap.dedent("""\
            ];
                        
            var data_arrows = [
            """))
        for i, j in arrows:
            file.write(f"[{i}, {i + 1}],\n")
        file.write("]\n")
        
        file.write("\nbullets_tilt_angle_deg = -30;\n")

    c.close()
conn.close()