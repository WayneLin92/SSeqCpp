#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <fstream>

/* Add the differentials from diagram1 to diagram2 */
void Migrate_ss(const Diagram& diagram1, Diagram& diagram2)
{
    auto flag = SSFlag::no_op;
    int count = 0;
    int1d diff;
    const size_t num_cw = diagram1.GetRings().size() + diagram1.GetModules().size();
    for (size_t iCw = 0; iCw < num_cw; ++iCw) {
        if (iCw < diagram1.GetRings().size()) {
            auto& name = diagram1.GetRings()[iCw].name;
            int iCw2 = diagram2.GetRingIndexByName(name);
            if (iCw2 != -1) {
                auto& nodes_ss1 = diagram1.GetRings()[iCw].nodes_ss;
                auto& nodes_ss2 = diagram2.GetRings()[iCw2].nodes_ss;

                int t_max2 = nodes_ss2.front().rbegin()->first.t;
                for (auto& [deg, _] : nodes_ss1.front()) {
                    if (deg.t > t_max2)
                        break;
                    const auto& sc1 = ut::GetRecentValue(nodes_ss1, deg);
                    for (size_t i = 0; i < sc1.levels.size(); ++i) {
                        if (sc1.levels[i] > LEVEL_MAX / 2) {
                            int r = LEVEL_MAX - sc1.levels[i];
                            if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                                diff = NULL_DIFF;
                            else
                                diff = sc1.diffs[i];
                            if (diff != NULL_DIFF) {
                                if (diagram2.IsNewDiff(nodes_ss2, deg, sc1.basis[i], diff, r)) {
                                    Logger::LogDiff(0, EnumReason::migrate, name, deg, sc1.basis[i], diff, r);
                                    count += diagram2.SetRingDiffGlobal(iCw2, deg, sc1.basis[i], diff, r, true, flag);
                                }
                            }
                            else {
                                if (diagram2.IsNewDiff(nodes_ss2, deg, sc1.basis[i], int1d{}, r - 1)) {
                                    Logger::LogDiff(0, EnumReason::migrate, name, deg, sc1.basis[i], int1d{}, r - 1);
                                    count += diagram2.SetRingDiffGlobal(iCw2, deg, sc1.basis[i], int1d{}, r - 1, true, flag);
                                }
                            }
                        }
                        else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                            int r = sc1.levels[i] + 1;
                            const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                            if (diagram2.IsNewDiff(nodes_ss2, deg_src, {}, sc1.basis[i], r)) {
                                Logger::LogDiffInv(0, EnumReason::migrate, name, deg_src, deg, {}, sc1.basis[i], r);
                                count += diagram2.SetRingDiffGlobal(iCw2, deg_src, {}, sc1.basis[i], r, true, flag);
                            }
                        }
                    }
                }
            }
        }
        else {
            size_t iMod = iCw - diagram1.GetRings().size();
            const auto& name = diagram1.GetModules()[iMod].name;
            int iMod2 = diagram2.GetModuleIndexByName(name);
            if (iMod2 != -1) {
                const auto& nodes_ss1 = diagram1.GetModules()[iMod].nodes_ss;
                const auto& nodes_ss2 = diagram2.GetModules()[iMod2].nodes_ss;
                int t_max2 = nodes_ss2.front().rbegin()->first.t;
                for (auto& [deg, _] : nodes_ss1.front()) {
                    if (deg.t > t_max2)
                        break;
                    const auto& sc1 = ut::GetRecentValue(nodes_ss1, deg);
                    for (size_t i = 0; i < sc1.levels.size(); ++i) {
                        if (sc1.levels[i] > LEVEL_MAX / 2) {
                            int r = LEVEL_MAX - sc1.levels[i];
                            if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                                diff = NULL_DIFF;
                            else
                                diff = sc1.diffs[i];
                            if (diff != NULL_DIFF) {
                                if (diagram2.IsNewDiff(nodes_ss2, deg, sc1.basis[i], diff, r)) {
                                    Logger::LogDiff(0, EnumReason::migrate, name, deg, sc1.basis[i], diff, r);
                                    count += diagram2.SetModuleDiffGlobal(iMod2, deg, sc1.basis[i], diff, r, true, flag);
                                }
                            }
                            else {
                                if (diagram2.IsNewDiff(nodes_ss2, deg, sc1.basis[i], int1d{}, r - 1)) {
                                    Logger::LogDiff(0, EnumReason::migrate, name, deg, sc1.basis[i], int1d{}, r - 1);
                                    count += diagram2.SetModuleDiffGlobal(iMod2, deg, sc1.basis[i], int1d{}, r - 1, true, flag);
                                }
                            }
                        }
                        else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                            int r = sc1.levels[i] + 1;
                            const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                            if (diagram2.IsNewDiff(nodes_ss2, deg_src, {}, sc1.basis[i], r)) {
                                Logger::LogDiffInv(0, EnumReason::migrate, name, deg_src, deg, {}, sc1.basis[i], r);
                                count += diagram2.SetModuleDiffGlobal(iMod2, deg_src, {}, sc1.basis[i], r, true, flag);
                            }
                        }
                    }
                }
            }
        }
    }
    Logger::LogSummary("Changed differentials", count);
}

/* Add the homotopy from diagram1 to diagram2
 * Currently we only support migration from small range to bigger range
 */
//void Migrate_htpy(const Diagram& diagram1, Diagram& diagram2)
//{
    // for (size_t iSS = 0; iSS < diagram1.GetAllBasisSs().size(); ++iSS) {
    //     if (iSS == 0) {
    //         auto& pi_gb1 = diagram1.GetRings().pi_gb;
    //         auto& pi_gb2 = diagram2.GetRings().pi_gb;
    //         int t_max2 = diagram2.GetRings().nodes_ss.front().rbegin()->first.t;
    //         if (pi_gb2.gen_degs().size() > 1) {  //// TODO: Support merge of homotopy information
    //             Logger::LogException(0, 0x43dbbfcaU, "Diagram2 should be reset first.\n");
    //             throw MyException(0x43dbbfcaU, "Diagram2 should be reset first.");
    //         }
    //         algZ::Poly1d rels = pi_gb1.data();
    //         uint32_t i = 0;
    //         for (auto& deg : pi_gb1.gen_degs()) {
    //             if (deg.t <= t_max2) {
    //                 pi_gb2.AddGen(deg);
    //                 if (deg.stem() % 2 == 1)
    //                     rels.push_back(algZ::Mon::two_x_square(i, deg.s));
    //             }
    //             ++i;
    //         }
    //         diagram2.GetRings().pi_gen_Einf = diagram1.GetRings().pi_gen_Einf;
    //         diagram2.GetRings().pi_gen_defs = diagram1.GetRings().pi_gen_defs;
    //         diagram2.GetRings().pi_gen_def_mons = diagram1.GetRings().pi_gen_def_mons;
    //         diagram2.AddPiRelsRing(std::move(rels));
    //     }
    //     else {
    //         auto& ssCof1 = diagram1.GetModules()[iSS - 1];
    //         auto& ssCof2 = diagram2.GetModules()[iSS - 1];
    //         auto& pi_gb1 = ssCof1.pi_gb;
    //         auto& pi_gb2 = ssCof2.pi_gb;
    //         int t_max2 = ssCof2.nodes_ss.front().rbegin()->first.t;
    //         if (pi_gb2.v_degs().size() > 1) {  //// TODO: Support merge of homotopy information
    //             Logger::LogException(0, 0x92cb2691, "Diagram2 should be reset first.\n");
    //             throw MyException(0x92cb2691, "Diagram2 should be reset first.");
    //         }
    //         for (auto& deg : pi_gb1.v_degs())
    //             if (deg.t <= t_max2)
    //                 pi_gb2.AddGen(deg);
    //         ssCof2.pi_gen_Einf = ssCof1.pi_gen_Einf;
    //         ssCof2.pi_gen_defs = ssCof1.pi_gen_defs;
    //         ssCof2.pi_gen_def_mons = ssCof1.pi_gen_def_mons;
    //         ssCof2.nodes_pi_qt = ssCof1.nodes_pi_qt;
    //         diagram2.AddPiRelsCof(iSS - 1, pi_gb1.data());
    //     }
    // }
    // int count_ss = 0, count_htpy = 0;
    // diagram2.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_htpy, 0);
//}

void ImportChuaD2(Diagram& diagram)
{
    myio::DbAdamsSS db_S0("chua_d2/S0_AdamsSS_t261.db");
    myio::DbAdamsSS db_Chua_d2("chua_d2/Chua_d2.db");

    std::map<AdamsDeg, int2d> d2_Chua;
    {
        myio::Statement stmt(db_Chua_d2, "select stem, Adams_filtration, \"index\", REPLACE(SUBSTR(target, 2, LENGTH(target) - 2), \" \", \"\") from Chua_d2 order by rowid;");
        while (stmt.step() == MYSQLITE_ROW) {
            int stem = stmt.column_int(0), s = stmt.column_int(1), index1 = stmt.column_int(2);
            AdamsDeg deg = {s, stem + s};
            if (index1 != (int)d2_Chua[deg].size())
                throw MyException(0x794fed62, "Not ordered by index");
            int1d target = myio::Deserialize<int1d>(stmt.column_str(3));
            int1d target1;
            for (int i = 0; i < (int)target.size(); ++i)
                if (target[i])
                    target1.push_back(i);
            d2_Chua[deg].push_back(std::move(target1));
        }
    }

    std::map<AdamsDeg, int2d> BFCC_to_res_id;
    {
        myio::Statement stmt(db_Chua_d2, "select stem, s, i, res_index from BFCC_to_res_id order by rowid;");
        while (stmt.step() == MYSQLITE_ROW) {
            int stem = stmt.column_int(0), s = stmt.column_int(1), index1 = stmt.column_int(2);
            AdamsDeg deg = {s, stem + s};
            if (index1 != (int)BFCC_to_res_id[deg].size())
                throw MyException(0x507f4941, "Not ordered by index");
            int1d res_id = myio::Deserialize<int1d>(stmt.column_str(3));
            int1d res_id1;
            for (int i = 0; i < (int)res_id.size(); ++i)
                if (res_id[i])
                    res_id1.push_back(i);
            BFCC_to_res_id[deg].push_back(std::move(res_id1));
        }
    }

    std::map<AdamsDeg, int2d> basis_to_res_id;
    {
        myio::Statement stmt(db_S0, "select s, t, repr from S0_AdamsE2_basis order by id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            AdamsDeg deg = {s, t};
            int1d repr = myio::Deserialize<int1d>(stmt.column_str(2));
            basis_to_res_id[deg].push_back(std::move(repr));
        }
    }
    for (auto& [deg, map] : basis_to_res_id) {
        int1d min_res_id1d;
        for (auto& map_i : map)
            if (!map_i.empty())
                min_res_id1d.push_back(*std::min_element(map_i.begin(), map_i.end()));
        int min_res_id = *std::min_element(min_res_id1d.begin(), min_res_id1d.end());
        for (auto& map_i : map)
            for (auto& j : map_i)
                j -= min_res_id;
    }

    std::map<AdamsDeg, int2d> res_id_to_basis;
    for (auto& [deg, map] : basis_to_res_id) {
        for (int i = 0; i < map.size(); ++i) {
            int2d image, kernel, g;
            lina::SetLinearMap(map, image, kernel, g);
            res_id_to_basis[deg].push_back(lina::GetImage(image, g, {i}));
        }
    }

    std::map<AdamsDeg, int2d> x, dx;
    for (auto& [deg, map] : d2_Chua) {
        AdamsDeg deg_dx = deg + AdamsDeg(2, 1);
        if (ut::has(BFCC_to_res_id, deg) && ut::has(res_id_to_basis, deg) && ut::has(BFCC_to_res_id, deg_dx) && ut::has(res_id_to_basis, deg_dx)) {
            for (int i = 0; i < (int)map.size(); ++i) {
                int1d x_res_id = BFCC_to_res_id[deg][i];
                int1d x_basis;
                for (int j : x_res_id)
                    x_basis = lina::add(x_basis, res_id_to_basis[deg][j]);
                x[deg].push_back(x_basis);

                int1d dx_res_id;
                for (int j : map[i])
                    dx_res_id = lina::add(dx_res_id, BFCC_to_res_id[deg_dx][j]);
                int1d dx_basis;
                for (int j : dx_res_id)
                    dx_basis = lina::add(dx_basis, res_id_to_basis[deg_dx][j]);
                dx[deg].push_back(dx_basis);
            }
        }
    }

    for (auto& [deg, x_d] : x) {
        for (size_t i = 0; i < x_d.size(); ++i)
            diagram.SetRingDiffGlobal(0, deg, x_d[i], dx[deg][i], 2, false, SSFlag::no_op);
    }
}

int main_import_chua_d2(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;

    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag_no_op = SSFlag::no_op;
    Diagram diagram(diagram_name, flag_no_op);

    try {
        ImportChuaD2(diagram);
        diagram.save(diagram_name, flag_no_op);
    }
    /*catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}

int main_migrate_ss(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name1, diagram_name2;

    myio::CmdArg1d args = {{"diagram1", &diagram_name1}, {"diagram2", &diagram_name2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag_no_op = SSFlag::no_op;
    Diagram diagram1(diagram_name1, flag_no_op, false);
    Diagram diagram2(diagram_name2, flag_no_op);

    try {
        Migrate_ss(diagram1, diagram2);
        diagram2.save(diagram_name2, flag_no_op);
    }
    /*catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}

int main_migrate_pi(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name1, diagram_name2;

    myio::CmdArg1d args = {{"diagram1", &diagram_name1}, {"diagram2", &diagram_name2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::pi | SSFlag::pi_def;
    Diagram diagram1(diagram_name1, flag);
    Diagram diagram2(diagram_name2, flag, false);

    try {
        //Migrate_htpy(diagram1, diagram2);
        diagram2.save(diagram_name2, flag);
    }
    /* catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}

int main_resetfrom(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;
    std::string diagram_name_from;

    myio::CmdArg1d args = {{"diagram", &diagram_name}, {"diagram_from", &diagram_name_from}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    std::vector<std::string> names, paths;
    std::vector<int> isRing;
    GetAllDbNames(diagram_name, names, paths, isRing);
    std::vector<std::string> names2, paths2;
    std::vector<int> isRing2;
    GetAllDbNames(diagram_name_from, names2, paths2, isRing2);

    for (size_t k = 0; k < names.size(); ++k) {
        std::ifstream src(paths2[k], std::ios::binary);
        std::ofstream dst(paths[k], std::ios::binary);
        dst << src.rdbuf();
    }

    return 0;
}

int main_truncate(int argc, char** argv, int index, const char* desc)
{
    std::string diagram_name;

    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    std::vector<std::string> names, paths;
    std::vector<int> isRing;
    GetAllDbNames(diagram_name, names, paths, isRing);
    for (size_t k = 0; k < names.size(); ++k) {
        DBSS db(paths[k]);
        auto nodes_ss = db.load_ss(names[k]);
        int t_max = nodes_ss.rbegin()->first.t;

        for (auto& [d, sc] : nodes_ss) {
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] > LEVEL_PERM) {
                    int r = LEVEL_MAX - sc.levels[i];
                    if (d.t + r - 1 > t_max && sc.diffs[i] != NULL_DIFF)
                        sc.diffs[i] = NULL_DIFF;
                }
            }
        }

        db.begin_transaction();
        db.update_ss(names[k], nodes_ss);
        db.end_transaction();
    }

    return 0;
}