#include "main.h"
#include "mylog.h"
#include <fstream>

/* Add the differentials from diagram1 to diagram2 */
void Migrate_ss(const Diagram& diagram1, Diagram& diagram2)
{
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
                    const auto& sc1 = diagram1.GetRecentSc(nodes_ss1, deg);
                    for (size_t i = 0; i < sc1.levels.size(); ++i) {
                        if (sc1.levels[i] > LEVEL_MAX / 2) {
                            int r = LEVEL_MAX - sc1.levels[i];
                            if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                                diff = NULL_DIFF;
                            else
                                diff = sc1.diffs[i];
                            if (diff != NULL_DIFF) {
                                int count = diagram2.SetRingDiffGlobal(iCw2, deg, sc1.basis[i], diff, r);
                                if (count > 0)
                                    Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], diff, r);
                            }
                            else {
                                int count = diagram2.SetRingDiffGlobal(iCw2, deg, sc1.basis[i], int1d{}, r - 1);
                                if (count > 0)
                                    Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], int1d{}, r - 1);
                            }
                        }
                        else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                            int r = sc1.levels[i] + 1;
                            const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                            int count = diagram2.SetRingDiffGlobal(iCw2, deg_src, {}, sc1.basis[i], r);
                            if (count > 0)
                                Logger::LogDiffInv(0, enumReason::migrate, name, deg, {}, sc1.basis[i], r);
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
                    const auto& sc1 = diagram1.GetRecentSc(nodes_ss1, deg);
                    for (size_t i = 0; i < sc1.levels.size(); ++i) {
                        if (sc1.levels[i] > LEVEL_MAX / 2) {
                            int r = LEVEL_MAX - sc1.levels[i];
                            if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                                diff = NULL_DIFF;
                            else
                                diff = sc1.diffs[i];
                            if (diff != NULL_DIFF) {
                                int count = diagram2.SetModuleDiffGlobal(iMod2, deg, sc1.basis[i], diff, r);
                                if (count > 0)
                                    Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], diff, r);
                            }
                            else {
                                int count = diagram2.SetModuleDiffGlobal(iMod2, deg, sc1.basis[i], int1d{}, r - 1);
                                if (count > 0)
                                    Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], int1d{}, r - 1);
                            }
                        }
                        else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                            int r = sc1.levels[i] + 1;
                            const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                            int count = diagram2.SetModuleDiffGlobal(iMod2, deg_src, {}, sc1.basis[i], r);
                            if (count > 0)
                                Logger::LogDiffInv(0, enumReason::migrate, name, deg, {}, sc1.basis[i], r);
                        }
                    }
                }
            }
        }
    }
}

/* Add the homotopy from diagram1 to diagram2
 * Currently we only support migration from small range to bigger range
 */
void Migrate_htpy(const Diagram& diagram1, Diagram& diagram2)
{
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
}

int main_migrate_ss(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name1, diagram_name2;

    myio::CmdArg1d args = {{"diagram1", &diagram_name1}, {"diagram2", &diagram_name2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag_no_op = DeduceFlag::no_op;
    Diagram diagram1(diagram_name1, flag_no_op);
    Diagram diagram2(diagram_name2, flag_no_op, false);

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

    DeduceFlag flag = DeduceFlag::homotopy | DeduceFlag::homotopy_def;
    Diagram diagram1(diagram_name1, flag);
    Diagram diagram2(diagram_name2, flag, false);

    try {
        Migrate_htpy(diagram1, diagram2);
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
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    std::vector<std::string> names, paths;
    std::vector<int> isRing;
    GetAllDbNames(diagram_name, names, paths, isRing);
    for (size_t k = 0; k < names.size(); ++k) {
        DBSS db(paths[k]);
        auto nodes_ss = db.load_basis_ss(names[k]);
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
        db.update_basis_ss(names[k], nodes_ss);
        db.end_transaction();
    }

    return 0;
}