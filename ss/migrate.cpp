#include "main.h"
#include "mylog.h"
#include <fstream>

/* Add the differentials from diagram1 to diagram2 */
void Migrate_ss(const Diagram& diagram1, Diagram& diagram2)
{
    int1d diff;
    for (size_t iSS = 0; iSS < diagram1.GetAllBasisSs().size(); ++iSS) {
        if (iSS == 0) {
            const auto& nodes_ss1 = diagram1.GetS0().nodes_ss;
            const auto& nodes_ss2 = diagram2.GetS0().nodes_ss;
            int t_max2 = nodes_ss2.front().rbegin()->first.t;
            for (auto& [deg, _] : nodes_ss1.front()) {
                if (deg.t > t_max2)
                    break;
                const auto& sc1 = diagram1.GetRecentStaircase(nodes_ss1, deg);
                for (size_t i = 0; i < sc1.levels.size(); ++i) {
                    if (sc1.levels[i] > LEVEL_MAX / 2) {
                        int r = LEVEL_MAX - sc1.levels[i];
                        if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                            diff = NULL_DIFF;
                        else
                            diff = sc1.diffs[i];
                        if (diff != NULL_DIFF) {
                            int count = diagram2.SetS0DiffGlobal(deg, sc1.basis[i], diff, r);
                            if (count > 0)
                                Logger::LogDiff(0, enumReason::migrate, "S0", deg, sc1.basis[i], diff, r);
                        }
                        else {
                            int count = diagram2.SetS0DiffGlobal(deg, sc1.basis[i], int1d{}, r - 1);
                            if (count > 0)
                                Logger::LogDiff(0, enumReason::migrate, "S0", deg, sc1.basis[i], int1d{}, r - 1);
                        }
                    }
                    else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                        int r = sc1.levels[i] + 1;
                        const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                        int count = diagram2.SetS0DiffGlobal(deg_src, {}, sc1.basis[i], r);
                        if (count > 0)
                            Logger::LogDiffInv(0, enumReason::migrate, "S0", deg, {}, sc1.basis[i], r);
                    }
                }
            }
        }
        else {
            size_t iCof = iSS - 1;
            const auto& nodes_ss1 = diagram1.GetCofs()[iCof].nodes_ss;
            const auto& nodes_ss2 = diagram2.GetCofs()[iCof].nodes_ss;
            const auto& name = diagram2.GetCofs()[iCof].name;
            int t_max2 = nodes_ss2.front().rbegin()->first.t;
            for (auto& [deg, _] : nodes_ss1.front()) {
                if (deg.t > t_max2)
                    break;
                const auto& sc1 = diagram1.GetRecentStaircase(nodes_ss1, deg);
                for (size_t i = 0; i < sc1.levels.size(); ++i) {
                    if (sc1.levels[i] > LEVEL_MAX / 2) {
                        int r = LEVEL_MAX - sc1.levels[i];
                        if (deg.t + r - 1 > t_max2 && sc1.diffs[i] != NULL_DIFF)
                            diff = NULL_DIFF;
                        else
                            diff = sc1.diffs[i];
                        if (diff != NULL_DIFF) {
                            int count = diagram2.SetCofDiffGlobal(iCof, deg, sc1.basis[i], diff, r);
                            if (count > 0)
                                Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], diff, r);
                        }
                        else {
                            int count = diagram2.SetCofDiffGlobal(iCof, deg, sc1.basis[i], int1d{}, r - 1);
                            if (count > 0)
                                Logger::LogDiff(0, enumReason::migrate, name, deg, sc1.basis[i], int1d{}, r - 1);
                        }
                    }
                    else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                        int r = sc1.levels[i] + 1;
                        const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                        int count = diagram2.SetCofDiffGlobal(iCof, deg_src, {}, sc1.basis[i], r);
                        if (count > 0)
                            Logger::LogDiffInv(0, enumReason::migrate, name, deg, {}, sc1.basis[i], r);
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
    for (size_t iSS = 0; iSS < diagram1.GetAllBasisSs().size(); ++iSS) {
        if (iSS == 0) {
            auto& pi_gb1 = diagram1.GetS0().pi_gb;
            auto& pi_gb2 = diagram2.GetS0().pi_gb;
            int t_max2 = diagram2.GetS0().nodes_ss.front().rbegin()->first.t;
            if (pi_gb2.gen_degs().size() > 1) {  //// TODO: Support merge of homotopy information
                Logger::LogException(0, 0x43dbbfcaU, "Diagram2 should be reset first.\n");
                throw MyException(0x43dbbfcaU, "Diagram2 should be reset first.");
            }
            algZ::Poly1d rels = pi_gb1.data();
            uint32_t i = 0;
            for (auto& deg : pi_gb1.gen_degs()) {
                if (deg.t <= t_max2) {
                    pi_gb2.AddGen(deg);
                    if (deg.stem() % 2 == 1)
                        rels.push_back(algZ::Mon::two_x_square(i, deg.s));
                }
                ++i;
            }
            diagram2.GetS0().pi_gen_Einf = diagram1.GetS0().pi_gen_Einf;
            diagram2.GetS0().pi_gen_defs = diagram1.GetS0().pi_gen_defs;
            diagram2.GetS0().pi_gen_def_mons = diagram1.GetS0().pi_gen_def_mons;
            diagram2.AddPiRelsS0(std::move(rels));
        }
        else {
            auto& ssCof1 = diagram1.GetCofs()[iSS - 1];
            auto& ssCof2 = diagram2.GetCofs()[iSS - 1];
            auto& pi_gb1 = ssCof1.pi_gb;
            auto& pi_gb2 = ssCof2.pi_gb;
            int t_max2 = ssCof2.nodes_ss.front().rbegin()->first.t;
            if (pi_gb2.v_degs().size() > 1) {  //// TODO: Support merge of homotopy information
                Logger::LogException(0, 0x92cb2691, "Diagram2 should be reset first.\n");
                throw MyException(0x92cb2691, "Diagram2 should be reset first.");
            }
            for (auto& deg : pi_gb1.v_degs())
                if (deg.t <= t_max2)
                    pi_gb2.AddGen(deg);
            ssCof2.pi_gen_Einf = ssCof1.pi_gen_Einf;
            ssCof2.pi_gen_defs = ssCof1.pi_gen_defs;
            ssCof2.pi_gen_def_mons = ssCof1.pi_gen_def_mons;
            ssCof2.pi_qt = ssCof1.pi_qt;
            diagram2.AddPiRelsCof(iSS - 1, pi_gb1.data());
        }
    }
    int count_ss = 0, count_htpy = 0;
    diagram2.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_htpy, 0);
}

int main_migrate_ss(int argc, char** argv, int index)
{
    std::string group1, group2;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss migrate_ss <diagram1> <diagram2>\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "diagram1", group1))
        return index;
    if (myio::load_arg(argc, argv, ++index, "diagram2", group2))
        return index;
    auto dbnames1 = GetDbNames(group1, false);
    auto dbnames2 = GetDbNames(group2);

    auto flag_no_op = DeduceFlag::no_op;
    Diagram diagram1(dbnames1, flag_no_op);
    Diagram diagram2(dbnames2, flag_no_op);

    try {
        Migrate_ss(diagram1, diagram2);
        diagram2.save(dbnames2, flag_no_op);
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

int main_migrate_htpy(int argc, char** argv, int index)
{
    std::string group1, group2;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss migrate_hpty <diagram1> <diagram2>\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "diagram1", group1))
        return index;
    if (myio::load_arg(argc, argv, ++index, "diagram2", group2))
        return index;
    auto dbnames1 = GetDbNames(group1, false);
    auto dbnames2 = GetDbNames(group2);

    DeduceFlag flag = DeduceFlag::homotopy | DeduceFlag::homotopy_def;
    Diagram diagram1(dbnames1, flag);
    Diagram diagram2(dbnames2, flag);

    try {
        Migrate_htpy(diagram1, diagram2);
        diagram2.save(dbnames2, flag);
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

int main_resetfrom(int argc, char** argv, int index)
{
    std::string selector;
    std::string selector_from;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Initialize the data from another diagram\n";
        std::cout << "Usage:\n  ss resetfrom <selector> <selector_from>\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";
        std::cout << "  selector_from = " << selector_from << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "selector", selector))
        return index;
    if (myio::load_arg(argc, argv, ++index, "selector_from", selector_from))
        return index;
    auto dbnames = GetDbNames(selector);
    auto dbnames_from = GetDbNames(selector_from, false);

    for (size_t k = 0; k < dbnames.size(); ++k) {
        std::ifstream src(dbnames_from[k], std::ios::binary);
        std::ofstream dst(dbnames[k], std::ios::binary);
        dst << src.rdbuf();
    }

    return 0;
}

int main_truncate(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss truncate [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    for (size_t k = 0; k < dbnames.size(); ++k) {
        DBSS db(dbnames[k]);
        auto table = GetE2TablePrefix(dbnames[k]);
        auto nodes_ss = db.load_basis_ss(table);
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
        db.update_basis_ss(table, nodes_ss);
        db.end_transaction();
    }

    return 0;
}