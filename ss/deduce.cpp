#include "algebras/linalg.h"
#include "main.h"
#include <set>

#ifndef MYDEPLOY
std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
#endif

/* Deduce zero differentials for degree reason
 * t_test is for DeduceDiffs() */
int Diagram::DeduceZeroDiffs()
{
    int old_count = 0, count_new_diffs = 0;
    while (true) {
        for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
            auto& basis_ss = *all_basis_ss_[k];
            int t_max = all_t_max_[k];
            for (auto& [d, basis_ss_d] : basis_ss.front()) {
                const Staircase& sc = GetRecentStaircase(basis_ss, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs_ind[i] == int1d{-1}) {
                        if (sc.levels[i] > kLevelPC) {
                            const int r = kLevelMax - sc.levels[i];
                            /* Find the first possible d_{r1} target for r1>=r */
                            int r1 = NextRTgt(basis_ss, t_max, d, r);
                            if (r != r1) {
                                SetDiffLeibnizV2(k, d, sc.basis_ind[i], {}, r);
                                ++count_new_diffs;
                            }
                        }
                        else if (sc.levels[i] < kLevelMax / 2) {
                            const int r = sc.levels[i];
                            int r1 = NextRSrc(basis_ss, d, r);
                            if (r != r1) {
                                SetDiffLeibnizV2(k, d - AdamsDeg(r, r - 1), {}, sc.basis_ind[i], r);
                                ++count_new_diffs;
                            }
                        }
                    }
                }
            }
        }
        if (old_count != count_new_diffs)
            old_count = count_new_diffs;
        else
            break;
    }
    return count_new_diffs;
}

/*
 * Stage 0: Fast nd cache, Try{Leibniz}, maxPoss = 4
 * Stage 1: Fast nd cache, Try{Leibniz}, maxPoss = 10
 * Stage 2: Full nd cache, Try{Leibniz}, maxPoss = 8
 * Stage 3: Fast nd cache, Try{Leibniz, {Stage<=0, epoch<=0, fast Leibniz}}, maxPoss = 4
 */
int Diagram::DeduceDiffs(int depth, int max_stage, Timer& timer)  // TODO: cache diff to add
{
    if (timer.timeout())
        return 0;
    int stage = depth == 0;
    int epoch = 0;
    size_t iSS = 0;

    int maxPoss = 4;
    int maxStem = depth == 0 ? (stage == 3 ? 127 : DEG_MAX) : 127;

    DeduceZeroDiffs();
    CacheNullDiffs(maxPoss, maxStem, stage == 2);

    std::string color, color_end = "\033[0m";  ////
    int old_count = 0, count = 0;
    while (iSS < all_basis_ss_.size()) {
        auto& basis_ss = *all_basis_ss_[iSS];
        auto& nds = *all_nd_[iSS];

        size_t index_nd = 0;
        while (index_nd < nds.back().size()) {
            if (depth == 0)
                std::cout << "stage=" << stage << " epoch=" << epoch << " iSS=" << iSS << " i=" << index_nd << '/' << nds.back().size() << "                     \r";

            const NullDiff nd = nds.back()[index_nd];
            int1d x, dx;
            AdamsDeg deg;
            int r;
            bool bNewDiff = false;
            /* Fixed source, find target. */
            if (nd.r > 0) {
                deg = nd.deg;
                r = nd.r;
                const AdamsDeg deg_tgt = nd.deg + AdamsDeg{r, r - 1};
                x = nd.x;

                int1d dx1;
                int count_pass = 0;
                unsigned i_max = 1 << nd.count;
                for (unsigned i = 0; i < i_max; ++i) {
                    const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
                    dx1.clear();
                    for (int j : two_expansion(i))
                        dx1 = lina::AddVectors(dx1, sc_tgt.basis_ind[(size_t)(nd.first + j)]);

                    AddNode(); /* Warning: the reallocation invalidates the reference sc_tgt above */
                    bool bException = false;
                    try {
                        SetDiffLeibnizV2(iSS, deg, x, dx1, r, depth == 1);
                        if (depth == 0 && stage == 3) {
                            // DeduceDiffs(depth + 1, 0, timer);
                            int count_ss1 = 0, count_homotopy1 = 0;
                            SyncHomotopy(count_ss1, count_homotopy1, depth + 1);
                            DeduceZeroExtensions(depth + 1);
                            DeduceExtensionsByExactness(depth + 1);
                        }
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        dx = std::move(dx1);
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                    }
                }
                if (count_pass == 0)
                    throw SSException(0x66d5bd9aU, "No compatible differentials");
                else if (count_pass == 1) {
                    ++count;
                    bNewDiff = true;
                }
                else
                    ++index_nd;
            }
            /* Fixed target, find source. */
            else {
                r = -nd.r;
                deg = nd.deg - AdamsDeg{r, r - 1};
                const Staircase& sc_src = GetRecentStaircase(basis_ss, deg);
                dx = nd.x;

                int1d x1;
                int count_pass = 0;
                unsigned i_max = 1 << nd.count;
                for (unsigned i = 0; i < i_max; ++i) {
                    x1.clear();
                    for (int j : two_expansion(i))
                        x1 = lina::AddVectors(x1, sc_src.basis_ind[(size_t)(nd.first + j)]);

                    AddNode(); /* Warning: the reallocation invalidates the reference sc_src above */
                    bool bException = false;
                    try {
                        SetDiffLeibnizV2(iSS, deg, x1, dx, r, depth == 1);
                        if (depth == 0 && stage == 3) {
                            // DeduceDiffs(depth + 1, 0, timer);
                            int count_ss1 = 0, count_homotopy1 = 0;
                            SyncHomotopy(count_ss1, count_homotopy1, depth + 1);
                            DeduceZeroExtensions(depth + 1);
                            DeduceExtensionsByExactness(depth + 1);
                        }
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        x = std::move(x1);
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                    }
                }
                if (count_pass == 0)
                    throw SSException(0x66d5bd9aU, "No compatible differentials");
                else if (count_pass == 1) {
                    ++count;
                    bNewDiff = true;
                }
                else
                    ++index_nd;
            }

            if (bNewDiff) {
                SetDiffLeibnizV2(iSS, deg, x, dx, r);
                if (depth == 0) {
                    if (nd.deg.stem() <= 127) {
                        if (x.empty() || dx.empty())
                            color = "\033[38;2;200;64;64m";
                        else
                            color = "\033[38;2;255;128;128m";
                    }
                    if (nd.r > 0)
                        std::cout << color + "iSS=" << iSS << "  " << nd.deg.StrAdams() << "  d_{" << r << '}' << x << '=' << dx << "                  " + color_end + '\n';
                    else
                        std::cout << color + "iSS=" << iSS << "  " << nd.deg.StrAdams() << "  " << dx << "=d_{" << r << '}' << x << "                  " + color_end + '\n';
                }
                DeduceZeroDiffs();
                CacheNullDiffs(maxPoss, maxStem, stage == 2);
            }
        }

        if (timer.timeout())
            break;
        if ((++iSS) == all_basis_ss_.size() && depth == 0) {
            if (old_count != count) {
                ++epoch;
                iSS = 0;
                old_count = count;
            }
            else if (stage < max_stage) {
                ++stage;
                epoch = 0;
                iSS = 0;
                if (stage == 1) {
                    maxPoss = 10;
                    CacheNullDiffs(maxPoss, maxStem, stage == 2);
                }
                else if (stage == 2) {
                    maxPoss = 8;
                    CacheNullDiffs(maxPoss, maxStem, stage == 2);
                }
                else if (stage == 3) {
                    maxPoss = 4;
                    maxStem = 126;
                    CacheNullDiffs(maxPoss, maxStem, stage == 2);
                }
            }
        }
    }
    return count;
}

int main_deduce_zero(int argc, char** argv, int index)
{
    std::string db_S0 = DB_S0;
    std::vector<std::string> dbnames = {
        DB_C2,
        DB_Ceta,
        DB_Cnu,
        DB_Csigma,
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce trivial differentials for degree reason\n";
        std::cout << "Usage:\n  ss deduce zero [db_S0] [db_Cofibs ...]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofibs", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);

    Diagram diagram(dbnames);

    int count = 0;
    try {
        count = diagram.DeduceZeroDiffs();

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetE2TablePrefix(dbnames[k]), diagram.GetChanges(k));
            db.end_transaction();
        }
    }
    /*catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    std::cout << "Changed differentials: " << count << '\n';
    return 0;
}

int main_deduce_diff(int argc, char** argv, int index)
{
    double stop_time = 600;
    std::string db_S0 = DB_S0;
    std::vector<std::string> dbnames = {
        DB_C2,
        DB_Ceta,
        DB_Cnu,
        DB_Csigma,
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce differentials by Leibniz rule\n";
        std::cout << "Usage:\n  ss deduce diff [stop_time] [db_S0] [db_Cofibs ...]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  stop_time = " << stop_time << "\n";
        std::cout << "  db_S0 = " << db_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "stop_time", stop_time))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofib", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);

    Diagram diagram(dbnames);

    int count = 0;
    try {
        Timer timer(stop_time);
        count = diagram.DeduceDiffs(0, 3, timer);
        diagram.ApplyChanges(1);

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetE2TablePrefix(dbnames[k]), (*diagram.GetAllBasisSs()[k])[1]);
            db.end_transaction();
        }
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    std::cout << "Changed differentials: " << count << '\n';
    std::cout << "Done" << std::endl;
    return 0;
}

/* This is for debugging */
int main_deduce_tmp(int argc, char** argv, int index)
{
    std::string db_S0 = DB_S0;
    std::vector<std::string> dbnames = {
        DB_C2,
        DB_Ceta,
        DB_Cnu,
        DB_Csigma,
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce tmp [db_S0] [db_Cofibs ...]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofib", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);

    Diagram diagram(dbnames);

    int count = 0;
    try {
        int count_ss = 0, count_homotopy = 0;
        /*diagram.SyncHomotopy(count_ss, count_homotopy);
        count_homotopy += diagram.DeduceZeroExtensions();
        count_homotopy += diagram.DeduceExtensionsByExactness();*/
        diagram.DeduceExtensions(count_ss, count_homotopy, 0);
        diagram.SimplifyPiRels();

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            auto pi_table = GetComplexName(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(pi_table + "_AdamsE2", (*diagram.GetAllBasisSs()[k])[1]);

            db.drop_and_create_pi_relations(pi_table);
            db.drop_and_create_pi_basis(pi_table);

            if (k == 0) {
                db.drop_and_create_pi_generators(pi_table);
                db.save_pi_generators(pi_table, diagram.GetS0().pi_gb.gen_degs(), diagram.GetS0().pi_gen_Einf);
                db.save_pi_gb(pi_table, diagram.GetS0().pi_gb.OutputForDatabase(), diagram.GetS0GbEinf());
                db.save_pi_basis(pi_table, diagram.GetS0().pi_basis, diagram.GetS0().pi_basis_Einf);
            }
            else {
                db.drop_and_create_pi_generators_mod(pi_table);
                auto& Cof = diagram.GetCofs()[k - 1];
                if (Cof.pi_f_top_cell.size() != 1)
                    throw MyException(0x925afecU, "Not on the initial node");
                db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_f_top_cell.front());
                db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), diagram.GetCofGbEinf(int(k - 1)));
                db.save_pi_basis_mod(pi_table, Cof.pi_basis, Cof.pi_basis_Einf);
            }

            db.end_transaction();
        }
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    std::cout << "Changed differentials: " << count << '\n';
    // bench::Counter::print();
    return 0;
}

/* generate the table of the spectral sequence */
int main_deduce(int argc, char** argv, int index)
{
    std::string cmd;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss deduce <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  zero: deduce trivial differentials for degree reason\n";
        std::cout << "  j: deduce trivial differentials because of image of J\n";
        std::cout << "  diff: deduce differentials by Leibniz rule\n";
        std::cout << "  migrate: migrate ss data from an old database to a brandnew one with bigger range.\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "zero")
        return main_deduce_zero(argc, argv, index);
    else if (cmd == "diff")
        return main_deduce_diff(argc, argv, index);
    else if (cmd == "migrate")
        return main_deduce_migrate(argc, argv, index);
    else if (cmd == "tmp")
        return main_deduce_tmp(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << cmd << '\n';
    return 0;
}
