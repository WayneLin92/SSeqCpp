#include "algebras/linalg.h"
#include "main.h"
#include <set>

#ifndef MYDEPLOY
std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
#endif

/* If n = 2^k1 + ... + 2^kn,
 * return the array k1, ..., kn. */
int1d two_expansion(unsigned n)
{
    int1d result;
    int k = 0;
    while (n > 0) {
        if (n & 1)
            result.push_back(k);
        n >>= 1;
        ++k;
    }
    return result;
}

/* Deduce zero differentials for degree reason
 * t_test is for DeduceDiffs() */
int Diagram::DeduceZeroDiffs()
{
    int count_new_diffs = 0;
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& basis_ss = *all_basis_ss_[k];
        int t_max = all_t_max_[k];
        for (auto& [d, basis_ss_d] : basis_ss.front()) {
            const Staircase& sc = GetRecentStaircase(basis_ss, d);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.diffs_ind[i] == int1d{-1} && sc.levels[i] > kLevelPC) {
                    int r = kLevelMax - sc.levels[i];
                    /* Find the first possible d_{r1} target for r1>=r */
                    int r1 = NextRTgt(basis_ss, t_max, d, r);
                    if (r != r1) {
                        SetDiffLeibnizV2(k, d, sc.basis_ind[i], {}, r);
                        ++count_new_diffs;
                    }
                }
            }
        }
    }
    return count_new_diffs;
}

// TODO: keep only one of (nd1->nd2), (nd2->nd1) in the deduction tree
int Diagram::DeduceDiffs(int r_max, int maxPoss, int top_depth, int depth, Timer& timer)
{
    if (timer.timeout())
        return 0;
    int old_count = 0, count = 0;

    DeduceZeroDiffs();
    CacheNullDiffs(maxPoss);

    size_t k = 0;
    while (k < all_basis_ss_.size()) {
        auto& basis_ss = *all_basis_ss_[k];
        auto& nd = *all_nd_[k];

        size_t index_nd = 0;
        while (index_nd < nd.back().size()) {
            if (depth == top_depth)
                std::cout << "k=" << k << " index_nd = " << index_nd << '/' << nd.back().size() << "                    \r ";
            const NullDiff ndi = nd.back()[index_nd];
            const Staircase& sc = GetRecentStaircase(basis_ss, ndi.deg);
            /* Fixed source, find target. */
            if (sc.levels[ndi.index] > kLevelMax / 2) {
                const int r = kLevelMax - sc.levels[ndi.index];
                const AdamsDeg deg_tgt = ndi.deg + AdamsDeg{r, r - 1};
                const int1d src = sc.basis_ind[ndi.index];

                int1d tgt, tgt_pass;
                int count_pass = 0;
                unsigned i_max = 1 << ndi.count;
                for (unsigned i = 0; i < i_max; ++i) {
                    tgt.clear();
                    if (i) {
                        const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
                        bool bBreak = false;
                        for (int j : two_expansion(i)) {
                            if (sc_tgt.levels[(size_t)(ndi.first + j)] == 5000) { /* skip permanent cycles */
                                bBreak = true;
                                break;
                            }
                            tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)(ndi.first + j)]);
                        }
                        if (bBreak)
                            continue;
                    }

                    AddNode(); /* Warning: the reallocation invalidates the reference sc above */
                    bool bException = false;
                    try {
                        int nChanges = SetDiffLeibnizV2(k, ndi.deg, src, tgt, r);
                        if (!tgt.empty() || r <= r_max) {
                            if (depth > 1) {
                                int nDeductions = DeduceDiffs(r_max, maxPoss, top_depth, 1, timer);
                                if (depth > 2 && nChanges >= 20 && nDeductions > 0) {
                                    DeduceDiffs(r_max, maxPoss, top_depth, depth - 1, timer);
                                }
                            }
                        }
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        tgt_pass = std::move(tgt);
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                    }
                }
                if (count_pass == 0)
                    throw SSException(0x66d5bd9aU, "No compatible differentials");
                else if (count_pass == 1) {
                    ++count;
                    int nChanges = SetDiffLeibnizV2(k, ndi.deg, src, tgt_pass, r);
                    DeduceZeroDiffs();  // TODO: Implement this after CacheNullDiffs
                    CacheNullDiffs(maxPoss);
                    if (depth == top_depth)
                        std::clog << "k=" << k << " (" << ndi.deg.t - ndi.deg.s << ',' << ndi.deg.s << ") d_{" << r << '}' << src << '=' << tgt_pass << "          \n";
                    if (timer.timeout()) {
                        k = 1000;
                        break;
                    }
                    if (depth > 1)
                        count += DeduceDiffs(r_max, maxPoss, top_depth == depth ? 1 : top_depth, 1, timer);
                }
                else
                    ++index_nd;
            }
            /* Fixed target, find source. */
            else {
                const int r = sc.levels[ndi.index];
                const AdamsDeg deg_src = ndi.deg - AdamsDeg{r, r - 1};
                const int1d tgt = sc.basis_ind[ndi.index];

                int1d src, src_pass;
                int count_pass = 0;
                unsigned i_max = 1 << ndi.count;
                for (unsigned i = 0; i < i_max; ++i) {
                    src.clear();
                    if (i) {
                        const Staircase& sc_src = GetRecentStaircase(basis_ss, deg_src);
                        for (int j : two_expansion(i))
                            src = lina::AddVectors(src, sc_src.basis_ind[(size_t)(ndi.first + j)]);
                    }

                    AddNode(); /* Warning: the reallocation invalidates the reference sc above */
                    bool bException = false;
                    try {
                        int nChanges = SetDiffLeibnizV2(k, deg_src, src, tgt, r);
                        if (!tgt.empty() || r <= r_max) {
                            if (depth > 1) {
                                int nDeductions = DeduceDiffs(r_max, maxPoss, top_depth, 1, timer);
                                if (depth > 2 && nChanges >= 20 && nDeductions > 0) {
                                    DeduceDiffs(r_max, maxPoss, top_depth, depth - 1, timer);
                                }
                            }
                        }
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        src_pass = std::move(src);
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                    }
                }
                if (count_pass == 0)
                    throw SSException(0x66d5bd9aU, "No compatible differentials");
                else if (count_pass == 1) {
                    ++count;
                    int nChanges = SetDiffLeibnizV2(k, deg_src, src_pass, tgt, r);
                    DeduceZeroDiffs();  // TODO: Implement this after CacheNullDiffs
                    CacheNullDiffs(maxPoss);
                    if (depth == top_depth)
                        std::clog << "k=" << k << " (" << ndi.deg.t - ndi.deg.s + 1 << ',' << ndi.deg.s - r << ") d_{" << r << '}' << src_pass << '=' << tgt << "          \n";
                    if (timer.timeout()) {
                        k = 1000;
                        break;
                    }
                    if (depth > 1)
                        count += DeduceDiffs(r_max, maxPoss, top_depth == depth ? 1 : top_depth, 1, timer);
                }
                else
                    ++index_nd;
            }
        }
        if ((++k) == all_basis_ss_.size()) {
            if (old_count != count) {
                k = 0;
                old_count = count;
            }
        }
    }
    return count;
}

/* Deduce zero differentials using the image of J
 * t_test is for DeduceDiffs() */
int Diagram::DeduceImageJ()
{
    int count_new_diffs = 0;
    int prev_t = -1;
    AdamsDeg1d arr_degs_image_of_J = {AdamsDeg(1, 2), AdamsDeg(2, 4), AdamsDeg(3, 6), AdamsDeg(1, 4), AdamsDeg(2, 5), AdamsDeg(2, 9), AdamsDeg(3, 10), AdamsDeg(4, 11), AdamsDeg(3, 11), AdamsDeg(4, 13)};
    std::set<AdamsDeg> degs_image_of_J;
    for (int i = 0; i < 30; ++i)
        for (size_t j = 0; j < arr_degs_image_of_J.size(); ++j)
            degs_image_of_J.insert(arr_degs_image_of_J[j] + AdamsDeg(4, 12) * i);
    for (auto& [d, basis_ss_d] : ssS0_.basis_ss.front()) {
        if (d.t != prev_t) {
            std::cout << d.t << " " << count_new_diffs << "     \r";
            prev_t = d.t;
        }
        const Staircase& sc = GetRecentStaircase(ssS0_.basis_ss, d);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.diffs_ind[i] == int1d{-1} && sc.levels[i] > kLevelPC) {
                int r = kLevelMax - sc.levels[i];
                AdamsDeg deg_tgt = d + AdamsDeg(r, r - 1);
                if (degs_image_of_J.find(deg_tgt) != degs_image_of_J.end()) {
                    SetS0DiffLeibnizV2(d, sc.basis_ind[i], {}, kLevelPC - 1);
                    ++count_new_diffs;
                }
            }
        }
    }
    return count_new_diffs;
}

int main_deduce_zero(int argc, char** argv, int index)
{
    std::string db_S0 = DB_DEFAULT;
    std::vector<std::string> dbnames = {
        "C2_AdamsSS_t200.db",
        "Ceta_AdamsSS_t200.db",
        "Cnu_AdamsSS_t200.db",
        "Csigma_AdamsSS_t200.db",
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

    bench::Timer timer;
    Diagram diagram(dbnames);

    int count = 0;
    try {
        count = diagram.DeduceZeroDiffs();

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetTablePrefix(dbnames[k]), diagram.GetChanges(k));
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

int main_deduce_j(int argc, char** argv, int index)
{
    std::string db_filename = DB_DEFAULT;
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce trivial differentials because of image of J\n";
        std::cout << "Usage:\n  ss deduce j [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-12)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    bench::Timer timer;
    DBSS db(db_filename);
    Diagram ss({db_filename});

    int count = 0;
    try {
        count = ss.DeduceImageJ();

        db.begin_transaction();
        db.update_basis_ss(table_prefix, ss.GetChanges(0));
        db.end_transaction();
    }
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }

    std::cout << "Changed differentials: " << count << '\n';
    return 0;
}

int main_deduce_diff(int argc, char** argv, int index)
{
    int r_max = 10, maxPoss = 10, depth = 1;
    double stop_time = 600;
    size_t kInput = 0;
    std::string db_S0 = DB_DEFAULT;
    std::vector<std::string> dbnames = {
        "C2_AdamsSS_t200.db",
        "Ceta_AdamsSS_t200.db",
        "Cnu_AdamsSS_t200.db",
        "Csigma_AdamsSS_t200.db",
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce differentials by Leibniz rule\n";
        std::cout << "Usage:\n  ss deduce diff [r_max] [maxPoss] [depth] [stop_time] [db_S0] [db_Cofibs ...]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  r_max = " << r_max << "\n";
        std::cout << "  maxPoss = " << maxPoss << "\n";
        std::cout << "  depth = " << depth << "\n";
        std::cout << "  stop_time = " << stop_time << "\n";
        std::cout << "  db_S0 = " << db_S0 << "\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "r_max", r_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "maxPoss", maxPoss))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "depth", depth))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stop_time", stop_time))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofib", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);

    bench::Timer timer;
    Diagram diagram(dbnames);

    int count = 0;
    try {
        Timer timer(stop_time);
        count = diagram.DeduceDiffs(r_max, maxPoss, depth, depth, timer);

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetTablePrefix(dbnames[k]), (*diagram.GetAllBasisSs()[k])[1]);
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
    std::cout << "Done" << std::endl;
    return 0;
}

/* This is for debugging */
int main_deduce_tmp(int argc, char** argv, int index)
{
    std::string db_filename = DB_DEFAULT;
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce tmp [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-19)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    // bench::Timer timer;
    DBSS db(db_filename);
    Diagram ss({db_filename});

    // int count = 0;
    try {
        ss.CacheNullDiffs(10);
        for (size_t i = 0; i < ss.GetS0().nd.back().size(); ++i)
            auto& nd = ss.GetS0().nd.back()[i];
        std::cout << "i=" << ss.GetFirstIndexOfFixedLevels(ss.GetS0().basis_ss, AdamsDeg(13, 47 + 13), 9994) << std::endl;
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

    // std::cout << "Changed differentials: " << count << '\n';
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
    else if (cmd == "j")
        return main_deduce_j(argc, argv, index);
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
