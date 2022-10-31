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
                            DeduceDiffs(depth + 1, 0, timer);
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
                            DeduceDiffs(depth + 1, 0, timer);
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
                    std::string color, color_end;
                    if (x.empty() || dx.empty()) {
                        if (nd.deg.stem() <= 127) {
                            color = "\033[38;2;156;64;64m";
                            color_end = "\033[0m";
                        }
                        else {
                            color = "\033[38;2;128;128;128m";
                            color_end = "\033[0m";
                        }
                    }
                    else if (nd.deg.stem() <= 127) {
                        color = "\033[38;2;255;128;128m";
                        color_end = "\033[0m";
                    }
                    std::cout << color + "iSS=" << iSS << " (" << nd.deg.stem() << ',' << nd.deg.s << ") d_{" << r << '}' << x << '=' << dx << "                  " + color_end + '\n';
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

int main_deduce_j(int argc, char** argv, int index)
{
    std::string db_filename = DB_S0;
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

    DBSS db(db_S0);
    Diagram diagram(dbnames);

    int count = 0;
    try {
        int count_ss = 0, count_homotopy = 0;
        diagram.SyncHomotopy(count_ss, count_homotopy);
        count_homotopy += diagram.DeduceZeroExtensions();
        diagram.SimplifyPiRels();

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            auto pi_table = GetComplexName(dbnames[k]);
            db.begin_transaction();
            //db.update_basis_ss(table, (*diagram.GetAllBasisSs()[k])[1]);

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
                db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_f_top_cell);
                db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), diagram.GetCofGbEinf(int(k - 1)));
                db.save_pi_basis_mod(pi_table, Cof.pi_basis, Cof.pi_basis_Einf);
            }

            db.end_transaction();
        }

        auto& gen_degs = diagram.GetS0().pi_gb.gen_degs();
        auto& gen_repr = diagram.GetS0().pi_gen_Einf;
        auto& data = diagram.GetS0().pi_gb.data();

        std::cout << "Generators\n";
        for (size_t i = 0; i < gen_degs.size(); ++i) {
            if (gen_degs[i].stem() <= 30)
                std::cout << i << ' ' << gen_degs[i].StrAdams() << ' ' << gen_repr[i] << '\n';
        }
        std::cout << "\nRelations\n";
        for (size_t i = 0; i < data.size(); ++i) {
            if (algZ::GetDeg(data[i].GetLead(), gen_degs).stem() <= 30)
                std::cout << data[i] << '\n';
        }
        std::cout << "\nBasis\n";
        for (auto& [d, basis_d] : diagram.GetS0().pi_basis) {
            if (d.stem() <= 30 && d.stem() > 0 && basis_d.size() > 0) {
                std::cout << d.StrAdams() << ' ';
                for (auto& m : basis_d)
                    std::cout << m << ' ';
                std::cout << '\n';
            }
        }
        auto& v_degs = diagram.GetCofs()[0].pi_gb.v_degs();
        auto& gen_repr1 = diagram.GetCofs()[0].pi_gen_Einf;
        auto& data1 = diagram.GetCofs()[0].pi_gb.data();

        std::cout << "Generators\n";
        for (size_t i = 0; i < v_degs.size(); ++i) {
            if (v_degs[i].stem() <= 30)
                std::cout << i << ' ' << v_degs[i].StrAdams() << "  detected by " << gen_repr1[i] << "  f(this)=" << diagram.GetCofs()[0].pi_f_top_cell[i] << '\n';
        }
        std::cout << "\nRelations\n";
        for (size_t i = 0; i < data1.size(); ++i) {
            if (algZ::GetDeg(data1[i].GetLead(), gen_degs, v_degs).stem() <= 30)
                std::cout << data1[i] << '\n';
        }
        std::cout << "\nBasis\n";
        for (auto& [d, basis_d] : diagram.GetCofs()[0].pi_basis) {
            if (d.stem() <= 30 && d.stem() > 0 && basis_d.size() > 0) {
                std::cout << d.StrAdams() << ' ';
                for (auto& m : basis_d)
                    std::cout << m << ' ';
                std::cout << '\n';
            }
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
    bench::Counter::print();
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
