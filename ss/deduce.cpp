#include "algebras/linalg.h"
#include "main.h"
#include <set>

#ifndef MYDEPLOY
std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
#endif

/* Deduce zero differentials for degree reason
 * t_test is for DeduceDiffs() */
int Diagram::DeduceTrivialDiffs()
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

int Diagram::TryDiff(size_t iSS, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, DeduceFlag flag)
{
    AddNode();
    bool bException = false;
    try {
        SetDiffLeibnizV2(iSS, deg_x, x, dx, r, flag & DeduceFlag::fast_try_diff);
        /*if (flag & DeduceFlag::set_diff)
            DeduceDiffs(depth + 1, 0);*/
        if (flag & DeduceFlag::homotopy) {
            int count_ss1 = 0, count_homotopy1 = 0;
            AdamsDeg deg_min = deg_x - AdamsDeg(0, 1);
            if (iSS > 0)
                deg_min = deg_min - ssCofs_[iSS - 1].deg_qt;
            SyncHomotopy(deg_min, count_ss1, count_homotopy1, depth + 1);
            DeduceTrivialExtensions(depth + 1);
            if (flag & DeduceFlag::check_exactness)
                DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth + 1);
        }
    }
    catch (SSException&) {
        //std::cout << iSS << ' ' << deg_x.StrAdams() << ' ' << x << ' ' << r << ' ' << dx << '\n';
        //std::cout << std::hex << e.id() << '\n';
        bException = true;
    }
    PopNode();

    if (bException)
        return 1;
    else
        return 0;
}

int Diagram::DeduceDiffs(size_t iSS, AdamsDeg deg, int depth, DeduceFlag flag)
{
    auto& basis_ss = *all_basis_ss_[iSS];
    auto& name = all_names_[iSS];

    int count = 0;
    std::string color, color_end = "\033[0m";
    NullDiff1d nds;
    CacheNullDiffs(iSS, deg, flag, nds);

    size_t index_nd = 0;
    while (index_nd < nds.size()) {
        const NullDiff nd = nds[index_nd];
        int1d x, dx;
        int r;
        bool bNewDiff = false;
        /* Fixed source, find target. */
        AdamsDeg deg_src;
        if (nd.r > 0) {
            r = nd.r;
            deg_src = deg;
            const AdamsDeg deg_tgt = deg_src + AdamsDeg{r, r - 1};
            x = nd.x;

            int1d dx1;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;
            for (unsigned i = 1; i < i_max; ++i) {
                const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
                dx1.clear();
                for (int j : two_expansion(i))
                    dx1 = lina::AddVectors(dx1, sc_tgt.basis_ind[(size_t)(nd.first + j)]);

                if (!TryDiff(iSS, deg_src, x, dx1, r, depth, flag)) {
                    ++count_pass;
                    if (count_pass > 1)
                        break;
                    dx = std::move(dx1);
                }
            }
            if (count_pass == 0) {
                dx.clear();
                bNewDiff = true;
            }
            else if (count_pass > 1)
                ++index_nd;
            else {
                dx1.clear();
                if (TryDiff(iSS, deg_src, x, dx1, r, depth, flag))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }
        /* Fixed target, find source. */
        else {
            r = -nd.r;
            deg_src = deg - AdamsDeg{r, r - 1};
            const Staircase& sc_src = GetRecentStaircase(basis_ss, deg_src);
            dx = nd.x;

            int1d x1;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;
            for (unsigned i = 1; i < i_max; ++i) {
                x1.clear();
                for (int j : two_expansion(i))
                    x1 = lina::AddVectors(x1, sc_src.basis_ind[(size_t)(nd.first + j)]);

                if (!TryDiff(iSS, deg_src, x1, dx, r, depth, flag)) {
                    ++count_pass;
                    if (count_pass > 1)
                        break;
                    x = std::move(x1);
                }
            }
            if (count_pass == 0) {
                x.clear();
                bNewDiff = true;
            }
            else if (count_pass > 1)
                ++index_nd;
            else {
                x1.clear();
                if (TryDiff(iSS, deg_src, x1, dx, r, depth, flag))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }

        if (bNewDiff) {
            ++count;
            SetDiffLeibnizV2(iSS, deg_src, x, dx, r);
            if (depth == 0) {
                if (deg.stem() <= 127) {
                    if (x.empty() || dx.empty())
                        color = "\033[38;2;255;150;150m";
                    else
                        color = "\033[38;2;255;100;100m";
                }
                else
                    color = "";
                if (nd.r > 0)
                    std::cout << color << name << "  " << deg.StrAdams() << "  d_{" << r << '}' << x << '=' << dx << "                  " + color_end + '\n';
                else
                    std::cout << color << name << "  " << deg.StrAdams() << "  " << dx << "=d_{" << r << '}' << x << "                  " + color_end + '\n';
            }
            int count_trivial = DeduceTrivialDiffs();
            count += count_trivial;
            CacheNullDiffs(iSS, deg, flag, nds);
            if (flag & DeduceFlag::homotopy) {
                int count_homotopy1 = 0;
                AdamsDeg deg_min = deg_src - AdamsDeg(0, 1);
                if (iSS > 0)
                    deg_min = deg_min - ssCofs_[iSS - 1].deg_qt;
                if (count_trivial)
                    deg_min.s = 0;
                SyncHomotopy(deg_min, count, count_homotopy1, depth);
                DeduceTrivialExtensions(depth);
                if (flag & DeduceFlag::check_exactness)
                    DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth);
            }
        }
    }
    return count;
}

int Diagram::DeduceDiffs(int stem_min, int stem_max, int depth, DeduceFlag flag)
{
    int count = 0;

    DeduceTrivialDiffs();
    if (flag & DeduceFlag::homotopy) {
        int count_homotopy1 = 0;
        SyncHomotopy(AdamsDeg(0, 0), count, count_homotopy1, depth + 1);
        DeduceTrivialExtensions(depth + 1);
        if (flag & DeduceFlag::check_exactness)
            DeduceExtensionsByExactness(0, stem_max_exactness_, depth + 1);
    }

    for (size_t iSS = 0; iSS < all_basis_ss_.size(); ++iSS) {
        // auto& basis_ss = *all_basis_ss_[iSS];
        auto& name = all_names_[iSS];
        auto& degs = [&]() {
            if (iSS == 0)
                return ssS0_.degs_basis_order_by_stem;
            else
                return ssCofs_[iSS - 1].degs_basis_order_by_stem;
        }();

        for (AdamsDeg deg : degs) {
            if (depth == 0)
                std::cout << name << "  deg=" << deg.StrAdams() << "                        \r";
            if (deg.stem() < stem_min)
                continue;
            else if (deg.stem() > stem_max)
                break;

            count += DeduceDiffs(iSS, deg, depth, flag);
        }
    }
    return count;
}

int main_deduce_diff(int argc, char** argv, int index)
{
    int stem_min = 0, stem_max = 261;
    std::string selector = "default";
    std::vector<std::string> strFlags;
    DeduceFlag flag = DeduceFlag::no_op;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce differentials by Leibniz rule\n";
        std::cout << "Usage:\n  ss deduce diff [stem_min] [stem_max] [selector] [flags...]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  stem_min = " << stem_min << "\n";
        std::cout << "  stem_max = " << stem_max << "\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "stem_min", stem_min))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    if (myio::load_args(argc, argv, ++index, "flags", strFlags))
        return index;
    auto dbnames = GetDbNames(selector);
    for (auto& f : strFlags) {
        if (f == "all_x")
            flag = flag | DeduceFlag::all_x;
        else if (f == "homotopy")
            flag = flag | DeduceFlag::homotopy;
        else if (f == "exact")
            flag = flag | DeduceFlag::check_exactness;
        else {
            std::cout << "Not a supported flag: " << f << '\n';
            return 100;
        }
    }

    Diagram diagram(dbnames);

    int count = 0;
    try {
        if (flag & DeduceFlag::homotopy) {
            int count_homotopy1 = 0;
            diagram.SyncHomotopy(AdamsDeg(0, 0), count, count_homotopy1, 0);
            diagram.DeduceTrivialExtensions(0);
            if (flag & DeduceFlag::check_exactness)
                diagram.DeduceExtensionsByExactness(0, diagram.stem_max_exactness_, 0);
        }
        count = diagram.DeduceDiffs(stem_min, stem_max, 0, flag);
        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetE2TablePrefix(dbnames[k]), (*diagram.GetAllBasisSs()[k])[1]);

            if (flag & DeduceFlag::homotopy) {
                auto pi_table = GetComplexName(dbnames[k]);
                db.drop_and_create_pi_relations(pi_table);
                db.drop_and_create_pi_basis(pi_table);

                if (k == 0) {
                    db.drop_and_create_pi_generators(pi_table);
                    db.save_pi_generators(pi_table, diagram.GetS0().pi_gb.gen_degs(), diagram.GetS0().pi_gen_Einf);
                    db.save_pi_gb(pi_table, diagram.GetS0().pi_gb.OutputForDatabase(), diagram.GetS0GbEinf());
                    db.save_pi_basis(pi_table, diagram.GetS0().pi_basis.front());
                }
                else {
                    db.drop_and_create_pi_generators_mod(pi_table);
                    auto& Cof = diagram.GetCofs()[k - 1];
                    if (Cof.pi_qt.size() != 1)
                        throw MyException(0x925afecU, "Not on the initial node");
                    db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_qt.front());
                    db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), diagram.GetCofGbEinf(int(k - 1)));
                    db.save_pi_basis_mod(pi_table, Cof.pi_basis.front());
                }
            }

            db.end_transaction();
        }
        myio::Logger::dout() << "Changed differentials: " << count << '\n';
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

    return 0;
}

/* This is for debugging */
int main_deduce_tmp(int argc, char** argv, int index)
{
    std::string selector = "default";
    DeduceFlag flag = DeduceFlag::no_op;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce tmp [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    Diagram diagram(dbnames);

    int count = 0;
    try {
        int count_ss = 0, count_homotopy = 0;
        diagram.DeduceExtensions(0, 30, count_ss, count_homotopy, 0, flag);
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
                db.save_pi_basis(pi_table, diagram.GetS0().pi_basis.front());
            }
            else {
                db.drop_and_create_pi_generators_mod(pi_table);
                auto& Cof = diagram.GetCofs()[k - 1];
                if (Cof.pi_qt.size() != 1)
                    throw MyException(0x925afecU, "Not on the initial node");
                db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_qt.front());
                db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), diagram.GetCofGbEinf(int(k - 1)));
                db.save_pi_basis_mod(pi_table, Cof.pi_basis.front());
            }

            db.end_transaction();
        }
        std::cout << "Changed differentials: " << count << '\n';
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

    // bench::Counter::print();
    return 0;
}

/* generate the table of the spectral sequence */
int main_deduce(int argc, char** argv, int index)
{
    std::string cmd;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss deduce_diffs <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  diff: deduce differentials\n";
        std::cout << "  ext: deduce extensions\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "diff")
        return main_deduce_diff(argc, argv, index);
    else if (cmd == "ext")
        return main_deduce_ext(argc, argv, index);
    else if (cmd == "tmp")
        return main_deduce_tmp(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << cmd << '\n';
    return 0;
}
