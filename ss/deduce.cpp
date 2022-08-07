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
int SS::DeduceZeroDiffs()
{
    int count_new_diffs = 0;
    int prev_t = -1;
    for (auto& [d, basis_ss_d] : basis_ss_.front()) {
        if (d.t != prev_t) {
            std::cout << d.t << " " << count_new_diffs << "     \r";
            prev_t = d.t;
        }
        const Staircase& sc = GetRecentStaircase(d);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.diffs_ind[i] == int1d{-1} && sc.levels[i] > kLevelPC) {
                int r = kLevelMax - sc.levels[i];
                /* Find the first possible d_{r1} target for r1>=r */
                int r1 = NextRTgt(d, r);  //
                if (r != r1) {
                    SetDiffLeibnizV2(d, sc.basis_ind[i], {}, r);
                    ++count_new_diffs;
                }
            }
        }
    }
    return count_new_diffs;
}

/* Deduce zero differentials using the image of J
 * t_test is for DeduceDiffs() */
int SS::DeduceImageJ()
{
    int count_new_diffs = 0;
    int prev_t = -1;
    AdamsDeg1d arr_degs_image_of_J = {AdamsDeg(1, 2), AdamsDeg(2, 4), AdamsDeg(3, 6), AdamsDeg(1, 4), AdamsDeg(2, 5), AdamsDeg(2, 9), AdamsDeg(3, 10), AdamsDeg(4, 11), AdamsDeg(3, 11), AdamsDeg(4, 13)};
    std::set<AdamsDeg> degs_image_of_J;
    for (int i = 0; i < 30; ++i)
        for (size_t j = 0; j < arr_degs_image_of_J.size(); ++j)
            degs_image_of_J.insert(arr_degs_image_of_J[j] + AdamsDeg(4, 12) * i);
    for (auto& [d, basis_ss_d] : basis_ss_.front()) {
        if (d.t != prev_t) {
            std::cout << d.t << " " << count_new_diffs << "     \r";
            prev_t = d.t;
        }
        const Staircase& sc = GetRecentStaircase(d);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.diffs_ind[i] == int1d{-1} && sc.levels[i] > kLevelPC) {
                int r = kLevelMax - sc.levels[i];
                AdamsDeg deg_tgt = d + AdamsDeg(r, r - 1);
                if (degs_image_of_J.find(deg_tgt) != degs_image_of_J.end()) {
                    SetDiffLeibnizV2(d, sc.basis_ind[i], {}, kLevelPC - 1);
                    ++count_new_diffs;
                }
            }
        }
    }
    return count_new_diffs;
}

// TODO: keep only one of (nd1->nd2), (nd2->nd1) in the deduction tree
int SS::DeduceDiffs(int r_max, int maxPoss, int top_depth, int depth, Timer& timer)
{
    if (timer.timeout())
        return 0;
    int old_count = 0, count = 0;

    DeduceZeroDiffs();
    CacheNullDiffs(maxPoss);
    size_t index_nd = 0;
    while (nd_.back().size() > 0 && (old_count != count || index_nd < nd_.back().size())) {
        if (index_nd >= nd_.back().size()) {
            index_nd = 0;
            old_count = count;
        }

        if (depth == top_depth)
            std::cout << "index_nd = " << index_nd << '/' << nd_.back().size() << "                    \r";
        const NullDiff nd = nd_.back()[index_nd];
        const Staircase& sc = GetRecentStaircase(nd.deg);
        if (sc.levels[nd.index] > kLevelMax / 2) {
            const int r = kLevelMax - sc.levels[nd.index];
            const AdamsDeg deg_tgt = nd.deg + AdamsDeg{r, r - 1};
            const int1d src = sc.basis_ind[nd.index];

            int1d tgt, tgt_pass;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;
            for (unsigned i = 0; i < i_max; ++i) {
                tgt.clear();
                if (i) {
                    const Staircase& sc_tgt = GetRecentStaircase(deg_tgt);
                    bool bBreak = false;
                    for (int j : two_expansion(i)) {
                        if (sc_tgt.levels[(size_t)(nd.first + j)] == 5000) { /* skip permanent cycles */
                            bBreak = true;
                            break;
                        }
                        tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)(nd.first + j)]);
                    }
                    if (bBreak)
                        continue;
                }

                AddNode(); /* Warning: the reallocation invalidates the reference sc above */
                bool bException = false;
                try {
                    int nChanges = SetDiffLeibnizV2(nd.deg, src, tgt, r);
                    if (!tgt.empty() || r <= r_max) {
                        if (depth > 1) {
                            int nDeductions = DeduceDiffs(r_max, maxPoss, top_depth, 1, timer);
                            if (depth > 2 && nChanges >= 20 && nDeductions > 0) {
                                if (timer.timeout())
                                    break;
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
                SetDiffLeibnizV2(nd.deg, src, tgt_pass, r);
                DeduceZeroDiffs(); //TODO: Implement this after CacheNullDiffs
                CacheNullDiffs(maxPoss);
                if (depth == top_depth)
                    std::clog << '(' << nd.deg.t - nd.deg.s << ',' << nd.deg.s << ") d_{" << r << '}' << src << '=' << tgt_pass << "          \n";
                if (timer.timeout())
                    break;
                if (depth > 1)
                    count += DeduceDiffs(r_max, maxPoss, top_depth == depth ? 1 : top_depth, 1, timer);
            }
            else
                ++index_nd;
        }
        else {
            const int r = sc.levels[nd.index];
            const AdamsDeg deg_src = nd.deg - AdamsDeg{r, r - 1};
            const int1d tgt = sc.basis_ind[nd.index];

            int1d src, src_pass;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;
            for (unsigned i = 0; i < i_max; ++i) {
                src.clear();
                if (i) {
                    const Staircase& sc_src = GetRecentStaircase(deg_src);
                    for (int j : two_expansion(i))
                        src = lina::AddVectors(src, sc_src.basis_ind[(size_t)(nd.first + j)]);
                }

                AddNode(); /* Warning: the reallocation invalidates the reference sc above */
                bool bException = false;
                try {
                    int nChanges = SetDiffLeibnizV2(deg_src, src, tgt, r);
                    if (!tgt.empty() || r <= r_max) {
                        if (depth > 1) {
                            int nDeductions = DeduceDiffs(r_max, maxPoss, top_depth, 1, timer);
                            if (depth > 2 && nChanges >= 20 && nDeductions > 0) {
                                if (timer.timeout())
                                    break;
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
                SetDiffLeibnizV2(deg_src, src_pass, tgt, r);
                DeduceZeroDiffs();  // TODO: Implement this after CacheNullDiffs
                CacheNullDiffs(maxPoss);
                if (depth == top_depth)
                    std::clog << '(' << nd.deg.t - nd.deg.s + 1 << ',' << nd.deg.s - r << ") d_{" << r << '}' << src_pass << '=' << tgt << "          \n";
                if (timer.timeout())
                    break;
                if (depth > 1)
                    count += DeduceDiffs(r_max, maxPoss, top_depth == depth ? 1 : top_depth, 1, timer);
            }
            else
                ++index_nd;
        }
    }
    return count;
}


// int DeduceDiffsByNat(const Groebner& gb, const std::map<AdamsDeg, Mon1d>& basis, Staircases1d& basis_ss, const Poly1d& gen_images, const std::vector<AdamsDeg>& gen_degs1, const Groebner& gb1, const std::map<AdamsDeg, Mon1d>& basis1,
//                     Staircases1d& basis_ss1, int t_try, int t_test, bool bForEt)
//{
//    bench::Timer timer;
//    int count_new_diffs = 0;
//    NullDiffs nullDiffs;
//    nullDiffs.InitNullDiffs(basis_ss, t_try, true);
//
//    /* List the degrees */
//    std::vector<AdamsDeg> degs;
//    int1d nums_tgts;
//    for (auto p = nullDiffs.null_diffs_.begin(); p != nullDiffs.null_diffs_.end(); ++p) {
//        int num_tgts = nullDiffs.null_diffs_.at(p->first).num_tgts;
//        if (num_tgts <= 10 && num_tgts >= 0) {
//            degs.push_back(p->first);
//            nums_tgts.push_back(num_tgts);
//        }
//    }
//    int1d indices_degs = ut::int_range((int)degs.size());
//    std::stable_sort(indices_degs.begin(), indices_degs.end(), [nums_tgts](int i1, int i2) { return nums_tgts[i1] < nums_tgts[i2]; });
//
//    int index_count = 0;
//    for (int index_degs : indices_degs) {
//        std::cout << ++index_count << '/' << indices_degs.size() << " " << count_new_diffs << "      \r";
//        const AdamsDeg& deg = degs[index_degs];
//        const Staircase& sc = GetRecentStaircase(basis_ss, deg);
//        const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
//        if (nd.index == -1 || nd.num_tgts > 10)  // skip if there are no null diffs or a null diff has too many possible targets
//            continue;
//        const int r = kLevelMax - sc.levels[nd.index];
//        const AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};
//        int1d src = sc.basis_ind[nd.index];
//        int1d tgt_pass;
//
//        Poly poly_src = src.empty() ? Poly{} : Indices2Poly(src, basis.at(deg));
//        Poly poly_src_image = poly_src ? subsMGb(poly_src, gb1, gen_images) : Poly();
//        AdamsDeg deg1 = poly_src_image.GetAdamsDeg(gen_degs1);
//        int1d src_image = poly_src_image ? Poly2Indices(poly_src_image, basis1.at(deg1)) : int1d{};
//
//        int count_pass = 0;
//        for (int i = 0; i < (1 << nd.num_tgts); ++i) {
//            /* Calculate tgt */
//            int1d tgt;
//            if (nd.num_tgts > 0) {
//                const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
//                for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
//                    tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
//            }
//
//            Poly poly_tgt = tgt.empty() ? Poly{} : Indices2Poly(tgt, basis.at(deg_tgt));
//            Poly poly_tgt_image = poly_tgt ? subsMGb(poly_tgt, gb1, gen_images) : Poly();
//            AdamsDeg deg1_tgt = poly_tgt_image.GetAdamsDeg(gen_degs1);
//            int1d tgt_image = poly_tgt_image ? Poly2Indices(poly_tgt_image, basis1.at(deg1_tgt)) : int1d{};
//
//            // std::cout << "Try " << i + 1 << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt;
//            // std::cout << "-> " << deg1 << " d_{" << r << '}' << src_image1 << '=' << tgt_image1;
//
//            if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
//                // std::cout << " Accepted\n";
//                tgt_pass = tgt;
//                ++count_pass;
//                break;
//            }
//
//            basis_ss.push_back({}); /* Warning: invalidate references sc above */
//            basis_ss1.push_back({});
//            bool bException = false;
//            try {
//                if (!src_image.empty() && poly_tgt_image && deg1_tgt != deg1 + AdamsDeg{r, r - 1})
//                    throw MyException(0xf434b447U, "degrees not compatible in E4t");
//                SetDiffV2(gb1, basis1, basis_ss1, deg1, src_image, tgt_image, r);
//                SetDiffV2(gb, basis, basis_ss, deg, src, tgt, r);
//            }
//            catch (SSException&) {
//                bException = true;
//            }
//            basis_ss.pop_back();
//            basis_ss1.pop_back();
//            if (!bException) {
//                // std::cout << " Success\n";
//                tgt_pass = std::move(tgt);
//                ++count_pass;
//                if (count_pass > 1)
//                    break;
//            }
//            /*else
//                std::cout << " Fail\n";*/
//        }
//        if (count_pass == 0)
//            throw SSException(0x5b3d7e35U, "no compatible differentials");
//        else if (count_pass == 1) {
//            basis_ss.push_back({});
//            SetDiffV2(gb, basis, basis_ss, deg, src, tgt_pass, r);
//            std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=' << tgt_pass << '\n';
//            /*if (count_new_diffs % 50 == 0) {
//                DeduceZeroDiffs(gb, basis, basis_ss, -1, false, -1);
//                if (bForEt)
//                    DeduceDiffsForEt(gb, basis, basis_ss, -1);
//            }*/
//            nullDiffs.InitNullDiffs(basis_ss, t_try, false);
//            ApplyRecentChanges(basis_ss);
//            ++count_new_diffs;
//
//            if (timer.Elapsed() > 3600.0 * 7) {  //
//                timer.SuppressPrint();
//                break;
//            }
//        }
//    }
//    return count_new_diffs;
//}
//
///* For a map SS1->SS2, add induced differentials from SS1 to SS2 */
// int DeduceDiffsByNatV2(const Groebner& gb, const std::map<AdamsDeg, Mon1d>& basis, Staircases1d& basis_ss, const Poly1d& gen_to_E4t, const std::vector<AdamsDeg>& gen_degs1, const Groebner& gb1, const std::map<AdamsDeg, Mon1d>& basis1,
//                        Staircases1d& basis_ss1)
//{
//     bench::Timer timer;
//
//     int index_count = 0;
//     for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
//         const AdamsDeg& deg = p->first;
//         std::cout << deg << "      \r";
//         const Staircase& sc = GetRecentStaircase(basis_ss, deg);
//         for (size_t i = 0; i < sc.levels.size(); ++i) {
//             int r = kLevelMax - sc.levels[i];
//             const AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};
//             const int1d& src = sc.basis_ind[i];
//             int1d tgt;
//             if (r < 5000) {
//                 if (sc.diffs_ind[i] == int1d{-1}) {
//                     tgt = int1d{};
//                     r -= 2;
//                     if (r < kLevelMin)
//                         continue;
//                 }
//                 else {
//                     tgt = sc.diffs_ind[i];
//                 }
//             }
//             else {
//                 if (sc.diffs_ind[i] == int1d{-1}) {
//                     tgt = int1d{};
//                     r = kLevelMax / 4;
//                 }
//                 else {
//                     continue;
//                 }
//             }
//
//             Poly poly_src = src.empty() ? Poly{} : Indices2Poly(src, basis.at(deg));
//             Poly poly_src_image = poly_src ? subsMGb(poly_src, gb1, gen_to_E4t) : Poly();
//             AdamsDeg deg1 = poly_src_image.GetAdamsDeg(gen_degs1);
//             int1d src_image = poly_src_image ? Poly2Indices(poly_src_image, basis1.at(deg1)) : int1d{};
//
//             Poly poly_tgt = tgt.empty() ? Poly{} : Indices2Poly(tgt, basis.at(deg_tgt));
//             Poly poly_tgt_image = poly_tgt ? subsMGb(poly_tgt, gb1, gen_to_E4t) : Poly();
//             AdamsDeg deg1_tgt = poly_tgt_image.GetAdamsDeg(gen_degs1);
//             int1d tgt_image = poly_tgt_image ? Poly2Indices(poly_tgt_image, basis1.at(deg1_tgt)) : int1d{};
//
//             if (deg1.IsNull())
//                 deg1 = deg1_tgt - AdamsDeg{r, r - 1};
//             int nNewDiffs = SetDiffV2(gb1, basis1, basis_ss1, deg1, src_image, tgt_image, r);
//             index_count += nNewDiffs;
//             if (nNewDiffs)
//                 std::cout << deg << ' ' << deg1 << ' ' << deg1_tgt << ' ' << " d_{" << r << '}' << src_image << '=' << tgt_image << '\n';
//         }
//     }
//     return index_count;
// }

int main_deduce_zero(int argc, char** argv, int index)
{
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce trivial differentials for degree reason\n";
        std::cout << "Usage:\n  ss deduce zero <db_filename> <table_prefix>\n\n";

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
    SSDB db(db_filename);
    SS ss = db.load_ss("AdamsE2");

    int count = 0;
    try {
        count = ss.DeduceZeroDiffs();

        db.begin_transaction();
        db.update_basis_ss(table_prefix, ss.GetChanges());
        db.end_transaction();
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
    std::string db_filename = "AdamsE2Export_t220.db";
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
    SSDB db(db_filename);
    SS ss = db.load_ss("AdamsE2");

    int count = 0;
    try {
        count = ss.DeduceImageJ();

        db.begin_transaction();
        db.update_basis_ss(table_prefix, ss.GetChanges());
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
    int r_max = 5, maxPoss = 10, depth = 2;
    double stop_time = 600;
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce differentials by Leibniz rule\n";
        std::cout << "Usage:\n  ss deduce diff [r_max] [maxPoss] [depth] [stop_time] [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  r_max = " << r_max << "\n";
        std::cout << "  maxPoss = " << maxPoss << "\n";
        std::cout << "  depth = " << depth << "\n";
        std::cout << "  stop_time = " << stop_time << "\n";
        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-19)" << std::endl;
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
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    bench::Timer timer;
    SSDB db(db_filename);
    SS ss = db.load_ss("AdamsE2");

    int count = 0;
    try {
        Timer timer(stop_time);
        count = ss.DeduceDiffs(r_max, maxPoss, depth, depth, timer);

        db.begin_transaction();
        db.update_basis_ss(table_prefix, ss.GetChanges());
        db.end_transaction();
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
    //bench::Counter::print();
    return 0;
}

/* This is for debugging */
int main_deduce_tmp(int argc, char** argv, int index)
{
    std::string db_filename = "AdamsE2Export_t220.db";
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

    //bench::Timer timer;
    SSDB db(db_filename);
    SS ss = db.load_ss("AdamsE2");

    //int count = 0;
    try {
        ss.CacheNullDiffs(10);
        for (size_t i = 0; i < ss.GetND().back().size(); ++i) {
            auto& nd = ss.GetND().back()[i];
            if (nd.deg.s == 9 && nd.deg.t == 70 + 9)
                std::cout << nd.first << ' ' << nd.count << '\n';
        }
        std::cout << "i=" << ss.GetFirstIndexOfFixedLevels(AdamsDeg(13, 47 + 13), 9994) << std::endl;
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

    //std::cout << "Changed differentials: " << count << '\n';
    //bench::Counter::print();
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
        std::cout << "  diff: deduce differentials by Leibniz rule\n\n";

        std::cout << "Version:\n  1.0 (2022-7-10)" << std::endl;
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
    else if (cmd == "tmp")
        return main_deduce_tmp(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << argv[1] << '\n';
    return 0;
}
