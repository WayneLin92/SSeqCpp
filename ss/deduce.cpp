#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <set>

/* Deduce zero differentials for degree reason
 * t_test is for DeduceDiffs() */
int Diagram::DeduceTrivialDiffs()
{
    int old_count = 0, count_new_diffs = 0;
    const size_t num_cw = rings_.size() + modules_.size();
    while (true) {
        for (size_t iCw = 0; iCw < num_cw; ++iCw) {
            auto& nodes_ss = iCw < rings_.size() ? rings_[iCw].nodes_ss : modules_[iCw - rings_.size()].nodes_ss;
            int t_max = iCw < rings_.size() ? rings_[iCw].t_max : modules_[iCw - rings_.size()].t_max;
            for (auto& [d, basis_ss_d] : nodes_ss.front()) {
                const Staircase& sc = GetRecentSc(nodes_ss, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs[i] == NULL_DIFF) {
                        if (sc.levels[i] > LEVEL_PERM) {
                            const int r = LEVEL_MAX - sc.levels[i];
                            /* Find the first possible d_{r1} target for r1>=r */
                            int r1 = NextRTgt(nodes_ss, t_max, d, r);
                            if (r != r1) {
                                SetCwDiffGlobal(iCw, d, sc.basis[i], {}, r);
                                ++count_new_diffs;
                            }
                        }
                        else if (sc.levels[i] < LEVEL_MAX / 2) {
                            const int r = sc.levels[i];
                            int r1 = NextRSrc(nodes_ss, d, r);
                            if (r != r1) {
                                SetCwDiffGlobal(iCw, d - AdamsDeg(r, r - 1), {}, sc.basis[i], r);
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

int Diagram::TryDiff(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, DeduceFlag flag)
{
    AddNode(flag);
    bool bException = false;
    try {
        std::string_view name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;
        Logger::LogDiff(depth + 1, enumReason::try_, name, deg_x, x, dx, r);
        SetCwDiffGlobal(iCw, deg_x, x, dx, r, flag & DeduceFlag::fast_try_diff);
        int count_trivial = DeduceTrivialDiffs();

        /*if (flag & DeduceFlag::set_diff)
            DeduceDiffs(depth + 1, 0);*/
        if (flag & DeduceFlag::homotopy) {
            int count_ss1 = 0, count_homotopy1 = 0;
            AdamsDeg deg_min = deg_x - AdamsDeg(0, 1);
            if (iCw >= rings_.size())
                deg_min = deg_min - modules_[iCw - rings_.size()].deg_qt;
            if (count_trivial) {
                deg_min.t = deg_min.stem();
                deg_min.s = 0;
            }
            SyncHomotopy(deg_min, count_ss1, count_homotopy1, depth + 1);
            DeduceTrivialExtensions(depth + 1);
            if (flag & DeduceFlag::homotopy_exact)
                DeduceExtensionsByExactness(deg_min.stem(), 100, depth + 1);
        }
    }
    catch (SSException&) {
        bException = true;
    }
    PopNode(flag);

    if (bException)
        return 1;
    else
        return 0;
}

int Diagram::DeduceDiffs(size_t iCw, AdamsDeg deg, int depth, DeduceFlag flag)
{
    auto& nodes_ss = iCw < rings_.size() ? rings_[iCw].nodes_ss : modules_[iCw - rings_.size()].nodes_ss;
    int t_max = iCw < rings_.size() ? rings_[iCw].t_max : modules_[iCw - rings_.size()].t_max;
    std::string_view name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;

    int count = 0;
    std::string color, color_end = "\033[0m";
    NullDiff1d nds;
    CacheNullDiffs(nodes_ss, t_max, deg, flag, nds);

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
                const Staircase& sc_tgt = GetRecentSc(nodes_ss, deg_tgt);
                dx1.clear();
                for (int j : two_expansion(i))
                    dx1 = lina::add(dx1, sc_tgt.basis[(size_t)(nd.first + j)]);

                if (!TryDiff(iCw, deg_src, x, dx1, r, depth, flag)) {
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
                if (TryDiff(iCw, deg_src, x, dx1, r, depth, flag))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }
        /* Fixed target, find source. */
        else {
            r = -nd.r;
            deg_src = deg - AdamsDeg{r, r - 1};
            const Staircase& sc_src = GetRecentSc(nodes_ss, deg_src);
            dx = nd.x;

            int1d x1;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;

            for (unsigned i = 1; i < i_max; ++i) {
                x1.clear();
                for (int j : two_expansion(i))
                    x1 = lina::add(x1, sc_src.basis[(size_t)(nd.first + j)]);

                if (!TryDiff(iCw, deg_src, x1, dx, r, depth, flag)) {
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
                if (TryDiff(iCw, deg_src, x1, dx, r, depth, flag))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }

        if (bNewDiff) {
            ++count;
            Logger::PrintDepth();
            if (nd.r > 0)
                Logger::LogDiff(depth, enumReason::deduce, name, deg, x, dx, r);
            else
                Logger::LogDiffInv(depth, enumReason::deduce, name, deg, x, dx, r);
            SetCwDiffGlobal(iCw, deg_src, x, dx, r);
            int count_trivial = DeduceTrivialDiffs();
            count += count_trivial;
            CacheNullDiffs(nodes_ss, t_max, deg, flag, nds);
            if (flag & DeduceFlag::homotopy) {
                int count_homotopy1 = 0;
                AdamsDeg deg_min = deg_src - AdamsDeg(0, 1);
                if (iCw >= rings_.size())
                    deg_min = deg_min - modules_[iCw - rings_.size()].deg_qt;
                if (count_trivial) {
                    deg_min.t = deg_min.stem();
                    deg_min.s = 0;
                }
                SyncHomotopy(deg_min, count, count_homotopy1, depth);
                DeduceTrivialExtensions(depth);
                if (flag & DeduceFlag::homotopy_exact)
                    DeduceExtensionsByExactness(deg_min.stem(), 100, depth);
            }
        }
        else
            Logger::ClearDepth();
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
        if (flag & DeduceFlag::homotopy_exact)
            DeduceExtensionsByExactness(0, 100, depth + 1);
    }

    const size_t num_cw = rings_.size() + modules_.size();
    for (size_t iCw = 0; iCw < num_cw; ++iCw) {
        std::string_view name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;
        if (name == "j")
            continue;
        auto& degs = iCw < rings_.size() ? rings_[iCw].degs_basis_order_by_stem : modules_[iCw - rings_.size()].degs_basis_order_by_stem;

        for (AdamsDeg deg : degs) {
            if (depth == 0)
                std::cout << name << "  deg=" << deg.StrAdams() << "                        \r";
            if (deg.stem() < stem_min)
                continue;
            else if (deg.stem() > stem_max)
                break;

            count += DeduceDiffs(iCw, deg, depth, flag);
        }
    }
    return count;
}

int main_deduce_diff(int argc, char** argv, int& index, const char* desc)
{
    int stem_min = 0, stem_max = 261;
    std::string diagram_name = "default";
    std::vector<std::string> strFlags;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"diagram", &diagram_name}, {"flags...", &strFlags}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    DeduceFlag flag = DeduceFlag::no_op;
    for (auto& f : strFlags) {
        if (f == "all_x")
            flag = flag | DeduceFlag::all_x;
        else if (f == "homotopy")
            flag = flag | DeduceFlag::homotopy;
        else if (f == "exact")
            flag = flag | DeduceFlag::homotopy_exact;
        else {
            std::cout << "Not a supported flag: " << f << '\n';
            return 100;
        }
    }

    Diagram diagram(diagram_name, flag);

    int count = 0;
    if (flag & DeduceFlag::homotopy) {
        int count_homotopy1 = 0;
        diagram.SyncHomotopy(AdamsDeg(0, 0), count, count_homotopy1, 0);
        diagram.DeduceTrivialExtensions(0);
        if (flag & DeduceFlag::homotopy_exact)
            diagram.DeduceExtensionsByExactness(0, 100, 0);
    }
    count = diagram.DeduceDiffs(stem_min, stem_max, 0, flag);
    if (flag & DeduceFlag::homotopy) {
        diagram.SimplifyPiRels();
    }
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

/* This is for debugging */
int main_deduce_tmp(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    DeduceFlag flag = DeduceFlag::no_op;
    Diagram diagram(diagram_name, flag);
    int count = 0;
    int count_ss = 0, count_homotopy = 0;
    diagram.DeduceExtensions(0, 30, count_ss, count_homotopy, 0, flag);
    diagram.SimplifyPiRels();
    diagram.save(diagram_name, flag);
    std::cout << "Changed differentials: " << count << '\n';

    return 0;
}

/* Generate the table of the spectral sequence */
int main_deduce(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"diff", "Deduce differentials", main_deduce_diff},
        {"ext", "Deduce extensions", main_deduce_ext},
        {"ext_def", "Define decomposables", main_deduce_ext_def},
        {"ext_def2", "Define indecomposables", main_deduce_ext_def2},
        {"ext_2tor", "Compute 2-torsion degrees of generators of rings", main_deduce_ext_2tor},  //// TODO: check all rings
        {"tmp", "Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting", main_deduce_tmp},
    };
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, "Make deductions on ss or homotopy", VERSION, subcmds))
        return error;

    return 0;
}
