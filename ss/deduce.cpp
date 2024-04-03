#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <set>

/* Deduce zero differentials for degree reason */
int Diagram::DeduceTrivialDiffs(SSFlag flag)
{
    int old_count = 0, count = 0;
    const size_t num_cw = rings_.size() + modules_.size();
    while (true) {
        for (size_t iCw = 0; iCw < num_cw; ++iCw) {
            auto& nodes_ss = iCw < rings_.size() ? rings_[iCw].nodes_ss : modules_[iCw - rings_.size()].nodes_ss;
            if (nodes_ss.size() > 2 && nodes_ss.back().empty())
                continue;
            auto& name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;
            int t_max = iCw < rings_.size() ? rings_[iCw].t_max : modules_[iCw - rings_.size()].t_max;
            for (auto& [d, _] : nodes_ss.front()) {
                const auto& sc = ut::GetRecentValue(nodes_ss, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs[i] == NULL_DIFF) {
                        if (sc.levels[i] > LEVEL_PERM) {
                            const int r = LEVEL_MAX - sc.levels[i];
                            /* Find the first possible d_{r1} target for r1>=r */
                            int r1 = NextRTgt(nodes_ss, t_max, d, r);
                            if (r != r1) {
                                Logger::LogDiff(int(nodes_ss.size() - 2), EnumReason::degree, name, d, sc.basis[i], {}, r1 - 1);
                                SetCwDiffGlobal(iCw, d, sc.basis[i], {}, r1 - 1, true, flag);
                                ++count;
                            }
                        }
                        else if (sc.levels[i] < LEVEL_MAX / 2) {
                            const int r = sc.levels[i];
                            int r1 = NextRSrc(nodes_ss, d, r);
                            if (r != r1) {
                                AdamsDeg d_src = d - AdamsDeg(r1 + 1, r1);
                                Logger::LogDiffInv(int(nodes_ss.size() - 2), EnumReason::degree, name, d_src, d, {}, sc.basis[i], r1 + 1);
                                SetCwDiffGlobal(iCw, d_src, {}, sc.basis[i], r1 + 1, true, flag);
                                ++count;
                            }
                        }
                    }
                }
            }
        }
        if (old_count != count)
            old_count = count;
        else
            break;
    }
    return count;
}

/* Deduce zero differentials for degree reason */
int Diagram::DeduceTrivialDiffsCofseq(SSFlag flag)
{
    int old_count = 0, count = 0;
    while (true) {
        for (auto& cofseq : cofseqs_) {
            for (size_t iCs = 0; iCs < 3; ++iCs) {
                auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
                if (nodes_cofseq.size() > 2 && nodes_cofseq.back().empty())
                    continue;
                size_t iCs1 = (iCs + 2) % 3;
                size_t iCs2 = (iCs + 1) % 3;
                const int stem_map1 = cofseq.degMap[iCs1].stem();
                auto& name = cofseq.name;
                for (auto& [d, _] : nodes_cofseq.front()) {
                    const auto& sc = ut::GetRecentValue(nodes_cofseq, d);
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        if (sc.diffs[i] == NULL_DIFF) {
                            if (sc.levels[i] > LEVEL_PERM) {
                                const int r = LEVEL_MAX - sc.levels[i];
                                const auto& map = maps_[cofseq.indexMap[iCs]];
                                if (r <= cofseq.degMap[iCs].s && d.t <= map->t_max) {
                                    int1d dx = Residue(map->map(sc.basis[i], d, *this), *cofseq.nodes_ss[iCs2], d + cofseq.degMap[iCs], LEVEL_PERM);
                                    SetDiffLeibnizCofseq(cofseq, iCs, d, sc.basis[i], dx, cofseq.degMap[iCs].s, flag);
                                    continue;
                                }
                                /* Find the first possible d_{r1} target for r1>=r */
                                int r1 = NextRTgtCofseq(cofseq, iCs, d, r);
                                if (r != r1) {
                                    if (r1 != R_PERM) {
                                        Logger::LogDiff(int(nodes_cofseq.size() - 2), EnumReason::degree, fmt::format("{}:{}", name, iCs), d, sc.basis[i], {}, r1 - 1);
                                        SetDiffLeibnizCofseq(cofseq, iCs, d, sc.basis[i], {}, r1 - 1, flag);
                                        ++count;
                                    }
                                    else {
                                        r1 = NextRSrcCofseq(cofseq, iCs, d, R_PERM);
                                        AdamsDeg d_src = d - AdamsDeg(r1 + 1, r1 + 1 + stem_map1);
                                        Logger::LogDiffInv(int(nodes_cofseq.size() - 2), EnumReason::degree, fmt::format("{}:{}", name, iCs), d_src, d, {}, sc.basis[i], r1 + 1);
                                        SetDiffLeibnizCofseq(cofseq, iCs1, d_src, {}, sc.basis[i], r1 + 1, flag);
                                        ++count;
                                    }
                                }
                            }
                            else if (sc.levels[i] <= LEVEL_PERM) {
                                const int r = sc.levels[i];
                                int r1 = NextRSrcCofseq(cofseq, iCs, d, r);
                                if (r1 != r) {
                                    int1d dx = sc.basis[i];
                                    AdamsDeg d_src = d - AdamsDeg(r1 + 1, r1 + 1 + stem_map1);
                                    Logger::LogDiffInv(int(nodes_cofseq.size() - 2), EnumReason::degree, fmt::format("{}:{}", name, iCs), d_src, d, {}, dx, r1 + 1);
                                    SetDiffLeibnizCofseq(cofseq, iCs1, d_src, {}, dx, r1 + 1, flag);
                                    ++count;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (old_count != count)
            old_count = count;
        else
            break;
    }
    return count;
}

int Diagram::DeduceManual(SSFlag flag)
{
    /* S0__Fphi = 0 in positive dimensions */
    int old_count = 0, count = 0;
    while (true) {
        for (auto& cofseq : cofseqs_) {
            if (cofseq.name != "Fphi__RP1_256__S0")
                continue;
            size_t iCs = 2;
            size_t iCs_prev = (iCs + 2) % 3;
            auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
            if (nodes_cofseq.size() > 2 && nodes_cofseq.back().empty())
                continue;
            const int stem_map_prev = cofseq.degMap[iCs_prev].stem();
            for (auto& [d, _] : nodes_cofseq.front()) {
                if (d.stem() == 0)
                    continue;
                const auto& sc = ut::GetRecentValue(nodes_cofseq, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs[i] == NULL_DIFF && sc.levels[i] > LEVEL_PERM) {
                        int r = NextRSrcCofseq(cofseq, iCs, d, R_PERM);
                        AdamsDeg d_src = d - AdamsDeg(r + 1, r + 1 + stem_map_prev);
                        Logger::LogDiffInv(int(nodes_cofseq.size() - 2), EnumReason::degree, maps_[cofseq.indexMap[iCs_prev]]->name, d_src, d, {}, sc.basis[i], r + 1);
                        SetDiffLeibnizCofseq(cofseq, iCs_prev, d_src, {}, sc.basis[i], r + 1, flag);
                        ++count;
                    }
                }
            }
        }
        if (old_count != count)
            old_count = count;
        else
            break;
    }
    return count;
}

int Diagram::TryDiff(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, SSFlag flag, bool tryY)
{
    AddNode(flag);
    bool bException = false;
    try {
        const std::string& name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;

        Logger::Checkpoint();
        Logger::LogDiff(depth + 1, tryY ? EnumReason::try1 : EnumReason::try2, name, deg_x, x, dx, r);
        SetCwDiffGlobal(iCw, deg_x, x, dx, r, true, flag);
        if (depth == 0 && (flag & SSFlag::depth_ss_ss)) {
            DeduceDiffs(deg_x.stem(), deg_x.stem(), depth + 1, flag);
            if (tryY && !dx.empty()) {
                auto& nodes_ss = iCw < rings_.size() ? rings_[iCw].nodes_ss : modules_[iCw - rings_.size()].nodes_ss;
                if (!IsNewDiff(nodes_ss, deg_x, x, {}, r)) {
                    Logger::LogSSSSException(depth + 1, 0x311a);
                    throw SSException(0x311a, "Equivalent to trivial differential");
                }
            }
        }
        if (depth == 0 && (flag & SSFlag::depth_ss_cofseq)) {
            const auto& ind_cofs = iCw < rings_.size() ? rings_[iCw].ind_cofs : modules_[iCw - rings_.size()].ind_cofs;
            for (auto& ind_cof : ind_cofs) {
                auto& cofseq = cofseqs_[ind_cof.iCof];
                auto iCs = (size_t)ind_cof.iCs;
                auto iCs_prev = (iCs + 2) % 3;
                size_t iCw_prev = cofseq.indexCw[iCs_prev] + (cofseq.isRing[iCs_prev] ? 0 : rings_.size());
                int stem_prev = deg_x.stem() - cofseq.degMap[iCs_prev].stem();
                DeduceDiffs({iCw_prev}, stem_prev - 1, stem_prev, depth + 1, flag);
                DeduceDiffsNbhdCofseq(cofseq, iCs, deg_x.stem(), depth + 1, flag);
            }
        }
    }
    catch (SSException&) {
        bException = true;
    }
    PopNode(flag);

    if (bException)
        return 1;
    else {
        Logger::RollBackToCheckpoint();
        return 0;
    }
}

int Diagram::DeduceDiffs(size_t iCw, AdamsDeg deg, int depth, SSFlag flag)
{
    auto& nodes_ss = iCw < rings_.size() ? rings_[iCw].nodes_ss : modules_[iCw - rings_.size()].nodes_ss;
    int t_max = iCw < rings_.size() ? rings_[iCw].t_max : modules_[iCw - rings_.size()].t_max;
    const std::string& name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;

    int count = 0;
    NullDiff1d nds;
    CacheNullDiffs(nodes_ss, t_max, deg, flag, nds);

    size_t index_nd = 0;
    while (index_nd < nds.size()) {
        const auto& nd = nds[index_nd];
        int1d x, dx;
        int r;
        bool bNewDiff = false;
        /* Fixed source, find target. */
        AdamsDeg deg_src;
        if (nd.r > 0) {
            r = nd.r;
            deg_src = deg;
            x = nd.x;

            if (nd.count == 0) {
                dx.clear();
                bNewDiff = true;
            }
            else {
                int1d dx1;
                int count_pass = 0;
                unsigned i_max = 1 << nd.count;
                const AdamsDeg deg_tgt = deg_src + AdamsDeg{r, r - 1};
                const auto& sc_tgt = ut::GetRecentValue(nodes_ss, deg_tgt);

                for (unsigned i = 1; i < i_max; ++i) {
                    dx1.clear();
                    for (int j : two_exp(i))
                        dx1 = lina::add(dx1, sc_tgt.basis[(size_t)(nd.first + j)]);

                    if (!TryDiff(iCw, deg_src, x, dx1, r, depth, flag, true)) {
                        ++count_pass;
                        if (!(flag & SSFlag::try_all) && count_pass > 1)
                            break;
                        dx = std::move(dx1);
                    }
                }
                if (count_pass == 0) {
                    dx.clear();
                    bNewDiff = true;
                }
                else if (count_pass == 1) {
                    dx1.clear();
                    if (TryDiff(iCw, deg_src, x, dx1, r, depth, flag, true))
                        bNewDiff = true;
                    else
                        ++index_nd;
                }
                else
                    ++index_nd;
            }
        }
        /* Fixed target, find source. */
        else {
            r = -nd.r;
            deg_src = deg - AdamsDeg{r, r - 1};
            dx = nd.x;

            if (nd.count == 0) {
                x.clear();
                bNewDiff = true;
            }
            else {
                int1d x1;
                int count_pass = 0;
                unsigned i_max = 1 << nd.count;
                const auto& sc_src = ut::GetRecentValue(nodes_ss, deg_src);

                for (unsigned i = 1; i < i_max; ++i) {
                    x1.clear();
                    for (int j : two_exp(i))
                        x1 = lina::add(x1, sc_src.basis[(size_t)(nd.first + j)]);

                    if (!TryDiff(iCw, deg_src, x1, dx, r, depth, flag, false)) {
                        ++count_pass;
                        if (!(flag & SSFlag::try_all) && count_pass > 1)
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
                    if (TryDiff(iCw, deg_src, x1, dx, r, depth, flag, false))
                        bNewDiff = true;
                    else
                        ++index_nd;
                }
            }
        }

        if (bNewDiff) {
            ++count;
            /*if (deg == AdamsDeg(45, 197 + 45) && name == "S0" && x == int1d{3} && r == 3) {
                 fmt::print("Interupted\n");
                 throw InteruptAndSaveException(0, "debug");
            }*/
            if (nd.r > 0)
                Logger::LogDiff(depth, nd.count > 0 ? EnumReason::deduce : EnumReason::degree, name, deg, x, dx, r);
            else
                Logger::LogDiffInv(depth, nd.count > 0 ? EnumReason::deduce : EnumReason::degree, name, deg_src, deg, x, dx, r);

            SetCwDiffGlobal(iCw, deg_src, x, dx, r, true, flag);
            // if (flag & SSFlag::cofseq)
            //     count += DeduceTrivialDiffsCofseq(flag);  ////
            // DeduceTrivialDiffs(flag);
            CacheNullDiffs(nodes_ss, t_max, deg, flag, nds);
        }
        else {
            if ((flag & SSFlag::xy) && nd.r > 0) {
                if (iCw < rings_.size())
                    count += SetRingDiffLeibnizV2(iCw, deg, nd.x, nd.r, flag);
                else
                    count += SetModuleDiffLeibnizV2(iCw - rings_.size(), deg, nd.x, nd.r, flag);
            }
        }
    }
    return count;
}

int Diagram::DeduceDiffs(const size_t1d& cws, int stem_min, int stem_max, int depth, SSFlag flag)
{
    int count = 0;
    if (depth == 0)
        DeduceTrivialDiffs(flag);

    for (size_t iCw : cws) {
        std::string_view name = iCw < rings_.size() ? rings_[iCw].name : modules_[iCw - rings_.size()].name;
        auto& degs = iCw < rings_.size() ? rings_[iCw].degs_basis_order_by_stem : modules_[iCw - rings_.size()].degs_basis_order_by_stem;
        for (AdamsDeg deg : degs) {
            if (depth == 0)
                fmt::print("{} deg={}                        \r", name, deg);
            if (!BelowS0VanishingLine(deg))
                continue;
            if (deg.stem() < stem_min)
                continue;
            if (deg.stem() > stem_max)
                break;

            count += DeduceDiffs(iCw, deg, depth, flag);
        }
        if (depth == 0)
            DeduceTrivialDiffs(flag);
    }

    return count;
}

int Diagram::DeduceDiffs(int stem_min, int stem_max, int depth, SSFlag flag)
{
    return DeduceDiffs(deduce_list_spectra_, stem_min, stem_max, depth, flag);
}

int Diagram::DeduceDiffsV2()
{
    int count = 0;
    for (size_t iRing = 0; iRing < rings_.size(); ++iRing) {
        auto& ring = rings_[iRing];
        auto& nodes_ss = ring.nodes_ss;
        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            const auto& sc = ut::GetRecentValue(nodes_ss, deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.diffs[i] == NULL_DIFF) {
                    if (sc.levels[i] > LEVEL_PERM) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        auto& x = sc.basis[i];
                        count += SetRingDiffLeibnizV2(iRing, deg, x, r, SSFlag::no_op);
                    }
                }
            }
        }
    }
    for (size_t iMod = 0; iMod < modules_.size(); ++iMod) {
        auto& mod = modules_[iMod];
        auto& name = mod.name;
        auto& nodes_ss = mod.nodes_ss;
        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            fmt::print("{} deg={}                        \r", name, deg);
            const auto& sc = ut::GetRecentValue(nodes_ss, deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.diffs[i] == NULL_DIFF) {
                    if (sc.levels[i] > LEVEL_PERM) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        auto& x = sc.basis[i];
                        count += SetModuleDiffLeibnizV2(iMod, deg, x, r, SSFlag::no_op);
                    }
                }
            }
        }
    }
    return count;
}

int main_deduce_diff(int argc, char** argv, int& index, const char* desc)
{
    int stem_min = 0, stem_max = 261;
    std::string diagram_name;
    std::map<std::string, std::vector<std::string>> options;

    myio::CmdArg1d args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {{"deduce/flags...", &options}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::no_op;
    if (ut::has(options, "flags")) {
        for (auto& f : options.at("flags")) {
            if (f == "all_x")
                flag = flag | SSFlag::all_x;
            else if (f == "xy")
                flag = flag | SSFlag::xy;
            else if (f == "cofseq")
                flag = flag | SSFlag::cofseq;
            else if (f == "ss_cofseq")
                flag = flag | SSFlag::depth_ss_cofseq;
            else if (f == "ss_ss")
                flag = flag | SSFlag::depth_ss_ss;
            else if (f == "pi")
                flag = flag | SSFlag::pi;
            else if (f == "try_all")
                flag = flag | SSFlag::try_all;
            else {
                std::cout << "Not a supported flag: " << f << '\n';
                return 100;
            }
        }
    }

    Diagram diagram(diagram_name, flag);
    if (ut::has(options, "deduce"))
        diagram.SetDeduceList(options.at("deduce"));

    int count = 0;
    try {
        count = diagram.DeduceDiffs(stem_min, stem_max, 0, flag);
    }
    catch (InteruptAndSaveException&) {
    }
    /*catch (SSException&) {
    }*/
    if (flag & SSFlag::pi) {
        diagram.SimplifyPiRels();
    }
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

int main_deduce_diff_v2(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;

    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::no_op;
    Diagram diagram(diagram_name, flag, false);
    int count = diagram.DeduceDiffsV2();
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

int main_deduce_cofseq(int argc, char** argv, int& index, const char* desc)
{
    int stem_min = 0, stem_max = 261;
    std::string diagram_name;
    std::vector<std::string> strFlags;

    myio::CmdArg1d args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {{"flags...", &strFlags}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::cofseq;
    for (auto& f : strFlags) {
        if (f == "all_x")
            flag = flag | SSFlag::all_x;
        else if (f == "xy")
            flag = flag | SSFlag::xy;
        else if (f == "pi")
            flag = flag | SSFlag::pi;
        else {
            std::cout << "Not a supported flag: " << f << '\n';
            return 100;
        }
    }

    Diagram diagram(diagram_name, flag);

    int count = 0;
    count = diagram.DeduceDiffsCofseq(stem_min, stem_max, 0, flag);
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

int main_deduce_manual(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;
    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::cofseq;
    Diagram diagram(diagram_name, flag);
    int count = diagram.DeduceManual(flag);
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

/* This is for debugging */
int main_deduce_test(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;

    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::cofseq;
    Diagram diagram(diagram_name, flag, true);
    int count = diagram.CommuteCofseq(flag);
    diagram.save(diagram_name, flag);
    Logger::LogSummary("Changed differentials", count);

    return 0;
}

int main_deduce_ext(int, char**, int&, const char*);
int main_deduce_ext_def(int, char**, int&, const char*);
int main_deduce_ext_def2(int, char**, int&, const char*);
int main_deduce_ext_2tor(int, char**, int&, const char*);

/* Deduce differentials and extensions */
int main_deduce(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"diff", "Deduce differentials in ss", main_deduce_diff},
        {"diff_v2", "Deduce d(xy) or d(f(x)) when dx is uncertain", main_deduce_diff_v2},
        {"cofseq", "Deduce differentials in cofseq", main_deduce_cofseq},
        {"ext", "Deduce extensions", main_deduce_ext},
        {"ext_def", "Define decomposables", main_deduce_ext_def},
        {"ext_def2", "Define indecomposables", main_deduce_ext_def2},
        {"ext_2tor", "Compute 2-torsion degrees of generators of rings", main_deduce_ext_2tor},  //// TODO: check all rings
        {"manual", "Deduce by hard-coded human knowledge", main_deduce_manual},
        {"test", "For debugging", main_deduce_test},
    };
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;

    return 0;
}
