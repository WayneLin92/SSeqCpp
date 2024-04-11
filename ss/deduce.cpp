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
        for (size_t jCw = 0; jCw < num_cw; ++jCw) {
            auto iCw = jCw < rings_.size() ? IndexRing(jCw) : IndexMod(jCw - rings_.size());
            auto& nodes_ss = GetSS(iCw);
            if (depth_ > 0 && nodes_ss.back().empty())
                continue;
            auto& name = GetCwName(iCw);
            int t_max = GetTMax(iCw);
            for (auto& [d, _] : nodes_ss.front()) {
                const auto& sc = ut::GetRecentValue(nodes_ss, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs[i] == NULL_DIFF) {
                        if (sc.levels[i] > LEVEL_PERM) {
                            const int r = LEVEL_MAX - sc.levels[i];
                            /* Find the first possible d_{r1} target for r1>=r */
                            int r1 = NextRTgt(nodes_ss, t_max, d, r);
                            if (r != r1) {
                                Logger::LogDiff(depth_, EnumReason::degree, name, d, sc.basis[i], {}, r1 - 1);
                                SetCwDiffGlobal(iCw, d, sc.basis[i], {}, r1 - 1, true, flag);
                                ++count;
                            }
                        }
                        else if (sc.levels[i] < LEVEL_MAX / 2) {
                            const int r = sc.levels[i];
                            int r1 = NextRSrc(nodes_ss, d, r);
                            if (r != r1) {
                                AdamsDeg d_src = d - AdamsDeg(r1 + 1, r1);
                                Logger::LogDiffInv(depth_, EnumReason::degree, name, d_src, d, {}, sc.basis[i], r1 + 1);
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
                size_t iCs_prev = (iCs + 2) % 3;
                size_t iCs_next = (iCs + 1) % 3;
                const int stem_map1 = cofseq.degMap[iCs_prev].stem();
                auto& name = cofseq.name;
                for (auto& [d, _] : nodes_cofseq.front()) {
                    const auto& sc = ut::GetRecentValue(nodes_cofseq, d);
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        if (sc.diffs[i] == NULL_DIFF) {
                            if (sc.levels[i] > LEVEL_PERM) {
                                const int r = LEVEL_MAX - sc.levels[i];
                                const auto& map = maps_[cofseq.indexMap[iCs]];
                                if (r <= cofseq.degMap[iCs].s && d.t <= map->t_max) {
                                    int1d dx = Residue(map->map(sc.basis[i], d, *this), *cofseq.nodes_ss[iCs_next], d + cofseq.degMap[iCs], LEVEL_PERM);
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
                                        SetDiffLeibnizCofseq(cofseq, iCs_prev, d_src, {}, sc.basis[i], r1 + 1, flag);
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
                                    SetDiffLeibnizCofseq(cofseq, iCs_prev, d_src, {}, dx, r1 + 1, flag);
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

/* Deduce zero differentials for degree reason */
int Diagram::DeduceDiffBySynthetic(SSFlag flag)
{
    int old_count = 0, count = 0;
    const size_t num_cw = rings_.size() + modules_.size();
    while (true) {
        for (size_t jCw = 0; jCw < num_cw; ++jCw) {
            auto iCw = jCw < rings_.size() ? IndexRing(jCw) : IndexMod(jCw - rings_.size());
            auto& nodes_ss = GetSS(iCw);
            if (depth_ > 0 && nodes_ss.back().empty())
                continue;
            auto& degs = iCw.isRing ? rings_[iCw.index].degs_basis_order_by_stem : modules_[iCw.index].degs_basis_order_by_stem;
            for (auto& d : degs) {
                const auto& sc = ut::GetRecentValue(nodes_ss, d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.levels[i] > LEVEL_PERM && sc.diffs[i] != NULL_DIFF) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        bool hasCross = GetCrossR(nodes_ss, d, GetTMax(iCw)) <= r;
                        count += SetCwDiffSynthetic(iCw, d, sc.basis[i], sc.diffs[i], r, hasCross, flag);
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
int Diagram::DeduceDiffBySyntheticCofseq(SSFlag flag)
{
    int old_count = 0, count = 0;
    while (true) {
        for (size_t iCof = 0; iCof < cofseqs_.size(); ++iCof) {
            auto& cofseq = cofseqs_[iCof];
            for (size_t iCs = 0; iCs < 3; ++iCs) {
                auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
                auto& map_name = maps_[cofseq.indexMap[iCs]]->name;
                size_t iCs_next = (iCs + 1) % 3;
                for (auto& [deg, _] : nodes_cofseq.front()) {
                    // fmt::print("deg={}\n", deg);
                    const auto& sc = ut::GetRecentValue(nodes_cofseq, deg);
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        if (sc.diffs[i] == NULL_DIFF && sc.levels[i] > LEVEL_MAX / 2) {
                            AdamsDeg deg_fx;
                            int1d fx;
                            if (GetSynImage(IndexCof{iCof, iCs}, deg, sc.basis[i], LEVEL_PERM, deg_fx, fx, -1, true) == 0) {
                                int r = deg_fx.s - deg.s;
                                if (r == maps_[cofseq.indexMap[iCs]]->deg.s)
                                    continue;
                                fx = Residue(fx, *cofseq.nodes_ss[iCs_next], deg_fx, LEVEL_PERM);
                                if (fx.size() && IsNewDiff(*cofseq.nodes_ss[iCs_next], deg_fx, fx, int1d{}, R_PERM - 1)) {
                                    Logger::LogDiff(depth_, EnumReason::synext, cofseq.nameCw[iCs_next], deg_fx, fx, int1d{}, R_PERM - 1);
                                    SetCwDiffGlobal(cofseq.indexCw[iCs_next], deg_fx, fx, int1d{}, R_PERM - 1, true, flag);
                                }
                                if (IsNewDiffCofseq(cofseq, iCs, deg, sc.basis[i], fx, r)) {
                                    Logger::LogDiff(depth_, EnumReason::synext, map_name, deg, sc.basis[i], fx, r);
                                    SetDiffLeibnizCofseq(cofseq, iCs, deg, sc.basis[i], fx, r, flag);
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

int Diagram::TryDiff(IndexCw iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag, bool tryY)
{
    AddNode(flag);
    bool bException = false;
    try {
        const std::string& name = GetCwName(iCw);

        Logger::Checkpoint();
        Logger::LogDiff(depth_, tryY ? EnumReason::try1 : EnumReason::try2, name, deg_x, x, dx, r);
        SetCwDiffGlobal(iCw, deg_x, x, dx, r, true, flag);
        if (depth_ == 1 && (flag & SSFlag::depth_ss_ss)) {
            auto& index_maps_prev = iCw.isRing ? rings_[iCw.index].ind_maps_prev : modules_[iCw.index].ind_maps_prev;
            DeduceDiffs(iCw, deg_x.stem(), deg_x.stem(), deg_x.s, deg_x.s + r - 2, depth_, flag);
            DeduceDiffs(iCw, deg_x.stem() + 1, deg_x.stem() + 1, deg_x.s - r, deg_x.s - 2, depth_, flag);
            for (size_t iMap : index_maps_prev) {
                auto& map = maps_[iMap];
                int stem_prev = deg_x.stem() - map->deg.stem();
                DeduceDiffs(map->from, stem_prev, stem_prev, deg_x.s, deg_x.s + r - 2, depth_, flag);
                DeduceDiffs(map->from, stem_prev + 1, stem_prev + 1, deg_x.s - r, deg_x.s - 2, depth_, flag);
            }
            if (tryY && !dx.empty() && !IsNewDiff(GetSS(iCw), deg_x, x, {}, r)) {
                Logger::LogSSSSException(depth_, 0x311a);
                throw SSException(0x311a, "Equivalent to trivial differential");
            }
        }
        if (depth_ == 1 && (flag & SSFlag::depth_ss_cofseq)) {
            const auto& ind_cofs = GetIndexCof(iCw);
            for (auto& ind_cof : ind_cofs) {
                auto& cofseq = cofseqs_[ind_cof.iCof];
                auto iCs = (size_t)ind_cof.iCs;
                auto iCs_prev = (iCs + 2) % 3;
                int stem_prev = deg_x.stem() - cofseq.degMap[iCs_prev].stem();
                DeduceDiffs(cofseq.indexCw[iCs_prev], stem_prev, stem_prev, deg_x.s, deg_x.s + r - 2, depth_, flag);
                DeduceDiffs(cofseq.indexCw[iCs_prev], stem_prev + 1, stem_prev + 1, deg_x.s - r, deg_x.s - 2, depth_, flag);
                DeduceDiffsNbhdCofseq(cofseq, iCs, deg_x.stem(), depth_, flag);
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

int Diagram::DeduceDiffs(IndexCw iCw, AdamsDeg deg, int depth, SSFlag flag)
{
    auto& nodes_ss = GetSS(iCw);
    int t_max = GetTMax(iCw);
    const auto& name = GetCwName(iCw);

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
                    for (int j : ut::two_exp(i))
                        dx1 = lina::add(dx1, sc_tgt.basis[(size_t)(nd.first + j)]);

                    if (!TryDiff(iCw, deg_src, x, dx1, r, flag, true)) {
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
                    if (TryDiff(iCw, deg_src, x, dx1, r, flag, true))
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
                    for (int j : ut::two_exp(i))
                        x1 = lina::add(x1, sc_src.basis[(size_t)(nd.first + j)]);

                    if (!TryDiff(iCw, deg_src, x1, dx, r, flag, false)) {
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
                    if (TryDiff(iCw, deg_src, x1, dx, r, flag, false))
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
            CacheNullDiffs(nodes_ss, t_max, deg, flag, nds);
        }
        else {
            if ((flag & SSFlag::xy) && nd.r > 0) {
                if (iCw.isRing)
                    count += SetRingDiffLeibnizV2(iCw.index, deg, nd.x, nd.r, flag);
                else
                    count += SetModuleDiffLeibnizV2(iCw.index, deg, nd.x, nd.r, flag);
            }
        }
    }
    return count;
}

int Diagram::DeduceDiffs(IndexCw& iCw, int stem_min, int stem_max, int s_min, int s_max, int depth, SSFlag flag)
{
    int count = 0;
    if (depth == 0)
        DeduceTrivialDiffs(flag);

    std::string_view name = GetCwName(iCw);
    auto& degs = iCw.isRing ? rings_[iCw.index].degs_basis_order_by_stem : modules_[iCw.index].degs_basis_order_by_stem;
    for (AdamsDeg deg : degs) {
        if (depth == 0)
            fmt::print("{} deg={}                        \r", name, deg);
        if (!BelowS0VanishingLine(deg))
            continue;
        if (deg.stem() < stem_min || deg.s < s_min || deg.s > s_max)
            continue;
        if (deg.stem() > stem_max)
            break;

        count += DeduceDiffs(iCw, deg, depth, flag);
    }
    if (depth == 0)
        DeduceTrivialDiffs(flag);

    return count;
}

int Diagram::DeduceDiffs(int stem_min, int stem_max, int depth, SSFlag flag)
{
    int count = 0;
    for (auto iCw : deduce_list_spectra_)
        count += DeduceDiffs(iCw, stem_min, stem_max, 0, R_PERM, depth, flag);
    return count;
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
            else if (f == "synthetic")
                flag = flag | SSFlag::synthetic | SSFlag::cofseq;
            else if (f == "ss_cofseq")
                flag = flag | SSFlag::depth_ss_cofseq | SSFlag::cofseq;
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
    Diagram diagram(diagram_name, flag);
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

/* This is for debugging */
int main_deduce_synthetic(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name;

    myio::CmdArg1d args = {{"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::synthetic | SSFlag::cofseq;
    Diagram diagram(diagram_name, flag, true);
    int count = diagram.DeduceDiffBySynthetic(flag);
    count += diagram.DeduceDiffBySyntheticCofseq(flag);
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

    SSFlag flag = SSFlag::synthetic | SSFlag::cofseq;
    Diagram diagram(diagram_name, flag, true);
    int count = diagram.DeduceDiffBySynthetic(flag);
    count += diagram.DeduceDiffBySyntheticCofseq(flag);
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
        {"synthetic", "Deduce differentials by synthetic method", main_deduce_synthetic},
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
