#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <set>

/* Deduce zero differentials for degree reason */
SSRet Category::DeduceTrivialCwDiffs(IndexUniv iCw, SSFlag flag)  //// TODO: skip if there is no change in the current depth
{
    SSRet rt;
    auto& nodes_ss = GetNodesSS(iCw);
    auto& name = GetCwName(iCw);
    int t_max = GetTMax(iCw);
    for (AdamsDeg d : GetSSDegs(iCw)) {
        const auto* sc = &nodes_ss.GetRecentSc(d);
        for (size_t i = 0; i < sc->levels.size(); ++i) {
            if (sc->diffs[i] == NULL_DIFF) {
                if (sc->levels[i] > LEVEL_PERM) {
                    const int r = LEVEL_MAX - sc->levels[i];
                    /* Find the first possible d_{r1} target for r1>=r */
                    int r1 = NextRTgt(nodes_ss, t_max, d, r);
                    if (r != r1) {
                        Logger::LogDiff(depth_, EnumReason::degree, name, d, r1 - 1, sc->basis[i], {}, "", flag);
                        if (rt += SetCwDiffGlobal(iCw, d, sc->basis[i], {}, r1 - 1, true, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                        sc = &nodes_ss.GetRecentSc(d);
                    }
                }
                else if (sc->levels[i] < LEVEL_MAX / 2) {
                    const int r = sc->levels[i];
                    int r1 = NextRSrc(nodes_ss, d, r);
                    if (r != r1) {
                        AdamsDeg d_src = d - AdamsDeg(r1 + 1, r1);
                        Logger::LogDiffInv(depth_, EnumReason::degree2, name, d, r1 + 1, {}, sc->basis[i], "", flag);
                        if (rt += SetCwDiffGlobal(iCw, d_src, {}, sc->basis[i], r1 + 1, true, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                        sc = &nodes_ss.GetRecentSc(d);
                    }
                }
            }
        }
    }

    return rt;
}

SSRet Category::DeduceTrivialCwDiffs(SSFlag flag)
{
    SSRet rt;
    while (true) {
        SSRet rt1;
        for (auto iCw = IndexRing(0); iCw.index < rings_.size(); ++iCw.index) {
            if (rt1 += DeduceTrivialCwDiffs(iCw, flag))
                return rt1;
        }
        for (auto iCw = IndexMod(0); iCw.index < modules_.size(); ++iCw.index) {
            if (rt1 += DeduceTrivialCwDiffs(iCw, flag))
                return rt1;
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

/* Deduce zero differentials for degree reason */
SSRet Category::DeduceTrivialCofDiffs(size_t iCs_index, SSFlag flag)  //// TODO: skip if there is no change in the current depth
{
    SSRet rt;
    auto& cofseq = cofseqs_[iCs_index];
    auto& name = cofseq.name;
    for (size_t iTri = 0; iTri < 3; ++iTri) {
        auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
        if (depth_ > 0 && nodes_cofseq.back().empty())
            continue;
        const size_t iTri_prev = PreviTri(iTri);
        const size_t iTri_next = NextiTri(iTri);
        const int stem_map1 = cofseq.degMap[iTri_prev].stem();
        for (AdamsDeg deg : nodes_cofseq.front().degs()) {
            const auto* sc = &nodes_cofseq.GetRecentSc(deg);
            for (size_t i = 0; i < sc->levels.size(); ++i) {
                int1d x = sc->basis[i];
                if (sc->diffs[i] != NULL_DIFF)
                    continue;
                if (sc->levels[i] > LEVEL_PERM) {
                    const int r = LEVEL_MAX - sc->levels[i];
                    const auto& map = maps_[cofseq.indexMap[iTri]];
                    /* Initialize by the graded map */
                    if (r <= cofseq.degMap[iTri].s && deg.t <= map->t_max) {
                        int1d dx = Residue(map->map(x, deg, *this), *cofseq.nodes_ss[iTri_next], deg + cofseq.degMap[iTri], LEVEL_PERM);
                        Logger::LogDiff(depth_, EnumReason::degree, fmt::format("{}:{}", name, iTri), deg, cofseq.degMap[iTri].s, x, dx, "", flag);
                        if (rt += SetDiffLeibnizCofseq(cofseq, iTri, deg, x, dx, cofseq.degMap[iTri].s, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                        sc = &nodes_cofseq.GetRecentSc(deg);
                        if (const auto iCw = cofseq.indexCw[iTri]; iCw.isRing()) { /* For newly added ring permenant cycle Apply Leibniz to all cofseqs */
                            if (rt += SetDiffLeibnizCofseq(iCw, deg, x, flag))
                                return rt;
                        }
                        sc = &nodes_cofseq.GetRecentSc(deg);
                        continue;
                    }
                    /* Find the first possible d_{r1} target for r1>=r */
                    int r1 = NextRTgtCofseq(cofseq, iTri, deg, r);
                    if (r != r1) {
                        if (r1 != R_PERM) {
                            Logger::LogDiff(depth_, EnumReason::degree, fmt::format("{}:{}", name, iTri), deg, r1 - 1, x, {}, "", flag);
                            if (rt += SetDiffLeibnizCofseq(cofseq, iTri, deg, x, {}, r1 - 1, flag)) {
                                if (flag & SSFlag::log_proof)
                                    rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                                return rt;
                            }
                            rt += SSRet::CHANGE();
                            sc = &nodes_cofseq.GetRecentSc(deg);
                        }
                        else {
                            r1 = NextRSrcCofseq(cofseq, iTri, deg, R_PERM);
                            AdamsDeg d_src = deg - AdamsDeg(r1 + 1, r1 + 1 + stem_map1);
                            Logger::LogDiffInv(depth_, EnumReason::degree2, fmt::format("{}:{}", name, iTri), deg, r1 + 1, {}, x, "", flag);
                            if (rt += SetDiffLeibnizCofseq(cofseq, iTri_prev, d_src, {}, x, r1 + 1, flag)) {
                                if (flag & SSFlag::log_proof)
                                    rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                                return rt;
                            }
                            rt += SSRet::CHANGE();
                            sc = &nodes_cofseq.GetRecentSc(deg);
                        }
                    }
                }
                else if (sc->levels[i] <= LEVEL_PERM) {
                    const int r = sc->levels[i];
                    int r1 = NextRSrcCofseq(cofseq, iTri, deg, r);
                    if (r1 != r) {
                        AdamsDeg d_src = deg - AdamsDeg(r1 + 1, r1 + 1 + stem_map1);
                        Logger::LogDiffInv(depth_, EnumReason::degree2, fmt::format("{}:{}", name, iTri), deg, r1 + 1, {}, x, "", flag);
                        if (rt += SetDiffLeibnizCofseq(cofseq, iTri_prev, d_src, {}, x, r1 + 1, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("For degree reason, {}", rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                        sc = &nodes_cofseq.GetRecentSc(deg);
                    }
                }
            }
        }
    }
    return rt;
}

SSRet Category::DeduceTrivialCofDiffs(SSFlag flag)  //// TODO: skip if there is no change in the current depth
{
    SSRet rt;
    while (true) {
        SSRet rt1;
        for (size_t index = 0; index < rings_.size(); ++index) {
            if (rt1 += DeduceTrivialCofDiffs(index, flag))
                return rt1;
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

SSRet Category::DeduceTrivialDiffs(SSFlag flag)
{
    SSRet rt;
    while (true) {
        SSRet rt1;
        for (auto iCw = IndexRing(0); iCw.index < (int)rings_.size(); ++iCw.index) {
            if (rt1 += DeduceTrivialCwDiffs(iCw, flag))
                return rt1;
        }
        for (auto iCw = IndexMod(0); iCw.index < (int)modules_.size(); ++iCw.index) {
            if (rt1 += DeduceTrivialCwDiffs(iCw, flag))
                return rt1;
        }
        if (flag & SSFlag::cofseq) {
            for (size_t index = 0; index < cofseqs_.size(); ++index) {
                if (rt1 += DeduceTrivialCofDiffs(index, flag))
                    return rt1;
            }
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

SSRet Category::DeduceManual(SSFlag flag)
{
    /* S0__Fphi = 0 in positive dimensions */
    SSRet rt;
    while (true) {
        SSRet rt1;
        for (auto& cofseq : cofseqs_) {
            if (cofseq.name != "Fphi__RP1_256__S0")
                continue;
            size_t iTri = 2;
            const size_t iTri_prev = PreviTri(iTri);
            auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
            if (depth_ > 0 && nodes_cofseq.back().empty())
                continue;
            const int stem_map_prev = cofseq.degMap[iTri_prev].stem();
            for (AdamsDeg d : nodes_cofseq.front().degs()) {
                if (d.stem() == 0)
                    continue;
                const auto& sc = nodes_cofseq.GetRecentSc(d);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.diffs[i] == NULL_DIFF && sc.levels[i] > LEVEL_PERM) {
                        int r = NextRSrcCofseq(cofseq, iTri, d, R_PERM);
                        AdamsDeg d_src = d - AdamsDeg(r + 1, r + 1 + stem_map_prev);
                        Logger::LogDiffInv(depth_, EnumReason::degree2, maps_[cofseq.indexMap[iTri_prev]]->name, d, r + 1, {}, sc.basis[i], "", flag);  //// TODO: add proof
                        if (rt1 += SetDiffLeibnizCofseq(cofseq, iTri_prev, d_src, {}, sc.basis[i], r + 1, flag)) {
                            return rt1;
                        }
                        rt1 += SSRet::CHANGE();
                    }
                }
            }
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

SSRet Category::DeduceDiffBySynthetic(SSFlag flag)
{
    SSRet rt;
    const size_t num_cw = rings_.size() + modules_.size();
    while (true) {
        SSRet rt1;
        for (size_t jCw = 0; jCw < num_cw; ++jCw) {
            auto iCw = jCw < rings_.size() ? IndexRing(jCw) : IndexMod(jCw - rings_.size());
            auto& nodes_ss = GetNodesSS(iCw);
            if (depth_ > 0 && nodes_ss.back().empty())
                continue;
            for (AdamsDeg deg : GetSSDegs(iCw)) {
                const auto& sc = nodes_ss.GetRecentSc(deg);
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    if (sc.levels[i] > LEVEL_PERM && sc.diffs[i] != NULL_DIFF) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        bool hasCross = GetCrossR(nodes_ss, deg, GetTMax(iCw), r) <= r;
                        if (rt1 += SetCwDiffSynthetic(iCw, deg, sc.basis[i], sc.diffs[i], r, hasCross, flag))
                            return rt1;
                    }
                }
            }
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

/* Deduce differentials by synthetic method */
SSRet Category::DeduceDiffBySyntheticCofseq(SSFlag flag)
{
    SSRet rt;
    while (true) {
        SSRet rt1;
        for (size_t iCof = 0; iCof < cofseqs_.size(); ++iCof) {
            auto& cofseq = cofseqs_[iCof];
            for (size_t iTri = 0; iTri < 3; ++iTri) {
                auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
                auto& map_name = maps_[cofseq.indexMap[iTri]]->name;
                size_t iTri_next = NextiTri(iTri);
                for (AdamsDeg deg : nodes_cofseq.front().degs()) {
                    const auto& sc = nodes_cofseq.GetRecentSc(deg);
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        if (sc.diffs[i] == NULL_DIFF && sc.levels[i] > LEVEL_MAX / 2) {
                            AdamsDeg deg_fx;
                            int1d fx;
                            if (GetSynImage(IndexCof(iCof, iTri), deg, sc.basis[i], LEVEL_PERM, deg_fx, fx, -1, 0) == 0) {  ////
                                int r = deg_fx.s - deg.s;
                                if (r == maps_[cofseq.indexMap[iTri]]->deg.s)
                                    continue;
                                fx = Residue(fx, *cofseq.nodes_ss[iTri_next], deg_fx, LEVEL_PERM);
                                if (fx.size() && IsNewDiff(*cofseq.nodes_ss[iTri_next], deg_fx, fx, int1d{}, R_PERM - 1)) {
                                    Logger::LogDiff(depth_, EnumReason::synext_p, cofseq.nameCw[iTri_next], deg_fx, R_PERM - 1, fx, int1d{}, "", flag);
                                    if (rt1 += SetCwDiffGlobal(cofseq.indexCw[iTri_next], deg_fx, fx, int1d{}, R_PERM - 1, true, flag)) {
                                        if (flag & SSFlag::log_proof)
                                            rt1.err_msg = fmt::format("Apply Mahowald trick to `{} {} [{}]` and map `{}`. {}", cofseq.nameCw[iTri], deg, myio::Serialize(sc.basis[i]), map_name, rt1.err_msg);
                                        return rt1;
                                    }
                                }
                                if (IsNewDiffCofseq(cofseq, iTri, deg, sc.basis[i], fx, r)) {
                                    Logger::LogDiff(depth_, EnumReason::synext, fmt::format("{}:{}", cofseq.name, iTri), deg, r, sc.basis[i], fx, "", flag);
                                    if (rt1 += SetDiffLeibnizCofseq(cofseq, iTri, deg, sc.basis[i], fx, r, flag)) {
                                        if (flag & SSFlag::log_proof)
                                            rt1.err_msg = fmt::format("Apply Mahowald trick to `{} {} [{}]` and map `{}`. {}", cofseq.nameCw[iTri], deg, myio::Serialize(sc.basis[i]), map_name, rt1.err_msg);
                                        return rt1;
                                    }
                                    rt1 += SSRet::CHANGE();
                                }
                            }
                        }
                    }
                }
            }
        }
        if (rt1.IsChanged())
            rt += rt1;
        else
            break;
    }
    return rt;
}

SSRet Category::TryDiff(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag, bool tryY)
{
    auto logid = Logger::GetCheckpoint();

    AddNode(flag);
    SSRet rt = [&]() {
        SSRet rt;
        const std::string& name = GetCwName(iCw);
        if (rt += SetCwDiffGlobal(iCw, deg_x, x, dx, r, true, flag)) {
            if (tryY)
                Logger::LogDiff(depth_, EnumReason::try1, name, deg_x, r, x, dx, rt.err_msg, flag);
            else
                Logger::LogDiffInv(depth_, EnumReason::try2, name, deg_x + AdamsDeg(r, r - 1), r, x, dx, rt.err_msg, flag);
            return rt;
        }

        if (depth_ == 1) {
            auto deg_dx = deg_x + AdamsDeg(r, r - 1);
            if (tryY)
                Logger::LogDiff(depth_, EnumReason::try1, name, deg_x, r, x, dx, "", flag);
            else
                Logger::LogDiffInv(depth_, EnumReason::try2, name, deg_x + AdamsDeg(r, r - 1), r, x, dx, "", flag);
            if (flag & SSFlag::deduce_zero) {
                if (tryY) {
                    if (dx.empty() && (rt += DeduceDiff4XDepth(iCw, deg_x, x, LEVEL_MAX - r, flag))) {
                        return rt;
                    }
                }
                else {
                    if (x.empty() && (rt += DeduceDiff4XDepth(iCw, deg_dx, dx, r, flag))) {
                        return rt;
                    }
                }
            }
            if (flag & SSFlag::deduce_pullback) {
                auto& index_maps_prev = iCw.isRing() ? rings_[iCw.index].ind_maps_prev : modules_[iCw.index].ind_maps_prev;
                for (size_t iMap : index_maps_prev) {
                    auto& map = maps_[iMap];
                    auto deg_x_prev = deg_x - map->deg;
                    auto& nodes_ss_prev = GetNodesSS(map->from);
                    if (!x.empty() && nodes_ss_prev.front().has(deg_x - map->deg) && (rt += DeduceDiffs4Deg(map->from, deg_x - map->deg, flag)))
                        return rt;
                    if (!dx.empty() && nodes_ss_prev.front().has(deg_dx - map->deg) && (rt += DeduceDiffs4Deg(map->from, deg_dx - map->deg, flag)))
                        return rt;
                }
            }
        }
        return rt;
    }();
    PopNode(flag);

    if (!rt) /* if no contradiction */
        Logger::RollBackToCheckpoint(logid);
    else if (depth_ == 0 && !(flag & SSFlag::no_exclusions))
        Logger::LogExclusion(logid + 1);
    return rt;
}

SSRet Category::DeduceDiff4Nd(IndexUniv iCw, AdamsDeg deg, const NullDiff& nd, SSFlag flag)
{
    SSRet rt;
    auto& nodes_ss = GetNodesSS(iCw);
    const auto& name = GetCwName(iCw);
    int2d* exclusions = nullptr;
    if (auto pE = nodes_ss.exclusions.has(deg, nd.x, nd.r))
        exclusions = &pE->dxs;
    else
        exclusions = nullptr;

    int1d x, dx;
    int r;
    bool bNewDiff = false;
    AdamsDeg deg_src;
    auto logid = Logger::GetCheckpoint();
    /* Fixed source, find target. */
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
            const auto& sc_tgt = nodes_ss.GetRecentSc(deg_tgt);

            for (unsigned i = 0; i < i_max; ++i) {
                dx1.clear();
                for (int j : ut::two_exp(i))
                    dx1 = lina::add(dx1, sc_tgt.basis[(size_t)(nd.first + j)]);

                if (exclusions && ut::has(*exclusions, dx1))
                    continue;
                if (TryDiff(iCw, deg_src, x, dx1, r, flag, true))
                    continue;

                ++count_pass;
                if (!(flag & SSFlag::try_all) && count_pass > 1)
                    break;
                dx = std::move(dx1);
            }
            if (count_pass == 0) {
                if (depth_ == 0)
                    throw RunTimeError(fmt::format("Contradiction. name={} deg_x={} x=[{}] r={}", name, deg, myio::Serialize(x), r));
                else
                    return SSRet::FAIL_SS();
            }
            else if (count_pass == 1)
                bNewDiff = true;
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
            const auto& sc_src = nodes_ss.GetRecentSc(deg_src);

            for (unsigned i = 0; i < i_max; ++i) {
                x1.clear();
                for (int j : ut::two_exp(i))
                    x1 = lina::add(x1, sc_src.basis[(size_t)(nd.first + j)]);

                if (exclusions && ut::has(*exclusions, x1))
                    continue;
                if (TryDiff(iCw, deg_src, x1, dx, r, flag, false))
                    continue;

                ++count_pass;
                if (!(flag & SSFlag::try_all) && count_pass > 1)
                    break;
                x = std::move(x1);
            }
            if (count_pass == 0) {
                if (depth_ == 0)
                    throw std::runtime_error(fmt::format("Contradiction. name={} deg_dx={} dx=[{}] r={}", name, deg, myio::Serialize(dx), r));
                else
                    return SSRet::FAIL_SS();
            }
            else if (count_pass == 1)
                bNewDiff = true;
        }
    }

    if (bNewDiff) {
        rt += SSRet::CHANGE_AT_X();
        if (depth_ == 0)
            Logger::RollBackToExclusionsCheckpoint(logid);
        // if (deg == AdamsDeg(2, 42 + 2) && name == "Cnu" && r == 3) {
        //      fmt::print("Interupted\n");
        //      //throw InteruptAndSaveException(0, "debug");
        // }
        if (nd.r > 0)
            Logger::LogDiff(depth_, nd.count > 0 ? EnumReason::deduce : EnumReason::degree, name, deg, r, x, dx, "", flag);
        else
            Logger::LogDiffInv(depth_, nd.count > 0 ? EnumReason::deduce2 : EnumReason::degree2, name, deg, r, x, dx, "", flag);

        if (rt += SetCwDiffGlobal(iCw, deg_src, x, dx, r, true, flag))
            return rt;
    }
    else {
        if ((flag & SSFlag::deduce_dxy) && nd.r > 0) {  //// TODO: treat the case nd.r < 0
            if (iCw.isRing())
                rt += SetRingDiffLeibnizV2(iCw.index, deg, nd.x, nd.r, flag);
            else
                rt += SetModuleDiffLeibnizV2(iCw.index, deg, nd.x, nd.r, flag);
            if (rt)
                return rt;
        }
    }
    return rt;
}

SSRet Category::DeduceDiff4X(IndexUniv iCw, AdamsDeg deg, int1d x, int level, SSFlag flag)
{
    auto& nodes_ss = GetNodesSS(iCw);
    auto t_max = GetTMax(iCw);
    NullDiff nd;
    auto [diff1, level1] = GetDiffAndLevel(nodes_ss, deg, x);
    if (level1 < level || (level1 == level && diff1 != NULL_DIFF))
        return SSRet::CHANGE_AT_X(); /* d_r(x) is already known */
    else if (level1 > level)
        return SSRet::NUL();

    nd.x = std::move(x);
    if (level > LEVEL_PERM) {
        int r = LEVEL_MAX - level;
        AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};

        auto [index, count] = CountPossDrTgt(nodes_ss, t_max, deg_tgt, r);
        if (count > deduce_count_max_)
            return SSRet::NUL();
        nd.r = r;
        nd.first = index;
        nd.count = count;
    }
    else if (level < LEVEL_MAX / 2) {
        int r = level;
        AdamsDeg deg_src = deg - AdamsDeg{r, r - 1};

        auto [index, count] = CountPossDrSrc(nodes_ss, deg_src, r);
        if (count > deduce_count_max_)
            return SSRet::NUL();
        nd.r = -r;
        nd.first = index;
        nd.count = count;
    }
    return DeduceDiff4Nd(iCw, deg, nd, flag);
}

SSRet Category::DeduceDiff4XDepth(IndexUniv iCw, AdamsDeg deg, int1d x, int original_level, SSFlag flag)  //// TODO: x may change to x - x1
{
    SSRet rt;
    if (depth_ <= 1 && (flag & SSFlag::deduce_zero)) {
        int prev_level = original_level, level_x;
        auto& nodes_ss = GetNodesSS(iCw);
        while (true) {
            level_x = GetLevel(nodes_ss, deg, x);
            if (level_x < prev_level && level_x != LEVEL_PERM) {
                if (rt += DeduceDiff4X(iCw, deg, x, level_x, flag))
                    return rt;
                prev_level = level_x;
            }
            else
                break;
        }
        if (level_x == LEVEL_PERM) {
            for (IndexUniv iCof : GetIndexCofs(iCw)) {
                x = Residue(std::move(x), nodes_ss, deg, LEVEL_PERM);
                if (x.empty())
                    break;
                if (rt += DeduceCofDiff4XDepth(iCof, deg, x, LEVEL_MAX + 1, flag))
                    return rt;
            }
        }
    }
    return rt;
}

SSRet Category::DeduceDiffs4Deg(IndexUniv iCw, AdamsDeg deg, SSFlag flag)
{
    SSRet rt;
    NullDiff1d nds;
    int t_max = GetTMax(iCw);
    CacheNullDiffs(GetNodesSS(iCw), t_max, deg, flag, nds);

    size_t index_nd = 0;
    while (index_nd < nds.size()) {
        auto rt1 = DeduceDiff4Nd(iCw, deg, nds[index_nd], flag);
        rt += rt1;
        if (rt)
            return rt;
        if (rt1.IsChangedAtX())
            CacheNullDiffs(GetNodesSS(iCw), t_max, deg, flag, nds);
        else
            ++index_nd;
    }
    return rt;
}

SSRet Category::DeduceDiffs(IndexUniv& iCw, int stem_min, int stem_max, SSFlag flag)
{
    SSRet rt;
    if (depth_ == 0)
        if (rt += DeduceTrivialCwDiffs(flag))
            return rt;

    std::string_view name = GetCwName(iCw);
    for (AdamsDeg deg : GetSSDegs(iCw)) {
        if (depth_ == 0)
            fmt::print("{} t={} deg={}                        \r", name, deg.t, deg);
        if (!BelowS0VanishingLine(deg))
            continue;
        if (deg.stem() < stem_min || deg.stem() > stem_max)
            continue;

        if (rt += DeduceDiffs4Deg(iCw, deg, flag))
            return rt;
    }
    if (depth_ == 0) {
        if (rt += DeduceTrivialCwDiffs(flag))
            return rt;
    }

    return rt;
}

SSRet Category::DeduceDiffs(int stem_min, int stem_max, int T, int id_thread, SSFlag flag)
{
    SSRet rt;
    for (IndexUniv iCw : deduce_list_spectra_) {
        if (T > 1 && iCw.index % T != id_thread)
            continue;
        if (rt += DeduceDiffs(iCw, stem_min, stem_max, flag))
            return rt;
    }
    return rt;
}

SSRet Category::DeduceDiffsV2(SSFlag flag)
{
    SSRet rt;
    for (size_t iRing = 0; iRing < rings_.size(); ++iRing) {
        auto& ring = rings_[iRing];
        auto& nodes_ss = ring.nodes_ss;
        for (AdamsDeg deg : ring.degs_ss) {
            Staircase sc = nodes_ss.GetRecentSc(deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.diffs[i] == NULL_DIFF) {
                    if (sc.levels[i] > LEVEL_PERM) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        auto& x = sc.basis[i];
                        if (rt += SetRingDiffLeibnizV2(iRing, deg, x, r, flag))
                            return rt;
                    }
                }
            }
        }
    }
    for (size_t iMod = 0; iMod < modules_.size(); ++iMod) {
        auto& mod = modules_[iMod];
        auto& name = mod.name;
        auto& nodes_ss = mod.nodes_ss;
        for (AdamsDeg deg : mod.degs_ss) {
            fmt::print("{} deg={}                        \r", name, deg);
            Staircase sc = nodes_ss.GetRecentSc(deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.diffs[i] == NULL_DIFF) {
                    if (sc.levels[i] > LEVEL_PERM) {
                        const int r = LEVEL_MAX - sc.levels[i];
                        auto& x = sc.basis[i];
                        if (rt += SetModuleDiffLeibnizV2(iMod, deg, x, r, flag))
                            return rt;
                    }
                }
            }
        }
    }
    return rt;
}

SSFlag GetFlags(const myio::string1d& str_flags)
{
    SSFlag flags = SSFlag::no_op;
    for (auto& f : str_flags) {
        if (f == "cofseq")
            flags = flags | SSFlag::cofseq;
        else if (f == "zero")
            flags = flags | SSFlag::deduce_zero;
        else if (f == "pullback")
            flags = flags | SSFlag::deduce_pullback;
        else if (f == "all_x")
            flags = flags | SSFlag::deduce_4_all_x;
        else if (f == "dxy")
            flags = flags | SSFlag::deduce_dxy;
        else if (f == "syn")
            flags = flags | SSFlag::synthetic | SSFlag::cofseq;
        else if (f == "log_proof")
            flags = flags | SSFlag::log_proof;
        else if (f == "log_deg")
            flags = flags | SSFlag::log_deg;
        else if (f == "log_nat")
            flags = flags | SSFlag::log_nat;
        else if (f == "try_all")
            flags = flags | SSFlag::try_all;
        else if (f == "pi")
            flags = flags | SSFlag::pi;
        else if (f == "no_save")
            flags = flags | SSFlag::no_save;
        else
            throw RunTimeError(fmt::format("Not a supported flag: {}", f));
    }
    return flags;
}

int main_deduce_diff(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    int stem_min = 0, stem_max = 261;
    int T = 1, id_thread = 0;
    myio::string1d str_flags, deduce_list;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"flags", &str_flags}, {"deduce_list", &deduce_list}, {"T", &T}, {"id_thread", &id_thread}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    SSFlag flags = GetFlags(str_flags) | SSFlag::log_proof;

    Category category(cat_name, "", flags);
    if (!deduce_list.empty())
        category.SetDeduceList(deduce_list);

    try {
        if (SSRet rt = category.DeduceDiffs(stem_min, stem_max, T, id_thread, flags))
            throw ErrorIdMsg(rt.code, "SSError");
    }
    catch (InteruptAndSaveException&) {
    }
    if (flags & SSFlag::pi) {
        category.SimplifyPiRels();
    }
    category.SaveNodes(cat_name, "", true, flags);
    category.PrintSummary();

    return 0;
}

int main_deduce_diff_v2(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::no_op | SSFlag::log_proof;
    Category category(cat_name, "", flag);
    if (auto rt = category.DeduceDiffsV2(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    category.SaveNodes(cat_name, "", true, flag);
    category.PrintSummary();

    return 0;
}

int main_deduce_cofseq(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    int stem_min = 0, stem_max = 261;
    int T = 1, id_thread = 0;
    myio::string1d str_flags, deduce_list;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"T", &T}, {"id_thread", &id_thread}, {"flags", &str_flags}, {"deduce_list", &deduce_list}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    SSFlag flags = GetFlags(str_flags) | SSFlag::cofseq | SSFlag::log_proof;

    Category category(cat_name, "", flags);
    if (!deduce_list.empty())
        category.SetCofDeduceList(deduce_list);

    if (auto rt = category.DeduceCofDiffs(stem_min, stem_max, T, id_thread, flags))
        throw ErrorIdMsg(rt.code, "SSError");
    category.SaveNodes(cat_name, "", true, flags);
    category.PrintSummary();

    return 0;
}

/* This is for debugging */
int main_deduce_synthetic(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::synthetic | SSFlag::cofseq | SSFlag::log_proof;
    Category category(cat_name, "", flag, true);
    if (auto rt = category.DeduceDiffBySynthetic(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    if (auto rt = category.DeduceDiffBySyntheticCofseq(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    category.SaveNodes(cat_name, "", true, flag);
    category.PrintSummary();

    return 0;
}

int main_deduce_manual(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::cofseq;
    Category category(cat_name, "", flag);
    if (auto rt = category.DeduceManual(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    category.SaveNodes(cat_name, "", true, flag);
    category.PrintSummary();

    return 0;
}

/* This is for debugging */
int main_deduce_test(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SSFlag flag = SSFlag::synthetic | SSFlag::cofseq | SSFlag::log_proof;
    Category category(cat_name, "", flag, true);
    // count += category.CommuteCofseq(flag);
    // count += category.DeduceTrivialCofDiffs(flag);

    // count += category.DeduceCofDiffs(0, 300, 0, flag);
    if (auto rt = category.DeduceDiffBySynthetic(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    if (auto rt = category.DeduceDiffBySyntheticCofseq(flag))
        throw ErrorIdMsg(rt.code, "SSError");
    category.PrintSummary();
    // category.save(cat_name, flag);
    return 0;
}

/* This is for debugging */
int main_deduce_auto(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    int stem_min = 0, stem_max = 261;
    int T = 1, id_thread = 0;
    myio::string1d str_flags, mode, deduce_list;
    int num = 10;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"T", &T}, {"id_thread", &id_thread}, {"flags", &str_flags}, {"mode", &mode}, {"deduce_list", &deduce_list}, {"num", &num}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    SSFlag flags = GetFlags(str_flags) | SSFlag::cofseq | SSFlag::log_proof;

    bool do_ss = false, do_cs = false;
    if (!mode.empty()) {
        for (auto& m : mode) {
            if (m == "ss")
                do_ss = true;
            else if (m == "cs")
                do_cs = true;
            else {
                std::cout << "Not a supported mode: " << m << '\n';
                return 101;
            }
        }
    }
    else {
        do_ss = true;
        do_cs = true;
    }
    if (do_cs)
        flags = flags | SSFlag::cofseq;

    Category category(cat_name, "", flags);
    int count = 0;
    try {
        while (count++ < num) {
            SSRet rt;
            auto id_checkpt = Logger::GetCheckpoint();
            if (do_ss)
                if (rt += category.DeduceDiffs(stem_min, stem_max, T, id_thread, flags); rt)
                    throw RunTimeError("SSError");
            if (do_cs)
                if (rt += category.DeduceCofDiffs(stem_min, stem_max, T, id_thread, flags); rt)
                    throw RunTimeError("SSError");
            if (!rt.IsChanged())
                break;
            if (!(flags & SSFlag::no_exclusions))
                category.LoadExclusions(id_checkpt + 1, flags);
        }
    }
    catch (InteruptAndSaveException&) {
        category.SaveNodes(cat_name, "debug", true, flags);
        return -37;
    }
    category.SaveNodes(cat_name, "", true, flags);
    fmt::print("round={}     \n", count);
    category.PrintSummary();

    return 0;
}

int main_deduce_ext(int, char**, int&, const char*);
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
        {"ext_2tor", "Compute 2-torsion degrees of generators of rings", main_deduce_ext_2tor},  //// TODO: check all rings
        {"manual", "Deduce by hard-coded human knowledge", main_deduce_manual},
        {"auto", "Deduce with a combination of methods", main_deduce_auto},
        {"test", "For debugging", main_deduce_test},
    };
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;

    return 0;
}
