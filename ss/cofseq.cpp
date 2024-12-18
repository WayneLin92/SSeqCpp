#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

bool IsPossTgtCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg, int r_max)
{
    const size_t iTri_prev = PreviTri(iTri);
    const auto& deg_map = cofseq.degMap[iTri_prev];
    if (deg.t - deg_map.t > cofseq.t_max[iTri_prev])
        return true;
    r_max = std::min(r_max, deg.s);
    if (r_max < 0)
        return false;
    const auto& nodes_cofseq_from = cofseq.nodes_cofseq[iTri_prev];
    const auto& nodes_ss_from = *cofseq.nodes_ss[iTri_prev];
    if (IsPossTgt(*cofseq.nodes_ss[iTri], deg))  // Fixed a bug
        return true;
    for (int r1 = deg_map.s; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, deg_map.stem() + r1};
        if (PossEinf(nodes_ss_from, d_src))
            if (PossMoreEinf(nodes_ss_from, d_src) || GetMaxLevelWithND(nodes_cofseq_from.GetRecentSc(d_src)) >= LEVEL_MAX - r1)
                return true;
    }
    return false;
}

size_t GetFirstIndexOfFixedLevelsCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg, int level_min)
{
    const size_t iTri_next = NextiTri(iTri);
    auto& nodes_ss_next = *cofseq.nodes_ss[iTri_next];
    const int stem_map = cofseq.degMap[iTri].stem();
    const auto& sc = cofseq.nodes_cofseq[iTri].GetRecentSc(deg);
    size_t result = sc.levels.size();
    for (size_t i = sc.levels.size(); i-- > 0;) {
        if (sc.diffs[i] == NULL_DIFF || sc.levels[i] < level_min)
            break;
        if (i == 0 || sc.levels[i] != sc.levels[i - 1]) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg_tgt = deg + AdamsDeg{r, r + stem_map};
            if (IsPossTgt(nodes_ss_next, deg_tgt) || IsPossTgtCofseq(cofseq, iTri_next, deg_tgt, r - 1))
                break;
            else
                result = i;
        }
    }
    return result;
}

int Category::GetFirstFixedLevelForPlotCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg)
{
    const size_t iTri_next = NextiTri(iTri);
    auto& nodes_ss2 = *cofseq.nodes_ss[iTri_next];
    const int stem_map = cofseq.degMap[iTri].stem();
    const auto& sc = cofseq.nodes_cofseq[iTri].GetRecentSc(deg);
    int result = LEVEL_MAX + 1;
    for (size_t i = sc.levels.size(); i-- > 0 && sc.levels[i] >= LEVEL_PERM;) {
        if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg2 = deg + AdamsDeg{r, r + stem_map};
            if (IsPossTgt(nodes_ss2, deg2) || IsPossTgtCofseq(cofseq, iTri_next, deg2, r - 1))
                break;
            else
                result = sc.levels[i];
        }
    }
    return result;
}

std::pair<int, int> Category::CountPossDrTgtCofseq(const CofSeq& cofseq, size_t iTri, const AdamsDeg& deg_tgt, int r, bool possMorePerm) const
{
    std::pair<int, int> result;
    if (cofseq.nodes_cofseq[iTri].front().has(deg_tgt)) {
        const auto& sc_tgt = cofseq.nodes_cofseq[iTri].GetRecentSc(deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        int end = possMorePerm ? (int)sc_tgt.basis.size() : (int)GetFirstIndexOfFixedLevelsCofseq(cofseq, iTri, deg_tgt, LEVEL_MAX - r + 1);
        result.second = end - result.first;
    }
    else if (deg_tgt.t > cofseq.t_max[iTri])
        result = {-1, 100000}; /* Infinitely many possibilities */
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Category::CountPossDrSrcCofseq(const CofSeq& cofseq, size_t iTri, const AdamsDeg& deg_src, int r, bool possMorePerm) const
{
    std::pair<int, int> result;
    if (cofseq.nodes_cofseq[iTri].front().has(deg_src)) {
        const auto& sc_src = cofseq.nodes_cofseq[iTri].GetRecentSc(deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        int end = possMorePerm ? (int)sc_src.basis.size() : (int)GetFirstIndexOfFixedLevelsCofseq(cofseq, iTri, deg_src, LEVEL_MAX - r + 1);
        result.second = end - result.first;
    }
    else if (deg_src.t > cofseq.t_max[iTri])
        result = {-1, 100000}; /* Infinitely many possibilities */
    else
        result = {-1, 0};
    return result;
}

int Category::NextRTgtCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg, int r) const
{
    size_t iTri_next = NextiTri(iTri);
    auto& nodes_ss_tgt = *cofseq.nodes_ss[iTri_next];
    int stem_map = cofseq.degMap[iTri].stem();
    int t_max2 = cofseq.t_max[iTri_next];
    for (int r1 = r; r1 <= R_PERM; ++r1) {
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 + stem_map};
        if (d_tgt.t > t_max2)
            return r1;
        if (r1 >= 20 && AboveJ(d_tgt) && BelowCokerJ(deg)) /* Image of J */
            return R_PERM;
        if (AboveS0Vanishing(d_tgt) && !PossEinf(nodes_ss_tgt, d_tgt))
            return R_PERM;
        if (PossMoreEinf(nodes_ss_tgt, d_tgt))
            return r1;
        auto [_, count] = CountPossDrTgtCofseq(cofseq, iTri_next, d_tgt, r1, false);
        if (count > 0)
            return r1;
    }
    return R_PERM;
}

int Category::NextRSrcCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg, int r) const
{
    int r_max = std::min(r, deg.s);
    const size_t iTri_prev = PreviTri(iTri);
    int t_max1 = cofseq.t_max[iTri_prev];
    auto& nodes_ss_src = *cofseq.nodes_ss[iTri_prev];
    int stem_map = cofseq.degMap[iTri_prev].stem();
    for (int r1 = r_max; r1 >= 0; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 + stem_map};
        if (d_src.t > t_max1)
            return r1;
        if (r1 >= 20 && AboveJ(deg) && BelowCokerJ(d_src)) /* Image of J */
            continue;
        if (PossMoreEinf(nodes_ss_src, d_src))
            return r1;
        auto [_, count] = CountPossDrSrcCofseq(cofseq, iTri_prev, d_src, r1, false);
        if (count > 0)
            return r1;
    }
    return -1;
}

void Category::CacheNullCofDiffs(const CofSeq& cofseq, size_t iTri, AdamsDeg deg, SSFlag flag, NullDiffCofseq1d& nds) const
{
    const int count_ss_max = depth_ == 0 ? 3 : 0;
    nds.clear();
    auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
    const size_t iTri_prev = PreviTri(iTri);
    const size_t iTri_next = NextiTri(iTri);
    auto& nodes_ss_prev = *cofseq.nodes_ss[iTri_prev];
    auto& nodes_ss_next = *cofseq.nodes_ss[iTri_next];
    int stem_map_prev = cofseq.degMap[iTri_prev].stem();
    int stem_map = cofseq.degMap[iTri].stem();
    const auto& sc = nodes_cofseq.GetRecentSc(deg);
    for (size_t i = 0; i < sc.diffs.size(); ++i) {
        if (sc.diffs[i] != NULL_DIFF)
            continue;

        NullDiffCofseq nd;
        if (sc.levels[i] > LEVEL_PERM) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg_tgt = deg + AdamsDeg{r, r + stem_map};
            if (deg_tgt.t > cofseq.t_max[iTri_next])
                continue;
            auto [index_ss, count_ss] = CountPossMorePerm(nodes_ss_next, deg_tgt);
            auto [index, count] = CountPossDrTgtCofseq(cofseq, iTri_next, deg_tgt, r, count_ss > 0);
            if (count + count_ss > deduce_count_max_ || count_ss > count_ss_max)
                continue;
            nd.r = r;
            nd.first = index;
            nd.count = count;
            nd.first_ss = index_ss;
            nd.count_ss = count_ss;
        }
        else if (sc.levels[i] < LEVEL_MAX / 2) {
            int r = sc.levels[i];
            AdamsDeg deg_src = deg - AdamsDeg{r, r + stem_map_prev};
            if (deg_src.t > cofseq.t_max[iTri_prev])
                continue;
            auto [index_ss, count_ss] = CountPossMorePerm(nodes_ss_prev, deg_src);
            auto [index, count] = CountPossDrSrcCofseq(cofseq, iTri_prev, deg_src, r, count_ss > 0);
            if (count + count_ss > deduce_count_max_ || count_ss > count_ss_max)
                continue;
            nd.r = -r - 1;
            nd.first = index;
            nd.count = count;
            nd.first_ss = index_ss;
            nd.count_ss = count_ss;
        }
        else
            continue;

        size_t j = i + 1;
        while (j < sc.levels.size() && sc.diffs[j] == NULL_DIFF && sc.levels[i] == sc.levels[j])
            ++j;
        if (flag & SSFlag::deduce_4_all_x) {
            const unsigned k_max = unsigned(1) << (j - i);
            for (unsigned k = 1; k < k_max; ++k) {
                nd.x.clear();
                for (int l : ut::two_exp(k))
                    nd.x = lina::add(nd.x, sc.basis[i + l]);
                nds.push_back(nd);
            }
        }
        else {
            for (size_t k = 0; k < j - i; ++k) {
                nd.x = sc.basis[i + k];
                nds.push_back(nd);
            }
        }
        i = j - 1;
    }
}

bool Category::IsNewDiffCofseq(const CofSeq& cofseq, size_t iTri, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const
{
    size_t iTri_next = NextiTri(iTri);
    int stem_map = cofseq.degMap[iTri].stem();
    auto& nodes_cofseq_dx = cofseq.nodes_cofseq[iTri_next];
    auto& nodes_ss_dx = *cofseq.nodes_ss[iTri_next];
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, stem_map + r};
    if (x.empty() || IsZeroOnLevel(cofseq.nodes_ss[iTri]->GetRecentSc(deg_x), x, LEVEL_PERM))
        return !dx.empty() && !IsZeroOnLevel(nodes_ss_dx.GetRecentSc(deg_dx), dx, LEVEL_PERM) && !IsZeroOnLevel(nodes_cofseq_dx.GetRecentSc(deg_dx), dx, r);
    int1d dx1 = GetDiff(cofseq.nodes_cofseq[iTri], deg_x, x, r);  //// TODO: optimize the allocation
    if (dx1 == NULL_DIFF)
        return true;
    auto diff = Residue(lina::add(dx, dx1), nodes_ss_dx, deg_dx, LEVEL_PERM);
    return !diff.empty() && !IsZeroOnLevel(nodes_cofseq_dx.GetRecentSc(deg_dx), diff, r);
}

SSRet Category::ReSetScCofseq(CofSeq& cofseq, size_t iTri, AdamsDeg deg, SSFlag flag)
{
    SSRet rt;
    const size_t iTri_next = NextiTri(iTri);
    const size_t iTri_prev = PreviTri(iTri);
    const auto& nodes_ss = GetNodesSS(cofseq.indexCw[iTri]);

    auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];

    const auto& sc_ss = nodes_ss.GetRecentSc(deg);
    size_t first_PC = GetFirstIndexOnLevel(sc_ss, LEVEL_PERM);
    auto& sc_old = nodes_cofseq.GetRecentSc(deg);
    Staircase sc = sc_old;
    Staircase sc_new;

    int stem_map = cofseq.degMap[iTri].stem();
    int stem_prev = cofseq.degMap[iTri_prev].stem();

    int2d images;
    int1d levels;
    bool changed = false;

    for (size_t i = 0; i < sc.basis.size(); ++i) {
        /* Reduce sc.basis by ss boundaries */
        for (size_t j = 0; j < first_PC; ++j) {
            if (std::binary_search(sc.basis[i].begin(), sc.basis[i].end(), sc_ss.basis[j].front())) {
                sc.basis[i] = lina::add(sc.basis[i], sc_ss.basis[j]);
                changed = true;
            }
        }

        for (size_t j = 0; j < sc_new.basis.size(); ++j) {
            if (std::binary_search(sc.basis[i].begin(), sc.basis[i].end(), sc_new.basis[j][0])) {
                sc.basis[i] = lina::add(sc.basis[i], sc_new.basis[j]);
                if (sc.levels[i] == sc_new.levels[j] && sc.diffs[i] != NULL_DIFF)
                    sc.diffs[i] = lina::add(sc.diffs[i], sc_new.diffs[j]);
            }
        }
        if (sc.basis[i].empty()) {
            if (!sc.diffs[i].empty() && sc.diffs[i] != NULL_DIFF) {
                images.push_back(std::move(sc.diffs[i]));
                levels.push_back(LEVEL_MAX - sc.levels[i]);
            }
        }
        else {
            sc_new.basis.push_back(std::move(sc.basis[i]));  //// Improve: modify inplace
            sc_new.diffs.push_back(std::move(sc.diffs[i]));
            sc_new.levels.push_back(sc.levels[i]);
        }
    }

    if (changed)
        nodes_cofseq.back()[deg] = std::move(sc_new);

    /* Add images and cycles */
    for (size_t i = 0; i < levels.size(); ++i) {
        if (levels[i] < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg + AdamsDeg{levels[i], stem_map + levels[i]};
            if (rt += SetImageScCofseq(cofseq, iTri_next, deg_image_new, images[i], NULL_DIFF, levels[i] - 1, flag))
                return rt;
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - levels[i];
            AdamsDeg deg_image_new = deg - AdamsDeg{r_image, stem_prev + r_image};
            if (rt += SetDiffScCofseq(cofseq, iTri_prev, deg_image_new, images[i], NULL_DIFF, r_image + 1, flag))
                return rt;
        }
    }
    return rt;
}

SSRet Category::SetDiffGlobalCofseq(CofSeq& cofseq, size_t iTri, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    SSRet rt;
    if (newCertain || IsNewDiffCofseq(cofseq, iTri, deg_x, x, dx, r)) {
        int stem_map = cofseq.degMap[iTri].stem();
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r + stem_map);
        if (x.empty()) {
            int r1 = r;
            r = NextRSrcCofseq(cofseq, NextiTri(iTri), deg_dx, r1 - 1) + 1;
            if (r != r1)
                deg_x = deg_dx - AdamsDeg(r, r + stem_map);
        }
        else if (dx.empty()) {
            r = NextRTgtCofseq(cofseq, iTri, deg_x, r + 1) - 1;
        }
        if (!x.empty() && r == R_PERM - 1) {
            int r1 = NextRSrcCofseq(cofseq, iTri, deg_x, R_PERM);

            const size_t iTri_prev = PreviTri(iTri);
            AdamsDeg d_src = deg_x - AdamsDeg(r1 + 1, r1 + 1 + cofseq.degMap[iTri_prev].stem());
            if (rt += SetDiffLeibnizCofseq(cofseq, iTri_prev, d_src, {}, x, r1 + 1, flag))
                return rt;
        }
        else if (rt += SetDiffLeibnizCofseq(cofseq, iTri, deg_x, x, dx, r, flag))
            return rt;
    }
    return rt;
}

void Category::SetUnivDiffGlobal(std::string& name, AdamsDeg deg, int r, const int1d& x, const int1d& dx, bool isCs, bool isDInv, SSFlag flag)
{
    AdamsDeg deg_x, deg_dx;
    if (!isCs) {
        auto iCw_x = GetIndexCwByName(name);
        if (!isDInv) {
            deg_x = deg;
            deg_dx = deg_x + AdamsDeg(r, r - 1);
        }
        else {
            deg_dx = deg;
            deg_x = deg_dx - AdamsDeg(r, r - 1);
        }
        if (auto rt = SetCwDiffGlobal(iCw_x, deg_x, x, dx, r, false, flag))
            throw RunTimeError("Failed to SetCwDiffGlobal()");
    }
    else {
        auto iTri_x = GetIndexCofByName(name);
        auto& cofseq = GetCofSeqs()[iTri_x.index];
        if (!isDInv) {
            int stem_map = cofseq.degMap[iTri_x.iTri].stem();
            deg_x = deg;
            deg_dx = deg_x + AdamsDeg(r, stem_map + r);
        }
        else {
            iTri_x = iTri_x.prev();
            int stem_map = cofseq.degMap[iTri_x.iTri].stem();
            deg_dx = deg;
            deg_x = deg_dx - AdamsDeg(r, stem_map + r);
        }
        if (x.size())
            if (auto rt = SetCwDiffGlobal(cofseq.indexCw[iTri_x.iTri], deg_x, x, int1d{}, R_PERM - 1, false, flag))
                throw RunTimeError("Failed to Set permanent cycle");
        if (dx.size())
            if (auto rt = SetCwDiffGlobal(cofseq.indexCw[NextiTri(iTri_x.iTri)], deg_dx, dx, int1d{}, R_PERM - 1, false, flag))
                throw RunTimeError("Failed to Set permanent cycle");

        if (auto rt = SetDiffGlobalCofseq(cofseq, iTri_x.iTri, deg_x, x, dx, r, false, flag))
            throw RunTimeError("Failed to SetCwDiffGlobal()");
    }
}

SSRet Category::TryCofDiff(IndexUniv iCof, AdamsDeg deg_x, AdamsDeg deg_dx, const int1d& x, const int1d& dx, const int1d& perm, int r, SSFlag flag, bool tryY)
{
    /*if (cofseq.name == "S0__C2__S0" && iCs == 2 && deg_dx == AdamsDeg(6, 104 + 6) && r == 0)
        fmt::print("debug\n");*/
    auto logid = Logger::GetCheckpoint();

    auto& cofseq = cofseqs_[iCof.index];
    const auto iTri = iCof.iTri;
    AddNode(flag);
    SSRet rt = [&]() {
        SSRet rt;
        if (!perm.empty()) {
            if (tryY) {
                rt += SetCwDiffGlobal(cofseq.indexCw[NextiTri(iTri)], deg_dx, perm, int1d{}, R_PERM - 1, true, flag);
                if (rt) {
                    Logger::LogDiff(depth_, EnumReason::try1, fmt::format("{}:{}", cofseq.name, iTri), deg_x, r, x, dx, rt.err_msg, flag);
                    return rt;
                }
            }
            else {
                rt += SetCwDiffGlobal(cofseq.indexCw[iTri], deg_x, perm, int1d{}, R_PERM - 1, true, flag);
                if (rt) {
                    Logger::LogDiffInv(depth_, EnumReason::try2, fmt::format("{}:{}", cofseq.name, NextiTri(iTri)), deg_dx, r, x, dx, rt.err_msg, flag);
                    return rt;
                }
            }
        }
        if (rt += SetDiffGlobalCofseq(cofseq, iTri, deg_x, x, dx, r, true, flag)) {
            if (tryY)
                Logger::LogDiff(depth_, EnumReason::try1, fmt::format("{}:{}", cofseq.name, iTri), deg_x, r, x, dx, rt.err_msg, flag);
            else
                Logger::LogDiffInv(depth_, EnumReason::try2, fmt::format("{}:{}", cofseq.name, NextiTri(iTri)), deg_dx, r, x, dx, rt.err_msg, flag);
            return rt;
        }
        if (depth_ == 1) {
            if (flag & SSFlag::deduce_zero) {
                if (tryY)
                    Logger::LogDiff(depth_, EnumReason::try1, fmt::format("{}:{}", cofseq.name, iTri), deg_x, r, x, dx, "", flag);
                else
                    Logger::LogDiffInv(depth_, EnumReason::try2, fmt::format("{}:{}", cofseq.name, NextiTri(iTri)), deg_dx, r, x, dx, "", flag);
                if (tryY) {
                    if (dx.empty() && (rt += DeduceCofDiff4XDepth(iCof, deg_x, x, LEVEL_MAX - r, flag))) {
                        return rt;
                    }
                }
                else {
                    if (x.empty() && (rt += DeduceCofDiff4XDepth(iCof.next(), deg_dx, dx, r, flag))) {
                        return rt;
                    }
                }
            }
        }

        return rt;
    }();
    PopNode(flag);

    if (!rt)
        Logger::RollBackToCheckpoint(logid);
    else if (depth_ == 0 && !(flag & SSFlag::no_exclusions))
        Logger::LogExclusion(logid + 1);
    return rt;
}

SSRet Category::DeduceCofDiffs4Deg(IndexUniv iCof, AdamsDeg deg, SSFlag flag)
{
    SSRet rt;
    auto& cofseq = cofseqs_[iCof.index];
    auto iTri = size_t(iCof.iTri);
    NullDiffCofseq1d nds;
    CacheNullCofDiffs(cofseq, iTri, deg, flag, nds);

    size_t index_nd = 0;
    while (index_nd < nds.size()) {
        auto rt1 = DeduceCofDiff4Nd(iCof, deg, nds[index_nd], flag);
        if (rt1)
            return rt1;
        rt += rt1;
        if (rt1.IsChangedAtX())
            CacheNullCofDiffs(cofseq, iTri, deg, flag, nds);
        else
            ++index_nd;
    }
    return rt;
}

SSRet Category::DeduceCofDiff4Nd(IndexUniv iCof, AdamsDeg deg, const NullDiffCofseq& nd, SSFlag flag)
{
#ifdef MYDEBUG
    MyException::Assert(iCof.isCofseq(), "Need iCof.isCof()");
#endif

    SSRet rt;
    auto& cofseq = cofseqs_[iCof.index];
    auto iTri = (size_t)iCof.iTri;
    const size_t iTri_prev = PreviTri(iTri);
    const size_t iTri_next = NextiTri(iTri);
    auto& nodes_ss_prev = *cofseq.nodes_ss[iTri_prev];
    auto& nodes_ss_next = *cofseq.nodes_ss[iTri_next];
    int stem_map_prev = cofseq.degMap[iTri_prev].stem();
    int stem_map = cofseq.degMap[iTri].stem();
    auto& nodes_cofseq_prev = cofseq.nodes_cofseq[iTri_prev];
    auto& nodes_cofseq_next = cofseq.nodes_cofseq[iTri_next];

    int2d* exclusions = nullptr;
    int1d x, dx, perm;

    if (auto pE = cofseq.nodes_cofseq[iTri].exclusions.has(deg, nd.x, nd.r))
        exclusions = &pE->dxs;
    else
        exclusions = nullptr;

    int r;
    bool bNewDiff = false;
    AdamsDeg deg_src, deg_tgt;
    size_t iTri_src;
    auto logid = Logger::GetCheckpoint();
    /* Fixed source, find target. */
    if (nd.r >= 0) {
        r = nd.r;
        iTri_src = iTri;
        deg_src = deg;
        deg_tgt = deg_src + AdamsDeg{r, r + stem_map};
        x = nd.x;
        if (nd.count + nd.count_ss == 0) {
            dx.clear();
            perm.clear();
            bNewDiff = true;
        }
        else {
            int1d dx1;
            int1d perm1;
            int count_pass = 0;
            unsigned i_max = 1 << (nd.count + nd.count_ss);
            const auto& sc_ss = nodes_ss_next.GetRecentSc(deg_tgt);

            for (unsigned i = 0; i < i_max; ++i) {
                dx1.clear();
                perm1.clear();
                const Staircase* sc_tgt = nd.count ? &nodes_cofseq_next.GetRecentSc(deg_tgt) : nullptr;
                for (int j : ut::two_exp(i)) {
                    if (j < nd.count)
                        dx1 = lina::add(dx1, sc_tgt->basis[(size_t)(nd.first + j)]);
                    else
                        perm1 = lina::add(perm1, sc_ss.basis[(size_t)(nd.first_ss + j - nd.count)]);
                }
                dx1 = lina::add(dx1, perm1);

                if (exclusions && ut::has(*exclusions, dx1))
                    continue;
                if (TryCofDiff(iCof, deg_src, deg_tgt, x, dx1, perm1, r, flag, true))
                    continue;

                ++count_pass;
                if (!(flag & SSFlag::try_all) && count_pass > 1)
                    break;
                dx = std::move(dx1);
                perm = std::move(perm1);
            }
            if (count_pass == 0) {
                if (depth_ == 0)
                    throw RunTimeError(fmt::format("Contradiction. name={} deg_x={} x=[{}] r={}", fmt::format("{}:{}", cofseq.name, iTri), deg, myio::Serialize(x), r));
                else
                    return SSRet::FAIL_SS();
            }
            else if (count_pass == 1)
                bNewDiff = true;
        }
    }
    /* Fixed target, find source. */
    else {
        r = -nd.r - 1;
        iTri_src = PreviTri(iTri);
        deg_src = deg - AdamsDeg{r, r + stem_map_prev};
        deg_tgt = deg;
        dx = nd.x;

        /*if (cofseq.name == "S0__C2__S0" && iCs == 2 && deg_tgt == AdamsDeg(6, 104 + 6) && r == 0)
            fmt::print("debug\n");*/
        if (nd.count + nd.count_ss == 0) {
            x.clear();
            perm.clear();
            bNewDiff = true;
        }
        else {
            int1d x1;
            int1d perm1;
            int count_pass = 0;
            unsigned i_max = 1 << (nd.count + nd.count_ss);
            const auto& sc_ss = nodes_ss_prev.GetRecentSc(deg_src);

            for (unsigned i = 0; i < i_max; ++i) {
                x1.clear();
                perm1.clear();
                const Staircase* sc_src = nd.count ? &nodes_cofseq_prev.GetRecentSc(deg_src) : nullptr;
                for (int j : ut::two_exp(i)) {
                    if (j < nd.count)
                        x1 = lina::add(x1, sc_src->basis[(size_t)(nd.first + j)]);
                    else
                        perm1 = lina::add(perm1, sc_ss.basis[(size_t)(nd.first_ss + j - nd.count)]);
                }
                x1 = lina::add(x1, perm1);

                if (exclusions && ut::has(*exclusions, x1))
                    continue;
                if (TryCofDiff(iCof.prev(), deg_src, deg, x1, dx, perm1, r, flag, false))
                    continue;

                ++count_pass;
                if (!(flag & SSFlag::try_all) && count_pass > 1)
                    break;
                x = std::move(x1);
                perm = std::move(perm1);
            }
            if (count_pass == 0) {
                if (depth_ == 0)
                    throw std::runtime_error(fmt::format("Contradiction. name={} deg_dx={} dx=[{}] r={}", fmt::format("{}:{}", cofseq.name, iTri_next), deg, myio::Serialize(dx), r));
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
        if (!perm.empty()) {
            if (nd.r >= 0) {
                Logger::LogDiff(depth_, EnumReason::cofseq_in, cofseq.nameCw[iTri_next], deg_tgt, R_PERM - 1, perm, int1d{}, "", flag);
                if (rt += SetCwDiffGlobal(cofseq.indexCw[iTri_next], deg_tgt, perm, int1d{}, R_PERM - 1, true, flag))
                    return rt;
            }
            else {
                Logger::LogDiff(depth_, EnumReason::cofseq_in, cofseq.nameCw[iTri_prev], deg_src, R_PERM - 1, perm, int1d{}, "", flag);
                if (rt += SetCwDiffGlobal(cofseq.indexCw[iTri_prev], deg_src, perm, int1d{}, R_PERM - 1, true, flag))
                    return rt;
            }
        }
        if (nd.r >= 0)
            Logger::LogDiff(depth_, EnumReason::dd_cof, fmt::format("{}:{}", cofseq.name, iTri), deg, r, x, dx, "", flag);
        else
            Logger::LogDiffInv(depth_, EnumReason::dd_cof2, fmt::format("{}:{}", cofseq.name, iTri), deg, r, x, dx, "", flag);

        if (rt += SetDiffGlobalCofseq(cofseq, iTri_src, deg_src, x, dx, r, true, flag))
            return rt;
        if (depth_ == 0) {
            if (rt += DeduceTrivialDiffs(flag))
                return rt;
        }
    }
    else {
        /*if ((flag & SSFlag::xy) && nd.direction > 0) {
            if (iUniv < rings_.size())
                count += SetRingDiffLeibnizV2(iUniv, deg, nd.x, nd.r);
            else
                count += SetModuleDiffLeibnizV2(iUniv - rings_.size(), deg, nd.x, nd.r);
        }*/
    }
    return rt;
}

SSRet Category::DeduceCofDiff4X(IndexUniv iCof, AdamsDeg deg, int1d x, int level, SSFlag flag)
{
    auto& cofseq = cofseqs_[iCof.index];
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCof.iTri];
    NullDiffCofseq nd;
    int level_x_ss = GetLevel(*cofseq.nodes_ss[iCof.iTri], deg, x);
    if (level_x_ss > LEVEL_PERM)
        return SSRet::NUL();
    if (level_x_ss < LEVEL_MAX / 2) /* It was a bug without this. */
        return SSRet::CHANGE_AT_X();
    auto [diff1, level1] = GetDiffAndLevel(nodes_cofseq, deg, x);
    if (level1 < level || (level1 == level && diff1 != NULL_DIFF))
        return SSRet::CHANGE_AT_X();
    else if (level1 > level)
        return SSRet::NUL();

    nd.x = std::move(x);
    const int count_ss_max = 3;
    if (level > LEVEL_PERM) {
        int stem_map = cofseq.degMap[iCof.iTri].stem();
        size_t iTri_next = size_t(iCof.iTri + 1) % 3;
        auto& nodes_ss_next = *cofseq.nodes_ss[iTri_next];
        int r = LEVEL_MAX - level;
        AdamsDeg deg_tgt = deg + AdamsDeg{r, r + stem_map};

        auto [index_ss, count_ss] = CountPossMorePerm(nodes_ss_next, deg_tgt);
        auto [index, count] = CountPossDrTgtCofseq(cofseq, iTri_next, deg_tgt, r, count_ss > 0);
        if (count + count_ss > deduce_count_max_ || count_ss > count_ss_max)
            return SSRet::NUL();
        nd.r = r;
        nd.first = index;
        nd.count = count;
        nd.first_ss = index_ss;
        nd.count_ss = count_ss;
    }
    else if (level < LEVEL_MAX / 2) {
        size_t iTri_prev = PreviTri(iCof.iTri);
        int stem_map_prev = cofseq.degMap[iTri_prev].stem();
        auto& nodes_ss_prev = *cofseq.nodes_ss[iTri_prev];
        int r = level;
        AdamsDeg deg_src = deg - AdamsDeg(r, stem_map_prev + r);
        auto [index_ss, count_ss] = CountPossMorePerm(nodes_ss_prev, deg_src);
        auto [index, count] = CountPossDrSrcCofseq(cofseq, iTri_prev, deg_src, r, count_ss > 0);
        if (count + count_ss > deduce_count_max_ || count_ss > count_ss_max)
            return SSRet::NUL();
        nd.r = -r - 1;
        nd.first = index;
        nd.count = count;
        nd.first_ss = index_ss;
        nd.count_ss = count_ss;
    }

    return DeduceCofDiff4Nd(iCof, deg, nd, flag);
}

SSRet Category::DeduceCofDiff4XDepth(IndexUniv iCof, AdamsDeg deg, int1d x, int original_level, SSFlag flag)
{
    SSRet rt;
    if (depth_ <= 1 && (flag & SSFlag::deduce_zero)) {
        int prev_level = original_level, level_ss;
        auto& nodes_cofseq = GetNodesSS(iCof);
        auto& nodes_ss = *cofseqs_[iCof.index].nodes_ss[iCof.iTri];
        while (true) {
            level_ss = GetLevel(nodes_ss, deg, x);
            if (level_ss == LEVEL_PERM) {
                int level_x = GetLevel(nodes_cofseq, deg, x);
                if (level_x < prev_level) {
                    if (level_x != -1) {
                        if (rt += DeduceCofDiff4X(iCof, deg, x, level_x, flag))
                            return rt;
                    }
                    prev_level = level_x;
                }
                else
                    break;
            }
            else
                break;
        }
        if (level_ss < LEVEL_PERM) {
            if (rt += DeduceDiff4XDepth(cofseqs_[iCof.index].indexCw[iCof.iTri], deg, x, LEVEL_PERM, flag))
                return rt;
        }
    }
    return rt;
}

SSRet Category::DeduceCofDiffs(int stem_min, int stem_max, int T, int id_thread, SSFlag flag)
{
    SSRet rt;
    if (depth_ == 0) {
        if (rt += DeduceTrivialCofDiffs(flag))
            return rt;
        if (rt += CommuteCofseq(flag))
            return rt;
        if (rt += DeduceTrivialCofDiffs(flag))
            return rt;
    }
    for (size_t iCof : deduce_list_cofseq_) {
        if (T > 1 && iCof % T != id_thread)
            continue;
        auto& cofseq = cofseqs_[iCof];
        for (size_t iTri = 0; iTri < 3; ++iTri) {
            for (auto deg : cofseq.nodes_cofseq[iTri].front().degs()) {
                if (deg.stem() >= stem_min && deg.stem() <= stem_max && BelowS0VanishingLine(deg)) {
                    if (depth_ == 0)
                        fmt::print("{}:{} t={} deg={}                        \r", cofseq.name, iTri, deg.t, deg);
                    if (rt += DeduceCofDiffs4Deg(IndexCof(iCof, iTri), deg, flag))
                        return rt;
                }
            }
        }
    }
    return rt;
}

SSRet Category::DeduceCofDiffsNbhd(IndexUniv iCof, int stem, SSFlag flag)
{
    SSRet rt;
    auto& cofseq = cofseqs_[iCof.index];
    size_t iTri = size_t(iCof.iTri);
    size_t iTri_prev = PreviTri(iCof.iTri);
    size_t iTri_next = NextiTri(iCof.iTri);
    std::array<size_t, 3> iCss = {iTri_prev, iTri, iTri_next};
    std::array<int, 3> stems = {stem - cofseq.degMap[iTri_prev].stem(), stem, stem + cofseq.degMap[iTri].stem()};

    for (size_t i = 0; i < 3; ++i) {
        size_t iTri_ = iCss[i];
        int stem_ = stems[i];
        for (AdamsDeg deg : cofseq.nodes_cofseq[iTri_].front().degs()) {
            if (deg.stem() == stem_ && (rt += DeduceCofDiffs4Deg(iCof, deg, flag)))
                return rt;
        }
    }
    return rt;
}

SSRet Category::SyncToCofseq(SSFlag flag)
{
    SSRet rt;
    size_t size_rings = rings_.size();
    const size_t num_cw = size_rings + modules_.size();
    AdamsDeg1d unsynced_degs;
    for (size_t j = 0; j < num_cw; ++j) {
        auto iCw = GetIndexCw(j);
        auto& nodes_ss = GetNodesSS(iCw);
        for (AdamsDeg deg : nodes_ss.unsynced_degs()) {
            for (auto& ind_cof : GetIndexCofs(iCw)) {
                auto& cofseq = cofseqs_[ind_cof.index];
                if (cofseq.nodes_cofseq[ind_cof.iTri].front().has(deg)) {
                    if (rt += ReSetScCofseq(cofseq, ind_cof.iTri, deg, flag))
                        return rt;
                }
            }
        }
    }
    return rt;
}

SSRet Category::SyncCofseq(SSFlag flag)
{
    SSRet rt;
    size_t size_rings = rings_.size();
    const size_t num_cw = size_rings + modules_.size();
    AdamsDeg1d unsynced_degs;
    for (size_t j = 0; j < num_cw; ++j) {
        auto iCw = GetIndexCw(j);
        auto& nodes_ss = GetNodesSS(iCw);
        const auto& ind_cofs = GetIndexCofs(iCw);

        /* CsIn */
        for (AdamsDeg deg : GetSSDegs(iCw)) {
            const auto& sc = nodes_ss.GetRecentSc(deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] == LEVEL_PERM) {
                    for (auto& ind_cof : ind_cofs) {
                        auto& cofseq = cofseqs_[ind_cof.index];
                        auto iTri = (size_t)ind_cof.iTri;
                        if (!cofseq.nodes_cofseq[iTri].front().has(deg))
                            cofseq.nodes_cofseq[iTri].front()[deg] = {};
                        if (!IsZeroOnLevel(cofseq.nodes_cofseq[ind_cof.iTri].GetRecentSc(deg), sc.basis[i], LEVEL_MAX)) {
                            if (rt += SetDiffScCofseq(cofseq, ind_cof.iTri, deg, sc.basis[i], NULL_DIFF, 0, flag))
                                return rt;
                        }
                    }
                }
            }
        }

        /* CsOut */
        for (auto& ind_cof : ind_cofs) {
            auto& cofseq = cofseqs_[ind_cof.index];
            for (AdamsDeg deg : nodes_ss.front().degs()) {
                if (cofseq.nodes_cofseq[ind_cof.iTri].front().has(deg)) {
                    if (rt += ReSetScCofseq(cofseq, ind_cof.iTri, deg, flag))
                        return rt;
                }
            }
        }
    }
    return rt;
}

/* Return the minimal length of the crossing differentials */
int GetCofseqCrossR(const SSNodes& nodes_cofseq, const SSNodes& nodes_ss, AdamsDeg deg, int t_max, int r_min, int result_min)
{
    int result = R_PERM;
    for (int s = 1; s <= R_PERM; ++s) {
        auto deg_x = deg + AdamsDeg{s, s};
        if (s + r_min > result)
            return result;
        if (deg_x.t > t_max || PossMoreEinf(nodes_ss, deg_x)) {
            if (s + r_min < result)
                result = std::max(result_min, s + r_min);
            return result;
        }
        if (AboveS0Vanishing(deg_x) && !nodes_cofseq.front().has(deg_x))
            return result;
        if (!nodes_cofseq.front().has(deg_x)) {
            if (AboveS0Vanishing(deg_x))
                return result;
            continue;
        }
        auto& sc = nodes_cofseq.GetRecentSc(deg_x);
        for (size_t i = sc.levels.size(); i-- > 0;) {
            if (sc.levels[i] < LEVEL_MAX / 2)
                break;
            int r1 = LEVEL_MAX - sc.levels[i];
            if (s + r1 >= result)
                break;
            if (sc.diffs[i] == NULL_DIFF) {
                result = std::max(result_min, s + r1);
                break;
            }
            if (s + r1 >= result_min) {
                result = s + r1;
                break;
            }
        }
    }
    return result;
}

void GetRAndDiff(const Staircase& sc, size_t i, int& r, int1d& diff)
{
    diff.clear();
    if (sc.levels[i] < LEVEL_MAX / 2)
        r = R_PERM - 1;
    else if (sc.diffs[i] != NULL_DIFF) {
        diff = sc.diffs[i];
        r = LEVEL_MAX - sc.levels[i];
    }
    else {
        r = LEVEL_MAX - sc.levels[i] - 1;
    }
}

void GetRAndDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x, int& r, int1d& diff)
{
    int level = -1;
    std::tie(diff, level) = GetDiffAndLevel(nodes_ss, deg_x, x);
    if (level < LEVEL_MAX / 2) {
        r = R_PERM - 1;
        diff.clear();
    }
    else {
        r = LEVEL_MAX - level;
        if (diff == NULL_DIFF) {
            diff.clear();
            --r;
        }
    }
}

SSRet Category::CommuteCofseq(size_t iComm, SSFlag flag)
{
    SSRet rt;
    auto& comm = comms_[iComm];
    auto& f0 = maps_[comm.f0];
    auto& f1 = maps_[comm.f1];
    auto& cofseq_f0 = cofseqs_[f0->iCof.index];
    auto& nodes_cof_f0_tgt = cofseq_f0.nodes_cofseq[((size_t)f0->iCof.iTri + 1) % 3];
    auto& cofseq_f1 = cofseqs_[f1->iCof.index];

    if (comm.g0 == -1) {
        for (AdamsDeg deg : nodes_cof_f0_tgt.front().degs()) {
            Staircase sc = nodes_cof_f0_tgt.GetRecentSc(deg);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] < LEVEL_MAX / 2) {
                    if (IsNewDiffCofseq(cofseq_f1, f1->iCof.iTri, deg, sc.basis[i], int1d{}, R_PERM - 1)) {
                        Logger::LogDiff(depth_, EnumReason::comm, fmt::format("{} => {}", comm.name, f1->name), deg, R_PERM - 1, sc.basis[i], int1d{}, "", flag); ////
                        if (rt += SetDiffGlobalCofseq(cofseq_f1, f1->iCof.iTri, deg, sc.basis[i], int1d{}, R_PERM - 1, true, flag))
                            return rt;
                        rt += SSRet::CHANGE();
                    }
                }
            }
        }
    }
    else {
        auto& g0 = maps_[comm.g0];
        auto& cofseq_g0 = cofseqs_[g0->iCof.index];
        int stem_f0 = cofseq_f0.degMap[f0->iCof.iTri].stem();
        auto& nodes_cof_f0_src = cofseq_f0.nodes_cofseq[f0->iCof.iTri];
        auto& nodes_cof_f1_src = cofseq_f1.nodes_cofseq[f1->iCof.iTri];
        auto& nodes_cof_g0_src = cofseqs_[g0->iCof.index].nodes_cofseq[g0->iCof.iTri];
        auto& nodes_ss_f0_src = *cofseq_f0.nodes_ss[f0->iCof.iTri];
        auto& nodes_ss_f1_src = *cofseq_f1.nodes_ss[f1->iCof.iTri];
        auto& nodes_ss_g0_src = *cofseqs_[g0->iCof.index].nodes_ss[g0->iCof.iTri];

        int r_f0 = -1, r_f1 = -1, r_g0 = -1;
        int1d f0x, f1f0x, g0x;
        if (comm.g1 == -1) {
            for (AdamsDeg deg_x : nodes_cof_f0_src.front().degs()) {
                Staircase sc_f0 = nodes_cof_f0_src.GetRecentSc(deg_x);
                for (size_t i = 0; i < sc_f0.levels.size(); ++i) {
                    GetRAndDiff(sc_f0, i, r_f0, f0x);
                    AdamsDeg deg_f0x = deg_x + AdamsDeg(r_f0, r_f0 + stem_f0);
                    f0x = Residue(f0x, nodes_ss_f1_src, deg_f0x, LEVEL_PERM);
                    int cross_f0 = GetCofseqCrossR(nodes_cof_f0_src, nodes_ss_f0_src, deg_x, cofseq_f0.t_max[f0->iCof.iTri], cofseq_f0.degMap[f0->iCof.iTri].s, 0);

                    /* f1(f0x) = g0x
                     * We choose {x} such that g0{x}={g0x}
                     */
                    if (!f0x.empty()) {
                        GetRAndDiff(nodes_cof_g0_src, deg_x, sc_f0.basis[i], r_g0, g0x);
                        if (r_f0 >= cross_f0) { /* When f0 has crossing, then g0 should not have crossing */
                            int cross_g0 = GetCofseqCrossR(nodes_cof_g0_src, nodes_ss_g0_src, deg_x, cofseq_g0.t_max[g0->iCof.iTri], cofseq_g0.degMap[g0->iCof.iTri].s, 0);
                            if (r_g0 >= cross_g0) {
                                r_g0 = cross_g0 - 1;
                                g0x.clear();
                            }
                        }
                        r_f1 = r_g0 - r_f0;
                        if (r_f1 >= 0 && IsNewDiffCofseq(cofseq_f1, f1->iCof.iTri, deg_f0x, f0x, g0x, r_f1)) {
                            Logger::LogDiff(depth_, EnumReason::comm, fmt::format("{} => {}", comm.name, f1->name), deg_f0x, r_f1, f0x, g0x, "", flag);
                            if (rt += SetDiffGlobalCofseq(cofseq_f1, f1->iCof.iTri, deg_f0x, f0x, g0x, r_f1, true, flag))
                                return rt;
                            rt += SSRet::CHANGE();
                        }
                    }

                    /* g0x = f1f0x
                     * We choose {x} such that f0{x}={f0x}
                     */
                    {
                        if (!f0x.empty()) {
                            GetRAndDiff(nodes_cof_f1_src, deg_f0x, f0x, r_f1, f1f0x);
                        }
                        else {
                            r_f1 = R_PERM - 1;
                            f1f0x.clear();
                        }

                        GetRAndDiff(nodes_cof_g0_src, deg_x, sc_f0.basis[i], r_g0, g0x);
                        int cross_g0 = GetCofseqCrossR(nodes_cof_g0_src, nodes_ss_g0_src, deg_x, cofseq_g0.t_max[g0->iCof.iTri], cofseq_g0.degMap[g0->iCof.iTri].s, 0);
                        if (r_g0 >= cross_g0) {
                            r_g0 = cross_g0 - 1;
                        }
                        int cross_f1_min = r_g0 - r_f0 + (g0x.empty() ? 1 : 0);
                        /*if (comm.name == "CW_eta_nu_sigma__S0" && deg_x == AdamsDeg(11, 105 + 11) && sc_f0.basis[i] == int1d{2}) {
                            fmt::print("r_g0={} cross_g0={}, cross_f1_min={}\n", r_g0, cross_g0, cross_f1_min);
                            fmt::print("debug\n");
                        }*/
                        int cross_f1 = GetCofseqCrossR(nodes_cof_f1_src, nodes_ss_f1_src, deg_f0x, cofseq_f1.t_max[f1->iCof.iTri], cofseq_f1.degMap[f1->iCof.iTri].s, cross_f1_min);
                        if (r_f1 >= cross_f1) {
                            r_f1 = cross_f1 - 1;
                            f1f0x.clear();
                        }

                        r_g0 = r_f0 + r_f1;
                        if (r_g0 > R_PERM - 1)
                            r_g0 = R_PERM - 1;

                        /*if (comm.name == "CW_eta_nu_sigma__S0" && deg_x == AdamsDeg(11, 105 + 11) && sc_f0.basis[i] == int1d{2}) {
                            fmt::print("r_f0={}, f0x={}, r_f1={}, f1f0x={}\n", r_f0, myio::Serialize(f0x), r_f1, myio::Serialize(f1f0x));
                            fmt::print("cross_f1={}\n", cross_f1);
                            fmt::print("debug\n");
                        }*/
                        if (r_g0 >= 0 && IsNewDiffCofseq(cofseq_g0, g0->iCof.iTri, deg_x, sc_f0.basis[i], f1f0x, r_g0)) {
                            Logger::LogDiff(depth_, EnumReason::comm, fmt::format("{} => {}", comm.name, g0->name), deg_x, r_g0, sc_f0.basis[i], f1f0x, "", flag);
                            if (rt += SetDiffGlobalCofseq(cofseq_g0, g0->iCof.iTri, deg_x, sc_f0.basis[i], f1f0x, r_g0, true, flag))
                                return rt;
                            rt += SSRet::CHANGE();
                        }
                    }
                }
            }
        }
        else {
            auto& g1 = maps_[comm.g1];
            auto& cofseq_g1 = cofseqs_[g1->iCof.index];
            auto& nodes_cof_g1_src = cofseqs_[g1->iCof.index].nodes_cofseq[g1->iCof.iTri];
            auto& nodes_ss_g1_src = *cofseqs_[g1->iCof.index].nodes_ss[g1->iCof.iTri];
            int stem_g0 = cofseq_g0.degMap[g0->iCof.iTri].stem();
            // auto& nodes_cof_g0_tgt = cofseq_g0.nodes_cofseq[((size_t)g0->ind_cof.iCs + 1) % 3];

            int r_g1 = -1;
            int1d g1g0x;
            for (AdamsDeg deg_x : nodes_cof_f0_src.front().degs()) {
                Staircase sc_f0 = nodes_cof_f0_src.GetRecentSc(deg_x);
                for (size_t i = 0; i < sc_f0.levels.size(); ++i) {
                    /* f1(f0x) = g1(g0x)
                     * We choose {x} such that g0{x}={g0x}
                     */
                    {
                        GetRAndDiff(sc_f0, i, r_f0, f0x);
                        AdamsDeg deg_f0x = deg_x + AdamsDeg(r_f0, r_f0 + stem_f0);
                        f0x = Residue(f0x, nodes_ss_f1_src, deg_f0x, LEVEL_PERM);
                        int cross_f0 = GetCofseqCrossR(nodes_cof_f0_src, nodes_ss_f0_src, deg_x, cofseq_f0.t_max[f0->iCof.iTri], cofseq_f0.degMap[f0->iCof.iTri].s, 0);

                        if (!f0x.empty()) {
                            GetRAndDiff(nodes_cof_g0_src, deg_x, sc_f0.basis[i], r_g0, g0x);
                            if (r_f0 >= cross_f0) { /* When f0 has crossing, then g0 should not have crossing */
                                int cross_g0 = GetCofseqCrossR(nodes_cof_g0_src, nodes_ss_g0_src, deg_x, cofseq_g0.t_max[g0->iCof.iTri], cofseq_g0.degMap[g0->iCof.iTri].s, 0);
                                if (r_g0 >= cross_g0) {
                                    r_g0 = cross_g0 - 1;
                                    g0x.clear();
                                }
                            }
                            AdamsDeg deg_g0x = deg_x + AdamsDeg(r_g0, r_g0 + stem_g0);
                            if (!g0x.empty()) {
                                GetRAndDiff(nodes_cof_g1_src, deg_g0x, g0x, r_g1, g1g0x);
                            }
                            else {
                                r_g1 = R_PERM - 1;
                                g1g0x.clear();
                            }

                            GetRAndDiff(nodes_cof_f1_src, deg_f0x, f0x, r_f1, f1f0x);
                            int cross_f1 = GetCofseqCrossR(nodes_cof_f1_src, nodes_ss_f1_src, deg_x, cofseq_f1.t_max[f1->iCof.iTri], cofseq_f1.degMap[f1->iCof.iTri].s, 0);
                            if (r_f1 >= cross_f1) {
                                r_f1 = cross_f1 - 1;
                                f1f0x.clear();
                            }
                            int cross_g1_min = r_f0 + r_f1 - r_g0 + (f1f0x.empty() ? 1 : 0);
                            int cross_g1 = GetCofseqCrossR(nodes_cof_g1_src, nodes_ss_g1_src, deg_g0x, cofseq_g1.t_max[g1->iCof.iTri], cofseq_g1.degMap[g1->iCof.iTri].s, cross_g1_min);
                            if (r_g1 >= cross_g1) {
                                r_g1 = cross_g1 - 1;
                                g1g0x.clear();
                            }

                            r_f1 = r_g0 + r_g1 - r_f0;
                            if (r_f1 > R_PERM - 1)
                                r_f1 = R_PERM - 1;

                            if (r_f1 >= 0 && IsNewDiffCofseq(cofseq_f1, f1->iCof.iTri, deg_f0x, f0x, g1g0x, r_f1)) {
                                Logger::LogDiff(depth_, EnumReason::comm, fmt::format("{} => {}", comm.name, f1->name), deg_f0x, r_f1, f0x, g1g0x, "", flag); ////
                                if (rt += SetDiffGlobalCofseq(cofseq_f1, f1->iCof.iTri, deg_f0x, f0x, g1g0x, r_f1, true, flag))
                                    return rt;
                                rt += SSRet::CHANGE();
                            }
                        }
                    }

                    /* g1(g0x) = f1(f0x)
                     * We choose {x} such that f0{x}={f0x}
                     */
                    {
                        GetRAndDiff(nodes_cof_g0_src, deg_x, sc_f0.basis[i], r_g0, g0x);
                        AdamsDeg deg_g0x = deg_x + AdamsDeg(r_g0, r_g0 + stem_g0);
                        g0x = Residue(g0x, nodes_ss_g1_src, deg_g0x, LEVEL_PERM);
                        int cross_g0 = GetCofseqCrossR(nodes_cof_g0_src, nodes_ss_g0_src, deg_x, cofseq_g0.t_max[g0->iCof.iTri], cofseq_g0.degMap[g0->iCof.iTri].s, 0);

                        if (!g0x.empty()) {
                            GetRAndDiff(sc_f0, i, r_f0, f0x);
                            if (r_g0 >= cross_g0) { /* When g0 has crossing, then f0 should not have crossing */
                                int cross_f0 = GetCofseqCrossR(nodes_cof_f0_src, nodes_ss_f0_src, deg_x, cofseq_f0.t_max[f0->iCof.iTri], cofseq_f0.degMap[f0->iCof.iTri].s, 0);
                                if (r_f0 >= cross_f0) {
                                    r_f0 = cross_f0 - 1;
                                    f0x.clear();
                                }
                            }
                            AdamsDeg deg_f0x = deg_x + AdamsDeg(r_f0, r_f0 + stem_f0);
                            if (!f0x.empty()) {
                                GetRAndDiff(nodes_cof_f1_src, deg_f0x, f0x, r_f1, f1f0x);
                            }
                            else {
                                r_f1 = R_PERM - 1;
                                f1f0x.clear();
                            }

                            GetRAndDiff(nodes_cof_g1_src, deg_g0x, g0x, r_g1, g1g0x);
                            int cross_g1 = GetCofseqCrossR(nodes_cof_g1_src, nodes_ss_g1_src, deg_x, cofseq_g1.t_max[g1->iCof.iTri], cofseq_g1.degMap[g1->iCof.iTri].s, 0);
                            if (r_g1 >= cross_g1) {
                                r_g1 = cross_g1 - 1;
                                g1g0x.clear();
                            }
                            int cross_f1_min = r_g0 + r_g1 - r_f0 + (g1g0x.empty() ? 1 : 0);
                            int cross_f1 = GetCofseqCrossR(nodes_cof_f1_src, nodes_ss_f1_src, deg_f0x, cofseq_f1.t_max[f1->iCof.iTri], cofseq_f1.degMap[f1->iCof.iTri].s, cross_f1_min);
                            if (r_f1 >= cross_f1) {
                                r_f1 = cross_f1 - 1;
                                f1f0x.clear();
                            }

                            r_g1 = r_f0 + r_f1 - r_g0;
                            if (r_g1 > R_PERM - 1)
                                r_g1 = R_PERM - 1;

                            if (r_g1 >= 0 && IsNewDiffCofseq(cofseq_g1, g1->iCof.iTri, deg_g0x, g0x, f1f0x, r_g1)) {
                                Logger::LogDiff(depth_, EnumReason::comm, fmt::format("{} => {}", comm.name, g1->name), deg_g0x, r_g1, g0x, f1f0x, "", flag); ////
                                if (rt += SetDiffGlobalCofseq(cofseq_g1, g1->iCof.iTri, deg_g0x, g0x, f1f0x, r_g1, true, flag))
                                    return rt;
                                rt += SSRet::CHANGE();
                            }
                        }
                    }
                }
            }
        }
    }

    return rt;
}

SSRet Category::CommuteCofseq(SSFlag flag)
{
    SSRet rt;
    for (size_t iComm = 0; iComm < comms_.size(); ++iComm) {
        if (rt += CommuteCofseq(iComm, flag))
            return rt;
    }
    return rt;
}

/* This is for debugging */
int main_test(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"category", &cat_name}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;
    Category category(cat_name, "", flag, false);
    fmt::print("{}\n", category.GetRingByName("S0").nodes_ss.GetSc4Display(AdamsDeg(17, 79 + 17)));

    return 0;
}
