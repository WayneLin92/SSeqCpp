#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

bool Diagram::IsPossTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r_max)
{
    const size_t iCs1 = (iCs + 2) % 3;
    const auto& deg_map = cofseq.degMap[iCs1];
    if (deg.t - deg_map.t > cofseq.t_max[iCs1])
        return true;
    r_max = std::min(r_max, deg.s);
    if (r_max < 0)
        return false;
    const auto& nodes_cofseq_from = cofseq.nodes_cofseq[iCs1];
    const auto& nodes_ss_from = *cofseq.nodes_ss[iCs1];
    for (int r1 = deg_map.s; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, deg_map.stem() + r1};
        if (PossEinf(nodes_ss_from, d_src))
            if (PossMoreEinf(nodes_ss_from, d_src) || GetMaxLevelWithND(ut::GetRecentValue(nodes_cofseq_from, d_src)) >= LEVEL_MAX - r1)
                return true;
    }
    return false;
}

size_t Diagram::GetFirstIndexOfFixedLevelsCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int level_min)
{
    const size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_ss2 = *cofseq.nodes_ss[iCs2];
    const int stem_map = cofseq.degMap[iCs].stem();
    auto& sc = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg);
    size_t result = sc.levels.size();
    for (size_t i = sc.levels.size(); i-- > 0;) {
        if (sc.diffs[i] == NULL_DIFF || sc.levels[i] < level_min)
            break;
        if (i == 0 || sc.levels[i] != sc.levels[i - 1]) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg2 = deg + AdamsDeg{r, r + stem_map};
            if (IsPossTgt(nodes_ss2, deg2) || IsPossTgtCofseq(cofseq, iCs2, deg2, r - 1))
                break;
            else
                result = i;
        }
    }
    return result;
}

int Diagram::GetFirstFixedLevelForPlotCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg)
{
    const size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_ss2 = *cofseq.nodes_ss[iCs2];
    const int stem_map = cofseq.degMap[iCs].stem();
    auto& sc = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg);
    int result = LEVEL_MAX + 1;
    for (size_t i = sc.levels.size(); i-- > 0 && sc.levels[i] >= LEVEL_PERM;) {
        if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg2 = deg + AdamsDeg{r, r + stem_map};
            if (IsPossTgt(nodes_ss2, deg2) || IsPossTgtCofseq(cofseq, iCs2, deg2, r - 1))
                break;
            else
                result = sc.levels[i];
        }
    }
    return result;
}

std::pair<int, int> Diagram::CountPossDrTgtCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_tgt, int r) const
{
    std::pair<int, int> result;
    if (ut::has(cofseq.nodes_cofseq[iCs].front(), deg_tgt)) {
        auto& sc_tgt = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevelsCofseq(cofseq, iCs, deg_tgt, LEVEL_MAX - r + 1) - result.first;
    }
    else if (deg_tgt.t > cofseq.t_max[(iCs + 1) % 3])
        result = {-1, 100000}; /* Infinitely many possibilities */
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Diagram::CountPossDrSrcCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_src, int r) const
{
    std::pair<int, int> result;
    if (ut::has(cofseq.nodes_cofseq[iCs].front(), deg_src)) {
        auto& sc_src = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        result.second = (int)GetFirstIndexOfFixedLevelsCofseq(cofseq, iCs, deg_src, LEVEL_MAX - r + 1) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Diagram::NextRTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const
{
    size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_ss_tgt = *cofseq.nodes_ss[iCs2];
    auto& nodes_cofseq_tgt = cofseq.nodes_cofseq[iCs2];
    int stem_map = cofseq.degMap[iCs].stem();
    int t_max2 = cofseq.t_max[iCs2];
    for (int r1 = r; r1 <= R_PERM; ++r1) {
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 + stem_map};
        if (d_tgt.t > t_max2)
            return r1;
        if (AboveS0Vanishing(d_tgt) && !PossEinf(nodes_ss_tgt, d_tgt))
            return R_PERM;
        if (PossMoreEinf(nodes_ss_tgt, d_tgt))
            return r1;
        auto [_, count] = CountPossDrTgtCofseq(cofseq, iCs2, d_tgt, r1);
        if (count > 0)
            return r1;
    }
    return R_PERM;
}

int Diagram::NextRSrcCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const
{
    int r_max = std::min(r, deg.s);
    size_t iCs1 = (iCs + 2) % 3;
    int t_max1 = cofseq.t_max[iCs1];
    auto& nodes_ss_src = *cofseq.nodes_ss[iCs1];
    int stem_map = cofseq.degMap[iCs1].stem();
    for (int r1 = r_max; r1 >= 0; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 + stem_map};
        if (d_src.t > t_max1)
            return r1;
        if (PossMoreEinf(nodes_ss_src, d_src))
            return r1;
        auto [_, count] = CountPossDrSrcCofseq(cofseq, iCs1, d_src, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

void Diagram::CacheNullDiffsCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, DeduceFlag flag, NullDiff1d& nds) const
{
    nds.clear();
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    size_t iCs1 = (iCs + 2) % 3;
    auto& nodes_ss_src = *cofseq.nodes_ss[iCs1];
    size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_ss_tgt = *cofseq.nodes_ss[iCs2];
    int stem_map1 = cofseq.degMap[iCs1].stem();
    int stem_map = cofseq.degMap[iCs].stem();
    const Staircase& sc = ut::GetRecentValue(nodes_cofseq, deg);
    for (size_t i = 0; i < sc.diffs.size(); ++i) {
        if (sc.diffs[i] != NULL_DIFF)
            continue;

        NullDiff nd;
        if (sc.levels[i] > LEVEL_PERM) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg_tgt = deg + AdamsDeg{r, r + stem_map};
            if (deg_tgt.t > cofseq.t_max[iCs2] || PossMoreEinf(nodes_ss_tgt, deg_tgt))
                continue;
            auto [index, count] = CountPossDrTgtCofseq(cofseq, iCs2, deg_tgt, r);
            nd.direction = 1;
            nd.r = r;
            nd.first = index;
            nd.count = count;
            if (nd.count > deduce_count_max_)
                continue;
        }
        else if (sc.levels[i] < LEVEL_MAX / 2) {
            int r = sc.levels[i];
            AdamsDeg deg_src = deg - AdamsDeg{r, r + stem_map1};
            if (deg_src.t > cofseq.t_max[iCs1] || PossMoreEinf(nodes_ss_src, deg_src))
                continue;
            auto [index, count] = CountPossDrSrcCofseq(cofseq, iCs1, deg_src, r);
            nd.direction = -1;
            nd.r = r;
            nd.first = index;
            nd.count = count;
            if (nd.count > deduce_count_max_)
                continue;
        }
        else
            continue;

        size_t j = i + 1;
        while (j < sc.levels.size() && sc.diffs[j] == NULL_DIFF && sc.levels[i] == sc.levels[j])
            ++j;
        if (flag & DeduceFlag::all_x) {
            const unsigned k_max = unsigned(1) << (j - i);
            for (unsigned k = 1; k < k_max; ++k) {
                nd.x.clear();
                for (int l : two_expansion(k))
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

bool Diagram::IsNewDiffCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const
{
    size_t iCs2 = (iCs + 1) % 3;
    int stem_map = cofseq.degMap[iCs].stem();
    auto& nodes_cofseq_dx = cofseq.nodes_cofseq[iCs2];
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r + stem_map};
    if (x.empty())
        return !dx.empty() && !IsZeroOnLevel(ut::GetRecentValue(nodes_cofseq_dx, deg_dx), dx, r);
    int1d dx1 = GetDiff(cofseq.nodes_cofseq[iCs], deg_x, x, r);  //// TODO: optimize the allocation
    if (dx1 == NULL_DIFF)
        return true;
    int1d diff = lina::add(dx, dx1);
    return !diff.empty() && !IsZeroOnLevel(ut::GetRecentValue(nodes_cofseq_dx, deg_dx), diff, r);
}

int Diagram::SetDiffGlobalCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, DeduceFlag flag)
{
    int count = 0;
    if (newCertain || IsNewDiffCofseq(cofseq, iCs, deg_x, x, dx, r)) {
        int stem_map = cofseq.degMap[iCs].stem();
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r + stem_map);
        if (x.empty()) {
            int r1 = r;
            r = NextRSrcCofseq(cofseq, (iCs + 1) % 3, deg_dx, r1 - 1) + 1;
            if (r != r1)
                deg_x = deg_dx - AdamsDeg(r, r + stem_map);
        }
        else if (dx.empty())
            r = NextRTgtCofseq(cofseq, iCs, deg_x, r + 1) - 1;
        count = SetDiffLeibnizCofseq(cofseq, iCs, deg_x, x, dx, r, flag);
    }
    return count;
}

int Diagram::TryDiffCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, DeduceFlag flag, bool tryY)
{
    AddNode(flag);
    bool bException = false;
    try {
        Logger::LogDiff(depth + 1, tryY ? EnumReason::try1 : EnumReason::try2, fmt::format("{}:{}", cofseq.name, iCs), deg_x, x, dx, r);
        SetDiffGlobalCofseq(cofseq, iCs, deg_x, x, dx, r, true, flag);
        //DeduceTrivialDiffsCofseq(flag);
        //DeduceTrivialDiffs(flag);
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

int Diagram::DeduceDiffsCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg, int depth, DeduceFlag flag)
{
    int count = 0;
    NullDiff1d nds;
    CacheNullDiffsCofseq(cofseq, iCs, deg, flag, nds);
    size_t iCs1 = (iCs + 2) % 3;
    auto& nodes_ss_src = *cofseq.nodes_ss[iCs1];
    size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_ss_tgt = *cofseq.nodes_ss[iCs2];
    int stem_map1 = cofseq.degMap[iCs1].stem();
    int stem_map = cofseq.degMap[iCs].stem();
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_cofseq1 = cofseq.nodes_cofseq[iCs1];
    auto& nodes_cofseq2 = cofseq.nodes_cofseq[iCs2];

    size_t index_nd = 0;
    while (index_nd < nds.size()) {
        const NullDiff nd = nds[index_nd];
        int1d x, dx;
        int r = nd.r;
        bool bNewDiff = false;
        /* Fixed source, find target. */
        AdamsDeg deg_src;
        size_t iCs_src;
        if (nd.direction > 0) {
            iCs_src = iCs;
            deg_src = deg;
            const AdamsDeg deg_tgt = deg_src + AdamsDeg{r, r + stem_map};
            x = nd.x;

            int1d dx1;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;
            for (unsigned i = 1; i < i_max; ++i) {
                const Staircase& sc_tgt = ut::GetRecentValue(nodes_cofseq2, deg_tgt);
                dx1.clear();
                for (int j : two_expansion(i))
                    dx1 = lina::add(dx1, sc_tgt.basis[(size_t)(nd.first + j)]);

                if (!TryDiffCofseq(cofseq, iCs_src, deg_src, x, dx1, r, depth, flag, true)) {
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
                if (TryDiffCofseq(cofseq, iCs_src, deg_src, x, dx1, r, depth, flag, true))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }
        /* Fixed target, find source. */
        else {
            iCs_src = iCs1;
            deg_src = deg - AdamsDeg{r, r + stem_map1};
            const Staircase& sc_src = ut::GetRecentValue(nodes_cofseq1, deg_src);
            dx = nd.x;

            int1d x1;
            int count_pass = 0;
            unsigned i_max = 1 << nd.count;

            for (unsigned i = 1; i < i_max; ++i) {
                x1.clear();
                for (int j : two_expansion(i))
                    x1 = lina::add(x1, sc_src.basis[(size_t)(nd.first + j)]);

                if (!TryDiffCofseq(cofseq, iCs_src, deg_src, x1, dx, r, depth, flag, false)) {
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
                if (TryDiffCofseq(cofseq, iCs_src, deg_src, x1, dx, r, depth, flag, false))
                    bNewDiff = true;
                else
                    ++index_nd;
            }
        }

        if (bNewDiff) {
            ++count;
            if (nd.direction > 0)
                Logger::LogDiff(depth, EnumReason::deduce, fmt::format("{}:{}", cofseq.name, iCs), deg, x, dx, r);
            else
                Logger::LogDiffInv(depth, EnumReason::deduce, fmt::format("{}:{}", cofseq.name, iCs), deg_src, deg, x, dx, r);
            SetDiffGlobalCofseq(cofseq, iCs_src, deg_src, x, dx, r, true, flag);
            int count_trivial = DeduceTrivialDiffsCofseq(flag);
            count += count_trivial;
            CacheNullDiffsCofseq(cofseq, iCs, deg, flag, nds);
        }
        else {
            /*if ((flag & DeduceFlag::xy) && nd.direction > 0) {
                if (iCw < rings_.size())
                    count += SetRingDiffLeibnizV2(iCw, deg, nd.x, nd.r);
                else
                    count += SetModuleDiffLeibnizV2(iCw - rings_.size(), deg, nd.x, nd.r);
            }*/
        }
    }
    return count;
}

int Diagram::DeduceDiffsCofseq(int stem_min, int stem_max, int depth, DeduceFlag flag)
{
    int count = 0;
    DeduceTrivialDiffsCofseq(flag);
    for (size_t iCof : deduce_list_cofseq_) {
        auto& cofseq = cofseqs_[iCof];
        for (size_t iCs = 0; iCs < 3; ++iCs) {
            auto degs = ut::get_keys(cofseq.nodes_cofseq[iCs].front());
            std::stable_sort(degs.begin(), degs.end(), [](AdamsDeg deg1, AdamsDeg deg2) { return deg1.stem() < deg2.stem(); });
            for (AdamsDeg deg : degs) {
                if (depth == 0)
                    fmt::print("{}:{} deg={}                        \r", cofseq.name, iCs, deg);
                if (!BelowS0VanishingLine(deg))
                    continue;
                if (deg.stem() < stem_min)
                    continue;
                if (deg.stem() > stem_max)
                    break;

                count += DeduceDiffsCofseq(cofseq, iCs, deg, depth, flag);
            }
        }
    }
    return count;
}

int Diagram::DeduceDiffsNbhdCofseq(CofSeq& cofseq, size_t iCs_, int stem, int depth, DeduceFlag flag)
{
    int count = 0;
    DeduceTrivialDiffsCofseq(flag);
    size_t iCs1 = (iCs_ + 2) % 3;
    size_t iCs2 = (iCs_ + 1) % 3;
    std::array<size_t, 3> iCss = {iCs1, iCs_, iCs2};
    std::array<int, 3> stems = {stem - cofseq.degMap[iCs1].stem(), stem, stem + cofseq.degMap[iCs_].stem()};

    for (size_t i = 0; i < 3; ++i) {
        size_t iCs = iCss[i];
        int stem = stems[i];
        auto degs = ut::get_keys(cofseq.nodes_cofseq[iCs].front());
        ut::RemoveIf(degs, [stem](AdamsDeg d) { return d.stem() != stem; });
        for (AdamsDeg deg : degs)
            count += DeduceDiffsCofseq(cofseq, iCs, deg, depth, flag);
    }
    return count;
}

void Diagram::SyncCofseq(DeduceFlag flag)
{
    size_t size_rings = rings_.size();
    const size_t num_cw = size_rings + modules_.size();
    for (size_t iCw = 0; iCw < num_cw; ++iCw) {
        auto& nodes_ss = iCw < size_rings ? rings_[iCw].nodes_ss : modules_[iCw - size_rings].nodes_ss;
        const auto& ind_cofs = iCw < size_rings ? rings_[iCw].ind_cofs : modules_[iCw - size_rings].ind_cofs;
        int t_max = iCw < size_rings ? rings_[iCw].t_max : modules_[iCw - size_rings].t_max;
        for (auto& [d, _] : nodes_ss.front()) {
            const Staircase& sc = ut::GetRecentValue(nodes_ss, d);
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] == LEVEL_PERM) {
                    for (auto& ind_cof : ind_cofs) {
                        auto& cofseq = cofseqs_[ind_cof.iCof];
                        auto iCs = (size_t)ind_cof.iCs;
                        if (!ut::has(cofseq.nodes_cofseq[iCs].front(), d))
                            cofseq.nodes_cofseq[iCs].front()[d] = {};
                        SetDiffScCofseq(cofseq, ind_cof.iCs, d, sc.basis[i], NULL_DIFF, 0, flag);
                    }
                }
            }
        }
        for (auto& [d, _] : nodes_ss.front()) {
            for (auto& ind_cof : ind_cofs) {
                auto& cofseq = cofseqs_[ind_cof.iCof];
                auto iCs = (size_t)ind_cof.iCs;
                if (ut::has(cofseq.nodes_cofseq[iCs].front(), d))
                    ReSetScCofseq(cofseq, ind_cof.iCs, d, flag);
            }
        }
    }
}
