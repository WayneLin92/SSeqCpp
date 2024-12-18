#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <fmt/ranges.h>

void pop_front(Staircase& sc)
{
    sc.basis.erase(sc.basis.begin());
    sc.diffs.erase(sc.diffs.begin());
    sc.levels.erase(sc.levels.begin());
}

bool IsPossTgt(const SSNodes& nodes_ss, AdamsDeg deg, int r_max)
{
    r_max = std::min(r_max, deg.s);
    for (int r1 = LEVEL_MIN; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (nodes_ss.front().has(d_src) && GetMaxLevelWithND(nodes_ss.GetRecentSc(d_src)) >= LEVEL_MAX - r1)
            return true;
    }
    return false;
}

size_t GetFirstIndexOfFixedLevels(const SSNodes& nodes_ss, AdamsDeg deg, int level_min)
{
    const auto& sc = nodes_ss.GetRecentSc(deg);
    size_t result = sc.levels.size();
    for (size_t i = sc.levels.size(); i-- > 0;) {
        if (sc.diffs[i] == NULL_DIFF || sc.levels[i] < level_min)
            break;
        if (i == 0 || sc.levels[i] != sc.levels[i - 1]) {
            int r = LEVEL_MAX - sc.levels[i];
            if (IsPossTgt(nodes_ss, deg + AdamsDeg{r, r - 1}, r - 1))
                break;
            else
                result = i;
        }
    }
    return result;
}

int Category::GetFirstFixedLevelForPlot(const SSNodes& nodes_ss, AdamsDeg deg)
{
    const auto& sc = nodes_ss.GetRecentSc(deg);
    int result = LEVEL_MAX - LEVEL_MIN;
    for (size_t i = sc.levels.size(); i-- > 0 && sc.levels[i] >= LEVEL_PERM;) {
        if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
            int r = LEVEL_MAX - sc.levels[i];
            if (IsPossTgt(nodes_ss, deg + AdamsDeg{r, r - 1}, r - 1))
                break;
            else
                result = sc.levels[i];
        }
    }
    return result;
}

std::pair<int, int> Category::CountPossDrTgt(const SSNodes& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const
{
    std::pair<int, int> result;
    if (nodes_ss.front().has(deg_tgt)) {
        const auto& sc_tgt = nodes_ss.GetRecentSc(deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_tgt, LEVEL_MAX - r) - result.first;
    }
    else if (deg_tgt.t > t_max)
        result = {-1, 100000}; /* Infinitely many possibilities */
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Category::CountPossDrSrc(const SSNodes& nodes_ss, const AdamsDeg& deg_src, int r) const
{
    std::pair<int, int> result;
    if (nodes_ss.front().has(deg_src)) {
        const auto& sc_src = nodes_ss.GetRecentSc(deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_src, LEVEL_MAX - r) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Category::NextRTgt(const SSNodes& nodes_ss, int t_max, AdamsDeg deg, int r) const
{
    for (int r1 = r; r1 <= R_PERM; ++r1) {
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 - 1};
        if (d_tgt.t > t_max)
            return r1;
        if (r1 >= 20 && AboveJ(d_tgt) && BelowCokerJ(deg)) /* Image of J */
            return R_PERM;
        if (AboveS0Vanishing(d_tgt) && !nodes_ss.front().has(d_tgt))
            return R_PERM;
        if (CountPossDrTgt(nodes_ss, t_max, d_tgt, r1).second > 0)
            return r1;
    }
    return R_PERM;
}

int Category::NextRSrc(const SSNodes& nodes_ss, AdamsDeg deg, int r) const
{
    int r_max = std::min(r, deg.s);
    for (int r1 = r_max; r1 >= LEVEL_MIN; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (r1 >= 20 && AboveJ(deg) && BelowCokerJ(d_src)) /* Image of J */
            continue;
        if (CountPossDrSrc(nodes_ss, d_src, r1).second > 0)
            return r1;
    }
    return -1;
}

void Category::CacheNullDiffs(const SSNodes& nodes_ss, int t_max, AdamsDeg deg, SSFlag flag, NullDiff1d& nds) const
{
    nds.clear();
    const auto& sc = nodes_ss.GetRecentSc(deg);
    for (size_t i = 0; i < sc.diffs.size(); ++i) {
        if (sc.diffs[i] != NULL_DIFF)
            continue;

        NullDiff nd;
        if (sc.levels[i] > LEVEL_PERM) {
            int r = LEVEL_MAX - sc.levels[i];
            AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};
            nd.r = r;

            auto [index, count] = CountPossDrTgt(nodes_ss, t_max, deg_tgt, r);
            nd.first = index;
            nd.count = count;
            if (nd.count > deduce_count_max_)
                continue;
        }
        else if (sc.levels[i] < LEVEL_MAX / 2) {
            int r = sc.levels[i];
            AdamsDeg deg_src = deg - AdamsDeg{r, r - 1};
            nd.r = -r;

            auto [index, count] = CountPossDrSrc(nodes_ss, deg_src, r);
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
        if (flag & SSFlag::deduce_4_all_x) {
            const unsigned k_max = unsigned(1) << (j - i);
            if (k_max <= 8) {
                for (unsigned k = 1; k < k_max; ++k) {
                    nd.x.clear();
                    for (int l : ut::two_exp(k))
                        nd.x = lina::add(nd.x, sc.basis[i + l]);
                    nds.push_back(nd);
                }
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

/* Return d_r(x) */
int1d GetDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, const int1d& x, int r)
{
    int1d result;
    if (x.empty())
        return result;
    const auto& sc = nodes_ss.GetRecentSc(deg_x);
    size_t first = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    size_t last = GetFirstIndexOfNullOnLevel(sc, LEVEL_MAX - r);
    /* Compute x mod [0,first) */
    int1d x1 = lina::Residue(sc.basis.begin(), sc.basis.begin() + first, x);

    /* If x is in [first,last) */
    if (lina::Residue(sc.basis.begin() + first, sc.basis.begin() + last, x1).empty())
        result = lina::GetImage(sc.basis.begin() + first, sc.basis.begin() + last, sc.diffs.begin() + first, x1);
    else
        result = NULL_DIFF;
    return result;
}

int GetLevel(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x)
{
    int level = -1;
#ifdef MYDEBUG
    MyException::Assert(!x.empty(), "GetDiffAndLevel() para: !x.empty()");
#endif
    const auto& sc = nodes_ss.GetRecentSc(deg_x);
    for (size_t i = 0; i < sc.levels.size(); ++i) {
        if (sc.levels[i] != level)
            level = sc.levels[i];
        if (std::binary_search(x.begin(), x.end(), sc.basis[i].front()))
            x = lina::add(x, sc.basis[i]);
        if (x.empty())
            break;
    }
#ifdef MYDEBUG
    MyException::Assert(x.empty(), "GetDiffAndLevel() end: x.empty()");
#endif
    return level;
}

std::pair<int1d, int> GetDiffAndLevel(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x)
{
    int1d dx;
    int level = -1;
    if (x.empty())
        return std::make_pair(dx, level);

    const auto& sc = nodes_ss.GetRecentSc(deg_x);
    for (size_t i = 0; i < sc.levels.size(); ++i) {
        if (sc.levels[i] != level) {
            level = sc.levels[i];
            dx.clear();
        }
        if (std::binary_search(x.begin(), x.end(), sc.basis[i].front())) {
            x = lina::add(x, sc.basis[i]);
            if (dx != NULL_DIFF) {
                if (sc.diffs[i] != NULL_DIFF)
                    dx = lina::add(dx, sc.diffs[i]);
                else
                    dx = NULL_DIFF;
            }
        }
        if (x.empty())
            break;
    }
#ifdef MYDEBUG
    ErrorIdMsg::Assert(x.empty(), "GetDiffAndLevel() end: x.empty()");
#endif
    return std::make_pair(dx, level);
}

bool Category::IsNewDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    if (x.empty())
        return !dx.empty() && !IsZeroOnLevel(nodes_ss.GetRecentSc(deg_dx), dx, r);
    int1d dx1 = GetDiff(nodes_ss, deg_x, x, r);  //// TODO: optimize the allocation
    if (dx1 == NULL_DIFF)
        return true;
    int1d diff = lina::add(dx, dx1);
    return !diff.empty() && !IsZeroOnLevel(nodes_ss.GetRecentSc(deg_dx), diff, r);
}

/* Return the minimal length of the crossing differentials */
int GetCrossR(const SSNodes& nodes_ss, AdamsDeg deg, int t_max, int Er)
{
    int result = R_PERM;
    for (int r = 1; r <= Er - 2; ++r) {
        auto deg_x = deg + AdamsDeg{r, r};
        if (r + LEVEL_MIN > result)
            return result;
        if (deg_x.t > t_max)
            return std::min(result, r + LEVEL_MIN);
        if (!nodes_ss.front().has(deg_x)) {
            if (AboveS0Vanishing(deg_x))
                return result;
            continue;
        }
        auto& sc = nodes_ss.GetRecentSc(deg_x);
        if (!sc.levels.empty() && sc.levels.back() > LEVEL_MAX / 2) {
            int r1 = LEVEL_MAX - sc.levels.back();
            if (r + r1 < result)
                result = r + r1;
        }
    }
    return result;
}

SSRet Category::SetDiffSc(IndexUniv iCw, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, SSFlag flag)
{
    // if (iCw == IndexRing(0) && deg_x == AdamsDeg(25, 77 + 25))
    //     fmt::print("debug r={}\n", r);  //////////////
    SSRet rt;
    int r_original = r;

    /* NULL_DIFF checking */
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    if (x_ == NULL_DIFF) {
        rt += SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r, flag);
        return rt;
    }

    /* Triangularize x */
    if (x_.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            rt += SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return rt;
    }
    auto& nodes_ss = GetNodesSS(iCw);
    const auto& sc = nodes_ss.GetRecentSc(deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    int1d x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, x_);
    if (x.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            rt += SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return rt;
    }

    int1d image_new;
    int level_image_new = -1;
    if (dx == NULL_DIFF) {
        /* If the target is uncertain, insert it to the end of level N-r. */
        size_t first_Nmrp1 = GetFirstIndexOnLevel(sc, LEVEL_MAX - r + 1);
        x = lina::Residue(sc.basis.begin() + first_Nmr, sc.basis.begin() + first_Nmrp1, x);
        if (!x.empty())
            UpdateStaircase(nodes_ss, deg_x, sc, first_Nmrp1, x, NULL_DIFF, LEVEL_MAX - r, image_new, level_image_new);
    }
    else if (dx.empty()) {
        /* If the target is zero, insert it to the end of level N-r-1 */
        ++r;
        UpdateStaircase(nodes_ss, deg_x, sc, first_Nmr, x, NULL_DIFF, LEVEL_MAX - r, image_new, level_image_new);
    }
    else {
        /* Otherwise insert it to the beginning of level N-r */
        UpdateStaircase(nodes_ss, deg_x, sc, first_Nmr, x, dx, LEVEL_MAX - r, image_new, level_image_new);
    }

    /* ss to cofseq */
    /*if (iUniv == IndexMod(0) && deg_x == AdamsDeg(5, 34 + 5) && r == R_PERM)
        fmt::print("debug\n");*/
    if (r == R_PERM && (flag & SSFlag::cofseq)) {
        for (auto& iCof : GetIndexCofs(iCw)) {
            auto& cofseq = cofseqs_[iCof.index];
            if (!cofseq.nodes_cofseq[iCof.iTri].front().has(deg_x))
                cofseq.nodes_cofseq[iCof.iTri].front()[deg_x] = {};
            if (rt += SetDiffScCofseq(cofseq, iCof.iTri, deg_x, x, NULL_DIFF, 0, flag))
                return rt;
        }
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_x + AdamsDeg{level_image_new, level_image_new - 1};
            if (rt += SetImageSc(iCw, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]` where the element on the left-hand side already supports d_{}.\n{}", GetCwName(iCw), deg_x, r_original, myio::Serialize(x_), myio::Serialize(dx), level_image_new, rt.err_msg);
                return rt;
            }
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, r_image - 1};
            if (rt += SetDiffSc(iCw, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]` where the element on the left-hand side is already hit by d_{}.\n{}", GetCwName(iCw), deg_x, r_original, myio::Serialize(x_), myio::Serialize(dx), level_image_new, rt.err_msg);
                return rt;
            }
        }
    }

    /* Set the image */
    if (!dx.empty() && dx != NULL_DIFF)
        if (rt += SetImageSc(iCw, deg_dx, dx, x, r, flag))
            return rt;
    return rt;
}

SSRet Category::SetImageSc(IndexUniv iCw, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, SSFlag flag)
{
    // if (iCw == IndexRing(0) && deg_dx == AdamsDeg(36, 84 + 36))
    //     fmt::print("debug r={}\n", r);  //////////////
    SSRet rt;
    const auto& name = GetCwName(iCw);
    auto& nodes_ss = GetNodesSS(iCw);
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, r - 1};

    /* If dx is in Im(d_{r-1}) then x is in Ker(d_r) */
    const auto& sc = nodes_ss.GetRecentSc(deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    int1d dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, dx_);
    if (dx.empty()) {
        if (x != NULL_DIFF)
            rt += SetDiffSc(iCw, deg_x, x, NULL_DIFF, r + 1, flag);
        return rt;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == NULL_DIFF) {
        /* If the source is uncertain, check if it can be hit and then insert it to the end of level r. */
        size_t first_rp1 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp1, std::move(dx));
        if (!dx.empty()) {
            if (!IsPossTgt(nodes_ss, deg_dx, r)) {
                /*if (deg_dx == AdamsDeg(6, 42 + 6) && name == "Cnu" && r == 2)
                    fmt::print("debug\n");*/
                if (depth_ == 0)
                    throw RunTimeError(fmt::format("Contradiction name={} deg_dx={} dx={} r={}", name, deg_dx, myio::Serialize(dx), r));
                rt.code = SSRet::FAIL_SS().code;
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("However, `{} {} [{}]` is not in B_{}.", name, deg_dx, myio::Serialize(dx), r);
                return rt;
            }
            UpdateStaircase(nodes_ss, deg_dx, sc, first_rp1, dx, x, r, image_new, level_image_new);
        }
    }
    else {
        /* Otherwise insert it to the beginning of level r */
        UpdateStaircase(nodes_ss, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_dx + AdamsDeg{level_image_new, level_image_new - 1};
            if (rt += SetImageSc(iCw, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{} {} [{}]=d_{}[{}]` where the element on the left-hand side already supports d_{}.\n{}", GetCwName(iCw), deg_dx, myio::Serialize(dx_), r, myio::Serialize(x), level_image_new, rt.err_msg);
                return rt;
            }
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, r_image - 1};
            if (rt += SetDiffSc(iCw, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{} {} [{}]=d_{}[{}]` where the element on the left-hand side is already hit by d_{}.\n{}", GetCwName(iCw), deg_dx, myio::Serialize(dx_), r, myio::Serialize(x), r_image, rt.err_msg);
                return rt;
            }
        }
    }

    return rt;
}

SSRet Category::SetDiffScCofseq(CofSeq& cofseq, size_t iTri, AdamsDeg deg_x, const int1d& x_, const int1d& dx_, int r, SSFlag flag)
{
    /*if (depth_ == 0 && cofseq.name == "S0__Cnu__S0" && iTri == 2 && deg_x == AdamsDeg(21, 111 + 21) && r >= 9 && x_ == int1d{0}) {
        fmt::print("debug\n");
    }*/
    SSRet rt;
    const auto iTri_prev = PreviTri(iTri);
    const auto iTri_next = NextiTri(iTri);
    auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
    auto& nodes_ss = *cofseq.nodes_ss[iTri];
    int stem_map = cofseq.degMap[iTri].stem();
    int stem_map_prev = cofseq.degMap[iTri_prev].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, stem_map + r};

    /* NULL_DIFF checking */
    if (x_ == NULL_DIFF) {
        rt += SetImageScCofseq(cofseq, iTri_next, deg_dx, dx_, NULL_DIFF, r, flag);
        return rt;
    }

    /* Triangularize x */
    if (x_.empty()) {
        if (dx_ != NULL_DIFF && !dx_.empty())
            rt += SetImageScCofseq(cofseq, iTri_next, deg_dx, dx_, NULL_DIFF, r - 1, flag);
        return rt;
    }
    int1d x = Residue(x_, nodes_ss, deg_x, LEVEL_PERM);
    int1d dx = dx_;
    if (dx_ != NULL_DIFF && !dx_.empty())
        dx = Residue(std::move(dx), *cofseq.nodes_ss[iTri_next], deg_dx, LEVEL_PERM);
    const auto& sc = nodes_cofseq.GetRecentSc(deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, std::move(x));
    if (x.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            rt += SetImageScCofseq(cofseq, iTri_next, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return rt;
    }

    int1d image_new;
    int level_image_new = -1;
    if (dx == NULL_DIFF) {
        /* If the target is uncertain, insert it to the end of level N-r. */
        size_t first_Nmrp1 = GetFirstIndexOnLevel(sc, LEVEL_MAX - r + 1);
        x = lina::Residue(sc.basis.begin() + first_Nmr, sc.basis.begin() + first_Nmrp1, std::move(x));
        if (!x.empty())
            UpdateStaircase(nodes_cofseq, deg_x, sc, first_Nmrp1, x, NULL_DIFF, LEVEL_MAX - r, image_new, level_image_new);
    }
    else if (dx.empty()) {
        /* If the target is zero, insert it to the end of level N-r-1 */
        UpdateStaircase(nodes_cofseq, deg_x, sc, first_Nmr, x, NULL_DIFF, LEVEL_MAX - r - 1, image_new, level_image_new);
    }
    else {
        /* Otherwise insert it to the beginning of level N-r */
        UpdateStaircase(nodes_cofseq, deg_x, sc, first_Nmr, x, dx, LEVEL_MAX - r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_x + AdamsDeg{level_image_new, stem_map + level_image_new};
            if (rt += SetImageScCofseq(cofseq, iTri_next, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{}:{} {} d_{}[{}]=[{}]` where the element on the left-hand side already supports d_{}.\n{}", cofseq.name, iTri, deg_x, r, myio::Serialize(x_), myio::Serialize(dx_), level_image_new, rt.err_msg);
                return rt;
            }
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, stem_map_prev + r_image};
            if (rt += SetDiffScCofseq(cofseq, iTri_prev, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{}:{} {} d_{}[{}]=[{}]` where the element on the left-hand side is already hit by d_{}.\n{}", cofseq.name, iTri, deg_x, r, myio::Serialize(x_), myio::Serialize(dx_), r_image, rt.err_msg);
                return rt;
            }
        }
    }

    /* Set image */
    if (dx != NULL_DIFF && !dx.empty())
        if (rt += SetImageScCofseq(cofseq, iTri_next, deg_dx, dx, x, r, flag))
            return rt;

    /*if (cofseq.name == "S0__Ceta__S0" && iCs == 1 && deg_x == AdamsDeg(29, 176)) {
        auto& sc = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg_x);
        fmt::print("sc at {} =\n{}\n", deg_x, sc);
        AdamsDeg deg_next = deg_x + AdamsDeg{0, cofseq.degMap[iCs].stem() + 0};
        fmt::print("sc_next at {} =\n{}\n", deg_next, ut::GetRecentValue(cofseq.nodes_cofseq[(iCs + 1) % 3], deg_next));
        std::cout << "debug\n";
    }*/
    return rt;
}

SSRet Category::SetImageScCofseq(CofSeq& cofseq, size_t iTri, AdamsDeg deg_dx, const int1d& dx_, const int1d& x_, int r, SSFlag flag)
{
    SSRet rt;
    const auto iTri_prev = PreviTri(iTri);
    auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
    auto& nodes_ss = *cofseq.nodes_ss[iTri];
    int stem_map = cofseq.degMap[iTri].stem();
    int stem_map_prev = cofseq.degMap[iTri_prev].stem();
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, stem_map_prev + r};

    int1d dx = Residue(dx_, nodes_ss, deg_dx, LEVEL_PERM);
    int1d x = x_;
    if (x_ != NULL_DIFF && !x_.empty())
        x = Residue(std::move(x), *cofseq.nodes_ss[iTri_prev], deg_x, LEVEL_PERM);
#ifdef MYDEBUG
    if (!nodes_cofseq.front().has(deg_dx)) {
        fmt::print("{}:{} deg_dx={} dx_={} x_={} r={}\n", cofseq.name, iTri, deg_dx, dx_, x_, r);
        fmt::print("sc_ss:\n{}\n", nodes_ss.GetRecentSc(deg_dx));
    }
#endif
    const auto& sc = nodes_cofseq.GetRecentSc(deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, std::move(dx));
    if (dx.empty()) {
        if (x != NULL_DIFF)
            rt += SetDiffScCofseq(cofseq, iTri_prev, deg_x, x, NULL_DIFF, r + 1, flag);
        if (rt && (flag & SSFlag::log_proof))
            rt.err_msg = fmt::format("Get `{}:{} {} [{}]=d_{}[{}]`.\n{}", cofseq.name, iTri, deg_dx, myio::Serialize(dx_), r, myio::Serialize(x_), rt.err_msg);
        return rt;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == NULL_DIFF) {
        /* If the source is uncertain, insert it to the end of level r. */
        size_t first_rp1 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp1, std::move(dx));
        if (!dx.empty()) {
            if (!IsPossTgtCofseq(cofseq, iTri, deg_dx, r)) {
                r = -1;
                first_rp1 = 0;
            }
            UpdateStaircase(nodes_cofseq, deg_dx, sc, first_rp1, dx, x, r, image_new, level_image_new);

            /* cofseq to ss */
            if (r == -1) {
                pop_front(nodes_cofseq.back()[deg_dx]);
                Logger::LogDiffInv(depth_, EnumReason::cofseq_out, cofseq.nameCw[iTri], deg_dx, R_PERM + 1, {}, dx, "", flag);
                if (depth_ == 0) {
                    if (rt += SetImageSc(cofseq.indexCw[iTri], deg_dx, dx, NULL_DIFF, R_PERM, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("For degree reason `{} {} [{}]` is an Adams boundary.\n{}", cofseq.nameCw[iTri], deg_dx, myio::Serialize(dx_), rt.err_msg);
                        return rt;
                    }
                }
                else {
                    // if (flag & SSFlag::deduce_zero) {
                    //     if (rt += DeduceDiff4XDepth(cofseq.indexCw[iTri], deg_dx, dx, LEVEL_PERM, flag))
                    //         return rt;
                    // }
                    // else {
                    if (rt += SetCwDiffGlobal(cofseq.indexCw[iTri], deg_dx - AdamsDeg(R_PERM, R_PERM - 1), {}, dx, R_PERM, true, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("For degree reason `{} {} [{}]` is an Adams boundary.\n{}", cofseq.nameCw[iTri], deg_dx, myio::Serialize(dx_), rt.err_msg);
                        return rt;
                    }
                    //}
                }
            }
        }
    }
    else {
        /* Otherwise insert it to the beginning of level r */  //// TODO: Improve. Insert it to the end of known level r
        UpdateStaircase(nodes_cofseq, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_dx + AdamsDeg{level_image_new, stem_map + level_image_new};
            if (rt += SetImageScCofseq(cofseq, NextiTri(iTri), deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{}:{} {} [{}]=d_{}[{}]` where the element on the left-hand side already supports d_{}.\n{}", cofseq.name, iTri, deg_dx, myio::Serialize(dx_), r, myio::Serialize(x_), level_image_new, rt.err_msg);
                return rt;
            }
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, stem_map_prev + r_image};
            if (rt += SetDiffScCofseq(cofseq, iTri_prev, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag)) {
                if (flag & SSFlag::log_proof)
                    rt.err_msg = fmt::format("Get `{}:{} {} [{}]=d_{}[{}]` where the element on the left-hand side is already hit by d_{}.\n{}", cofseq.name, iTri, deg_dx, myio::Serialize(dx_), r, myio::Serialize(x_), r_image, rt.err_msg);
                return rt;
            }
        }
    }
    return rt;
}

SSRet Category::SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag)  ////
{
    SSRet rt;
    auto& ring = rings_[iRing];
    Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x)), poly_zero;
    Poly poly_drx = !dx.empty() ? Indices2Poly(dx, ring.basis.at(deg_x + AdamsDeg(r, r - 1))) : poly_zero;
    Poly poly_a, poly_ax, poly_da, poly_dax, poly_tmp1, poly_tmp2;
    int1d ax, dax;
    Mod mod_y, mod_xy, mod_dy, mod_dxy, mod_tmp1, mod_tmp2;
    int1d xy, dxy;

    /* maximum r_zero such that d_{r_zero}(x)=0 */
    int r_zero = r + 1 == R_PERM ? R_PERM : (dx.empty() ? r : r - 1);
    {
        auto& nodes_ss = ring.nodes_ss;
        auto& degs_ss = ring.degs_ss;
        auto& basis = ring.basis;
        int t_max = ring.t_max;
        for (AdamsDeg deg_a : degs_ss) {
            /*if (deg_a + deg_x + AdamsDeg(r, r - 1) == AdamsDeg(49, 131 + 49))
                fmt::print("debug\n");*/
            const AdamsDeg deg_ax = deg_x + deg_a;
            const auto& sc_a = nodes_ss.GetRecentSc(deg_a);
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                /* ax=d_R[?], x is d_R-cycle, a is d_R-boudnary */
                if (sc_a.levels[i] >= r_min && sc_a.levels[i] <= r_zero && sc_a.diffs[i] == NULL_DIFF) {
                    const int R = sc_a.levels[i];
                    Indices2AlgP(sc_a.basis[i], basis.at(deg_a), poly_a);
                    mulP(poly_a, poly_x, poly_ax); /* poly_ax = poly_a * poly_x */
                    ring.gb.ReduceP(poly_ax, poly_tmp1, poly_tmp2);
                    if (poly_ax) {
                        Alg2IndicesP(poly_ax, basis.at(deg_ax), ax);
                        if (rt += SetImageSc(IndexRing(iRing), deg_ax, ax, NULL_DIFF, R, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[]`. Apply the Leibniz rule with `{} {} [{}]=d_{}[?]` and get `{} {} [{}]=d_{}[?]`.\n{}", ring.name, deg_x, R, myio::Serialize(x), ring.name, deg_a, myio::Serialize(sc_a.basis[i]), R,
                                                         ring.name, deg_ax, myio::Serialize(ax), R, rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                    }
                }
                else if (sc_a.levels[i] >= r_min || sc_a.diffs[i] == NULL_DIFF) {  //// TODO: deal with dx that was originally not in Z_{level_a - 1}. This is actually rare.
                    const int r_a = LEVEL_MAX - sc_a.levels[i] - (sc_a.diffs[i] == NULL_DIFF ? 1 : 0);
                    if (r_a < r_min)
                        break;
                    const int R = std::min(r, r_a);
                    AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                    AdamsDeg deg_dax = deg_ax + AdamsDeg(R, R - 1);

                    const Poly& poly_dx = (R == r) ? poly_drx : poly_zero;
                    Indices2AlgP(sc_a.basis[i], basis.at(deg_a), poly_a);
                    mulP(poly_a, poly_x, poly_ax); /* poly_ax = poly_a * poly_x */
                    ring.gb.ReduceP(poly_ax, poly_tmp1, poly_tmp2);
                    Alg2IndicesP(poly_ax, basis, deg_ax, ax);

                    if (R == LEVEL_MAX - sc_a.levels[i])
                        Indices2AlgP(sc_a.diffs[i], basis.at(deg_a + AdamsDeg(R, R - 1)), poly_da);
                    else
                        poly_da.data.clear();
                    /* poly_dax = poly_a * poly_dx + poly_da * poly_x */
                    mulP(poly_a, poly_dx, poly_dax);
                    mulP(poly_da, poly_x, poly_tmp1);
                    poly_dax.iaddP(poly_tmp1, poly_tmp2);
                    ring.gb.ReduceP(poly_dax, poly_tmp1, poly_tmp2);
                    if (poly_dax) {
                        if (deg_dax.t <= t_max)
                            Alg2IndicesP(poly_dax, basis.at(deg_dax), dax);
                        else
                            dax = NULL_DIFF;
                    }
                    else
                        dax.clear();

                    if (!ax.empty() || !dax.empty()) {
                        if (rt += SetDiffSc(IndexRing(iRing), deg_ax, ax, dax, R, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]` and get `{} {} d_{}[{}]=[{}]`.\n{}", ring.name, deg_x, R, myio::Serialize(x), myio::Serialize(dx), ring.name, deg_a, R,
                                                         myio::Serialize(sc_a.basis[i]), poly_da ? myio::Serialize(sc_a.diffs[i]) : std::string{}, ring.name, deg_ax, R, myio::Serialize(ax), myio::Serialize(dax), rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                    }
                }
            }
        }
    }
    for (size_t iMod : ring.ind_mods) {
        auto& mod = modules_[iMod];
        auto& nodes_ss = mod.nodes_ss;
        auto& degs_ss = mod.degs_ss;
        auto& basis = mod.basis;
        int t_max = mod.t_max;
        for (AdamsDeg deg_y : degs_ss) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const auto& sc_y = nodes_ss.GetRecentSc(deg_y);
            if (deg_xy.t > t_max)
                break;
            for (size_t i = 0; i < sc_y.levels.size(); ++i) {
                if (sc_y.levels[i] < LEVEL_MAX / 2 && sc_y.diffs[i] == NULL_DIFF && sc_y.levels[i] >= r_min && sc_y.levels[i] <= r_zero) { /* xy=d_R[?] */
                    const int R = sc_y.levels[i];
                    Indices2AlgP(sc_y.basis[i], basis.at(deg_y), mod_y);
                    mulP(poly_x, mod_y, mod_xy); /* mod_xy = poly_x * mod_y */
                    mod.gb.ReduceP(mod_xy, poly_tmp1, mod_tmp1, mod_tmp2);
                    if (mod_xy) {
                        Alg2IndicesP(mod_xy, basis.at(deg_xy), xy);
                        if (rt += SetImageSc(IndexMod(iMod), deg_xy, xy, NULL_DIFF, R, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[]`. Apply the Leibniz rule with `{} {} [{}]=d_{}[?]` and get `{} {} [{}]=d_{}[?]`.\n{}", ring.name, deg_x, R, myio::Serialize(x), mod.name, deg_y, myio::Serialize(sc_y.basis[i]), R,
                                                         mod.name, deg_xy, myio::Serialize(xy), R, rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                    }
                }
                else if (sc_y.levels[i] >= r_min || sc_y.diffs[i] == NULL_DIFF) {
                    const int r_y = LEVEL_MAX - sc_y.levels[i] - (sc_y.diffs[i] == NULL_DIFF ? 1 : 0);
                    if (r_y < r_min)
                        break;
                    const int R = std::min(r, r_y);
                    AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                    AdamsDeg deg_dxy = deg_xy + AdamsDeg(R, R - 1);

                    const Poly& poly_dx = (R == r) ? poly_drx : poly_zero;
                    Indices2AlgP(sc_y.basis[i], basis.at(deg_y), mod_y);

                    mulP(poly_x, mod_y, mod_xy); /* mod_xy = poly_x * mod_y */
                    mod.gb.ReduceP(mod_xy, poly_tmp1, mod_tmp1, mod_tmp2);
                    Alg2IndicesP(mod_xy, basis, deg_x + deg_y, xy);

                    if (R == LEVEL_MAX - sc_y.levels[i])
                        Indices2AlgP(sc_y.diffs[i], basis.at(deg_y + AdamsDeg(R, R - 1)), mod_dy);
                    else
                        mod_dy.data.clear();
                    /* mod_dxy = poly_y * poly_dx + poly_dy * poly_x */
                    mulP(poly_dx, mod_y, mod_dxy);
                    mulP(poly_x, mod_dy, mod_tmp1);
                    mod_dxy.iaddP(mod_tmp1, mod_tmp2);
                    mod.gb.ReduceP(mod_dxy, poly_tmp1, mod_tmp1, mod_tmp2);
                    if (mod_dxy) {
                        if (deg_dxy.t <= t_max)
                            Alg2IndicesP(mod_dxy, basis.at(deg_dxy), dxy);
                        else
                            dxy = NULL_DIFF;
                    }
                    else
                        dxy.clear();

                    if (!xy.empty() || !dxy.empty()) {
                        if (rt += SetDiffSc(IndexMod(iMod), deg_xy, xy, dxy, R, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]` and get `{} {} d_{}[{}]=[{}]`.\n{}", ring.name, deg_x, R, myio::Serialize(x), myio::Serialize(dx), mod.name, deg_y, R,
                                                         myio::Serialize(sc_y.basis[i]), mod_dy ? myio::Serialize(sc_y.diffs[i]) : std::string{}, mod.name, deg_xy, R, myio::Serialize(xy), myio::Serialize(dxy), rt.err_msg);
                            return rt;
                        }
                        rt += SSRet::CHANGE();
                    }
                }
            }
        }
    }
    return rt;
}

SSRet Category::SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag)  ////
{
    SSRet rt;
    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& basis = mod.basis;
    auto& gb = mod.gb;
    const int t_max = mod.t_max;
    Mod poly_x = Indices2Mod(x, basis.at(deg_x)), poly_zero;
    Mod poly_drx = !dx.empty() ? Indices2Mod(dx, basis.at(deg_x + AdamsDeg(r, r - 1))) : poly_zero;
    int r_zero = r + 1 == R_PERM ? R_PERM : (dx.empty() ? r : r - 1);
    Poly poly_a, poly_da, poly_tmp1;
    Mod mod_ax, mod_dax, mod_tmp1, mod_tmp2;
    int1d ax, dax;

    for (AdamsDeg deg_a : ring.degs_ss) {
        const auto& sc_a = ring.nodes_ss.GetRecentSc(deg_a);
        AdamsDeg deg_ax = deg_x + deg_a;
        if (deg_ax.t > t_max)
            break;
        for (size_t i = 0; i < sc_a.levels.size(); ++i) {
            /*if (deg_a == AdamsDeg(16, 56 + 16) && sc_a.basis[i] == int1d{0})
                fmt::print("debug\n");*/
            if (sc_a.levels[i] < LEVEL_MAX / 2 && sc_a.diffs[i] == NULL_DIFF && sc_a.levels[i] >= r_min && sc_a.levels[i] <= r_zero) { /* ax=d_R[?] */
                const int R = sc_a.levels[i];
                Indices2AlgP(sc_a.basis[i], ring.basis.at(deg_a), poly_a);
                mulP(poly_a, poly_x, mod_ax); /* poly_ax = poly_a * poly_x */
                mod.gb.ReduceP(mod_ax, poly_tmp1, mod_tmp1, mod_tmp2);
                if (mod_ax) {
                    Alg2IndicesP(mod_ax, basis.at(deg_ax), ax);
                    if (rt += SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, R, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Apply the Leibniz rule with `{} {} [{}]=d_{}[?]` and get `{} {} [{}]=d_{}[?]`.\n{}", mod.name, deg_x, R, myio::Serialize(x), myio::Serialize(dx), ring.name, deg_a,
                                                     myio::Serialize(sc_a.basis[i]), R, mod.name, deg_ax, myio::Serialize(ax), R, rt.err_msg);
                        return rt;
                    }
                    rt += SSRet::CHANGE();
                }
            }
            else if (sc_a.levels[i] >= r_min || sc_a.diffs[i] == NULL_DIFF) {
                const int r_a = LEVEL_MAX - sc_a.levels[i] - (sc_a.diffs[i] == NULL_DIFF ? 1 : 0);
                if (r_a < r_min)
                    break;
                int R = std::min(r, r_a);
                AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                AdamsDeg deg_dax = deg_ax + AdamsDeg(R, R - 1);

                const Mod& poly_dx = (R == r) ? poly_drx : poly_zero;
                Indices2AlgP(sc_a.basis[i], ring.basis.at(deg_a), poly_a);
                mulP(poly_a, poly_x, mod_ax); /* mod_ax = poly_a * mod_x */
                gb.ReduceP(mod_ax, poly_tmp1, mod_tmp1, mod_tmp2);
                Alg2IndicesP(mod_ax, basis, deg_x + deg_a, ax);

                if (R == LEVEL_MAX - sc_a.levels[i])
                    Indices2AlgP(sc_a.diffs[i], ring.basis.at(deg_a + AdamsDeg(R, R - 1)), poly_da);
                else
                    poly_da.data.clear();
                /* mod_dax = poly_a * mod_dx + poly_da * mod_x */
                mulP(poly_a, poly_dx, mod_dax);
                mulP(poly_da, poly_x, mod_tmp1);
                mod_dax.iaddP(mod_tmp1, mod_tmp2);
                gb.ReduceP(mod_dax, poly_tmp1, mod_tmp1, mod_tmp2);
                if (mod_dax) {
                    if (deg_dax.t <= t_max)
                        Alg2IndicesP(mod_dax, basis.at(deg_dax), dax);
                    else
                        dax = NULL_DIFF;
                }
                else
                    dax.clear();

                if (!ax.empty() || !dax.empty()) {
                    if (rt += SetDiffSc(IndexMod(iMod), deg_ax, ax, dax, R, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]` and get `{} {} d_{}[{}]=[{}]`.\n{}", mod.name, deg_x, R, myio::Serialize(x), myio::Serialize(dx), ring.name, deg_a, R,
                                                     myio::Serialize(sc_a.basis[i]), poly_da ? myio::Serialize(sc_a.diffs[i]) : "", mod.name, deg_ax, R, myio::Serialize(ax), myio::Serialize(dax), rt.err_msg);
                        return rt;
                    }
                    rt += SSRet::CHANGE();
                }
            }
        }
    }
    return rt;
}

SSRet Category::SetRingDiffLeibnizV2(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    SSRet rt;
    auto& ring = rings_[iRing];
    const AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    if (!ring.nodes_ss.front().has(deg_dx))
        return rt;
    const Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    const auto* sc_dx = &ring.nodes_ss.GetRecentSc(deg_dx);
    auto [first_dx, count_dx] = CountPossDrTgt(ring.nodes_ss, ring.t_max, deg_dx, r);
    if (count_dx == 0 || count_dx > deduce_count_max_)
        return rt;
    int2d* exclusions = nullptr;
    if (auto pE = ring.nodes_ss.exclusions.has(deg_x, x, r))
        exclusions = &pE->dxs;
    else
        exclusions = nullptr;

    unsigned j_max = 1 << count_dx;
    int1d dx;
    bool printed_dx = false;
    std::string proof;
    {
        auto& nodes_ss = ring.nodes_ss;
        const auto& basis = ring.basis;
        const int t_max = ring.t_max;
        for (AdamsDeg deg_y : ring.degs_ss) {
            const AdamsDeg deg_xy = deg_x + deg_y;
            const AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
            if (deg_dxy.t > t_max)
                break;
            const auto* sc_y = &nodes_ss.GetRecentSc(deg_y);
            for (size_t i = 0; i < sc_y->levels.size(); ++i) { /* Loop over y */
                const int r_y = LEVEL_MAX - sc_y->levels[i] - (sc_y->diffs[i] == NULL_DIFF ? 1 : 0);
                if (r_y < r)
                    break;

                Poly poly_y = Indices2Poly(sc_y->basis[i], basis.at(deg_y));
                Poly poly_xy = ring.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = Poly2Indices(poly_xy, basis.at(deg_xy));

                    bool ydx_constant = true;
                    int1d ydx_first_value = NULL_DIFF;
                    for (unsigned j = 0; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : ut::two_exp(j))
                            dx = lina::add(dx, sc_dx->basis[(size_t)(first_dx + k)]);
                        if (exclusions && ut::has(*exclusions, dx))
                            continue;
                        if (j == 0) {
                            ydx_first_value = {};
                            continue;
                        }
                        Poly poly_dx = Indices2Poly(dx, basis.at(deg_dx));
                        Poly poly_ydx = ring.gb.Reduce(poly_dx * poly_y);
                        int1d ydx = poly_ydx ? Residue(Poly2Indices(poly_ydx, basis.at(deg_dxy)), nodes_ss, deg_dxy, r) : int1d{};
                        if (ydx_first_value == NULL_DIFF)
                            ydx_first_value = std::move(ydx);
                        else if (ydx_first_value != ydx) {
                            ydx_constant = false;
                            break;
                        }
                    }
                    if (ydx_constant) {
                        Poly poly_dy = (r == LEVEL_MAX - sc_y->levels[i]) ? Indices2Poly(sc_y->diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
                        Poly poly_xdy = ring.gb.Reduce(poly_x * poly_dy);
                        int1d xdy = poly_xdy ? Poly2Indices(poly_xdy, basis.at(deg_dxy)) : int1d{};
                        int1d dxy = lina::add(xdy, ydx_first_value);

                        if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                            rt += SSRet::CHANGE();
                            if (!printed_dx) {
                                printed_dx = true;
                                Logger::LogNullDiff(depth_, ring.name, deg_x, r, x, flag);
                            }
                            if (flag & SSFlag::log_proof) {
                                proof.clear();
                                proof = fmt::format("Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]`.", ring.name, deg_y, r, myio::Serialize(sc_y->basis[i]), poly_dy ? myio::Serialize(sc_y->diffs[i]) : "");
                            }
                            Logger::LogDiff(depth_, EnumReason::deduce_xy, ring.name, deg_xy, r, xy, dxy, proof, flag);
                            if (rt += SetRingDiffGlobal(iRing, deg_xy, xy, dxy, r, true, flag))
                                return rt;
                            sc_dx = &ring.nodes_ss.GetRecentSc(deg_dx);
                            sc_y = &nodes_ss.GetRecentSc(deg_y);
                        }
                    }
                }
            }
        }
    }
    for (size_t iMod : ring.ind_mods) {
        auto& mod = modules_[iMod];
        auto& nodes_ss = mod.nodes_ss;
        auto& basis = mod.basis;
        int t_max = mod.t_max;
        for (AdamsDeg deg_y : mod.degs_ss) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
            if (deg_dxy.t > t_max)
                break;
            const auto* sc_y = &nodes_ss.GetRecentSc(deg_y);
            for (size_t i = 0; i < sc_y->levels.size(); ++i) { /* Loop over y */
                const int r_y = LEVEL_MAX - sc_y->levels[i] - (sc_y->diffs[i] == NULL_DIFF ? 1 : 0);
                if (r_y < r)
                    break;

                Mod poly_y = Indices2Mod(sc_y->basis[i], basis.at(deg_y));
                Mod poly_xy = mod.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    bool ydx_constant = true;
                    int1d ydx_first_value = NULL_DIFF;
                    for (unsigned j = 0; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : ut::two_exp(j))
                            dx = lina::add(dx, sc_dx->basis[(size_t)(first_dx + k)]);
                        if (exclusions && ut::has(*exclusions, dx))
                            continue;
                        if (j == 0) {
                            ydx_first_value = {};
                            continue;
                        }

                        Poly poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
                        Mod poly_ydx = mod.gb.Reduce(poly_dx * poly_y);
                        int1d ydx = poly_ydx ? Residue(Mod2Indices(poly_ydx, basis.at(deg_dxy)), nodes_ss, deg_dxy, r) : int1d{};
                        if (ydx_first_value == NULL_DIFF)
                            ydx_first_value = std::move(ydx);
                        else if (ydx_first_value != ydx) {
                            ydx_constant = false;
                            break;
                        }
                    }
                    if (ydx_constant) {
                        Mod poly_dy = (r == LEVEL_MAX - sc_y->levels[i]) ? Indices2Mod(sc_y->diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Mod();
                        Mod poly_xdy = mod.gb.Reduce(poly_x * poly_dy);
                        int1d xdy = poly_xdy ? Mod2Indices(poly_xdy, basis.at(deg_dxy)) : int1d{};
                        int1d dxy = lina::add(xdy, ydx_first_value);

                        if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                            rt += SSRet::CHANGE();
                            if (!printed_dx) {
                                printed_dx = true;
                                Logger::LogNullDiff(depth_, ring.name, deg_x, r, x, flag);
                            }
                            if (flag & SSFlag::log_proof) {
                                proof.clear();
                                proof = fmt::format("Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]`.", mod.name, deg_y, r, myio::Serialize(sc_y->basis[i]), poly_dy ? myio::Serialize(sc_y->diffs[i]) : "");
                            }
                            Logger::LogDiff(depth_, EnumReason::deduce_xy, mod.name, deg_xy, r, xy, dxy, proof, flag);
                            if (rt += SetModuleDiffGlobal(iMod, deg_xy, xy, dxy, r, true, flag))
                                return rt;
                            sc_dx = &ring.nodes_ss.GetRecentSc(deg_dx);
                            sc_y = &nodes_ss.GetRecentSc(deg_y);
                        }
                    }
                }
            }
        }
    }
    for (size_t iMap : ring.ind_maps) { /* Loop over map */
        auto map = (MapRing2Ring*)maps_[iMap].get();
        if (deg_dx.t > map->t_max)
            continue;
        auto fx = map->map(x, deg_x, *this);
        if (!fx.empty()) {
            bool fdx_constant = true;
            int1d fdx_first_value = NULL_DIFF;
            for (unsigned j = 0; j < j_max; ++j) { /* Loop over dx */
                dx.clear();
                for (int k : ut::two_exp(j))
                    dx = lina::add(dx, sc_dx->basis[(size_t)(first_dx + k)]);
                if (exclusions && ut::has(*exclusions, dx))
                    continue;
                if (j == 0) {
                    fdx_first_value = {};
                    continue;
                }

                int1d fdx = map->map(dx, deg_dx, *this);
                fdx = fdx.size() ? Residue(std::move(fdx), rings_[map->to.index].nodes_ss, deg_dx, r) : int1d{};
                if (fdx_first_value == NULL_DIFF)
                    fdx_first_value = std::move(fdx);
                else if (fdx_first_value != fdx) {
                    fdx_constant = false;
                    break;
                }
            }
            if (fdx_constant && IsNewDiff(rings_[map->to.index].nodes_ss, deg_x, fx, fdx_first_value, r)) {
                rt += SSRet::CHANGE();
                if (!printed_dx) {
                    printed_dx = true;
                    Logger::LogNullDiff(depth_, ring.name, deg_x, r, x, flag);
                }
                if (flag & SSFlag::log_proof) {
                    proof.clear();
                    proof = fmt::format("Consider the map `{}`.", map->name);
                }
                Logger::LogDiff(depth_, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, rings_[map->to.index].name), deg_x, r, fx, fdx_first_value, proof, flag);
                if (rt += SetRingDiffGlobal(map->to.index, deg_x, fx, fdx_first_value, r, true, flag))
                    return rt;
                sc_dx = &ring.nodes_ss.GetRecentSc(deg_dx);
            }
        }
    }

    return rt;
}

SSRet Category::SetModuleDiffLeibnizV2(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    SSRet rt;
    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& nodes_ss = mod.nodes_ss;
    const AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    if (!nodes_ss.front().has(deg_dx))
        return rt;
    auto& basis = mod.basis;
    auto& gb = mod.gb;
    const int t_max = mod.t_max;
    int depth = depth_;
    const Mod poly_x = Indices2Mod(x, basis.at(deg_x));
    const auto* sc_dx = &nodes_ss.GetRecentSc(deg_dx);
    auto [first_dx, count_dx] = CountPossDrTgt(nodes_ss, t_max, deg_dx, r);
    if (count_dx == 0 || count_dx > deduce_count_max_)
        return rt;
    int2d* exclusions = nullptr;
    if (auto pE = nodes_ss.exclusions.has(deg_x, x, r))
        exclusions = &pE->dxs;
    else
        exclusions = nullptr;

    unsigned j_max = 1 << count_dx;
    int1d dx;
    bool printed_dx = false;
    std::string proof;

    for (AdamsDeg deg_y : ring.degs_ss) {
        AdamsDeg deg_xy = deg_x + deg_y;
        AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
        if (deg_dxy.t > t_max)
            break;
        const auto* sc_y = &ring.nodes_ss.GetRecentSc(deg_y);
        for (size_t i = 0; i < sc_y->levels.size(); ++i) { /* Loop over y */
            const int r_y = LEVEL_MAX - sc_y->levels[i] - (sc_y->diffs[i] == NULL_DIFF ? 1 : 0);
            if (r_y < r)
                break;

            Poly poly_y = Indices2Poly(sc_y->basis[i], ring.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            if (poly_xy) {
                int1d xy = Mod2Indices(poly_xy, basis.at(deg_x + deg_y));

                bool ydx_constant = true;
                int1d ydx_first_value = NULL_DIFF;
                for (unsigned j = 0; j < j_max; ++j) { /* Loop over dx */
                    dx.clear();
                    for (int k : ut::two_exp(j))
                        dx = lina::add(dx, sc_dx->basis[(size_t)(first_dx + k)]);
                    if (exclusions && ut::has(*exclusions, dx))
                        continue;
                    if (j == 0) {
                        ydx_first_value = {};
                        continue;
                    }

                    Mod poly_dx = Indices2Mod(dx, basis.at(deg_dx));
                    Mod poly_ydx = mod.gb.Reduce(poly_y * poly_dx);
                    int1d ydx = poly_ydx ? Residue(Mod2Indices(poly_ydx, basis.at(deg_dxy)), nodes_ss, deg_dxy, r) : int1d{};
                    if (ydx_first_value == NULL_DIFF)
                        ydx_first_value = std::move(ydx);
                    else if (ydx_first_value != ydx) {
                        ydx_constant = false;
                        break;
                    }
                }
                if (ydx_constant) {
                    Poly poly_dy = (r == LEVEL_MAX - sc_y->levels[i]) ? Indices2Poly(sc_y->diffs[i], ring.basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
                    Mod poly_xdy = gb.Reduce(poly_dy * poly_x);
                    int1d xdy = poly_xdy ? Mod2Indices(poly_xdy, basis.at(deg_dxy)) : int1d{};
                    int1d dxy = lina::add(xdy, ydx_first_value);

                    if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                        if (!printed_dx) {
                            printed_dx = true;
                            Logger::LogNullDiff(depth, mod.name, deg_x, r, x, flag);
                        }
                        if (flag & SSFlag::log_proof) {
                            proof.clear();
                            proof = fmt::format("Apply the Leibniz rule with `{} {} d_{}[{}]=[{}]`.", ring.name, deg_y, r, myio::Serialize(sc_y->basis[i]), poly_dy ? myio::Serialize(sc_y->diffs[i]) : "");
                        }
                        Logger::LogDiff(depth, EnumReason::deduce_xy, mod.name, deg_xy, r, xy, dxy, proof, flag);
                        if (rt += SetModuleDiffGlobal(iMod, deg_xy, xy, dxy, r, true, flag))
                            return rt;
                        sc_dx = &nodes_ss.GetRecentSc(deg_dx);
                        sc_y = &ring.nodes_ss.GetRecentSc(deg_y);
                    }
                }
            }
        }
    }
    for (size_t iMap : mod.ind_maps) { /* Loop over map */
        auto& map = maps_[iMap];
        if (deg_dx.t > map->t_max)
            continue;
        auto fx = map->map(x, deg_x, *this);
        if (!fx.empty()) {
            auto& ss_to = GetNodesSS(map->to);
            auto& name_to = GetCwName(map->to);
            AdamsDeg deg_fx = deg_x + map->deg;
            AdamsDeg deg_fdx = deg_fx + AdamsDeg(r, r - 1);
            bool fdx_constant = true;
            int1d fdx_first_value = NULL_DIFF;
            for (unsigned j = 0; j < j_max; ++j) { /* Loop over dx */
                dx.clear();
                for (int k : ut::two_exp(j))
                    dx = lina::add(dx, sc_dx->basis[(size_t)(first_dx + k)]);
                if (exclusions && ut::has(*exclusions, dx))
                    continue;
                if (j == 0) {
                    fdx_first_value = {};
                    continue;
                }

                int1d fdx = map->map(dx, deg_dx, *this);
                fdx = fdx.size() ? Residue(std::move(fdx), ss_to, deg_fdx, r) : int1d{};
                if (fdx_first_value == NULL_DIFF)
                    fdx_first_value = std::move(fdx);
                else if (fdx_first_value != fdx) {
                    fdx_constant = false;
                    break;
                }
            }
            if (fdx_constant && IsNewDiff(ss_to, deg_fx, fx, fdx_first_value, r)) {
                if (!printed_dx) {
                    printed_dx = true;
                    Logger::LogNullDiff(depth, mod.name, deg_x, r, x, flag);
                }
                if (flag & SSFlag::log_proof) {
                    proof.clear();
                    proof = fmt::format("Consider the map `{}`.", map->name);
                }
                Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, name_to), deg_fx, r, fx, fdx_first_value, proof, flag);
                if (rt += SetCwDiffGlobal(map->to, deg_fx, fx, fdx_first_value, r, true, flag))
                    return rt;
                sc_dx = &nodes_ss.GetRecentSc(deg_dx);
            }
        }
    }

    return rt;
}

SSRet Category::SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    SSRet rt;
    auto& ring = rings_[iRing];

    const int r_original = r;
    r = NextRSrc(ring.nodes_ss, deg_x, r);
    if (r == -1) {
        if (depth_ == 0)
            throw RunTimeError(fmt::format("Contradiction. name={} deg_x={} x=[{}] r={}", ring.name, deg_x, myio::Serialize(x), r));
        else {
            rt.code = SSRet::FAIL_SS().code;
            if (flag & SSFlag::log_proof)
                rt.err_msg = fmt::format("However, `{} {} [{}]` is not in B_{}.", ring.name, deg_x, myio::Serialize(x), r_original);
            return rt;
        }
    }

    Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    {
        auto& nodes_ss = ring.nodes_ss;
        auto& basis = ring.basis;
        int t_max = ring.t_max;

        for (AdamsDeg deg_a : ring.degs_ss) {
            const auto& sc_a = nodes_ss.GetRecentSc(deg_a);
            AdamsDeg deg_ax = deg_x + deg_a;
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                if (sc_a.levels[i] >= LEVEL_MAX - r)
                    break;
                Poly poly_a = Indices2Poly(sc_a.basis[i], basis.at(deg_a));
                Poly poly_ax = ring.gb.Reduce(poly_x * poly_a);
                if (poly_ax) {
                    int1d ax = Poly2Indices(poly_ax, basis.at(deg_ax));  // TODO: consider moving allocations out of the loop
                    if (rt += SetImageSc(IndexRing(iRing), deg_ax, ax, NULL_DIFF, r, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} [{}]=d_{}[?]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[]` and get `{} {} [{}]=d_{}[?]`.\n{}", ring.name, deg_x, myio::Serialize(x), r, ring.name, deg_a, r, myio::Serialize(sc_a.basis[i]),
                                                     ring.name, deg_ax, myio::Serialize(ax), r, rt.err_msg);
                        return rt;
                    }
                    rt += SSRet::CHANGE();
                }
            }
        }
    }
    for (size_t iMod : ring.ind_mods) {
        auto& mod = modules_[iMod];
        auto& nodes_ss = mod.nodes_ss;
        int t_max = mod.t_max;

        for (AdamsDeg deg_a : mod.degs_ss) {
            const auto& sc_a = nodes_ss.GetRecentSc(deg_a);
            AdamsDeg deg_ax = deg_x + deg_a;
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                if (sc_a.levels[i] >= LEVEL_MAX - r)
                    break;
                Mod poly_a = Indices2Mod(sc_a.basis[i], mod.basis.at(deg_a));
                Mod poly_ax = mod.gb.Reduce(poly_x * poly_a);
                if (poly_ax) {
                    int1d ax = Mod2Indices(poly_ax, mod.basis.at(deg_ax));
                    if (rt += SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, r, flag)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} [{}]=d_{}[?]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[]` and get `{} {} [{}]=d_{}[?]`.\n{}", ring.name, deg_x, myio::Serialize(x), r, mod.name, deg_a, r, myio::Serialize(sc_a.basis[i]),
                                                     mod.name, deg_ax, myio::Serialize(ax), r, rt.err_msg);
                        return rt;
                    }
                    rt += SSRet::CHANGE();
                }
            }
        }
    }
    return rt;
}

SSRet Category::SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    SSRet rt;
    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    const int r_original = r;
    r = NextRSrc(nodes_ss, deg_x, r);
    if (r == -1) {
        if (depth_ == 0)
            throw RunTimeError(fmt::format("Contradiction. name={} deg_x={} x=[{}] r={}", mod.name, deg_x, myio::Serialize(x), r_original));
        else {
            rt.code = SSRet::FAIL_SS().code;
            if (flag & SSFlag::log_proof)
                rt.err_msg = fmt::format("However, `{} {} [{}]` is not in B_{}.", mod.name, deg_x, myio::Serialize(x), r_original);
            return rt;
        }
    }

    Mod poly_x = Indices2Mod(x, mod.basis.at(deg_x));
    for (AdamsDeg deg_a : ring.degs_ss) {
        const auto& sc_a = ring.nodes_ss.GetRecentSc(deg_a);
        AdamsDeg deg_ax = deg_x + deg_a;
        if (deg_ax.t > t_max)
            break;
        for (size_t i = 0; i < sc_a.levels.size(); ++i) {
            if (sc_a.levels[i] >= LEVEL_MAX - r)
                break;
            Poly poly_a = Indices2Poly(sc_a.basis[i], ring.basis.at(deg_a));
            Mod poly_ax = mod.gb.Reduce(poly_a * poly_x);
            if (poly_ax) {
                int1d ax = Mod2Indices(poly_ax, mod.basis.at(deg_ax));
                if (rt += SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, r, flag)) {
                    if (flag & SSFlag::log_proof)
                        rt.err_msg = fmt::format("Get `{} {} [{}]=d_{}[?]`. Apply the Leibniz rule with `{} {} d_{}[{}]=[]` and get `{} {} [{}]=d_{}[?]`.\n{}", mod.name, deg_x, myio::Serialize(x), r, ring.name, deg_a, r, myio::Serialize(sc_a.basis[i]), mod.name,
                                                 deg_ax, myio::Serialize(ax), r, rt.err_msg);
                    return rt;
                }
            }
        }
    }
    return rt;
}

SSRet Category::SetDiffLeibnizCofseq(CofSeq& cofseq, size_t iTri, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag)
{
    SSRet rt;
    size_t iTri_next = NextiTri(iTri);
    auto iCw = cofseq.indexCw[iTri];
    auto iCw_next = cofseq.indexCw[iTri_next];
    auto& ring = iCw.isRing() ? rings_[iCw.index] : rings_[modules_[iCw.index].iRing];  //// Assuming that the rings are the same
    using BasisVariant = std::variant<BasisMon*, BasisMMod*>;
    BasisVariant basis1 = iCw.isRing() ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw.index].basis);
    BasisVariant basis2 = iCw_next.isRing() ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw_next.index].basis);
    Poly poly_x, poly_dx;
    Mod mod_x, mod_dx;
    if (!x.empty()) {
        if (iCw.isRing())
            poly_x = Indices2Poly(x, std::get<0>(basis1)->at(deg_x));
        else
            mod_x = Indices2Mod(x, std::get<1>(basis1)->at(deg_x));
    }
    int stem_map = cofseq.degMap[iTri].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r + stem_map);
    if (!dx.empty()) {
        if (iCw_next.isRing())
            poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
        else
            mod_dx = Indices2Mod(dx, std::get<1>(basis2)->at(deg_dx));
    }
    int t_max = cofseq.t_max[iTri];
    int t_max2 = cofseq.t_max[iTri_next];

    for (AdamsDeg deg_a : ring.degs_ss) {
        const AdamsDeg deg_ax = deg_x + deg_a;
        AdamsDeg deg_adx = deg_ax + AdamsDeg(r, r + stem_map);
        const auto* sc_a = &ring.nodes_ss.GetRecentSc(deg_a);

        size_t first_PC = GetFirstIndexOnLevel(*sc_a, LEVEL_PERM);
        size_t last_PC = GetFirstIndexOnLevel(*sc_a, LEVEL_PERM + 1);
        for (size_t i = first_PC; i < last_PC; ++i) {
            Poly poly_a = Indices2Poly(sc_a->basis[i], ring.basis.at(deg_a));

            /* Compute ax */
            int1d ax;
            if (x.empty())
                ;
            else if (deg_ax.t > t_max)
                ax = NULL_DIFF;
            else {
                if (iCw.isRing()) {
                    Poly poly_ax = ring.gb.Reduce(poly_a * poly_x);
                    ax = poly_ax ? Poly2Indices(poly_ax, ring.basis.at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), ring.nodes_ss, deg_ax, LEVEL_PERM);
                }
                else {
                    Mod mod_ax = modules_[iCw.index].gb.Reduce(poly_a * mod_x);
                    ax = mod_ax ? Mod2Indices(mod_ax, std::get<1>(basis1)->at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), *cofseq.nodes_ss[iTri], deg_ax, LEVEL_PERM);
                }
            }

            /* Compute adx */
            int1d adx;
            if (dx.empty())
                ;
            else if (deg_adx.t > t_max2)
                adx = NULL_DIFF;
            else if (iCw_next.isRing()) {
                Poly poly_adx = ring.gb.Reduce(poly_a * poly_dx);
                if (poly_adx) {
                    adx = Poly2Indices(poly_adx, ring.basis.at(deg_adx));
                    adx = Residue(std::move(adx), ring.nodes_ss, deg_adx, LEVEL_PERM);
                }
            }
            else {
                Mod poly_adx = modules_[iCw_next.index].gb.Reduce(poly_a * mod_dx);
                if (poly_adx) {
                    adx = Mod2Indices(poly_adx, std::get<1>(basis2)->at(deg_adx));
                    adx = Residue(std::move(adx), *cofseq.nodes_ss[iTri_next], deg_adx, LEVEL_PERM);
                }
            }

            if ((!ax.empty() && ax != NULL_DIFF) || (!adx.empty() && adx != NULL_DIFF)) {
                if (rt += SetDiffScCofseq(cofseq, iTri, deg_ax, ax, adx, r, flag)) {
                    if (flag & SSFlag::log_proof)
                        rt.err_msg = fmt::format("Get `{}:{} {} d_{}[{}]=[{}]`. Apply the Leibniz rule with Permanent cycle `{} {} [{}]` and get `{}:{} {} d_{}[{}]=[{}]`.\n{}", cofseq.name, iTri, deg_x, r, myio::Serialize(x), myio::Serialize(dx), ring.name,
                                                 deg_a, myio::Serialize(sc_a->basis[i]), cofseq.name, iTri, deg_ax, r, myio::Serialize(ax), myio::Serialize(adx), rt.err_msg);
                    return rt;
                }
                sc_a = &ring.nodes_ss.GetRecentSc(deg_a);
            }
        }
    }
    if (rt += SyncToCofseq(flag))
        return rt;
    return rt;
}

SSRet Category::SetDiffLeibnizCofseq(IndexUniv iRing, AdamsDeg deg_a, const int1d& a, SSFlag flag)
{
    SSRet rt;
    Poly poly_a = Indices2Poly(a, rings_[iRing.index].basis.at(deg_a));
    Poly poly_x, poly_dx;
    Mod mod_x, mod_dx;
    int1d dx, ax, adx;
    for (auto& cofseq : cofseqs_) {
        IndexUniv iCw0 = cofseq.indexCw[0];
        size_t iRing1 = iCw0.isRing() ? iCw0.index : modules_[iCw0.index].iRing;
        if (iRing1 != iRing.index)
            continue;

        auto& ring = rings_[iRing.index];
        for (size_t iTri = 0; iTri < 3; ++iTri) {
            size_t iTri_prev = PreviTri(iTri);
            size_t iTri_next = NextiTri(iTri);
            int t_max = cofseq.t_max[iTri];
            int t_max2 = cofseq.t_max[iTri_next];
            auto iCw = cofseq.indexCw[iTri];
            auto iCw_prev = cofseq.indexCw[iTri_prev];
            auto iCw_next = cofseq.indexCw[iTri_next];
            int stem_map = cofseq.degMap[iTri].stem();
            int stem_map_prev = cofseq.degMap[iTri_prev].stem();
            using BasisVariant = std::variant<BasisMon*, BasisMMod*>;
            BasisVariant basis1 = iCw.isRing() ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw.index].basis);
            BasisVariant basis2 = iCw_next.isRing() ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw_next.index].basis);

            auto& nodes_cofseq = cofseq.nodes_cofseq[iTri];
            for (auto deg_x : nodes_cofseq.front().degs()) {
                const auto* sc = &nodes_cofseq.GetRecentSc(deg_x);
                const AdamsDeg deg_ax = deg_x + deg_a;
                for (size_t i = 0; i < sc->levels.size(); ++i) {
                    if (iCw.isRing())
                        poly_x = Indices2Poly(sc->basis[i], std::get<0>(basis1)->at(deg_x));
                    else
                        mod_x = Indices2Mod(sc->basis[i], std::get<1>(basis1)->at(deg_x));

                    /* Compute ax */
                    if (deg_ax.t > t_max)  // TODO: abstraction of multiplications
                        ax = NULL_DIFF;
                    else if (iCw.isRing()) {
                        Poly poly_ax = ring.gb.Reduce(poly_a * poly_x);
                        ax = poly_ax ? Poly2Indices(poly_ax, ring.basis.at(deg_ax)) : int1d{};
                        ax = Residue(std::move(ax), ring.nodes_ss, deg_ax, LEVEL_PERM);
                    }
                    else {
                        Mod mod_ax = modules_[iCw.index].gb.Reduce(poly_a * mod_x);
                        ax = mod_ax ? Mod2Indices(mod_ax, std::get<1>(basis1)->at(deg_ax)) : int1d{};
                        ax = Residue(std::move(ax), *cofseq.nodes_ss[iTri], deg_ax, LEVEL_PERM);
                    }

                    if (sc->levels[i] > LEVEL_MAX / 2) {
                        int r = LEVEL_MAX - sc->levels[i] - int(sc->diffs[i] == NULL_DIFF);
                        auto deg_dx = deg_x + AdamsDeg(r, stem_map + r);
                        auto deg_adx = deg_a + deg_dx;
                        dx = sc->diffs[i] == NULL_DIFF ? int1d{} : sc->diffs[i];

                        /* Compute adx */
                        if (dx.empty())
                            adx.clear();
                        else if (deg_adx.t > t_max2)
                            adx = NULL_DIFF;
                        else if (iCw_next.isRing()) {
                            poly_dx = dx.empty() ? Poly{} : Indices2Poly(dx, ring.basis.at(deg_dx));
                            Poly poly_adx = ring.gb.Reduce(poly_a * poly_dx);
                            if (poly_adx) {
                                adx = Poly2Indices(poly_adx, ring.basis.at(deg_adx));
                                adx = Residue(std::move(adx), ring.nodes_ss, deg_adx, LEVEL_PERM);
                            }
                            else
                                adx.clear();
                        }
                        else {
                            mod_dx = dx.empty() ? Mod{} : Indices2Mod(dx, std::get<1>(basis2)->at(deg_dx));
                            Mod poly_adx = modules_[iCw_next.index].gb.Reduce(poly_a * mod_dx);
                            if (poly_adx) {
                                adx = Mod2Indices(poly_adx, std::get<1>(basis2)->at(deg_adx));
                                adx = Residue(std::move(adx), *cofseq.nodes_ss[iTri_next], deg_adx, LEVEL_PERM);
                            }
                            else
                                adx.clear();
                        }

                        if ((!ax.empty() && ax != NULL_DIFF) || (!adx.empty() && adx != NULL_DIFF)) {
                            if (rt += SetDiffScCofseq(cofseq, iTri, deg_ax, ax, adx, r, flag)) {
                                if (flag & SSFlag::log_proof)
                                    rt.err_msg = fmt::format("Get Permanent cycle `{} {} [{}]`. Apply the Leibniz rule with `{}:{} {} d_{}[{}]=[{}]` and get `{}:{} {} d_{}[{}]=[{}]`.\n{}", ring.name, deg_a, myio::Serialize(a), cofseq.name, iTri, deg_x, r,
                                                             myio::Serialize(sc->basis[i]), myio::Serialize(dx), cofseq.name, iTri, deg_ax, r, myio::Serialize(ax), myio::Serialize(adx), rt.err_msg);
                                return rt;
                            }
                        }
                    }
                    else if (sc->diffs[i] == NULL_DIFF && (!ax.empty() && ax != NULL_DIFF)) {
                        int r = sc->levels[i] + 1;
                        AdamsDeg d_src = deg_ax - AdamsDeg(r, r + stem_map_prev);
                        if (rt += SetDiffScCofseq(cofseq, iTri_prev, d_src, {}, ax, r, flag))
                            return rt;
                    }
                }
            }
        }
    }
    return rt;
}

Staircase GetSc4Display(const SSNodes& nodes_this, const SSNodes& nodes_src, AdamsDeg deg, int stem_map)
{
    Staircase result;
    if (!nodes_this.front().has(deg))
        return result;
    auto& sc = nodes_this.GetRecentSc(deg);
    for (size_t i = 0; i < sc.levels.size();) {
        if (sc.levels[i] < LEVEL_PERM && sc.diffs[i] != NULL_DIFF) {
            int r = sc.levels[i];
            auto deg_src = deg - AdamsDeg(r, stem_map + r);
            auto& sc_src = nodes_src.GetRecentSc(deg_src);
            int level = sc.levels[i];
            size_t j = 0;
            for (; j < sc_src.levels.size(); ++j) {
                if (sc_src.levels[j] == LEVEL_MAX - r && sc_src.diffs[j] != NULL_DIFF)
                    break;
            }
            while (i < sc.levels.size() && sc.levels[i] == level && sc.diffs[i] != NULL_DIFF) {
                result.basis.push_back(sc_src.diffs[j]);
                result.diffs.push_back(sc_src.basis[j]);
                result.levels.push_back(LEVEL_MAX - sc_src.levels[j]);
                ++i, ++j;
            }
        }
        else {
            result.basis.push_back(sc.basis[i]);
            result.diffs.push_back(sc.diffs[i]);
            result.levels.push_back(sc.levels[i]);
            ++i;
        }
    }
    return result;
}