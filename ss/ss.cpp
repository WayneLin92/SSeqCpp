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

bool Diagram::IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg, int r_max)
{
    r_max = std::min(r_max, deg.s);
    for (int r1 = LEVEL_MIN; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (ut::has(nodes_ss.front(), d_src) && GetMaxLevelWithND(ut::GetRecentValue(nodes_ss, d_src)) >= LEVEL_MAX - r1)
            return true;
    }
    return false;
}

size_t Diagram::GetFirstIndexOfFixedLevels(const Staircases1d& nodes_ss, AdamsDeg deg, int level_min)
{
    const auto& sc = ut::GetRecentValue(nodes_ss, deg);
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

int Diagram::GetFirstFixedLevelForPlot(const Staircases1d& nodes_ss, AdamsDeg deg)
{
    const auto& sc = ut::GetRecentValue(nodes_ss, deg);
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

std::pair<int, int> Diagram::CountPossDrTgt(const Staircases1d& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const
{
    std::pair<int, int> result;
    if (ut::has(nodes_ss.front(), deg_tgt)) {
        const auto& sc_tgt = ut::GetRecentValue(nodes_ss, deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_tgt, LEVEL_MAX - r) - result.first;
    }
    else if (deg_tgt.t > t_max)
        result = {-1, 100000}; /* Infinitely many possibilities */
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Diagram::CountPossDrSrc(const Staircases1d& nodes_ss, const AdamsDeg& deg_src, int r) const
{
    std::pair<int, int> result;
    if (ut::has(nodes_ss.front(), deg_src)) {
        const auto& sc_src = ut::GetRecentValue(nodes_ss, deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_src, LEVEL_MAX - r) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Diagram::NextRTgt(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r) const
{
    for (int r1 = r; r1 <= R_PERM; ++r1) {
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 - 1};
        if (d_tgt.t > t_max)
            return r1;
        if (r1 >= 20 && AboveJ(d_tgt) && BelowCokerJ(deg)) /* Image of J */
            return R_PERM;
        if (AboveS0Vanishing(d_tgt) && !ut::has(nodes_ss.front(), d_tgt))
            return R_PERM;
        if (CountPossDrTgt(nodes_ss, t_max, d_tgt, r1).second > 0)
            return r1;
    }
    return R_PERM;
}

int Diagram::NextRSrc(const Staircases1d& nodes_ss, AdamsDeg deg, int r) const
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

void Diagram::CacheNullDiffs(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, SSFlag flag, NullDiff1d& nds) const
{
    nds.clear();
    const auto& sc = ut::GetRecentValue(nodes_ss, deg);
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
        if (flag & SSFlag::all_x) {
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

/* Return d_r(x) */
int1d Diagram::GetDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, int r) const
{
    int1d result;
    if (x.empty())
        return result;
    const auto& sc = ut::GetRecentValue(nodes_ss, deg_x);
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

int1d Diagram::GetLevelAndDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, int1d x, int& level) const
{
    int1d dx;
    level = -1;

#ifdef MYDEBUG
    MyException::Assert(!x.empty(), "GetLevelAndDiff() para: !x.empty()");
#endif

    const auto& sc = ut::GetRecentValue(nodes_ss, deg_x);
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
    MyException::Assert(x.empty(), "GetLevelAndDiff() end: x.empty()");
#endif
    return dx;
}

bool Diagram::IsNewDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    if (x.empty())
        return !dx.empty() && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dx), dx, r);
    int1d dx1 = GetDiff(nodes_ss, deg_x, x, r);  //// TODO: optimize the allocation
    if (dx1 == NULL_DIFF)
        return true;
    int1d diff = lina::add(dx, dx1);
    return !diff.empty() && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dx), diff, r);
}

/* Return the minimal length of the crossing differentials */
int Diagram::GetCrossR(const Staircases1d& nodes_ss, AdamsDeg deg, int t_max) const
{
    int result = R_PERM;
    for (int r = 1; r <= R_PERM; ++r) {
        auto deg_x = deg + AdamsDeg{r, r};
        if (r + LEVEL_MIN > result)
            return result;
        if (deg_x.t > t_max)
            return std::min(result, r + LEVEL_MIN);
        if (!ut::has(nodes_ss.front(), deg_x)) {
            if (AboveS0Vanishing(deg_x))
                return result;
            continue;
        }
        auto& sc = ut::GetRecentValue(nodes_ss, deg_x);
        if (!sc.levels.empty() && sc.levels.back() > LEVEL_MAX / 2) {
            int r1 = LEVEL_MAX - sc.levels.back();
            if (r + r1 < result)
                result = r + r1;
        }
    }
    return result;
}

void Diagram::SetDiffSc(IndexCw iCw, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, SSFlag flag)
{
    /*if (deg_x == AdamsDeg(10, 135) && r == 2 && (x_ == int1d{0} || x_ == int1d{1} || x_ == int1d{0, 1})) {
        auto& name = !(iCw & FLAG_MOD) ? rings_[iCw].name : modules_[iCw ^ FLAG_MOD].name;
        if (name == "Ctheta5")
            std::cout << "debug\n";
    }*/

    /* NULL_DIFF checking */
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    if (x_ == NULL_DIFF) {
        SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r, flag);
        return;
    }

    /* Triangularize x */
    if (x_.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return;
    }
    auto& nodes_ss = GetSS(iCw);
    const auto& sc = ut::GetRecentValue(nodes_ss, deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    int1d x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, x_);
    if (x.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            SetImageSc(iCw, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return;
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
    if ((flag & SSFlag::cofseq) && r == R_PERM) {
        const auto& ind_cofs = GetIndexCof(iCw);
        for (auto& ind_cof : ind_cofs) {
            auto& cofseq = cofseqs_[ind_cof.iCof];
            auto iCs = (size_t)ind_cof.iCs;
            if (!ut::has(cofseq.nodes_cofseq[iCs].front(), deg_x))
                cofseq.nodes_cofseq[iCs].front()[deg_x] = {};
            SetDiffScCofseq(cofseq, ind_cof.iCs, deg_x, x, NULL_DIFF, 0, flag);
        }
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_x + AdamsDeg{level_image_new, level_image_new - 1};
            SetImageSc(iCw, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, r_image - 1};
            SetDiffSc(iCw, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }

    /* Set the image */
    if (!dx.empty() && dx != NULL_DIFF)
        SetImageSc(iCw, deg_dx, dx, x, r, flag);
}

void Diagram::SetImageSc(IndexCw iCw, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, SSFlag flag)
{
    const auto& name = GetCwName(iCw);
    auto& nodes_ss = GetSS(iCw);
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, r - 1};

    // if (nodes_ss.size() == 2 && deg_dx == AdamsDeg(16, 162 + 16) && r == 2 && dx_ == int1d{6} && name == "Joker") ////
    //     throw InteruptAndSaveException(0, "bug");

    /* If dx is in Im(d_{r-1}) then x is in Ker(d_r) */
    const auto& sc = ut::GetRecentValue(nodes_ss, deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    int1d dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, dx_);
    if (dx.empty()) {
        if (x != NULL_DIFF)
            SetDiffSc(iCw, deg_x, x, NULL_DIFF, r + 1, flag);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == NULL_DIFF) {
        /* If the source is uncertain, check if it can be hit and then insert it to the end of level r. */
        size_t first_rp1 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp1, std::move(dx));
        if (!dx.empty()) {
            if (!IsPossTgt(nodes_ss, deg_dx, r)) {
                if (nodes_ss.size() == 2)
                    fmt::print("Contradiction\n");
                Logger::LogSSException(depth_, name, deg_dx, dx, r, 0x7U, deg_leibniz_, a_leibniz_);
                throw SSException(0x7U, "No source for the image.");
            }
            UpdateStaircase(nodes_ss, deg_dx, sc, first_rp1, dx, x, r, image_new, level_image_new);
        }
    }
    else {
        /* Otherwise insert it to the beginning of level r */
        UpdateStaircase(nodes_ss, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new);
    }

    /* ss to cofseq */
    if (flag & SSFlag::cofseq) {
        const auto& ind_cofs = GetIndexCof(iCw);
        for (auto& ind_cof : ind_cofs) {
            auto& cofseq = cofseqs_[ind_cof.iCof];
            auto iCs = (size_t)ind_cof.iCs;
            if (ut::has(cofseq.nodes_cofseq[iCs].front(), deg_dx))
                ReSetScCofseq(cofseq, ind_cof.iCs, deg_dx, flag);
        }
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg_dx + AdamsDeg{level_image_new, level_image_new - 1};
            SetImageSc(iCw, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, r_image - 1};
            SetDiffSc(iCw, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }
}

void Diagram::SetDiffScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x_, const int1d& dx_, int r, SSFlag flag)
{
    /*if (cofseq.name == "S0__Ceta__S0" && iCs == 1 && deg_x == AdamsDeg(29, 176)) {
        auto& sc = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg_x);
        fmt::print("adding d_{}{}={}\n", r, x_, dx_);
        fmt::print("sc at {} =\n{}\n", deg_x, sc);
        AdamsDeg deg_next = deg_x + AdamsDeg{0, cofseq.degMap[iCs].stem() + 0};
        fmt::print("sc_next at {} =\n{}\n", deg_next, ut::GetRecentValue(cofseq.nodes_cofseq[(iCs + 1) % 3], deg_next));
        std::cout << "debug\n";
    }*/
    const size_t iCs_prev = (iCs + 2) % 3;
    const size_t iCs_next = (iCs + 1) % 3;
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_ss = *cofseq.nodes_ss[iCs];
    int stem_map = cofseq.degMap[iCs].stem();
    int stem_map_prev = cofseq.degMap[iCs_prev].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, stem_map + r};

    /* NULL_DIFF checking */
    if (x_ == NULL_DIFF) {
        SetImageScCofseq(cofseq, iCs_next, deg_dx, dx_, NULL_DIFF, r, flag);
        return;
    }

    /* Triangularize x */
    if (x_.empty()) {
        if (dx_ != NULL_DIFF && !dx_.empty())
            SetImageScCofseq(cofseq, iCs_next, deg_dx, dx_, NULL_DIFF, r - 1, flag);
        return;
    }
    int1d x = Residue(x_, nodes_ss, deg_x, LEVEL_PERM);
    int1d dx = dx_;
    if (dx_ != NULL_DIFF && !dx_.empty())
        dx = Residue(std::move(dx), *cofseq.nodes_ss[iCs_next], deg_dx, LEVEL_PERM);
    const auto& sc = ut::GetRecentValue(nodes_cofseq, deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, std::move(x));
    if (x.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            SetImageScCofseq(cofseq, iCs_next, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return;
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
            SetImageScCofseq(cofseq, iCs_next, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, stem_map_prev + r_image};
            SetDiffScCofseq(cofseq, iCs_prev, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }

    /* Set image */
    if (dx != NULL_DIFF && !dx.empty())
        SetImageScCofseq(cofseq, iCs_next, deg_dx, dx, x, r, flag);

    /*if (cofseq.name == "S0__Ceta__S0" && iCs == 1 && deg_x == AdamsDeg(29, 176)) {
        auto& sc = ut::GetRecentValue(cofseq.nodes_cofseq[iCs], deg_x);
        fmt::print("sc at {} =\n{}\n", deg_x, sc);
        AdamsDeg deg_next = deg_x + AdamsDeg{0, cofseq.degMap[iCs].stem() + 0};
        fmt::print("sc_next at {} =\n{}\n", deg_next, ut::GetRecentValue(cofseq.nodes_cofseq[(iCs + 1) % 3], deg_next));
        std::cout << "debug\n";
    }*/
}

void Diagram::SetImageScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_dx, const int1d& dx_, const int1d& x_, int r, SSFlag flag)
{
    const size_t iCs_prev = (iCs + 2) % 3;
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_ss = *cofseq.nodes_ss[iCs];
    int stem_map = cofseq.degMap[iCs].stem();
    int stem_map_prev = cofseq.degMap[iCs_prev].stem();
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, stem_map_prev + r};

    int1d dx = Residue(dx_, nodes_ss, deg_dx, LEVEL_PERM);
    int1d x = x_;
    if (x_ != NULL_DIFF && !x_.empty())
        x = Residue(std::move(x), *cofseq.nodes_ss[iCs_prev], deg_x, LEVEL_PERM);
    const auto& sc = ut::GetRecentValue(nodes_cofseq, deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, std::move(dx));
    if (dx.empty()) {
        if (x != NULL_DIFF)
            SetDiffScCofseq(cofseq, iCs_prev, deg_x, x, NULL_DIFF, r + 1, flag);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == NULL_DIFF) {
        /* If the source is uncertain, insert it to the end of level r. */
        size_t first_rp1 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp1, std::move(dx));
        if (!dx.empty()) {
            if (!IsPossTgtCofseq(cofseq, iCs, deg_dx, r)) {
                r = -1;
                first_rp1 = 0;
            }
            UpdateStaircase(nodes_cofseq, deg_dx, sc, first_rp1, dx, x, r, image_new, level_image_new);

            /* cofseq to ss */
            if (r == -1) {
                pop_front(nodes_cofseq.back().at(deg_dx));
                Logger::LogDiffInv(int(nodes_cofseq.size() - 2), EnumReason::cofseq_b, cofseq.nameCw[iCs], deg_dx - AdamsDeg(R_PERM + 1, R_PERM), deg_dx, {}, dx, R_PERM + 1);
                a_leibniz_ = nullptr;
                SetImageSc(cofseq.indexCw[iCs], deg_dx, dx, NULL_DIFF, R_PERM, flag);
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
            SetImageScCofseq(cofseq, (iCs + 1) % 3, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, stem_map_prev + r_image};
            SetDiffScCofseq(cofseq, iCs_prev, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }
}

int Diagram::SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag)  ////
{
    int count = 0;
    auto& ring = rings_[iRing];
    Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x)), poly_zero;
    Poly poly_drx = !dx.empty() ? Indices2Poly(dx, ring.basis.at(deg_x + AdamsDeg(r, r - 1))) : poly_zero;
    Poly poly_a, poly_ax, poly_da, poly_dax, poly_tmp1, poly_tmp2;
    int1d ax, dax;
    Mod mod_y, mod_xy, mod_dy, mod_dxy, mod_tmp1, mod_tmp2;
    int1d xy, dxy;

    int r_zero = dx.empty() ? r : r - 1;
    {
        auto& nodes_ss = ring.nodes_ss;
        auto& basis = ring.basis;
        int t_max = ring.t_max;
        for (auto& [deg_a, _] : nodes_ss.front()) {
            const AdamsDeg deg_ax = deg_x + deg_a;
            const auto& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                if (sc_a.levels[i] < LEVEL_MAX / 2 && sc_a.diffs[i] == NULL_DIFF && sc_a.levels[i] >= r_min && sc_a.levels[i] <= r_zero) { /* ax=d_R[?] */
                    const int R = sc_a.levels[i];
                    Indices2AlgP(sc_a.basis[i], basis.at(deg_a), poly_a);
                    mulP(poly_a, poly_x, poly_ax); /* poly_ax = poly_a * poly_x */
                    ring.gb.ReduceP(poly_ax, poly_tmp1, poly_tmp2);
                    if (poly_ax) {
                        Alg2IndicesP(poly_ax, basis.at(deg_ax), ax);
                        deg_leibniz_ = deg_a;
                        a_leibniz_ = &sc_a.basis[i];
                        SetImageSc(IndexRing(iRing), deg_ax, ax, NULL_DIFF, R, flag);
                        a_leibniz_ = nullptr;
                        ++count;
                    }
                }
                else if (sc_a.levels[i] > LEVEL_MAX / 2) {
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
                        deg_leibniz_ = deg_a;
                        a_leibniz_ = &sc_a.basis[i];
                        SetDiffSc(IndexRing(iRing), deg_ax, ax, dax, R, flag);
                        a_leibniz_ = nullptr;
                        ++count;
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
        for (auto& [deg_y, _] : nodes_ss.front()) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const auto& sc_y = ut::GetRecentValue(nodes_ss, deg_y);
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
                        deg_leibniz_ = deg_y;
                        a_leibniz_ = &sc_y.basis[i];
                        SetImageSc(IndexMod(iMod), deg_xy, xy, NULL_DIFF, R, flag);
                        a_leibniz_ = nullptr;
                        ++count;
                    }
                }
                else if (sc_y.levels[i] > LEVEL_MAX / 2) {
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
                        deg_leibniz_ = deg_y;
                        a_leibniz_ = &sc_y.basis[i];
                        SetDiffSc(IndexMod(iMod), deg_xy, xy, dxy, R, flag);
                        a_leibniz_ = nullptr;
                        ++count;
                    }
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag)  ////
{
    int count = 0;

    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& basis = mod.basis;
    auto& gb = mod.gb;
    const int t_max = mod.t_max;
    Mod poly_x = Indices2Mod(x, basis.at(deg_x)), poly_zero;
    Mod poly_drx = !dx.empty() ? Indices2Mod(dx, basis.at(deg_x + AdamsDeg(r, r - 1))) : poly_zero;
    int r_zero = dx.empty() ? r : r - 1;
    Poly poly_a, poly_da, poly_tmp1;
    Mod mod_ax, mod_dax, mod_tmp1, mod_tmp2;
    int1d ax, dax;

    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const auto& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);
        AdamsDeg deg_ax = deg_x + deg_a;
        if (deg_ax.t > t_max)
            break;
        for (size_t i = 0; i < sc_a.levels.size(); ++i) {
            if (sc_a.levels[i] < LEVEL_MAX / 2 && sc_a.diffs[i] == NULL_DIFF && sc_a.levels[i] >= r_min && sc_a.levels[i] <= r_zero) { /* ax=d_R[?] */
                const int R = sc_a.levels[i];
                Indices2AlgP(sc_a.basis[i], ring.basis.at(deg_a), poly_a);
                mulP(poly_a, poly_x, mod_ax); /* poly_ax = poly_a * poly_x */
                mod.gb.ReduceP(mod_ax, poly_tmp1, mod_tmp1, mod_tmp2);
                if (mod_ax) {
                    Alg2IndicesP(mod_ax, basis.at(deg_ax), ax);
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, R, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
            else if (sc_a.levels[i] > LEVEL_MAX / 2) {
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
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetDiffSc(IndexMod(iMod), deg_ax, ax, dax, R, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetRingDiffLeibnizV2(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    const AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    if (!ut::has(ring.nodes_ss.front(), deg_dx))
        return count;
    int depth = int(ring.nodes_ss.size() - 2);
    const Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    const auto& sc_dx = ut::GetRecentValue(ring.nodes_ss, deg_dx);
    auto [first_dx, count_dx] = CountPossDrTgt(ring.nodes_ss, ring.t_max, deg_dx, r);
    if (count_dx == 0 || count_dx > deduce_count_max_)
        return 0;
    unsigned j_max = 1 << count_dx;
    int1d dx;
    bool printed_dx = false;

    if (auto& nodes_ss = ring.nodes_ss; ut::has(nodes_ss.front(), deg_dx)) {
        const auto& basis = ring.basis;
        const int t_max = ring.t_max;
        for (auto& [deg_y, _] : nodes_ss.front()) {
            const AdamsDeg deg_xy = deg_x + deg_y;
            const AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
            if (deg_dxy.t > t_max)
                break;
            const auto& sc_y = ut::GetRecentValue(nodes_ss, deg_y);
            for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
                const int r_y = LEVEL_MAX - sc_y.levels[i] - (sc_y.diffs[i] == NULL_DIFF ? 1 : 0);
                if (r_y < r)
                    break;

                Poly poly_y = Indices2Poly(sc_y.basis[i], basis.at(deg_y));
                Poly poly_xy = ring.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = Poly2Indices(poly_xy, basis.at(deg_xy));

                    bool ydx_always_zero = true;
                    for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : ut::two_exp(j))
                            dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                        Poly poly_dx = Indices2Poly(dx, basis.at(deg_dx));
                        Poly poly_ydx = ring.gb.Reduce(poly_dx * poly_y);
                        if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Poly2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                            ydx_always_zero = false;
                            break;
                        }
                    }
                    if (ydx_always_zero) {
                        Poly poly_dy = (r == LEVEL_MAX - sc_y.levels[i]) ? Indices2Poly(sc_y.diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
                        Poly poly_dxy = ring.gb.Reduce(poly_x * poly_dy);
                        int1d dxy = poly_dxy ? Poly2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};

                        if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                            if (!printed_dx) {
                                printed_dx = true;
                                Logger::LogNullDiff(depth, ring.name, deg_x, x, r);
                            }
                            Logger::LogDiff(depth, EnumReason::deduce_xy, ring.name, deg_xy, xy, dxy, r);
                            count += SetRingDiffGlobal(iRing, deg_xy, xy, dxy, r, true, flag);
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
        for (auto& [deg_y, _] : nodes_ss.front()) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
            if (deg_dxy.t > t_max)
                break;
            const auto& sc_y = ut::GetRecentValue(nodes_ss, deg_y);
            for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
                const int r_y = LEVEL_MAX - sc_y.levels[i] - (sc_y.diffs[i] == NULL_DIFF ? 1 : 0);
                if (r_y < r)
                    break;

                Mod poly_y = Indices2Mod(sc_y.basis[i], basis.at(deg_y));
                Mod poly_xy = mod.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    bool ydx_always_zero = true;
                    for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : ut::two_exp(j))
                            dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                        Poly poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
                        Mod poly_ydx = mod.gb.Reduce(poly_dx * poly_y);
                        if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Mod2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                            ydx_always_zero = false;
                            break;
                        }
                    }
                    if (ydx_always_zero) {
                        Mod poly_dy = (r == LEVEL_MAX - sc_y.levels[i]) ? Indices2Mod(sc_y.diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Mod();
                        Mod poly_dxy = mod.gb.Reduce(poly_x * poly_dy);
                        int1d dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};

                        if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                            if (!printed_dx) {
                                printed_dx = true;
                                Logger::LogNullDiff(depth, ring.name, deg_x, x, r);
                            }
                            Logger::LogDiff(depth, EnumReason::deduce_xy, mod.name, deg_xy, xy, dxy, r);
                            count += SetModuleDiffGlobal(iMod, deg_xy, xy, dxy, r, true, flag);
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
            bool fdx_always_zero = true;
            for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                dx.clear();
                for (int k : ut::two_exp(j))
                    dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                int1d fdx = map->map(dx, deg_dx, *this);
                if (!fdx.empty() && !IsZeroOnLevel(ut::GetRecentValue(rings_[map->to.index].nodes_ss, deg_dx), fdx, r)) {
                    fdx_always_zero = false;
                    break;
                }
            }
            if (fdx_always_zero && IsNewDiff(rings_[map->to.index].nodes_ss, deg_x, fx, {}, r)) {
                if (!printed_dx) {
                    printed_dx = true;
                    Logger::LogNullDiff(depth, ring.name, deg_x, x, r);
                }
                Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, rings_[map->to.index].name), deg_x, fx, {}, r);
                count += SetRingDiffGlobal(map->to.index, deg_x, fx, {}, r, true, flag);
            }
        }
    }

    return count;
}

int Diagram::SetModuleDiffLeibnizV2(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    int count = 0;

    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& nodes_ss = mod.nodes_ss;
    const AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    if (!ut::has(nodes_ss.front(), deg_dx))
        return count;
    auto& basis = mod.basis;
    auto& gb = mod.gb;
    const int t_max = mod.t_max;
    int depth = depth_;
    const Mod poly_x = Indices2Mod(x, basis.at(deg_x));
    const auto& sc_dx = ut::GetRecentValue(nodes_ss, deg_dx);
    auto [first_dx, count_dx] = CountPossDrTgt(nodes_ss, t_max, deg_dx, r);
    if (count_dx == 0 || count_dx > deduce_count_max_)
        return 0;
    unsigned j_max = 1 << count_dx;
    int1d dx;
    bool printed_dx = false;

    for (auto& [deg_y, _] : ring.nodes_ss.front()) {
        AdamsDeg deg_xy = deg_x + deg_y;
        AdamsDeg deg_dxy = deg_xy + AdamsDeg(r, r - 1);
        if (deg_dxy.t > t_max)
            break;
        const auto& sc_y = ut::GetRecentValue(ring.nodes_ss, deg_y);
        for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
            const int r_y = LEVEL_MAX - sc_y.levels[i] - (sc_y.diffs[i] == NULL_DIFF ? 1 : 0);
            if (r_y < r)
                break;

            Poly poly_y = Indices2Poly(sc_y.basis[i], ring.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            if (poly_xy) {
                int1d xy = Mod2Indices(poly_xy, basis.at(deg_x + deg_y));

                bool ydx_always_zero = true;
                for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                    dx.clear();
                    for (int k : ut::two_exp(j))
                        dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                    Mod poly_dx = Indices2Mod(dx, basis.at(deg_dx));
                    Mod poly_ydx = mod.gb.Reduce(poly_y * poly_dx);
                    if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Mod2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                        ydx_always_zero = false;
                        break;
                    }
                }
                if (ydx_always_zero) {
                    Poly poly_dy = (r == LEVEL_MAX - sc_y.levels[i]) ? Indices2Poly(sc_y.diffs[i], ring.basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
                    Mod poly_dxy = gb.Reduce(poly_dy * poly_x);
                    int1d dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};

                    if (IsNewDiff(nodes_ss, deg_xy, xy, dxy, r)) {
                        if (!printed_dx) {
                            printed_dx = true;
                            Logger::LogNullDiff(depth, mod.name, deg_x, x, r);
                        }
                        Logger::LogDiff(depth, EnumReason::deduce_xy, mod.name, deg_xy, xy, dxy, r);
                        count += SetModuleDiffGlobal(iMod, deg_xy, xy, dxy, r, true, flag);
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
            auto& ss_to = GetSS(map->to);
            auto& name_to = GetCwName(map->to);
            AdamsDeg deg_fx = deg_x + map->deg;
            AdamsDeg deg_fdx = deg_fx + AdamsDeg(r, r - 1);
            bool fdx_always_zero = true;
            for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                dx.clear();
                for (int k : ut::two_exp(j))
                    dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                int1d fdx = map->map(dx, deg_dx, *this);
                if (!fdx.empty() && !IsZeroOnLevel(ut::GetRecentValue(ss_to, deg_fdx), fdx, r)) {
                    fdx_always_zero = false;
                    break;
                }
            }
            if (fdx_always_zero && IsNewDiff(ss_to, deg_fx, fx, {}, r)) {
                if (!printed_dx) {
                    printed_dx = true;
                    Logger::LogNullDiff(depth, mod.name, deg_x, x, r);
                }
                Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, name_to), deg_fx, fx, {}, r);
                count += SetCwDiffGlobal(map->to, deg_fx, fx, {}, r, true, flag);
            }
        }
    }

    return count;
}

int Diagram::SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];

    const int r_original = r;
    r = NextRSrc(ring.nodes_ss, deg_x, r);
    if (r == -1) {
        a_leibniz_ = nullptr;
        Logger::LogSSException(int(ring.nodes_ss.size() - 2), ring.name, deg_x, x, r_original, 0x51274f1dU, deg_leibniz_, a_leibniz_);
        throw SSException(0x51274f1dU, "No source for the image.");
    }

    Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    {
        auto& nodes_ss = ring.nodes_ss;
        auto& basis = ring.basis;
        int t_max = ring.t_max;

        for (auto& [deg_a, _] : nodes_ss.front()) {
            const auto& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
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
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetImageSc(IndexRing(iRing), deg_ax, ax, NULL_DIFF, r, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    for (size_t iMod : ring.ind_mods) {
        auto& mod = modules_[iMod];
        auto& nodes_ss = mod.nodes_ss;
        int t_max = mod.t_max;

        for (auto& [deg_a, _] : nodes_ss.front()) {
            const auto& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
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
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, r, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag)
{
    int count = 0;
    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    const int r_original = r;
    r = NextRSrc(nodes_ss, deg_x, r);
    if (r == -1) {
        a_leibniz_ = nullptr;
        Logger::LogSSException(depth_, mod.name, deg_x, x, r_original, 0xda298807U, deg_leibniz_, a_leibniz_);
        throw SSException(0xda298807U, "No source for the image.");
    }

    Mod poly_x = Indices2Mod(x, mod.basis.at(deg_x));
    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const auto& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);
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
                deg_leibniz_ = deg_a;
                a_leibniz_ = &sc_a.basis[i];
                SetImageSc(IndexMod(iMod), deg_ax, ax, NULL_DIFF, r, flag);
                a_leibniz_ = nullptr;
                ++count;
            }
        }
    }
    return count;
}

int Diagram::SetDiffLeibnizCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag)
{
    int count = 0;
    size_t iCs_next = (iCs + 1) % 3;
    auto iCw = cofseq.indexCw[iCs];
    auto iCw_next = cofseq.indexCw[iCs_next];
    auto& ring = iCw.isRing ? rings_[iCw.index] : rings_[modules_[iCw.index].iRing]; //// Assuming that the rings are the same
    using BasisVariant = std::variant<std::map<AdamsDeg, Mon1d>*, std::map<AdamsDeg, MMod1d>*>;
    BasisVariant basis1 = iCw.isRing ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw.index].basis);
    BasisVariant basis2 = iCw_next.isRing ? BasisVariant(&ring.basis) : BasisVariant(&modules_[iCw_next.index].basis);
    Poly poly_x, poly_dx;
    Mod mod_x, mod_dx;
    if (!x.empty()) {
        if (iCw.isRing)
            poly_x = Indices2Poly(x, std::get<0>(basis1)->at(deg_x));
        else
            mod_x = Indices2Mod(x, std::get<1>(basis1)->at(deg_x));
    }
    int stem_map = cofseq.degMap[iCs].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r + stem_map);
    if (!dx.empty()) {
        if (iCw_next.isRing)
            poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
        else
            mod_dx = Indices2Mod(dx, std::get<1>(basis2)->at(deg_dx));
    }
    int t_max = cofseq.t_max[iCs];
    int t_max2 = cofseq.t_max[iCs_next];
    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const AdamsDeg deg_ax = deg_x + deg_a;
        AdamsDeg deg_adx = deg_ax + AdamsDeg(r, r + stem_map);
        const auto& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);

        size_t first_PC = GetFirstIndexOnLevel(sc_a, LEVEL_PERM);
        size_t last_PC = GetFirstIndexOnLevel(sc_a, LEVEL_PERM + 1);
        for (size_t i = first_PC; i < last_PC; ++i) {
            Poly poly_a = Indices2Poly(sc_a.basis[i], ring.basis.at(deg_a));

            int1d ax;
            if (x.empty())
                ;
            else if (deg_ax.t > t_max)
                ax = NULL_DIFF;
            else {
                if (iCw.isRing) {
                    Poly poly_ax = ring.gb.Reduce(poly_a * poly_x);
                    ax = poly_ax ? Poly2Indices(poly_ax, ring.basis.at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), ring.nodes_ss, deg_ax, LEVEL_PERM);
                }
                else {
                    Mod mod_ax = modules_[iCw.index].gb.Reduce(poly_a * mod_x);
                    ax = mod_ax ? Mod2Indices(mod_ax, std::get<1>(basis1)->at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), *cofseq.nodes_ss[iCs], deg_ax, LEVEL_PERM);
                }
            }

            int1d adx;
            if (dx.empty())
                ;
            else if (deg_adx.t > t_max2)
                adx = NULL_DIFF;
            else if (iCw_next.isRing) {
                Poly poly_adx = ring.gb.Reduce(poly_a * poly_dx);
                if (poly_adx) {
                    adx = Poly2Indices(poly_adx, ring.basis.at(deg_adx));
                    adx = Residue(std::move(adx), ring.nodes_ss, deg_adx, LEVEL_PERM);
                }
            }
            else {
                Mod poly_adx = modules_[iCw_next.index].gb.Reduce(poly_a * mod_dx);
                if (poly_adx)
                    adx = Residue(Mod2Indices(poly_adx, std::get<1>(basis2)->at(deg_adx)), *cofseq.nodes_ss[iCs_next], deg_adx, LEVEL_PERM);
            }

            if ((!ax.empty() && ax != NULL_DIFF) || (!adx.empty() && adx != NULL_DIFF)) {
                SetDiffScCofseq(cofseq, iCs, deg_ax, ax, adx, r, flag);
                ++count;
            }
        }
    }
    return count;
}
