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
    if (r_max < LEVEL_MIN)
        return false;
    for (int r1 = LEVEL_MIN; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (ut::has(nodes_ss.front(), d_src))
            if (GetMaxLevelWithND(ut::GetRecentValue(nodes_ss, d_src)) >= LEVEL_MAX - r1)
                return true;
    }
    return false;
}

size_t Diagram::GetFirstIndexOfFixedLevels(const Staircases1d& nodes_ss, AdamsDeg deg, int level_min)
{
    auto& sc = ut::GetRecentValue(nodes_ss, deg);
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
    auto& sc = ut::GetRecentValue(nodes_ss, deg);
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
        auto& sc_tgt = ut::GetRecentValue(nodes_ss, deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_tgt, LEVEL_MAX - r + 1) - result.first;
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
        auto& sc_src = ut::GetRecentValue(nodes_ss, deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_src, LEVEL_MAX - r + 1) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Diagram::NextRTgt(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r) const
{
    for (int r1 = r; r1 <= R_PERM; ++r1) {
        if (deg.t + r1 - 1 > t_max)
            return r1;
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 - 1};
        if (AboveS0Vanishing(d_tgt) && !ut::has(nodes_ss.front(), d_tgt))
            return R_PERM;
        auto [_, count] = CountPossDrTgt(nodes_ss, t_max, d_tgt, r1);
        if (count > 0)
            return r1;
    }
    return R_PERM;
}

int Diagram::NextRSrc(const Staircases1d& nodes_ss, AdamsDeg deg, int r) const
{
    int r_max = std::min(r, deg.s);
    for (int r1 = r_max; r1 >= LEVEL_MIN; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        auto [_, count] = CountPossDrSrc(nodes_ss, d_src, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

void Diagram::CacheNullDiffs(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, DeduceFlag flag, NullDiff1d& nds) const
{
    nds.clear();
    const Staircase& sc = ut::GetRecentValue(nodes_ss, deg);
    for (size_t i = 0; i < sc.diffs.size(); ++i) {
        if (sc.diffs[i] != NULL_DIFF)
            continue;

        NullDiff nd;
        if (sc.levels[i] > /*LEVEL_PERM*/ LEVEL_MAX - 4) { ////
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

/* Return d_r(x) */
int1d Diagram::GetDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, int r) const
{
    if (x.empty())
        return int1d{};
    const Staircase& sc = ut::GetRecentValue(nodes_ss, deg_x);
    size_t first = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    size_t last = GetFirstIndexOfNullOnLevel(sc, LEVEL_MAX - r);
    /* Compute x mod [0,first) */
    int1d x1 = lina::Residue(sc.basis.begin(), sc.basis.begin() + first, x);

    /* If x is in [first,last) */
    if (lina::Residue(sc.basis.begin() + first, sc.basis.begin() + last, x1).empty())
        return lina::GetImage(sc.basis.begin() + first, sc.basis.begin() + last, sc.diffs.begin() + first, sc.diffs.begin() + last, x1);
    else
        return NULL_DIFF;
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

void Diagram::SetDiffSc(size_t iCw, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, DeduceFlag flag)
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
    auto& nodes_ss = !(iCw & FLAG_MOD) ? rings_[iCw].nodes_ss : modules_[iCw ^ FLAG_MOD].nodes_ss;
    const Staircase& sc = ut::GetRecentValue(nodes_ss, deg_x);
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
    if ((flag & DeduceFlag::cofseq) && r == R_PERM) {
        const auto& ind_cofs = !(iCw & FLAG_MOD) ? rings_[iCw].ind_cofs : modules_[iCw ^ FLAG_MOD].ind_cofs;
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

void Diagram::SetImageSc(size_t iCw, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, DeduceFlag flag)
{
    const auto& name = !(iCw & FLAG_MOD) ? rings_[iCw].name : modules_[iCw ^ FLAG_MOD].name;
    auto& nodes_ss = !(iCw & FLAG_MOD) ? rings_[iCw].nodes_ss : modules_[iCw ^ FLAG_MOD].nodes_ss;
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, r - 1};

    if (nodes_ss.size() == 2 && deg_dx == AdamsDeg(16, 162 + 16) && r == 2 && dx_ == int1d{6} && name == "Joker") ////
        throw InteruptAndSaveException(0, "bug");

    /* If dx is in Im(d_{r-1}) then x is in Ker(d_r) */
    const Staircase& sc = ut::GetRecentValue(nodes_ss, deg_dx);
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
                    std::cout << "debug\n";
                Logger::LogSSException(int(nodes_ss.size() - 2), name, deg_dx, dx, r, 0x7U, deg_leibniz_, a_leibniz_);
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
    if (flag & DeduceFlag::cofseq) {
        const auto& ind_cofs = !(iCw & FLAG_MOD) ? rings_[iCw].ind_cofs : modules_[iCw ^ FLAG_MOD].ind_cofs;
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

void Diagram::SetDiffScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, DeduceFlag flag)
{
    const size_t iCs1 = (iCs + 2) % 3;
    const size_t iCs2 = (iCs + 1) % 3;
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_ss = *cofseq.nodes_ss[iCs];
    int stem_map = cofseq.degMap[iCs].stem();
    int stem_map1 = cofseq.degMap[iCs1].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, stem_map + r};

    /* NULL_DIFF checking */
    if (x_ == NULL_DIFF) {
        SetImageScCofseq(cofseq, iCs2, deg_dx, dx, NULL_DIFF, r, flag);
        return;
    }

    /* Triangularize x */
    if (x_.empty()) {  //// TODO: Remove this
        if (dx != NULL_DIFF && !dx.empty())
            SetImageScCofseq(cofseq, iCs2, deg_dx, dx, NULL_DIFF, r - 1, flag);
        return;
    }
    int1d x = Residue(x_, nodes_ss, deg_x, LEVEL_PERM);
    const Staircase& sc = ut::GetRecentValue(nodes_cofseq, deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, std::move(x));
    if (x.empty()) {
        if (dx != NULL_DIFF && !dx.empty())
            SetImageScCofseq(cofseq, iCs2, deg_dx, dx, NULL_DIFF, r - 1, flag);
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
            SetImageScCofseq(cofseq, iCs2, deg_image_new, image_new, NULL_DIFF, level_image_new - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, stem_map1 + r_image};
            SetDiffScCofseq(cofseq, iCs1, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }

    /* Set image */
    if (dx != NULL_DIFF && !dx.empty())
        SetImageScCofseq(cofseq, iCs2, deg_dx, dx, x, r, flag);
}

void Diagram::SetImageScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, DeduceFlag flag)
{
    const size_t iCs1 = (iCs + 2) % 3;
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_ss = *cofseq.nodes_ss[iCs];
    auto& name = cofseq.nameCw[iCs];
    int stem_map = cofseq.degMap[iCs].stem();
    int stem_map1 = cofseq.degMap[iCs1].stem();
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, stem_map1 + r};

    int1d dx = Residue(dx_, nodes_ss, deg_dx, LEVEL_PERM);
    const Staircase& sc = ut::GetRecentValue(nodes_cofseq, deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, std::move(dx));
    if (dx.empty()) {
        if (x != NULL_DIFF)
            SetDiffScCofseq(cofseq, iCs1, deg_x, x, NULL_DIFF, r + 1, flag);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == NULL_DIFF) {
        /* If the source is uncertain, insert it to the end of level r. */
        size_t first_rp1 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp1, std::move(dx));
        if (!dx.empty()) {
            UpdateStaircase(nodes_cofseq, deg_dx, sc, first_rp1, dx, x, r, image_new, level_image_new);

            /* cofseq to ss */
            if (r == -1) {
                Staircase& sc1 = ut::GetRecentValue(nodes_cofseq, deg_dx);
                pop_front(sc1);
                Logger::LogDiffInv(int(nodes_cofseq.size() - 2), EnumReason::cofseq_b, cofseq.nameCw[iCs], deg_dx - AdamsDeg(R_PERM + 1, R_PERM), deg_dx, {}, dx, R_PERM + 1);
                size_t iCw = cofseq.isRing[iCs] ? cofseq.indexCw[iCs] : cofseq.indexCw[iCs] | FLAG_MOD;
                a_leibniz_ = nullptr;
                SetImageSc(iCw, deg_dx, dx, NULL_DIFF, R_PERM, flag);
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
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, stem_map1 + r_image};
            SetDiffScCofseq(cofseq, iCs1, deg_image_new, image_new, NULL_DIFF, r_image + 1, flag);
        }
    }
}

void Diagram::ReSetScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg, DeduceFlag flag)
{
    const auto& nodes_ss = cofseq.isRing[iCs] ? rings_[cofseq.indexCw[iCs]].nodes_ss : modules_[cofseq.indexCw[iCs]].nodes_ss;
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];

    auto& sc_ss = ut::GetRecentValue(nodes_ss, deg);
    Staircase sc = ut::GetRecentValue(nodes_cofseq, deg);

    const size_t iCs1 = (iCs + 2) % 3;
    int stem_map = cofseq.degMap[iCs].stem();
    int stem_map1 = cofseq.degMap[iCs1].stem();

    Staircase sc1;
    int2d images;
    int1d levels;
    for (size_t i = 0; i < sc.basis.size(); ++i) {
        size_t first_PC = GetFirstIndexOnLevel(sc_ss, LEVEL_PERM);
        sc.basis[i] = lina::Residue(sc_ss.basis.begin(), sc_ss.basis.begin() + first_PC, std::move(sc.basis[i]));
        for (size_t j = 0; j < sc1.basis.size(); ++j) {
            if (std::binary_search(sc.basis[i].begin(), sc.basis[i].end(), sc1.basis[j][0])) {
                sc.basis[i] = lina::add(sc.basis[i], sc1.basis[j]);
                if (sc.levels[i] == sc1.levels[j] && sc.diffs[i] != NULL_DIFF)
                    sc.diffs[i] = lina::add(sc.diffs[i], sc1.diffs[j]);
            }
        }
        if (sc.basis[i].empty()) {
            if (!sc.diffs[i].empty() && sc.diffs[i] != NULL_DIFF) {
                images.push_back(std::move(sc.diffs[i]));
                levels.push_back(LEVEL_MAX - sc.levels[i]);
            }
        }
        else {
            sc1.basis.push_back(std::move(sc.basis[i]));
            sc1.diffs.push_back(std::move(sc.diffs[i]));
            sc1.levels.push_back(sc.levels[i]);
        }
    }

    nodes_cofseq.back()[deg] = std::move(sc1);

    /* Add images and cycles */
    for (size_t i = 0; i < levels.size(); ++i) {
        if (levels[i] < LEVEL_MAX / 2) {
            /* Set an image */
            AdamsDeg deg_image_new = deg + AdamsDeg{levels[i], stem_map + levels[i]};
            SetImageScCofseq(cofseq, (iCs + 1) % 3, deg_image_new, images[i], NULL_DIFF, levels[i] - 1, flag);
        }
        else {
            /* Set a cycle */
            int r_image = LEVEL_MAX - levels[i];
            AdamsDeg deg_image_new = deg - AdamsDeg{r_image, stem_map1 + r_image};
            SetDiffScCofseq(cofseq, iCs1, deg_image_new, images[i], NULL_DIFF, r_image + 1, flag);
        }
    }

    return;
}

int Diagram::SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, DeduceFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    {
        auto& nodes_ss = ring.nodes_ss;
        auto& basis = ring.basis;
        int t_max = ring.t_max;
        for (auto& [deg_a, _] : nodes_ss.front()) {
            const AdamsDeg deg_ax = deg_x + deg_a;
            const Staircase& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                if (sc_a.levels[i] > LEVEL_MAX - r_min)
                    break;
                const int r_a = LEVEL_MAX - sc_a.levels[i];
                const int R = std::min(r, r_a);
                if (R == r_a && sc_a.diffs[i] == NULL_DIFF)
                    continue;
                AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                AdamsDeg deg_dxy = deg_ax + AdamsDeg(R, R - 1);

                Poly poly_dx = (R == r && !dx.empty()) ? Indices2Poly(dx, ring.basis.at(deg_dx)) : Poly();  //// TODO: precompute
                Poly poly_a = Indices2Poly(sc_a.basis[i], basis.at(deg_a));
                Poly poly_ax = ring.gb.Reduce(poly_x * poly_a);
                int1d ax = poly_ax ? Poly2Indices(poly_ax, basis.at(deg_ax)) : int1d{};

                Poly poly_da = (R == r_a) ? Indices2Poly(sc_a.diffs[i], basis.at(deg_a + AdamsDeg(R, R - 1))) : Poly();
                Poly poly_dax = ring.gb.Reduce(poly_x * poly_da + poly_dx * poly_a);
                int1d dax;
                if (poly_dax) {
                    if (deg_dxy.t <= t_max)
                        dax = Poly2Indices(poly_dax, basis.at(deg_dxy));
                    else
                        dax = NULL_DIFF;
                }

                if (!ax.empty() || !dax.empty()) {
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetDiffSc(iRing, deg_ax, ax, dax, R, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    for (size_t iMod : ring.ind_mods) {
        auto& mod = modules_[iMod];
        auto& nodes_ss = mod.nodes_ss;
        auto& basis = mod.basis;
        int t_max = mod.t_max;
        for (auto& [deg_a, _] : nodes_ss.front()) {
            AdamsDeg deg_ax = deg_x + deg_a;
            const Staircase& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
            if (deg_ax.t > t_max)
                break;
            for (size_t i = 0; i < sc_a.levels.size(); ++i) {
                if (sc_a.levels[i] > LEVEL_MAX - r_min)
                    break;
                const int r_a = LEVEL_MAX - sc_a.levels[i];
                const int R = std::min(r, r_a);
                if (R == r_a && sc_a.diffs[i] == NULL_DIFF)
                    continue;
                AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                AdamsDeg deg_dax = deg_ax + AdamsDeg(R, R - 1);

                Poly poly_dx = (R == r && !dx.empty()) ? Indices2Poly(dx, ring.basis.at(deg_dx)) : Poly();

                Mod poly_a = Indices2Mod(sc_a.basis[i], basis.at(deg_a));
                Mod poly_ax = mod.gb.Reduce(poly_x * poly_a);
                int1d ax = poly_ax ? Mod2Indices(poly_ax, basis.at(deg_x + deg_a)) : int1d{};

                Mod poly_da = (R == r_a) ? Indices2Mod(sc_a.diffs[i], basis.at(deg_a + AdamsDeg(R, R - 1))) : Mod();
                Mod poly_dax = mod.gb.Reduce(poly_x * poly_da + poly_dx * poly_a);
                int1d dax;
                if (poly_dax) {
                    if (deg_dax.t <= t_max)
                        dax = Mod2Indices(poly_dax, basis.at(deg_dax));
                    else
                        dax = NULL_DIFF;
                }

                if (!ax.empty() || !dax.empty()) {
                    deg_leibniz_ = deg_a;
                    a_leibniz_ = &sc_a.basis[i];
                    SetDiffSc(iMod | FLAG_MOD, deg_ax, ax, dax, R, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, DeduceFlag flag)
{
    int count = 0;

    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& nodes_ss = mod.nodes_ss;
    auto& basis = mod.basis;
    auto& gb = mod.gb;
    const int t_max = mod.t_max;
    Mod poly_x = Indices2Mod(x, basis.at(deg_x));

    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const Staircase& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);
        AdamsDeg deg_ax = deg_x + deg_a;
        if (deg_ax.t > t_max)
            break;
        for (size_t i = 0; i < sc_a.levels.size(); ++i) {
            if (sc_a.levels[i] > LEVEL_MAX - r_min)
                break;
            int r_a = LEVEL_MAX - sc_a.levels[i];
            int R = std::min(r, r_a);
            if (R == r_a && sc_a.diffs[i] == NULL_DIFF)
                continue;
            AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
            AdamsDeg deg_dax = deg_ax + AdamsDeg(R, R - 1);

            Mod poly_dx = (R == r && !dx.empty()) ? Indices2Mod(dx, basis.at(deg_dx)) : Mod();

            Poly poly_a = Indices2Poly(sc_a.basis[i], ring.basis.at(deg_a));
            Mod poly_ax = gb.Reduce(poly_a * poly_x);
            int1d ax = poly_ax ? Mod2Indices(poly_ax, basis.at(deg_x + deg_a)) : int1d{};

            Poly poly_da = (R == r_a) ? Indices2Poly(sc_a.diffs[i], ring.basis.at(deg_a + AdamsDeg(R, R - 1))) : Poly();
            Mod poly_dax = gb.Reduce(poly_da * poly_x + poly_a * poly_dx);
            int1d dax;
            if (poly_dax) {
                if (deg_dax.t <= t_max)
                    dax = Mod2Indices(poly_dax, basis.at(deg_dax));
                else
                    dax = NULL_DIFF;
            }

            if (!ax.empty() || !dax.empty()) {
                deg_leibniz_ = deg_a;
                a_leibniz_ = &sc_a.basis[i];
                SetDiffSc(iMod | FLAG_MOD, deg_ax, ax, dax, R, flag);
                a_leibniz_ = nullptr;
                ++count;
            }
        }
    }
    return count;
}

int Diagram::SetRingDiffLeibnizV2(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, DeduceFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    const AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    if (!ut::has(ring.nodes_ss.front(), deg_dx))
        return count;
    int depth = int(ring.nodes_ss.size() - 2);
    const Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
    const Staircase& sc_dx = ut::GetRecentValue(ring.nodes_ss, deg_dx);
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
            const Staircase& sc_y = ut::GetRecentValue(nodes_ss, deg_y);
            for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
                if (sc_y.levels[i] > LEVEL_MAX - r)
                    break;
                const int r_y = LEVEL_MAX - sc_y.levels[i];
                if (r == r_y && sc_y.diffs[i] == NULL_DIFF)
                    continue;

                Poly poly_y = Indices2Poly(sc_y.basis[i], basis.at(deg_y));
                Poly poly_xy = ring.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = Poly2Indices(poly_xy, basis.at(deg_xy));

                    bool ydx_always_zero = true;
                    for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : two_expansion(j))
                            dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                        Poly poly_dx = Indices2Poly(dx, basis.at(deg_dx));
                        Poly poly_ydx = ring.gb.Reduce(poly_dx * poly_y);
                        if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Poly2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                            ydx_always_zero = false;
                            break;
                        }
                    }
                    if (ydx_always_zero) {
                        Poly poly_dy = (r == r_y) ? Indices2Poly(sc_y.diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
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
            const Staircase& sc_y = ut::GetRecentValue(nodes_ss, deg_y);
            for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
                if (sc_y.levels[i] > LEVEL_MAX - r)
                    break;
                const int r_y = LEVEL_MAX - sc_y.levels[i];
                if (r == r_y && sc_y.diffs[i] == NULL_DIFF)
                    continue;

                Mod poly_y = Indices2Mod(sc_y.basis[i], basis.at(deg_y));
                Mod poly_xy = mod.gb.Reduce(poly_x * poly_y);
                if (poly_xy) {
                    int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    bool ydx_always_zero = true;
                    for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                        dx.clear();
                        for (int k : two_expansion(j))
                            dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                        Poly poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
                        Mod poly_ydx = mod.gb.Reduce(poly_dx * poly_y);
                        if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Mod2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                            ydx_always_zero = false;
                            break;
                        }
                    }
                    if (ydx_always_zero) {
                        Mod poly_dy = (r == r_y) ? Indices2Mod(sc_y.diffs[i], basis.at(deg_y + AdamsDeg(r, r - 1))) : Mod();
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
                for (int k : two_expansion(j))
                    dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                int1d fdx = map->map(dx, deg_dx, *this);
                if (!fdx.empty() && !IsZeroOnLevel(ut::GetRecentValue(rings_[map->to].nodes_ss, deg_dx), fdx, r)) {
                    fdx_always_zero = false;
                    break;
                }
            }
            if (fdx_always_zero && IsNewDiff(rings_[map->to].nodes_ss, deg_x, fx, {}, r)) {
                if (!printed_dx) {
                    printed_dx = true;
                    Logger::LogNullDiff(depth, ring.name, deg_x, x, r);
                }
                Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, rings_[map->to].name), deg_x, fx, {}, r);
                count += SetRingDiffGlobal(map->to, deg_x, fx, {}, r, true, flag);
            }
        }
    }

    return count;
}

int Diagram::SetModuleDiffLeibnizV2(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, DeduceFlag flag)
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
    int depth = int(nodes_ss.size() - 2);
    const Mod poly_x = Indices2Mod(x, basis.at(deg_x));
    const Staircase& sc_dx = ut::GetRecentValue(nodes_ss, deg_dx);
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
        const Staircase& sc_y = ut::GetRecentValue(ring.nodes_ss, deg_y);
        for (size_t i = 0; i < sc_y.levels.size(); ++i) { /* Loop over y */
            if (sc_y.levels[i] > LEVEL_MAX - r)
                break;
            const int r_y = LEVEL_MAX - sc_y.levels[i];
            if (r == r_y && sc_y.diffs[i] == NULL_DIFF)
                continue;

            Poly poly_y = Indices2Poly(sc_y.basis[i], ring.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            if (poly_xy) {
                int1d xy = Mod2Indices(poly_xy, basis.at(deg_x + deg_y));

                bool ydx_always_zero = true;
                for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                    dx.clear();
                    for (int k : two_expansion(j))
                        dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                    Mod poly_dx = Indices2Mod(dx, basis.at(deg_dx));
                    Mod poly_ydx = mod.gb.Reduce(poly_y * poly_dx);
                    if (poly_ydx && !IsZeroOnLevel(ut::GetRecentValue(nodes_ss, deg_dxy), Mod2Indices(poly_ydx, basis.at(deg_dxy)), r)) {
                        ydx_always_zero = false;
                        break;
                    }
                }
                if (ydx_always_zero) {
                    Poly poly_dy = (r == r_y) ? Indices2Poly(sc_y.diffs[i], ring.basis.at(deg_y + AdamsDeg(r, r - 1))) : Poly();
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

        size_t to;
        if (map->IsToRing(to)) {
            auto fx = map->map(x, deg_x, *this);
            if (!fx.empty()) {
                AdamsDeg deg_fx = deg_x + map->deg;
                AdamsDeg deg_fdx = deg_fx + AdamsDeg(r, r - 1);
                bool fdx_always_zero = true;
                for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                    dx.clear();
                    for (int k : two_expansion(j))
                        dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                    int1d fdx = map->map(dx, deg_dx, *this);
                    if (!fdx.empty() && !IsZeroOnLevel(ut::GetRecentValue(rings_[to].nodes_ss, deg_fdx), fdx, r)) {
                        fdx_always_zero = false;
                        break;
                    }
                }
                if (fdx_always_zero && IsNewDiff(rings_[to].nodes_ss, deg_fx, fx, {}, r)) {
                    if (!printed_dx) {
                        printed_dx = true;
                        Logger::LogNullDiff(depth, mod.name, deg_x, x, r);
                    }
                    Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, rings_[to].name), deg_fx, fx, {}, r);
                    count += SetRingDiffGlobal(to, deg_fx, fx, {}, r, true, flag);
                }
            }
        }
        else {
            auto fx = map->map(x, deg_x, *this);
            if (!fx.empty()) {
                AdamsDeg deg_fx = deg_x + map->deg;
                AdamsDeg deg_fdx = deg_fx + AdamsDeg(r, r - 1);
                bool fdx_always_zero = true;
                for (unsigned j = 1; j < j_max; ++j) { /* Loop over dx */
                    dx.clear();
                    for (int k : two_expansion(j))
                        dx = lina::add(dx, sc_dx.basis[(size_t)(first_dx + k)]);
                    int1d fdx = map->map(dx, deg_dx, *this);
                    if (!fdx.empty() && !IsZeroOnLevel(ut::GetRecentValue(modules_[to].nodes_ss, deg_fdx), fdx, r)) {
                        fdx_always_zero = false;
                        break;
                    }
                }
                if (fdx_always_zero && IsNewDiff(modules_[to].nodes_ss, deg_fx, fx, {}, r)) {
                    if (!printed_dx) {
                        printed_dx = true;
                        Logger::LogNullDiff(depth, mod.name, deg_x, x, r);
                    }
                    Logger::LogDiff(depth, EnumReason::deduce_xy, fmt::format("({}) {}", map->display, modules_[to].name), deg_fx, fx, {}, r);
                    count += SetModuleDiffGlobal(to, deg_fx, fx, {}, r, true, flag);
                }
            }
        }
    }

    return count;
}

int Diagram::SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, DeduceFlag flag)
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
            const Staircase& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
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
                    SetImageSc(iRing, deg_ax, ax, NULL_DIFF, r, flag);
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
            const Staircase& sc_a = ut::GetRecentValue(nodes_ss, deg_a);
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
                    SetImageSc(iMod | FLAG_MOD, deg_ax, ax, NULL_DIFF, r, flag);
                    a_leibniz_ = nullptr;
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, DeduceFlag flag)
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
        Logger::LogSSException(int(nodes_ss.size() - 2), mod.name, deg_x, x, r_original, 0xda298807U, deg_leibniz_, a_leibniz_);
        throw SSException(0xda298807U, "No source for the image.");
    }

    Mod poly_x = Indices2Mod(x, mod.basis.at(deg_x));
    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const Staircase& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);
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
                SetImageSc(iMod | FLAG_MOD, deg_ax, ax, NULL_DIFF, r, flag);
                a_leibniz_ = nullptr;
                ++count;
            }
        }
    }
    return count;
}

int Diagram::SetDiffLeibnizCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, DeduceFlag flag)
{
    int count = 0;
    size_t iCs2 = (iCs + 1) % 3;
    bool isRing1 = cofseq.isRing[iCs], isRing2 = cofseq.isRing[iCs2];
    auto& ring = isRing1 ? rings_[cofseq.indexCw[iCs]] : rings_[modules_[cofseq.indexCw[iCs]].iRing];
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    using BasisVariant = std::variant<std::map<AdamsDeg, Mon1d>*, std::map<AdamsDeg, MMod1d>*>;
    BasisVariant basis1 = isRing1 ? BasisVariant(&ring.basis) : BasisVariant(&modules_[cofseq.indexCw[iCs]].basis);
    BasisVariant basis2 = isRing2 ? BasisVariant(&ring.basis) : BasisVariant(&modules_[cofseq.indexCw[iCs2]].basis);
    Poly poly_x, poly_dx;
    Mod mod_x, mod_dx;
    if (!x.empty()) {
        if (isRing1)
            poly_x = Indices2Poly(x, std::get<0>(basis1)->at(deg_x));
        else
            mod_x = Indices2Mod(x, std::get<1>(basis1)->at(deg_x));
    }
    int stem_map = cofseq.degMap[iCs].stem();
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r + stem_map);
    if (!dx.empty()) {
        if (isRing2)
            poly_dx = Indices2Poly(dx, ring.basis.at(deg_dx));
        else
            mod_dx = Indices2Mod(dx, std::get<1>(basis2)->at(deg_dx));
    }
    int t_max = cofseq.t_max[iCs];
    int t_max2 = cofseq.t_max[iCs2];
    for (auto& [deg_a, _] : ring.nodes_ss.front()) {
        const AdamsDeg deg_ax = deg_x + deg_a;
        AdamsDeg deg_adx = deg_ax + AdamsDeg(r, r + stem_map);
        const Staircase& sc_a = ut::GetRecentValue(ring.nodes_ss, deg_a);

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
                if (isRing1) {
                    Poly poly_ax = ring.gb.Reduce(poly_a * poly_x);
                    ax = poly_ax ? Poly2Indices(poly_ax, ring.basis.at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), ring.nodes_ss, deg_ax, LEVEL_PERM);
                }
                else {
                    Mod mod_ax = modules_[cofseq.indexCw[iCs]].gb.Reduce(poly_a * mod_x);
                    ax = mod_ax ? Mod2Indices(mod_ax, std::get<1>(basis1)->at(deg_ax)) : int1d{};
                    ax = Residue(std::move(ax), *cofseq.nodes_ss[iCs], deg_ax, LEVEL_PERM);
                }
            }

            int1d adx;
            if (dx.empty())
                ;
            else if (deg_adx.t > t_max2)
                adx = NULL_DIFF;
            else if (isRing2) {
                Poly poly_adx = ring.gb.Reduce(poly_a * poly_dx);
                if (poly_adx) {
                    adx = Poly2Indices(poly_adx, ring.basis.at(deg_adx));
                    adx = Residue(std::move(adx), ring.nodes_ss, deg_adx, LEVEL_PERM);
                }
            }
            else {
                Mod poly_adx = modules_[cofseq.indexCw[iCs2]].gb.Reduce(poly_a * mod_dx);
                if (poly_adx)
                    adx = Residue(Mod2Indices(poly_adx, std::get<1>(basis2)->at(deg_adx)), *cofseq.nodes_ss[iCs2], deg_adx, LEVEL_PERM);
            }

            if ((!ax.empty() && ax != NULL_DIFF) || (!adx.empty() && adx != NULL_DIFF)) {
                SetDiffScCofseq(cofseq, iCs, deg_ax, ax, adx, r, flag);
                ++count;
            }
        }
    }
    return count;
}
