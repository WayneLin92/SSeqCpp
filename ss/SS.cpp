#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <fmt/ranges.h>

std::string GetComplexName(const std::string& db)
{
    auto p = db.find_last_of("/\\");
    std::string db1 = db;
    if (p != std::string::npos)
        db1 = db.substr(p + 1);
    if (db1.size() < 2)
        throw MyException(0x68360e1dU, "The file name of the database is too short");
    std::string table;
    if (db1[1] == '0')
        table = "S0";
    else if (db1[1] == '2')
        table = "C2";
    else if (db1[1] == 'e')
        table = "Ceta";
    else if (db1[1] == 'n')
        table = "Cnu";
    else if (db1[1] == 's')
        table = "Csigma";
    else
        throw MyException(0x331d5e90U, db1 + " is not recognized as a db name");
    return table;
}

std::string GetE2TablePrefix(const std::string& db)
{
    return GetComplexName(db) + "_AdamsE2";
}

int GetTopCellT(const std::string& db)
{
    auto p = db.find_last_of("/\\");
    std::string db1 = db;
    if (p != std::string::npos)
        db1 = db.substr(p + 1);
    if (db.size() < 2)
        throw MyException(0xeedeca68U, "The file name of the database is too short");
    std::string table;
    if (db1[1] == '2')
        return 1;
    else if (db1[1] == 'e')
        return 2;
    else if (db1[1] == 'n')
        return 4;
    else if (db1[1] == 's')
        return 8;
    else
        throw MyException(0x8cafd31fU, db + " is not recognized as a db name");
    return -10000;
}

bool Diagram::IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg, int r_max)
{
    r_max = std::min(r_max, deg.s);
    for (int r1 = kLevelMin; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (nodes_ss.front().find(d_src) != nodes_ss.front().end())
            if (GetMaxLevelWithNull(GetRecentStaircase(nodes_ss, d_src)) >= LEVEL_MAX - r1)
                return true;
    }
    return false;
}

bool Diagram::IsPossSrc(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r_min)
{
    int r_max = (deg.t - deg.s * 3 + 2) / 2;
    for (int r = r_min; r <= r_max; ++r) {
        AdamsDeg d_tgt = deg + AdamsDeg{r, r - 1};
        if (d_tgt.t > t_max)
            return true;
        if (nodes_ss.front().find(d_tgt) != nodes_ss.front().end()) {
            if (GetMaxLevelWithNull(GetRecentStaircase(nodes_ss, d_tgt)) >= r)
                return true;
        }
    }
    return false;
}

int Diagram::GetFirstFixedLevelForPlot(const Staircases1d& nodes_ss, AdamsDeg deg)
{
    auto& sc = GetRecentStaircase(nodes_ss, deg);
    int result = LEVEL_MAX - kLevelMin;
    for (size_t i = sc.levels.size(); i-- > 0 && sc.levels[i] >= kLevelPC;) {
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

size_t Diagram::GetFirstIndexOfFixedLevels(const Staircases1d& nodes_ss, AdamsDeg deg, int level_min)
{
    auto& sc = GetRecentStaircase(nodes_ss, deg);
    size_t result = sc.levels.size();
    for (size_t i = sc.levels.size(); i-- > 0;) {
        if (sc.diffs[i] == int1d{-1} || sc.levels[i] < level_min)
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

std::pair<int, int> Diagram::CountPossDrTgt(const Staircases1d& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const
{
    std::pair<int, int> result;
    if (nodes_ss.front().find(deg_tgt) != nodes_ss.front().end()) {
        const Staircase& sc_tgt = GetRecentStaircase(nodes_ss, deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_tgt, LEVEL_MAX - r + 1) - result.first;
    }
    else if (deg_tgt.t > t_max)
        result = {-1, 100000};
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Diagram::CountPossDrSrc(const Staircases1d& nodes_ss, const AdamsDeg& deg_src, int r) const
{
    std::pair<int, int> result;
    if (nodes_ss.front().find(deg_src) != nodes_ss.front().end()) {
        const Staircase& sc_src = GetRecentStaircase(nodes_ss, deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, LEVEL_MAX - r);
        result.second = (int)GetFirstIndexOfFixedLevels(nodes_ss, deg_src, LEVEL_MAX - r + 1) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Diagram::NextRTgt(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r) const
{
    AdamsDeg deg_tgt;
    int count = 0, index = -1;
    const Staircase& sc = GetRecentStaircase(nodes_ss, deg);
    int r_max = (deg.t - deg.s * 3 + 2) / 2;
    for (int r1 = r; r1 <= r_max; ++r1) {
        if (deg.t + r1 - 1 > t_max)
            return r1;
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 - 1};
        auto [index, count] = CountPossDrTgt(nodes_ss, t_max, d_tgt, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

int Diagram::NextRSrc(const Staircases1d& nodes_ss, AdamsDeg deg, int r) const
{
    AdamsDeg deg_src;
    int count = 0, index = -1;
    int r_max = std::min(r, deg.s - 1);
    const Staircase& sc = GetRecentStaircase(nodes_ss, deg);
    for (int r1 = r_max; r1 >= kLevelMin; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        auto [index, count] = CountPossDrSrc(nodes_ss, d_src, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

void Diagram::CacheNullDiffs(size_t iSS, AdamsDeg deg, DeduceFlag flag, NullDiff1d& nds)
{
    auto& nodes_ss = *all_basis_ss_[iSS];
    int t_max = all_t_max_[iSS];

    nds.clear();
    const Staircase& sc = GetRecentStaircase(nodes_ss, deg);
    for (size_t i = 0; i < sc.diffs.size(); ++i) {
        if (sc.diffs[i] == int1d{-1}) {
            NullDiff nd;
            if (sc.levels[i] > kLevelPC) {
                int r = LEVEL_MAX - sc.levels[i];
                AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};
                auto [index, count] = CountPossDrTgt(nodes_ss, t_max, deg_tgt, r);
                nd.r = r;
                nd.first = index;
                nd.count = count;
                if (nd.count > 10)
                    continue;
            }
            else if (sc.levels[i] < LEVEL_MAX / 2) {
                int r = sc.levels[i];
                AdamsDeg deg_src = deg - AdamsDeg{r, r - 1};
                auto [index, count] = CountPossDrSrc(nodes_ss, deg_src, r);
                nd.r = -r;
                nd.first = index;
                nd.count = count;
                if (nd.count > 10)
                    continue;
            }
            else
                continue;

            size_t j = i + 1;
            while (j < sc.levels.size() && sc.diffs[j] == int1d{-1} && sc.levels[i] == sc.levels[j])
                ++j;
            if (!(flag & DeduceFlag::all_x)) {
                for (size_t k = 0; k < j - i; ++k) {
                    nd.x = sc.basis[i + k];
                    nds.push_back(nd);
                }
            }
            else {
                const unsigned k_max = unsigned(1) << (j - i);
                for (unsigned k = 1; k < k_max; ++k) {
                    nd.x.clear();
                    for (int l : two_expansion(k))
                        nd.x = lina::AddVectors(nd.x, sc.basis[i + l]);
                    nds.push_back(nd);
                }
            }
            i = j - 1;
        }
    }
}

/* Return d_r(x) */
int1d Diagram::GetDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, int r) const
{
    if (x.empty())
        return int1d{};
    const Staircase& sc = GetRecentStaircase(nodes_ss, deg_x);
    size_t first = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    size_t last = GetFirstIndexOfNullOnLevel(sc, LEVEL_MAX - r);
    /* Compute x mod [0,first) */
    int1d x1 = lina::Residue(sc.basis.begin(), sc.basis.begin() + first, x);

    /* If x is in [first,last) */
    if (lina::Residue(sc.basis.begin() + first, sc.basis.begin() + last, x1).empty())
        return lina::GetImage(sc.basis.begin() + first, sc.basis.begin() + last, sc.diffs.begin() + first, sc.diffs.begin() + last, x1);
    else
        return int1d{-1};
}

bool Diagram::IsNewDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    int1d dx1 = GetDiff(nodes_ss, deg_x, x, r);
    if (dx1 == int1d{-1})
        return true;
    if (nodes_ss.front().find(deg_dx) != nodes_ss.front().end()) {
        const Staircase& sc = GetRecentStaircase(nodes_ss, deg_dx);
        size_t first_r = GetFirstIndexOnLevel(sc, r);
        /* Check if dx-dx1 is trivial mod Im(d_{r-1}) */
        return !lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, lina::AddVectors(dx, dx1)).empty();
    }
    return false;
}  // TODO: replace int1d{-1}

void Diagram::SetDiffSc(std::string_view name, Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};

    /* If x is zero then dx is in Im(d_{r-1}) */
    if (x_.empty()) {
        if (dx != int1d{-1} && !dx.empty())
            SetImageSc(name, nodes_ss, deg_dx, dx, {-1}, r - 1);
        return;
    }

    const Staircase& sc = GetRecentStaircase(nodes_ss, deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, LEVEL_MAX - r);
    int1d x = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_Nmr, x_);
    if (x.empty()) {
        /* If x is in Ker(d_r) then dx is in Im(d_{r-1}) */
        if (dx != int1d{-1} && !dx.empty())
            SetImageSc(name, nodes_ss, deg_dx, dx, {-1}, r - 1);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (dx == int1d{-1}) {
        /* If the target is unknown, insert it to the end of level N-r. */
        size_t first_Nmrp1 = GetFirstIndexOnLevel(sc, LEVEL_MAX - r + 1);
        x = lina::Residue(sc.basis.begin() + first_Nmr, sc.basis.begin() + first_Nmrp1, x);
        if (!x.empty())
            UpdateStaircase(nodes_ss, deg_x, sc, first_Nmrp1, x, {-1}, LEVEL_MAX - r, image_new, level_image_new);
    }
    else if (dx.empty()) {
        /* If the target is zero, insert it to the end of level N-r-1 */
        UpdateStaircase(nodes_ss, deg_x, sc, first_Nmr, x, {-1}, LEVEL_MAX - r - 1, image_new, level_image_new);
    }
    else {
        /* Otherwise insert it to the beginning of level N-r */
        UpdateStaircase(nodes_ss, deg_x, sc, first_Nmr, x, dx, LEVEL_MAX - r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Add a d_{r1-1} image */
            AdamsDeg deg_image_new = deg_x + AdamsDeg{level_image_new, level_image_new - 1};
            SetImageSc(name, nodes_ss, deg_image_new, image_new, {-1}, level_image_new - 1);
        }
        else {
            /* Add a d_{r1} cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, r_image - 1};
            SetDiffSc(name, nodes_ss, deg_image_new, image_new, {}, r_image);
        }
    }

    /* Add image */
    if (dx != int1d{-1} && !dx.empty())
        SetImageSc(name, nodes_ss, deg_dx, dx, x, r);
}

void Diagram::SetImageSc(std::string_view name, Staircases1d& nodes_ss, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r)
{
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, r - 1};
    if (deg_x.s < 0) {
        Logger::LogException(int(nodes_ss.size() - 2), 0x7dc5fa8cU, "No source for the image. {} deg_dx={}, dx={}, r={}\n", name, deg_dx, dx_, r);
        throw SSException(0x7dc5fa8cU, "No source for the image.");
    }

    /* If dx is in Im(d_{r-1}) then x is in Ker(d_r) */
    const Staircase& sc = GetRecentStaircase(nodes_ss, deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    int1d dx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_r, dx_);
    if (dx.empty()) {
        if (x != int1d{-1} && !x.empty())
            SetDiffSc(name, nodes_ss, deg_x, x, {-1}, r + 1);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == int1d{-1}) {
        /* If the source is unknown, check if it can be hit and then insert it to the end of level r. */
        size_t first_rp2 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis.begin() + first_r, sc.basis.begin() + first_rp2, dx);
        if (!dx.empty()) {
            if (!IsPossTgt(nodes_ss, deg_dx, r)) {
                Logger::LogException(int(nodes_ss.size() - 2), 0x75989376U, "No source for the image. {} deg_dx={}, dx={}, r={}\n", name, deg_dx, dx, r);
                throw SSException(0x75989376U, "No source for the image.");
            }
            UpdateStaircase(nodes_ss, deg_dx, sc, first_rp2, dx, x, r, image_new, level_image_new);
        }
    }
    else {
        /* Otherwise insert it to the beginning of level r */
        UpdateStaircase(nodes_ss, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < LEVEL_MAX / 2) {
            /* Add a d_{r1-1} image */
            AdamsDeg deg_image_new = deg_dx + AdamsDeg{level_image_new, level_image_new - 1};
            SetImageSc(name, nodes_ss, deg_image_new, image_new, {-1}, level_image_new - 1);
        }
        else {
            /* Add a d_r1 cycle */
            int r_image = LEVEL_MAX - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, r_image - 1};
            SetDiffSc(name, nodes_ss, deg_image_new, image_new, {-1}, r_image + 1);
        }
    }
}

int Diagram::SetS0DiffLeibniz(AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry)
{
#ifndef NDEBUG
    if (dx == int1d{-1})
        throw MyException(0xbd067da0U, "dx should not be null.");
#endif
    int count = 0;
    Poly poly_x = Indices2Poly(x, ssS0_.basis.at(deg_x));

    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& nodes_ss = *all_basis_ss_[k];
        int t_max = all_t_max_[k];
        for (auto& [deg_y, basis_ss_d_original] : nodes_ss.front()) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const Staircase& basis_ss_d = GetRecentStaircase(nodes_ss, deg_y);
            if (deg_xy.t > t_max)
                break;
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                if (basis_ss_d.levels[i] > LEVEL_MAX - r_min)
                    break;
                const int r1 = LEVEL_MAX - basis_ss_d.levels[i];
                const int R = std::min(r, r1);
                if (R == r1 && basis_ss_d.diffs[i] == int1d{-1})
                    continue;
                AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                AdamsDeg deg_dxy = deg_xy + AdamsDeg(R, R - 1);

                if (bFastTry) {
                    auto& degs_changed = nodes_ss.back();
                    const auto degs_changed_end = degs_changed.end();
                    if (degs_changed.find(deg_y) == degs_changed_end && degs_changed.find(deg_xy) == degs_changed_end && degs_changed.find(deg_dxy) == degs_changed_end)
                        continue;
                }

                Poly poly_dx = (R == r && !dx.empty()) ? Indices2Poly(dx, ssS0_.basis.at(deg_dx)) : Poly();
                int1d xy, dxy;
                if (k == 0) {
                    auto& basis = ssS0_.basis;
                    auto& gb = ssS0_.gb;
                    Poly poly_y = Indices2Poly(basis_ss_d.basis[i], basis.at(deg_y));
                    Poly poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Poly2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    if (3 * deg_dxy.s <= deg_dxy.t + 3) {
                        if (deg_dxy.t > t_max)
                            dxy = int1d{-1};
                        else {
                            Poly poly_dy = (R == r1) ? Indices2Poly(basis_ss_d.diffs[i], basis.at(deg_y + AdamsDeg(R, R - 1))) : Poly();
                            Poly poly_dxy = gb.Reduce(poly_x * poly_dy + poly_dx * poly_y);
                            dxy = poly_dxy ? Poly2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                        }
                    }
                }
                else {
                    auto& basis = ssCofs_[k - 1].basis;
                    auto& gb = ssCofs_[k - 1].gb;
                    Mod poly_y = Indices2Mod(basis_ss_d.basis[i], basis.at(deg_y));
                    Mod poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    if (3 * deg_dxy.s <= deg_dxy.t + 3) {
                        if (deg_dxy.t > t_max)
                            dxy = int1d{-1};
                        else {
                            Mod poly_dy = (R == r1) ? Indices2Mod(basis_ss_d.diffs[i], basis.at(deg_y + AdamsDeg(R, R - 1))) : Mod();
                            Mod poly_dxy = gb.Reduce(poly_x * poly_dy + poly_dx * poly_y);
                            dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                        }
                    }
                }

                if (!xy.empty() || !dxy.empty()) {
                    SetDiffSc("S0", nodes_ss, deg_xy, xy, dxy, R);
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetCofDiffLeibniz(size_t iCof, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry)
{
#ifndef NDEBUG
    if (dx == int1d{-1})
        throw MyException(0x2f052b5aU, "dx should not be null.");
#endif
    int count = 0;

    auto& ssCof = ssCofs_[iCof];
    auto& nodes_ss = ssCof.nodes_ss;
    auto& basis = ssCof.basis;
    auto& gb = ssCof.gb;
    const int t_max = ssCof.t_max;
    Mod poly_x = Indices2Mod(x, basis.at(deg_x));

    for (auto& [deg_y, basis_ss_d_original] : ssS0_.nodes_ss.front()) {
        const Staircase& basis_ss_d = GetRecentStaircase(ssS0_.nodes_ss, deg_y);
        AdamsDeg deg_xy = deg_x + deg_y;
        if (deg_xy.t > t_max)
            break;
        for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
            if (basis_ss_d.levels[i] > LEVEL_MAX - r_min)
                break;
            int r1 = LEVEL_MAX - basis_ss_d.levels[i];
            int R = std::min(r, r1);
            if (R == r1 && basis_ss_d.diffs[i] == int1d{-1})
                continue;
            AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
            AdamsDeg deg_dxy = deg_xy + AdamsDeg(R, R - 1);

            if (bFastTry) {
                auto& degs_changed = nodes_ss.back();
                const auto degs_changed_end = degs_changed.end();
                if (degs_changed.find(deg_y) == degs_changed_end && degs_changed.find(deg_xy) == degs_changed_end && degs_changed.find(deg_dxy) == degs_changed_end)
                    continue;
            }

            Mod poly_dx = (R == r && !dx.empty()) ? Indices2Mod(dx, basis.at(deg_dx)) : Mod();

            Poly poly_y = Indices2Poly(basis_ss_d.basis[i], ssS0_.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

            int1d dxy;
            if (BelowS0VanishingLine(deg_dxy)) {
                if (deg_dxy.t > t_max)
                    dxy = int1d{-1};
                else {
                    Poly poly_dy = (R == r1) ? Indices2Poly(basis_ss_d.diffs[i], ssS0_.basis.at(deg_y + AdamsDeg(R, R - 1))) : Poly();
                    Mod poly_dxy = gb.Reduce(poly_dy * poly_x + poly_y * poly_dx);
                    dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                }
            }

            if (!xy.empty() || !dxy.empty()) {
                SetDiffSc(ssCof.name, nodes_ss, deg_xy, xy, dxy, R);
                ++count;
            }
        }
    }
    return count;
}

int Diagram::SetS0ImageLeibniz(AdamsDeg deg_x, const int1d& x, int r)
{
    int count = 0;

    const int r_original = r;
    r = NextRSrc(ssS0_.nodes_ss, deg_x, r);
    if (r == -1) {
        Logger::LogException(int(ssS0_.nodes_ss.size() - 2), 0x51274f1dU, "No source for the image. S0 deg_dx={}, dx={}, r={}\n", deg_x, x, r_original);
        throw SSException(0x51274f1dU, "No source for the image.");
    }

    Poly poly_x = Indices2Poly(x, ssS0_.basis.at(deg_x));
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& nodes_ss = *all_basis_ss_[k];
        int t_max = all_t_max_[k];

        for (auto& [deg_y, _] : nodes_ss.front()) {
            const Staircase& sc = GetRecentStaircase(nodes_ss, deg_y);
            AdamsDeg deg_xy = deg_x + deg_y;
            if (deg_xy.t > t_max)
                break;
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] >= LEVEL_MAX - r)
                    break;

                int1d xy;
                if (k == 0) {
                    auto& basis = ssS0_.basis;
                    auto& gb = ssS0_.gb;
                    Poly poly_y = Indices2Poly(sc.basis[i], basis.at(deg_y));
                    Poly poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Poly2Indices(poly_xy, basis.at(deg_xy)) : int1d();
                }
                else {
                    auto& basis = ssCofs_[k - 1].basis;
                    auto& gb = ssCofs_[k - 1].gb;
                    Mod poly_y = Indices2Mod(sc.basis[i], basis.at(deg_y));
                    Mod poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_xy)) : int1d();
                }

                if (!xy.empty()) {
                    SetImageSc("S0", nodes_ss, deg_xy, xy, {-1}, r);
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetCofImageLeibniz(size_t iCof, AdamsDeg deg_x, const int1d& x, int r)
{
    int count = 0;

    auto& ssCof = ssCofs_[iCof];
    auto& nodes_ss = ssCof.nodes_ss;
    auto& basis = ssCof.basis;
    auto& gb = ssCof.gb;
    int t_max = ssCof.t_max;

    const int r_original = r;
    r = NextRSrc(nodes_ss, deg_x, r);
    if (r == -1) {
        Logger::LogException(int(nodes_ss.size() - 2), 0xda298807U, "No source for the image. {} deg_dx={}, dx={}, r={}\n", ssCof.name, deg_x, x, r_original);
        throw SSException(0xda298807U, "No source for the image.");
    }

    Mod poly_x = Indices2Mod(x, basis.at(deg_x));

    for (auto& [deg_y, _] : ssS0_.nodes_ss.front()) {
        const Staircase& sc = GetRecentStaircase(ssS0_.nodes_ss, deg_y);
        AdamsDeg deg_xy = deg_x + deg_y;
        if (deg_xy.t > t_max)
            break;
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.levels[i] >= LEVEL_MAX - r)
                break;

            Poly poly_y = Indices2Poly(sc.basis[i], ssS0_.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_xy)) : int1d();

            if (!xy.empty()) {
                SetImageSc(ssCof.name, nodes_ss, deg_xy, xy, {-1}, r);
                ++count;
            }
        }
    }
    return count;
}
