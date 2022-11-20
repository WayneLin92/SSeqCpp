#include "algebras/linalg.h"
#include "main.h"

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

size_t GetFirstIndexOnLevel(const Staircase& sc, int level)
{
    return std::lower_bound(sc.levels.begin(), sc.levels.end(), level) - sc.levels.begin();
}

size_t GetFirstIndexOfNullOnLevel(const Staircase& sc, int level)
{
    int1d::const_iterator first = std::lower_bound(sc.levels.begin(), sc.levels.end(), level);
    int1d::const_iterator last = std::lower_bound(first, sc.levels.end(), level + 1);

    int1d::const_iterator it;
    ptrdiff_t count, step;
    count = std::distance(first, last);
    while (count > 0) {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if (*(sc.diffs_ind.begin() + (it - sc.levels.begin())) != int1d{-1}) {
            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    return first - sc.levels.begin();
}

/* Return -1 if not found */
int GetMaxLevelWithNull(const Staircase& sc)
{
    size_t i = sc.levels.size();
    while (i-- > 0) {
        if (sc.diffs_ind[i] == int1d{-1})
            return sc.levels[i];
    }
    return -1;
}

/* Return if x is in the vector space <level */
bool IsZeroOnLevel(const Staircase& sc, const int1d& x, int level)
{
    size_t first_l = GetFirstIndexOnLevel(sc, level);
    return lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_l, x).empty();
}

/* Add x, dx, level and triangularize.
 * Output the image of a differential that should be moved to the next level */
void triangularize(Staircase& sc, size_t i_insert, int1d x, int1d dx, int level, int1d& image, int& level_image)
{
    level_image = -1;

    size_t i = i_insert;
    while (!x.empty()) {
        std::swap(x, sc.basis_ind[i]);
        std::swap(dx, sc.diffs_ind[i]);
        std::swap(level, sc.levels[i]);

        ++i;
        for (size_t j = i_insert; j < i; ++j) {
            if (std::binary_search(x.begin(), x.end(), sc.basis_ind[j][0])) {
                x = lina::AddVectors(x, sc.basis_ind[j]);
                if (level == sc.levels[j] && dx != int1d{-1})
                    dx = lina::AddVectors(dx, sc.diffs_ind[j]);
            }
        }
    }
    if (dx != int1d{-1} && !dx.empty()) {
        image = std::move(dx);
        level_image = kLevelMax - level;
    }

    /* Triangularize the rest */
    for (; i < sc.basis_ind.size(); ++i) {
        for (size_t j = i_insert; j < i; ++j) {
            if (std::binary_search(sc.basis_ind[i].begin(), sc.basis_ind[i].end(), sc.basis_ind[j][0])) {
                sc.basis_ind[i] = lina::AddVectors(sc.basis_ind[i], sc.basis_ind[j]);
                if (sc.levels[i] == sc.levels[j] && sc.diffs_ind[i] != int1d{-1})
                    sc.diffs_ind[i] = lina::AddVectors(sc.diffs_ind[i], sc.diffs_ind[j]);
            }
        }
#ifndef NDEBUG
        if (sc.basis_ind[i].empty())
            throw MyException(0xfe35902dU, "BUG: triangularize()");
#endif
    }
}

Diagram::Diagram(const std::vector<std::string>& dbnames)
{
    {
        DBSS db(dbnames[0]);
        std::string complexName = "S0";
        std::string table_S0 = "S0_AdamsE2";

        ssS0_.basis = db.load_basis(table_S0);
        ssS0_.t_max = ssS0_.basis.rbegin()->first.t;
        ssS0_.basis_ss = {db.load_basis_ss(table_S0), {}};
        ssS0_.nd = {{}, {}};
        ssS0_.gb = Groebner(ssS0_.t_max, {}, db.load_gb(table_S0, DEG_MAX));
        ssS0_.pi_gen_Einf = db.get_column_from_str<Poly>(complexName + "_pi_generators", "Einf", "", myio::Deserialize<Poly>);
        ssS0_.pi_gb = algZ::Groebner(ssS0_.t_max, db.load_pi_gen_adamsdegs(complexName), db.load_pi_gb(complexName, DEG_MAX), true);

        all_basis_ss_.push_back(&ssS0_.basis_ss);
        all_nd_.push_back(&ssS0_.nd);
        all_t_max_.push_back(ssS0_.t_max);
    }

    for (size_t i = 1; i < dbnames.size(); ++i) {
        DBSS dbCof(dbnames[i]);
        std::string table_CW = GetE2TablePrefix(dbnames[i]);
        std::string complexName = GetComplexName(dbnames[i]);
        SSMod ssCof;
        ssCof.basis = dbCof.load_basis_mod(table_CW);
        ssCof.t_max = ssCof.basis.rbegin()->first.t;
        ssCof.basis_ss = {dbCof.load_basis_ss(table_CW), {}};
        ssCof.nd = {{}, {}};
        Mod1d xs = dbCof.load_gb_mod(table_CW, DEG_MAX);
        ssCof.gb = GroebnerMod(&ssS0_.gb, ssCof.t_max, {}, std::move(xs));
        ssCof.pi_gen_Einf = dbCof.get_column_from_str<Mod>(complexName + "_pi_generators", "Einf", "", myio::Deserialize<Mod>);
        ssCof.pi_gb = algZ::GroebnerMod(&ssS0_.pi_gb, ssCof.t_max, dbCof.load_pi_gen_adamsdegs(complexName), dbCof.load_pi_gb_mod(complexName, DEG_MAX), true);

        ssCof.f_top_cell = dbCof.get_column_from_str<Poly>(table_CW + "_generators", "to_S0", "", myio::Deserialize<Poly>);
        ssCof.deg_f_top_cell = AdamsDeg(0, GetTopCellT(dbnames[i]));
        ssCof.pi_f_top_cell = {dbCof.get_column_from_str<algZ::Poly>(complexName + "_pi_generators", "to_S0", "", myio::Deserialize<algZ::Poly>)};

        ssCofs_.push_back(std::move(ssCof));
    }
    for (size_t i = 0; i < ssCofs_.size(); ++i) {
        all_basis_ss_.push_back(&ssCofs_[i].basis_ss);
        all_nd_.push_back(&ssCofs_[i].nd);
        all_t_max_.push_back(ssCofs_[i].t_max);
    }
}

const Staircase& Diagram::GetRecentStaircase(const Staircases1d& basis_ss, AdamsDeg deg)
{
    for (auto p = basis_ss.rbegin(); p != basis_ss.rend(); ++p)
        if (p->find(deg) != p->end())
            return p->at(deg);
    throw MyException(0x553989e0U, "RecentStaircase not found. deg=" + deg.Str());
}
    /* Add a node */
void Diagram::AddNode()
{
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        all_basis_ss_[k]->push_back({});
        all_nd_[k]->push_back({});
    }

    ssS0_.pi_nodes_gen.push_back(ssS0_.pi_gb.gen_degs().size());
    ssS0_.pi_nodes_rel.push_back(ssS0_.pi_gb.data().size());
    ssS0_.pi_nodes_gen_2tor_degs.push_back(ssS0_.pi_gb.gen_2tor_degs());
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        ssCofs_[iCof].pi_nodes_gen.push_back(ssCofs_[iCof].pi_gb.v_degs().size());
        ssCofs_[iCof].pi_nodes_rel.push_back(ssCofs_[iCof].pi_gb.data().size());
        ssCofs_[iCof].pi_f_top_cell.push_back(ssCofs_[iCof].pi_f_top_cell.back());
    }
}

/* Pop the lastest node */
void Diagram::PopNode()
{
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        all_basis_ss_[k]->pop_back();
        all_nd_[k]->pop_back();
    }
    ssS0_.pi_gb.Pop(ssS0_.pi_nodes_gen.back(), ssS0_.pi_nodes_rel.back());
    ssS0_.pi_gb.set_gen_2tor_degs(std::move(ssS0_.pi_nodes_gen_2tor_degs.back()));
    ssS0_.pi_gen_Einf.resize(ssS0_.pi_nodes_gen.back());
    ssS0_.pi_nodes_gen.pop_back();
    ssS0_.pi_nodes_rel.pop_back();
    ssS0_.pi_nodes_gen_2tor_degs.pop_back();
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        ssCofs_[iCof].pi_gb.Pop(ssCofs_[iCof].pi_nodes_gen.back(), ssCofs_[iCof].pi_nodes_rel.back());
        ssCofs_[iCof].pi_gen_Einf.resize(ssCofs_[iCof].pi_nodes_gen.back());
        ssCofs_[iCof].pi_f_top_cell.pop_back();
        ssCofs_[iCof].pi_nodes_gen.pop_back();
        ssCofs_[iCof].pi_nodes_rel.pop_back();
    }
}

void Diagram::ApplyChanges(size_t index) // TODO: remove this function
{
    if (index == 0)
        throw MyException(0xc1b36735U, "The original basis_ss should not be changed");
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& basis_ss = *all_basis_ss_[k];
        for (auto p = basis_ss[index].begin(); p != basis_ss[index].end(); ++p) {
            auto& sc = GetRecentStaircase(basis_ss, p->first);
            if (&(p->second) != &sc)
                p->second = std::move(sc);
        }
        basis_ss.resize(index + 1);
        all_nd_[k]->reserve(index + 1);
    }
}

void Diagram::ApplyRecentChanges(std::vector<std::set<AdamsDeg>>& degs)  // TODO: remove this function
{
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& basis_ss = *all_basis_ss_[k];
        if (basis_ss.size() <= 2)
            throw MyException(0xc1b36735U, "ApplyRecentChanges() requires at least two new records");
        size_t index_before_last = basis_ss.size() - 2;
        degs[k].clear();
        for (auto p = basis_ss.back().begin(); p != basis_ss.back().end(); ++p) {
            basis_ss[index_before_last][p->first] = std::move(p->second);
            degs[k].insert(p->first);
        }
    }
    PopNode();
}

/**
 * Apply the change of the staircase to the current history
 */
void Diagram::UpdateStaircase(Staircases1d& basis_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, int1d x, int1d dx, int level, int1d& image, int& level_image)
{
    if (basis_ss.back().find(deg) == basis_ss.back().end())
        basis_ss.back()[deg] = sc_i;
    triangularize(basis_ss.back()[deg], i_insert, std::move(x), std::move(dx), level, image, level_image);
}

void Diagram::CacheNullDiffs(int maxPoss, int maxStem, bool bFull)
{
    for (size_t iSS = 0; iSS < all_basis_ss_.size(); ++iSS) {
        auto& basis_ss = *all_basis_ss_[iSS];
        int t_max = all_t_max_[iSS];
        all_nd_[iSS]->back().clear();
        for (auto& [deg, _] : basis_ss.front()) {
            if (deg.stem() > maxStem)
                continue;
            const Staircase& sc = GetRecentStaircase(basis_ss, deg);
            for (size_t i = 0; i < sc.diffs_ind.size(); ++i) {
                if (sc.diffs_ind[i] == int1d{-1}) {
                    NullDiff nd;
                    nd.deg = deg;
                    if (sc.levels[i] > kLevelPC) {
                        int r = kLevelMax - sc.levels[i];
                        AdamsDeg deg_tgt = deg + AdamsDeg{r, r - 1};
                        auto [index, count] = CountPossDrTgt(basis_ss, t_max, deg_tgt, r);
                        nd.r = r;
                        nd.first = index;
                        nd.count = count;
                        if (nd.count > maxPoss)
                            continue;
                    }
                    else if (sc.levels[i] < kLevelMax / 2) {
                        int r = sc.levels[i];
                        AdamsDeg deg_src = deg - AdamsDeg{r, r - 1};
                        auto [index, count] = CountPossDrSrc(basis_ss, deg_src, r);
                        nd.r = -r;
                        nd.first = index;
                        nd.count = count;
                        if (nd.count > maxPoss)
                            continue;
                    }
                    else
                        continue;

                    size_t j = i + 1;
                    while (j < sc.levels.size() && sc.diffs_ind[j] == int1d{-1} && sc.levels[i] == sc.levels[j])
                        ++j;
                    if (!bFull) {
                        for (size_t k = 0; k < j - i; ++k) {
                            nd.x = sc.basis_ind[i + k];
                            all_nd_[iSS]->back().push_back(nd);
                        }
                    }
                    else {
                        const unsigned k_max = unsigned(1) << (j - i);
                        for (unsigned k = 1; k < k_max; ++k) {
                            nd.x.clear();
                            for (int l : two_expansion(k))
                                nd.x = lina::AddVectors(nd.x, sc.basis_ind[i + l]);
                            all_nd_[iSS]->back().push_back(nd);
                        }
                    }
                    i = j - 1;
                }
            }
        }
        std::stable_sort(all_nd_[iSS]->back().begin(), all_nd_[iSS]->back().end(), [&](const NullDiff& nd1, const NullDiff& nd2) { return nd1.count < nd2.count; });
    }
}

bool Diagram::IsPossTgt(const Staircases1d& basis_ss, AdamsDeg deg, int r_max)
{
    r_max = std::min(r_max, deg.s - 1);
    for (int r1 = kLevelMin; r1 <= r_max; ++r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        if (basis_ss.front().find(d_src) != basis_ss.front().end())
            if (GetMaxLevelWithNull(GetRecentStaircase(basis_ss, d_src)) >= kLevelMax - r1)
                return true;
    }
    return false;
}

bool Diagram::IsPossSrc(const Staircases1d& basis_ss, int t_max, AdamsDeg deg, int r_min)
{
    int r_max = (deg.t - deg.s * 3 + 2) / 2;
    for (int r = r_min; r <= r_max; ++r) {
        AdamsDeg d_tgt = deg + AdamsDeg{r, r - 1};
        if (d_tgt.t > t_max)
            return true;
        if (basis_ss.front().find(d_tgt) != basis_ss.front().end()) {
            if (GetMaxLevelWithNull(GetRecentStaircase(basis_ss, d_tgt)) >= r)
                return true;
        }
    }
    return false;
}

int Diagram::GetFirstFixedLevelForPlot(const Staircases1d& basis_ss, AdamsDeg deg)
{
    auto& sc = GetRecentStaircase(basis_ss, deg);
    int result = kLevelMax - kLevelMin;
    for (size_t i = sc.levels.size(); i-- > 0 && sc.levels[i] >= kLevelPC;) {
        if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
            int r = kLevelMax - sc.levels[i];
            if (IsPossTgt(basis_ss, deg + AdamsDeg{r, r - 1}, r - 1))
                break;
            else
                result = sc.levels[i];
        }
    }
    return result;
}

size_t Diagram::GetFirstIndexOfFixedLevels(const Staircases1d& basis_ss, AdamsDeg deg, int level_min)
{
    auto& sc = GetRecentStaircase(basis_ss, deg);
    size_t result = sc.levels.size();
    for (size_t i = sc.levels.size(); i-- > 0;) {
        if (sc.diffs_ind[i] == int1d{-1} || sc.levels[i] < level_min)
            break;
        if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
            int r = kLevelMax - sc.levels[i];
            if (IsPossTgt(basis_ss, deg + AdamsDeg{r, r - 1}, r - 1))
                break;
            else
                result = i;
        }
    }
    return result;
}

std::pair<int, int> Diagram::CountPossDrTgt(const Staircases1d& basis_ss, int t_max, const AdamsDeg& deg_tgt, int r) const
{
    std::pair<int, int> result;
    if (basis_ss.front().find(deg_tgt) != basis_ss.front().end()) {
        const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
        result.first = (int)GetFirstIndexOnLevel(sc_tgt, r);
        result.second = (int)GetFirstIndexOfFixedLevels(basis_ss, deg_tgt, kLevelMax - r + 1) - result.first;
    }
    else if (deg_tgt.t > t_max)
        result = {-1, 10086};
    else
        result = {-1, 0};
    return result;
}

std::pair<int, int> Diagram::CountPossDrSrc(const Staircases1d& basis_ss, const AdamsDeg& deg_src, int r) const
{
    std::pair<int, int> result;
    if (basis_ss.front().find(deg_src) != basis_ss.front().end()) {
        const Staircase& sc_src = GetRecentStaircase(basis_ss, deg_src);
        result.first = (int)GetFirstIndexOnLevel(sc_src, kLevelMax - r);
        result.second = (int)GetFirstIndexOfFixedLevels(basis_ss, deg_src, kLevelMax - r + 1) - result.first;
    }
    else
        result = {-1, 0};
    return result;
}

int Diagram::NextRTgt(const Staircases1d& basis_ss, int t_max, AdamsDeg deg, int r) const
{
    AdamsDeg deg_tgt;
    int count = 0, index = -1;
    const Staircase& sc = GetRecentStaircase(basis_ss, deg);
    int r_max = (deg.t - deg.s * 3 + 2) / 2;
    for (int r1 = r; r1 <= r_max; ++r1) {
        if (deg.t + r1 - 1 > t_max)
            return r1;
        AdamsDeg d_tgt = deg + AdamsDeg{r1, r1 - 1};
        auto [index, count] = CountPossDrTgt(basis_ss, t_max, d_tgt, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

int Diagram::NextRSrc(const Staircases1d& basis_ss, AdamsDeg deg, int r) const
{
    AdamsDeg deg_src;
    int count = 0, index = -1;
    int r_max = std::min(r, deg.s - 1);
    const Staircase& sc = GetRecentStaircase(basis_ss, deg);
    for (int r1 = r_max; r1 >= kLevelMin; --r1) {
        AdamsDeg d_src = deg - AdamsDeg{r1, r1 - 1};
        auto [index, count] = CountPossDrSrc(basis_ss, d_src, r1);
        if (count > 0)
            return r1;
    }
    return -1;
}

/* Return d_r(x) */
int1d Diagram::GetDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int r) const
{
    if (x.empty())
        return int1d{};
    const Staircase& sc = GetRecentStaircase(basis_ss, deg_x);
    size_t first = GetFirstIndexOnLevel(sc, kLevelMax - r);
    size_t last = GetFirstIndexOfNullOnLevel(sc, kLevelMax - r);
    /* Compute x mod [0,first) */
    x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first, x);

    /* If x is in [first,last) */
    if (lina::Residue(sc.basis_ind.begin() + first, sc.basis_ind.begin() + last, x).empty())
        return lina::GetImage(sc.basis_ind.begin() + first, sc.basis_ind.begin() + last, sc.diffs_ind.begin() + first, sc.diffs_ind.begin() + last, x);
    else
        return int1d{-1};
}

bool Diagram::IsNewDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int1d dx, int r) const
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    int1d dx1 = GetDiff(basis_ss, deg_x, x, r);
    if (dx1 == int1d{-1})
        return true;
    if (basis_ss.front().find(deg_dx) != basis_ss.front().end()) {
        const Staircase& sc = GetRecentStaircase(basis_ss, deg_dx);
        size_t first_r = GetFirstIndexOnLevel(sc, r);
        /* Check if dx-dx1 is trivial mod Im(d_{r-1}) */
        return !lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, lina::AddVectors(dx, dx1)).empty();
    }
    return false;
}  // TODO: replace int1d{-1}

void Diagram::SetDiff(Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int1d dx, int r)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};

    /* If x is zero then dx is in Im(d_{r-1}) */
    if (x.empty()) {
        if (dx != int1d{-1} && !dx.empty())
            SetImage(basis_ss, deg_dx, std::move(dx), {-1}, r - 1);
        return;
    }

    const Staircase& sc = GetRecentStaircase(basis_ss, deg_x);
    size_t first_Nmr = GetFirstIndexOnLevel(sc, kLevelMax - r);
    x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_Nmr, x);
    if (x.empty()) {
        /* If x is in Ker(d_r) then dx is in Im(d_{r-1}) */
        if (dx != int1d{-1} && !dx.empty())
            SetImage(basis_ss, deg_dx, std::move(dx), {-1}, r - 1);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (dx == int1d{-1}) {
        /* If the target is unknown, insert it to the end of level N-r. */
        size_t first_Nmrp1 = GetFirstIndexOnLevel(sc, kLevelMax - r + 1);
        x = lina::Residue(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp1, x);
        if (!x.empty())
            UpdateStaircase(basis_ss, deg_x, sc, first_Nmrp1, x, {-1}, kLevelMax - r, image_new, level_image_new);
    }
    else if (dx.empty()) {
        /* If the target is zero, insert it to the end of level N-r-1 */
        UpdateStaircase(basis_ss, deg_x, sc, first_Nmr, x, {-1}, kLevelMax - r - 1, image_new, level_image_new);
    }
    else {
        /* Otherwise insert it to the beginning of level N-r */
        UpdateStaircase(basis_ss, deg_x, sc, first_Nmr, x, dx, kLevelMax - r, image_new, level_image_new);
    }

    if (level_image_new != -1) {
        if (level_image_new < kLevelMax / 2) {
            /* Add a d_{r1-1} image */
            AdamsDeg deg_image_new = deg_x + AdamsDeg{level_image_new, level_image_new - 1};
            SetImage(basis_ss, deg_image_new, std::move(image_new), {-1}, level_image_new - 1);
        }
        else {
            /* Add a d_{r1} cycle */
            int r_image = kLevelMax - level_image_new;
            AdamsDeg deg_image_new = deg_x - AdamsDeg{r_image, r_image - 1};
            SetDiff(basis_ss, deg_image_new, std::move(image_new), {}, r_image);
        }
    }

    /* Add image */
    if (dx != int1d{-1} && !dx.empty())
        SetImage(basis_ss, deg_dx, std::move(dx), std::move(x), r);
}

void Diagram::SetImage(Staircases1d& basis_ss, AdamsDeg deg_dx, int1d dx, int1d x, int r)
{
    AdamsDeg deg_x = deg_dx - AdamsDeg{r, r - 1};
    if (deg_x.s < 0)
        throw SSException(0x7dc5fa8cU, "7dc5fa8cU: No source for the image. deg_dx=" + deg_dx.StrAdams() + " r=" + std::to_string(r) + " dx=" + myio::Serialize(dx));

    /* If dx is in Im(d_{r-1}) then x is in Ker(d_r) */
    const Staircase& sc = GetRecentStaircase(basis_ss, deg_dx);
    size_t first_r = GetFirstIndexOnLevel(sc, r);
    dx = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, dx);
    if (dx.empty()) {
        if (x != int1d{-1} && !x.empty())
            SetDiff(basis_ss, deg_x, std::move(x), {-1}, r + 1);
        return;
    }

    int1d image_new;
    int level_image_new = -1;
    if (x == int1d{-1}) {
        /* If the source is unknown, check if it can be hit and then insert it to the end of level r. */
        size_t first_rp2 = GetFirstIndexOnLevel(sc, r + 1);
        dx = lina::Residue(sc.basis_ind.begin() + first_r, sc.basis_ind.begin() + first_rp2, dx);
        if (!dx.empty()) {
            if (!IsPossTgt(basis_ss, deg_dx, r))
                throw SSException(0x75989376U, "75989376U: No source for the image. deg_dx=" + deg_dx.StrAdams() + " r=" + std::to_string(r) + " dx=" + myio::Serialize(dx));
            UpdateStaircase(basis_ss, deg_dx, sc, first_rp2, dx, x, r, image_new, level_image_new);
        }
    }
    else {
        /* Otherwise insert it to the beginning of level r */
        UpdateStaircase(basis_ss, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new);  //
    }

    if (level_image_new != -1) {
        if (level_image_new < kLevelMax / 2) {
            /* Add a d_{r1-1} image */
            AdamsDeg deg_image_new = deg_dx + AdamsDeg{level_image_new, level_image_new - 1};
            SetImage(basis_ss, deg_image_new, std::move(image_new), {-1}, level_image_new - 1);
        }
        else {
            /* Add a d_r1 cycle */
            int r_image = kLevelMax - level_image_new;
            AdamsDeg deg_image_new = deg_dx - AdamsDeg{r_image, r_image - 1};
            SetDiff(basis_ss, deg_image_new, std::move(image_new), {-1}, r_image + 1);
        }
    }
}

int Diagram::SetS0DiffLeibniz(AdamsDeg deg_x, int1d x, int1d dx, int r, int r_min, bool bFastTry)
{
#ifndef NDEBUG
    if (dx == int1d{-1})
        throw MyException(0xbd067da0U, "dx should not be null.");
#endif
    int count = 0;
    Poly poly_x = Indices2Poly(x, ssS0_.basis.at(deg_x));

    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& basis_ss = *all_basis_ss_[k];
        int t_max = all_t_max_[k];
        for (auto& [deg_y, basis_ss_d_original] : basis_ss.front()) {
            AdamsDeg deg_xy = deg_x + deg_y;
            const Staircase& basis_ss_d = GetRecentStaircase(basis_ss, deg_y);
            if (deg_xy.t > t_max)
                break;
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                if (basis_ss_d.levels[i] > kLevelMax - r_min)
                    break;
                const int r1 = kLevelMax - basis_ss_d.levels[i];
                const int R = std::min(r, r1);
                if (R == r1 && basis_ss_d.diffs_ind[i] == int1d{-1})
                    continue;
                AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
                AdamsDeg deg_dxy = deg_xy + AdamsDeg(R, R - 1);

                if (bFastTry) {
                    auto& degs_changed = basis_ss.back();
                    const auto degs_changed_end = degs_changed.end();
                    if (degs_changed.find(deg_y) == degs_changed_end && degs_changed.find(deg_xy) == degs_changed_end && degs_changed.find(deg_dxy) == degs_changed_end)
                        continue;
                }

                Poly poly_dx = (R == r && !dx.empty()) ? Indices2Poly(dx, ssS0_.basis.at(deg_dx)) : Poly();
                int1d xy, dxy;
                if (k == 0) {
                    auto& basis = ssS0_.basis;
                    auto& gb = ssS0_.gb;
                    Poly poly_y = Indices2Poly(basis_ss_d.basis_ind[i], basis.at(deg_y));
                    Poly poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Poly2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    if (3 * deg_dxy.s <= deg_dxy.t + 3) {
                        if (deg_dxy.t > t_max)
                            dxy = int1d{-1};
                        else {
                            Poly poly_dy = (R == r1) ? Indices2Poly(basis_ss_d.diffs_ind[i], basis.at(deg_y + AdamsDeg(R, R - 1))) : Poly();
                            Poly poly_dxy = gb.Reduce(poly_x * poly_dy + poly_dx * poly_y);
                            dxy = poly_dxy ? Poly2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                        }
                    }
                }
                else {
                    auto& basis = ssCofs_[k - 1].basis;
                    auto& gb = ssCofs_[k - 1].gb;
                    Mod poly_y = Indices2Mod(basis_ss_d.basis_ind[i], basis.at(deg_y));
                    Mod poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

                    if (3 * deg_dxy.s <= deg_dxy.t + 3) {
                        if (deg_dxy.t > t_max)
                            dxy = int1d{-1};
                        else {
                            Mod poly_dy = (R == r1) ? Indices2Mod(basis_ss_d.diffs_ind[i], basis.at(deg_y + AdamsDeg(R, R - 1))) : Mod();
                            Mod poly_dxy = gb.Reduce(poly_x * poly_dy + poly_dx * poly_y);
                            dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                        }
                    }
                }

                if (!xy.empty() || !dxy.empty()) {
                    SetDiff(basis_ss, deg_xy, std::move(xy), std::move(dxy), R);
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetCofDiffLeibniz(size_t iCof, AdamsDeg deg_x, int1d x, int1d dx, int r, int r_min, bool bFastTry)
{
#ifndef NDEBUG
    if (dx == int1d{-1})
        throw MyException(0x2f052b5aU, "dx should not be null.");
#endif
    int count = 0;

    auto& basis_ss = ssCofs_[iCof].basis_ss;
    auto& basis = ssCofs_[iCof].basis;
    auto& gb = ssCofs_[iCof].gb;
    int t_max = ssCofs_[iCof].t_max;
    Mod poly_x = Indices2Mod(x, basis.at(deg_x));

    for (auto& [deg_y, basis_ss_d_original] : ssS0_.basis_ss.front()) {
        const Staircase& basis_ss_d = GetRecentStaircase(ssS0_.basis_ss, deg_y);
        AdamsDeg deg_xy = deg_x + deg_y;
        if (deg_xy.t > t_max)
            break;
        for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
            if (basis_ss_d.levels[i] > kLevelMax - r_min)
                break;
            int r1 = kLevelMax - basis_ss_d.levels[i];
            int R = std::min(r, r1);
            if (R == r1 && basis_ss_d.diffs_ind[i] == int1d{-1})
                continue;
            AdamsDeg deg_dx = deg_x + AdamsDeg(R, R - 1);
            AdamsDeg deg_dxy = deg_xy + AdamsDeg(R, R - 1);

            if (bFastTry) {
                auto& degs_changed = basis_ss.back();
                const auto degs_changed_end = degs_changed.end();
                if (degs_changed.find(deg_y) == degs_changed_end && degs_changed.find(deg_xy) == degs_changed_end && degs_changed.find(deg_dxy) == degs_changed_end)
                    continue;
            }

            Mod poly_dx = (R == r && !dx.empty()) ? Indices2Mod(dx, basis.at(deg_dx)) : Mod();

            Poly poly_y = Indices2Poly(basis_ss_d.basis_ind[i], ssS0_.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_x + deg_y)) : int1d{};

            int1d dxy;
            if (BelowS0VanishingLine(deg_dxy)) {
                if (deg_dxy.t > t_max)
                    dxy = int1d{-1};
                else {
                    Poly poly_dy = (R == r1) ? Indices2Poly(basis_ss_d.diffs_ind[i], ssS0_.basis.at(deg_y + AdamsDeg(R, R - 1))) : Poly();
                    Mod poly_dxy = gb.Reduce(poly_dy * poly_x + poly_y * poly_dx);
                    dxy = poly_dxy ? Mod2Indices(poly_dxy, basis.at(deg_dxy)) : int1d{};
                }
            }

            if (!xy.empty() || !dxy.empty()) {
                SetDiff(basis_ss, deg_xy, std::move(xy), std::move(dxy), R);
                ++count;
            }
        }
    }
    return count;
}

int Diagram::SetS0DiffLeibnizV2(AdamsDeg deg_x, int1d x, int1d dx, int r, bool bFastTry)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    int result = 0;

    auto& basis = ssS0_.basis;
    auto& basis_ss = ssS0_.basis_ss;
    int t_max = ssS0_.t_max;

    if (x.empty()) {
        if (!dx.empty() && !IsZeroOnLevel(GetRecentStaircase(basis_ss, deg_dx), dx, r))
            result += SetS0ImageLeibniz(deg_dx, dx, r - 1);
    }
    else if (IsNewDiff(basis_ss, deg_x, x, dx, r)) {
        int r_min = kLevelMin;
        while (r_min < r && !IsNewDiff(basis_ss, deg_x, x, {}, r_min))  // TODO: improve this
            ++r_min;
        if (dx.empty()) {
            int r_max = NextRTgt(basis_ss, t_max, deg_x, r + 1);
            if (r_max == -1)
                r = kRPC - 1;
            else
                r = r_max - 1;
        }
        result += SetS0DiffLeibniz(deg_x, x, dx, r, r_min, bFastTry);
    }
    return result;
}

int Diagram::SetCofDiffLeibnizV2(size_t iCof, AdamsDeg deg_x, int1d x, int1d dx, int r, bool bFastTry)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int count = 0;

    auto& basis = ssCofs_[iCof].basis;
    auto& basis_ss = ssCofs_[iCof].basis_ss;
    int t_max = ssCofs_[iCof].t_max;

    if (x.empty()) {
        if (!dx.empty() && !IsZeroOnLevel(GetRecentStaircase(basis_ss, deg_dx), dx, r)) {
            count += SetCofImageLeibniz(iCof, deg_dx, dx, r - 1);

            Mod mod_dx = Indices2Mod(dx, basis.at(deg_dx));
            Poly poly_dx_S0 = ssS0_.gb.Reduce(subsMod(mod_dx, ssCofs_[iCof].f_top_cell));
            if (poly_dx_S0) {
                AdamsDeg deg_dx_S0 = deg_dx - ssCofs_[iCof].deg_f_top_cell;
                int1d dx_S0 = Poly2Indices(poly_dx_S0, ssS0_.basis.at(deg_dx_S0));
                count += SetS0ImageLeibniz(deg_dx_S0, dx_S0, r - 1);
            }
        }
    }
    else if (IsNewDiff(basis_ss, deg_x, x, dx, r)) {
        int r_min = kLevelMin;
        while (r_min < r && !IsNewDiff(basis_ss, deg_x, x, {}, r_min))  // TODO: improve this
            ++r_min;
        if (dx.empty()) {
            int r_max = NextRTgt(basis_ss, t_max, deg_x, r + 1);
            if (r_max == -1)
                r = kRPC - 1;
            else
                r = r_max - 1;
        }
        count += SetCofDiffLeibniz(iCof, deg_x, x, dx, r, r_min, bFastTry);

        Mod mod_x = Indices2Mod(x, basis.at(deg_x));
        Poly poly_x_S0 = ssS0_.gb.Reduce(subsMod(mod_x, ssCofs_[iCof].f_top_cell));
        Mod mod_dx = !dx.empty() ? Indices2Mod(dx, basis.at(deg_dx)) : Mod();
        Poly poly_dx_S0 = ssS0_.gb.Reduce(subsMod(mod_dx, ssCofs_[iCof].f_top_cell));
        if (poly_x_S0 || poly_dx_S0) {
            AdamsDeg deg_x_S0 = deg_x - ssCofs_[iCof].deg_f_top_cell;
            int1d x_S0 = poly_x_S0 ? Poly2Indices(poly_x_S0, ssS0_.basis.at(deg_x_S0)) : int1d{};
            AdamsDeg deg_dx_S0 = deg_dx - ssCofs_[iCof].deg_f_top_cell;
            int1d dx_S0 = poly_dx_S0 ? Poly2Indices(poly_dx_S0, ssS0_.basis.at(deg_dx_S0)) : int1d{};
            count += SetS0DiffLeibnizV2(deg_x_S0, x_S0, dx_S0, r, bFastTry);
        }
    }
    return count;
}

int Diagram::SetDiffLeibnizV2(size_t index, AdamsDeg deg_x, int1d x, int1d dx, int r, bool bFastTry)
{
    if (index == 0)
        return SetS0DiffLeibnizV2(deg_x, x, dx, r, bFastTry);
    else
        return SetCofDiffLeibnizV2(index - 1, deg_x, x, dx, r, bFastTry);
}

int Diagram::SetS0ImageLeibniz(AdamsDeg deg_x, int1d x, int r)
{
    int count = 0;

    r = NextRSrc(ssS0_.basis_ss, deg_x, r);
    if (r == -1)
        throw SSException(0x51274f1dU, "bef9931bU: No source for the image. deg_dx=" + deg_x.StrAdams() + " r=" + std::to_string(r) + " dx=" + myio::Serialize(x));

    Poly poly_x = Indices2Poly(x, ssS0_.basis.at(deg_x));
    for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
        auto& basis_ss = *all_basis_ss_[k];
        int t_max = all_t_max_[k];

        for (auto& [deg_y, _] : basis_ss.front()) {
            const Staircase& sc = GetRecentStaircase(basis_ss, deg_y);
            AdamsDeg deg_xy = deg_x + deg_y;
            if (deg_xy.t > t_max)
                break;
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] >= kLevelMax - r)
                    break;

                int1d xy;
                if (k == 0) {
                    auto& basis = ssS0_.basis;
                    auto& gb = ssS0_.gb;
                    Poly poly_y = Indices2Poly(sc.basis_ind[i], basis.at(deg_y));
                    Poly poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Poly2Indices(poly_xy, basis.at(deg_xy)) : int1d();
                }
                else {
                    auto& basis = ssCofs_[k - 1].basis;
                    auto& gb = ssCofs_[k - 1].gb;
                    Mod poly_y = Indices2Mod(sc.basis_ind[i], basis.at(deg_y));
                    Mod poly_xy = gb.Reduce(poly_x * poly_y);
                    xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_xy)) : int1d();
                }

                if (!xy.empty()) {
                    SetImage(basis_ss, deg_xy, std::move(xy), {-1}, r);
                    ++count;
                }
            }
        }
    }
    return count;
}

int Diagram::SetCofImageLeibniz(size_t iCof, AdamsDeg deg_x, int1d x, int r)
{
    int count = 0;

    auto& basis_ss = ssCofs_[iCof].basis_ss;
    auto& basis = ssCofs_[iCof].basis;
    auto& gb = ssCofs_[iCof].gb;
    int t_max = ssCofs_[iCof].t_max;

    r = NextRSrc(basis_ss, deg_x, r);
    if (r == -1)
        throw SSException(0xda298807U, "bef9931bU: No source for the image. deg_dx=" + deg_x.StrAdams() + " dx=" + myio::Serialize(x));

    Mod poly_x = Indices2Mod(x, basis.at(deg_x));

    for (auto& [deg_y, _] : ssS0_.basis_ss.front()) {
        const Staircase& sc = GetRecentStaircase(ssS0_.basis_ss, deg_y);
        AdamsDeg deg_xy = deg_x + deg_y;
        if (deg_xy.t > t_max)
            break;
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.levels[i] >= kLevelMax - r)
                break;

            Poly poly_y = Indices2Poly(sc.basis_ind[i], ssS0_.basis.at(deg_y));
            Mod poly_xy = gb.Reduce(poly_y * poly_x);
            int1d xy = poly_xy ? Mod2Indices(poly_xy, basis.at(deg_xy)) : int1d();

            if (!xy.empty()) {
                SetImage(basis_ss, deg_xy, std::move(xy), {-1}, r);
                ++count;
            }
        }
    }
    return count;
}
