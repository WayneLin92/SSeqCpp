#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

Poly Proj(const algZ::Mon& mon, const Poly1d& map)
{
    Poly tmp;
    Poly result = pow(map[0], mon.c());
    result.imulP(subs(mon.m(), map), tmp);
    return result;
}

Mod Proj(const algZ::MMod& mon, const Poly1d& map, const Mod1d& map_v)
{
    return Proj(mon.m, map) * map_v[mon.v];
}

/* This version checks if W is a subspace of V */
int2d QuotientSpace(const int2d& spaceV, const int2d& spaceW)
{
    int2d quotient;
    size_t dimQuo = spaceV.size() - spaceW.size();
    for (size_t i = 0; i < spaceV.size(); i++) {
        auto v1 = lina::Residue(quotient, lina::Residue(spaceW, spaceV[i]));
        if (!v1.empty())
            quotient.push_back(std::move(v1));
    }
    if (quotient.size() != dimQuo)
        throw MyException(0x72715463U, "W is not a subspace of V!\n");
    return quotient;
}

/* This is for use in GetKernel */
int NextSExtV2(const ut::map_seq2d<int, 0>& basis_ss_possEinf, int t_max, int stem, int s_min, const ut::map_seq<int, 0>& num_leads)
{
    int s_max = (stem + 3) / 2;
    for (int s = s_min; s <= s_max; ++s) {
        if (stem + s > t_max || basis_ss_possEinf(stem, s) > num_leads[s])
            return s;
    }
    return algZ::FIL_MAX + 1;
}

template <typename T>
int ExtendRel(const ut::map_seq2d<int, 0>& basis_ss_possEinf, int stem, int t_max, const T& rel, T& rel_extended)
{
    auto& m = rel.data.back();
    if (m.IsUnKnown()) {
        int s = algZ::NextO(basis_ss_possEinf, t_max, stem, m.fil());
        if (s > m.fil()) {
            rel_extended = rel;
            if (s <= algZ::FIL_MAX)
                rel_extended.data.back() = T::O(s);
            else
                rel_extended.data.pop_back();
            return 1;
        }
    }
    return 0;
}

template <typename T>
int ExtendRelV2(const ut::map_seq2d<int, 0>& basis_ss_possEinf, int stem, int t_max, const T& rel, T& rel_extended, ut::map_seq<int, 0>& num_leads)
{
    auto& m = rel.data.back();
    if (m.IsUnKnown()) {
        int s = NextSExtV2(basis_ss_possEinf, t_max, stem, m.fil(), num_leads);
        if (s > m.fil()) {
            rel_extended = rel;
            if (s <= algZ::FIL_MAX)
                rel_extended.data.back() = T::O(s);
            else
                rel_extended.data.pop_back();
            return 1;
        }
    }
    return 0;
}

/* Get a Z2-basis of the kernel.
 * fx should all be certain.
 */
template <typename T1, typename T2, typename GB1, typename GB2>
void GetKernel(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, std::vector<T1>& kernel)
{
    /* Sort by certainty */
    auto indices = ut::size_t_range(x.size());
    std::stable_sort(indices.begin(), indices.end(), [&x](size_t i, size_t j) { return x[i].UnknownFil() > x[j].UnknownFil(); });

    std::vector<T2> image;
    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    /* f(g[i]) = image[i] */
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (tgt.data[index] == image[j].GetLead()) {
                        int c = tgt.data[index].c() - image[j].GetLead().c();
                        tgt.isubmulP(algZ::Mon::twoTo(c), image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.Reduce(std::move(tgt));
                        src.isubmulP(algZ::Mon::twoTo(c), g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.Reduce(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
        else {
            size_t index = 0;
            while (index < src.data.size() && !src.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < kernel.size(); j++) {
                    if (src.data[index] == kernel[j].GetLead()) {
                        int c = src.data[index].c() - kernel[j].GetLead().c();
                        src.isubmulP(algZ::Mon::twoTo(c), kernel[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.Reduce(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
            if (algZ::IsValidRel(src))
                kernel.push_back(std::move(src));
        }
    }
}

/* Get a Z2-basis of the image */
template <typename T1, typename T2, typename GB1, typename GB2>
void GetImage(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, std::vector<T2>& image, int& O1, int& O2, T1& gO)
{
    /* Sort by certainty */
    auto indices = ut::size_t_range(fx.size());
    std::stable_sort(indices.begin(), indices.end(), [&fx](size_t i, size_t j) { return fx[i].UnknownFil() > fx[j].UnknownFil(); });

    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    /* f(g[i]) = image[i] */
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size() && !tgt.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (tgt.data[index] == image[j].GetLead()) {
                        int c = tgt.data[index].c() - image[j].GetLead().c();
                        tgt.isubmulP(algZ::Mon::twoTo(c), image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.ReduceV2(std::move(tgt));
                        src.isubmulP(algZ::Mon::twoTo(c), g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.ReduceV2(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            if (!tgt.GetLead().IsUnKnown()) {
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
            else {
                int O = tgt.UnknownFil();
                if (O < O1) {
                    O2 = O1;
                    O1 = O;
                    gO = std::move(src);
                }
                else if (O == O1) {
                    O2 = O1;
                    gO = T1::O(0);
                }
                else if (O < O2)
                    O2 = O;
            }
        }
    }
}

/* Group the leads by s and return the numbers of leads for each s */
template <typename T1, typename T2, typename GB1, typename GB2>
void GetImageLeads(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, ut::map_seq<int, 0>& num_leads, ut::map_seq<int, algZ::FIL_MAX + 1>& src_Os)
{
    /* Sort by certainty */
    auto indices = ut::size_t_range(fx.size());
    std::stable_sort(indices.begin(), indices.end(), [&fx](size_t i, size_t j) { return fx[i].UnknownFil() > fx[j].UnknownFil(); });

    /* f(g[i]) = image[i] */
    std::vector<T2> image;
    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size() && !tgt.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (tgt.data[index] == image[j].GetLead()) {
                        tgt.isubP(image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.ReduceV2(std::move(tgt));
                        src.isubP(g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.ReduceV2(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            if (!tgt.GetLead().IsUnKnown()) {
                ++num_leads[tgt.GetLead().fil()];
                src_Os[tgt.GetLead().fil()] = std::min(src_Os[tgt.GetLead().fil()], src.GetLead().fil());
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
        }
    }
}

template <typename T, typename GB>
T Residue(const std::vector<T>& space, T x, const GB& gb)
{
    T tmp{};
    x = gb.ReduceV2(std::move(x));
    size_t index = 0;
    while (index < x.data.size() && !x.data[index].IsUnKnown()) {
        size_t d_index = 1;
        for (size_t j = 0; j < space.size(); j++) {
            if (x.data[index] == space[j].GetLead()) {
                int c = x.data[index].c() - space[j].GetLead().c();
                x.isubmulP(algZ::Mon::twoTo(c), space[j], tmp, gb.gen_2tor_degs());
                x = gb.ReduceV2(std::move(x));
                d_index = 0;
                break;
            }
        }
        index += d_index;
    }
    return x;
}

template <typename T>
void CountByS(std::vector<T>& cont, ut::map_seq<int, 0>& num_leads)
{
    for (auto& x : cont)
        ++num_leads[x.GetLead().fil()];
}

const PiBase* Diagram::GetRecentPiBasis(const PiBasis1d& nodes_pi_basis, AdamsDeg deg)
{
    for (auto p = nodes_pi_basis.rbegin(); p != nodes_pi_basis.rend(); ++p)
        if (p->find(deg) != p->end())
            return &p->at(deg);
    return nullptr;
}

const PiBaseMod* Diagram::GetRecentPiBasis(const PiBasisMod1d& nodes_pi_basis, AdamsDeg deg)
{
    for (auto p = nodes_pi_basis.rbegin(); p != nodes_pi_basis.rend(); ++p)
        if (p->find(deg) != p->end())
            return &p->at(deg);
    return nullptr;
}

int Diagram::PossEinf(const Staircases1d& nodes_ss, AdamsDeg deg)
{
    if (ut::has(nodes_ss.front(), deg)) {
        const Staircase& sc = ut::GetRecentValue(nodes_ss, deg);
        size_t i_start_perm = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
        size_t i_stable = GetFirstIndexOfFixedLevels(nodes_ss, deg, LEVEL_PERM + 1);
        return int(i_stable - i_start_perm);
    }
    return 0;
}

void Diagram::UpdatePossEinf(const Staircases1d& nodes_ss, ut::map_seq2d<int, 0>& basis_ss_possEinf)
{
    for (auto& [deg, _] : nodes_ss.front())
        basis_ss_possEinf(deg.stem(), deg.s) = PossEinf(nodes_ss, deg);
}

int Diagram::PossMoreEinf(const Staircases1d& nodes_ss, AdamsDeg deg)  //// TODO: improve
{
    if (ut::has(nodes_ss.front(), deg)) {
        const Staircase& sc = ut::GetRecentValue(nodes_ss, deg);
        size_t i_end_perm = GetFirstIndexOnLevel(sc, LEVEL_PERM + 1);
        size_t i_stable = GetFirstIndexOfFixedLevels(nodes_ss, deg, LEVEL_PERM + 1);
        return int(i_stable - i_end_perm);
    }
    return 0;
}

void Diagram::PossMoreEinfFirstS_Ring(size_t iRing, int1d& O1s, int1d& O2s, int1d& isSingle) const
{
    auto& ring = rings_[iRing];
    for (int i = 0; i <= ring.t_max; ++i) {
        O1s.push_back(ring.t_max - i + 1);
        O2s.push_back(ring.t_max - i + 1);
        isSingle.push_back(0);
    }
    for (auto& [deg, _] : ring.nodes_ss.front()) {
        if (int num = PossMoreEinf(ring.nodes_ss, deg)) {
            if (deg.s < O1s[deg.stem()]) {
                O2s[deg.stem()] = O1s[deg.stem()];
                O1s[deg.stem()] = deg.s;
                if (num == 1)
                    isSingle[deg.stem()] = 1;
            }
            else if (deg.s == O1s[deg.stem()]) {
                O2s[deg.stem()] = O1s[deg.stem()];
                isSingle[deg.stem()] = 0;
            }
            else if (deg.s < O2s[deg.stem()])
                O2s[deg.stem()] = deg.s;
        }
    }
}

void Diagram::PossMoreEinfFirstS_Mod(size_t iMod, int1d& O1s, int1d& O2s, int1d& isSingle) const
{
    auto& Mod = modules_[iMod];
    for (int i = 0; i <= Mod.t_max; ++i) {
        O1s.push_back(Mod.t_max - i + 1);
        O2s.push_back(Mod.t_max - i + 1);
        isSingle.push_back(0);
    }
    for (auto& [deg, _] : Mod.nodes_ss.front()) {
        if (int num = PossMoreEinf(Mod.nodes_ss, deg)) {
            if (deg.s < O1s[deg.stem()]) {
                O2s[deg.stem()] = O1s[deg.stem()];
                O1s[deg.stem()] = deg.s;
                if (num == 1)
                    isSingle[deg.stem()] = 1;
            }
            else if (deg.s == O1s[deg.stem()]) {
                O2s[deg.stem()] = O1s[deg.stem()];
                isSingle[deg.stem()] = 0;
            }
            else if (deg.s < O2s[deg.stem()])
                O2s[deg.stem()] = deg.s;
        }
    }
}

int Diagram::ExtendRelRing(size_t iRing, int stem, const algZ::Poly& rel, algZ::Poly& rel_extended) const
{
    return ExtendRel(rings_[iRing].basis_ss_possEinf, stem, rings_[iRing].t_max, rel, rel_extended);
}

int Diagram::ExtendRelMod(size_t iMod, int stem, const algZ::Mod& rel, algZ::Mod& rel_extended) const
{
    return ExtendRel(modules_[iMod].basis_ss_possEinf, stem, modules_[iMod].t_max, rel, rel_extended);
}

int Diagram::ExtendRelRing(size_t iRing, int stem, algZ::Poly& rel) const
{
    if (rel) {
        algZ::Poly rel1;
        if (ExtendRelRing(iRing, stem, rel, rel1)) {
            rel = std::move(rel1);
            return 1;
        }
    }
    return 0;
}

int Diagram::ExtendRelMod(size_t iMod, int stem, algZ::Mod& rel) const
{
    if (rel) {
        algZ::Mod rel1;
        if (ExtendRelMod(iMod, stem, rel, rel1)) {
            rel = std::move(rel1);
            return 1;
        }
    }
    return 0;
}

int Diagram::ExtendRelRingV2(size_t iRing, int stem, algZ::Poly& rel, ut::map_seq<int, 0>& num_leads) const
{
    if (rel) {
        algZ::Poly rel1;
        if (ExtendRelV2(rings_[iRing].basis_ss_possEinf, stem, rings_[iRing].t_max, rel, rel1, num_leads)) {
            rel = std::move(rel1);
            return 1;
        }
    }
    return 0;
}

int Diagram::ExtendRelCofV2(size_t iCof, int stem, algZ::Mod& rel, ut::map_seq<int, 0>& num_leads) const
{
    if (rel) {
        algZ::Mod rel1;
        if (ExtendRelV2(modules_[iCof].basis_ss_possEinf, stem, modules_[iCof].t_max, rel, rel1, num_leads)) {
            rel = std::move(rel1);
            return 1;
        }
    }
    return 0;
}

int2d Diagram::GetRingGbEinf(size_t iRing, AdamsDeg deg) const
{
    auto& ring = rings_[iRing];
    auto& gb = ring.gb;
    auto& pi_gb = ring.pi_gb;
    auto& pi_gen_Einf = ring.pi_gen_Einf;

    int2d result;
    Poly tmp;
    algZ::Poly1d gb_Einf = pi_gb.RelsLF(deg);
    for (auto& f_Einf : gb_Einf) {
        Poly LF;
        for (auto& m : f_Einf.data)
            LF.iaddP(Proj(m, pi_gen_Einf), tmp);
        LF = gb.Reduce(std::move(LF));
        if (LF)
            result.push_back(Poly2Indices(LF, ring.basis.at(deg)));
        else
            result.push_back(int1d{});
    }

    return result;
}

int2d Diagram::GetModuleGbEinf(size_t iMod, AdamsDeg deg) const
{
    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& gb = mod.gb;
    auto& pi_gb = mod.pi_gb;
    auto& pi_gen_Einf = mod.pi_gen_Einf;

    int2d result;
    Mod tmp;
    algZ::Mod1d gb_Einf = pi_gb.RelsLF(deg);
    for (auto& f_Einf : gb_Einf) {
        Mod LF;
        for (auto& m : f_Einf.data)
            LF.iaddP(Proj(m, ring.pi_gen_Einf, pi_gen_Einf), tmp);
        LF = gb.Reduce(std::move(LF));
        if (LF)
            result.push_back(Mod2Indices(LF, modules_[iMod].basis.at(deg)));
        else
            result.push_back(int1d{});
    }

    return result;
}

std::map<AdamsDeg, int2d> Diagram::GetRingGbEinf(size_t iRing) const
{
    std::map<AdamsDeg, int2d> result;
    for (auto& [deg, _] : rings_[iRing].pi_gb.leads_group_by_deg())
        result[deg] = GetRingGbEinf(iRing, deg);
    return result;
}

std::map<AdamsDeg, int2d> Diagram::GetModuleGbEinf(size_t iMod) const
{
    std::map<AdamsDeg, int2d> result;
    for (auto& [deg, _] : modules_[iMod].pi_gb.leads_group_by_deg())
        result[deg] = GetModuleGbEinf(iMod, deg);
    return result;
}

void Diagram::SetPermanentCycle(int depth, size_t iMod, AdamsDeg deg_x) //// TODO: correct this
{
    /*auto& nodes_ss = modules_[iMod].nodes_ss;
    if (nodes_ss.front().find(deg_x) != nodes_ss.front().end()) {
        const Staircase& sc = ut::GetRecentValue(nodes_ss, deg_x);
        size_t i_end_perm = GetFirstIndexOnLevel(sc, LEVEL_PERM + 1);
        size_t i_stable = GetFirstIndexOfFixedLevels(nodes_ss, deg_x, LEVEL_PERM + 1);
        if (i_stable - i_end_perm == 1) {
            int r = LEVEL_MAX - LEVEL_PERM;
            if (i_end_perm > 0)
                r = std::min(r, LEVEL_MAX - sc.levels[i_end_perm - 1]);
            if (r != LEVEL_MAX - sc.levels[i_end_perm]) {
                int1d x = {sc.basis[i_end_perm]};
                Logger::LogDiff(depth, enumReason::exact_hq, modules_[iMod].name, deg_x, x, {}, r - 1);
                SetModuleDiffGlobal(iMod, deg_x, x, {}, r - 1, true);
            }
        }
        else if (i_stable - i_end_perm > 1) {
            Logger::LogException(0, 0xc625fffU, "More than one possible permanent cycles. {} deg={}.\n", modules_[iMod].name, deg_x);
            throw MyException(0xc625fffU, "More than one possible permanent cycles");
        }
    }
    else {
        Logger::LogException(0, 0xccd5e2c4U, "Permanent cycle not found. {} deg={}.\n", modules_[iMod].name, deg_x);
        throw MyException(0xccd5e2c4U, "Permanent cycle not found");
    }*/
}

void Diagram::AddPiRelsRing(size_t iRing, algZ::Poly1d rels)
{
    auto& ring = rings_[iRing];
    ring.pi_gb.AddRels(std::move(rels), ring.t_max, ring.basis_ss_possEinf);
    for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
        modules_[iCof].pi_gb.AddRels({}, modules_[iCof].t_max, modules_[iCof].basis_ss_possEinf);
    ////TODO: add pi_maps
}

void Diagram::AddPiRelsCof(size_t iMod, algZ::Mod1d rels)
{
    auto& mod = modules_[iMod];
    mod.pi_gb.AddRels(std::move(rels), mod.t_max, mod.basis_ss_possEinf);
    ////TODO: add pi_maps
    /*algZ::Poly1d relsS0;
    for (auto& rel : rels) {
        auto relS0 = algZ::subs(rel, mod.nodes_pi_qt.back(), mod.pi_gb.v_degs());
        if (algZ::IsValidRel(relS0))
            relsS0.push_back(std::move(relS0));
    }
    AddPiRelsRing(std::move(relsS0));*/
}

void Diagram::AddPiRelsByNat(size_t iMap)
{
    /*algZ::Poly1d relsS0;
    for (auto& rel : modules_[iCof].pi_gb.data()) {
        auto relS0 = algZ::subs(rel, modules_[iCof].nodes_pi_qt.back(), modules_[iCof].pi_gb.v_degs());
        if (algZ::IsValidRel(relS0)) {
            relsS0.push_back(std::move(relS0));
        }
    }
    AddPiRelsRing(std::move(relsS0));*/
}

algZ::Mon1d Diagram::GenBasis(const algZ::Groebner& gb, AdamsDeg deg, const PiBasis1d& nodes_pi_basis)
{
    algZ::Mon1d result;
    for (size_t gen_id = 0; gen_id < gb.gen_degs().size(); ++gen_id) {
        AdamsDeg d1 = deg - gb.gen_degs()[gen_id];
        if (auto p = GetRecentPiBasis(nodes_pi_basis, d1)) {
            for (auto& m : p->nodes_pi_basis) {
                if (!m || (int)gen_id >= m.backg()) {
                    algZ::Mon mon = algZ::mul_unsigned(m, gb.Gen((uint32_t)gen_id));
                    if (gb.IsNewBaseByLastGen(mon, (uint32_t)gen_id))
                        result.push_back(std::move(mon));
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

algZ::MMod1d Diagram::GenBasis(const algZ::GroebnerMod& gb, AdamsDeg deg, const PiBasis1d& nodes_pi_basis)
{
    algZ::MMod1d result;
    auto& v_degs = gb.v_degs();
    for (size_t v = 0; v < v_degs.size(); ++v) {
        AdamsDeg d1 = deg - v_degs[v];
        if (auto p = GetRecentPiBasis(nodes_pi_basis, d1)) {
            for (auto& m : p->nodes_pi_basis) {
                if (gb.IsNewBaseByV(m, (uint32_t)v))
                    result.emplace_back(m, (int)v, v_degs[v].s);
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

void Diagram::SyncS0Homotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
{
    // int t_max = rings_.t_max;
    // auto& gb = rings_.gb;
    // auto& basis = rings_.basis;
    // auto& nodes_ss = rings_.nodes_ss;
    // auto& pi_gb = rings_.pi_gb;
    // auto& nodes_pi_basis = rings_.nodes_pi_basis;
    // auto& pi_gen_Einf = rings_.pi_gen_Einf;

    // int count_ss_old = count_ss;

    // Poly tmp;
    // for (int t = 1; t <= t_max; ++t) {
    //     for (int s = 0; s <= t; ++s) {
    //         AdamsDeg deg(s, t);
    //         if (deg.stem() < deg_min.stem() || deg.s < deg_min.s)
    //             continue;
    //         algZ::Mon1d pi_basis_d = GenBasis(pi_gb, deg, nodes_pi_basis);
    //         if (basis.find(deg) != basis.end()) {
    //             /* Add new boundaries to ss */
    //             int2d boundaries = GetRingGbEinf(deg);
    //             const int r = deg.s + 1;
    //             const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
    //             for (int1d& boundary : boundaries) {
    //                 if (!boundaries.empty()) {
    //                     const int count1 = SetRingDiffGlobal(deg_src, {}, boundary, r);
    //                     if (count1 > 0) {
    //                         Logger::LogDiffBoun(depth, enumReason::htpy2ss, "S0", deg, boundary);
    //                         count_ss += count1;
    //                     }
    //                 }
    //             }

    //            /* Construct the projection map */
    //            auto& basis_d = basis.at(deg);
    //            int2d projs;
    //            {
    //                auto& sc = ut::GetRecentValue(nodes_ss, deg);
    //                size_t first_PC = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
    //                for (auto& m : pi_basis_d) {
    //                    int1d proj = Poly2Indices(gb.Reduce(Proj(m, pi_gen_Einf)), basis_d);
    //                    proj = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_PC, proj);
    //                    projs.push_back(proj);
    //                }
    //            }
    //            int2d image, kernel, g;
    //            lina::SetLinearMap(projs, image, kernel, g);

    //            /* Compute pi_basis */
    //            int1d non_leads = lina::add(ut::int_range((int)pi_basis_d.size()), lina::GetLeads(kernel));
    //            algZ::Mon1d pi_basis_d_v2;
    //            int2d pi_basis_d_Einf;
    //            for (int i : non_leads) {
    //                pi_basis_d_v2.push_back(pi_basis_d[i]);
    //                pi_basis_d_Einf.push_back(projs[i]);
    //            }

    //            /* Add new relations in homotopy */
    //            algZ::Poly1d pi_rels;
    //            for (auto& k : kernel) {
    //                pi_rels.push_back(algZ::Indices2Poly(k, pi_basis_d) + algZ::Mon::O(deg.s + 1));
    //                Logger::LogHtpyRel(depth, enumReason::ss2htpy, "S0", deg, pi_rels.back());
    //            }
    //            pi_gb.AddRels(pi_rels, t_max, rings_.basis_ss_possEinf);
    //            count_homotopy += (int)pi_rels.size();

    //            /* Add new generators in homotopy */
    //            int2d Einf;
    //            {
    //                auto& sc = ut::GetRecentValue(nodes_ss, deg);
    //                size_t first_PC = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
    //                size_t last_PC = GetFirstIndexOnLevel(sc, LEVEL_PERM + 1);
    //                for (size_t i = first_PC; i < last_PC; ++i)
    //                    Einf.push_back(sc.basis[i]);
    //            }
    //            int2d new_generators = QuotientSpace(Einf, image); /* image is supposed to be a subspace of Einf */
    //            size_t gen_size_old = pi_gb.gen_degs().size();
    //            for (size_t i = 0; i < new_generators.size(); ++i) {
    //                pi_gb.AddGen(deg);
    //                pi_gen_Einf.push_back(Indices2Poly(new_generators[i], basis_d));
    //                pi_basis_d_v2.push_back(algZ::Mon::Gen(uint32_t(gen_size_old + i), 1, deg.s, deg.stem() % 2 == 0));
    //                pi_basis_d_Einf.push_back(new_generators[i]);
    //                Logger::LogHtpyGen(depth, enumReason::ss2htpy, "S0", deg, pi_gb.gen_degs().size() - 1, pi_gen_Einf.back());
    //            }
    //            count_homotopy += (int)new_generators.size();

    //            if (deg.stem() % 2 == 1) {
    //                algZ::Poly1d pi_trivial_rels;
    //                for (size_t i = 0; i < new_generators.size(); ++i)
    //                    pi_trivial_rels.push_back(algZ::Mon::two_x_square(uint32_t(gen_size_old + i), deg.s));
    //                pi_gb.AddRels(pi_trivial_rels, t_max, rings_.basis_ss_possEinf);
    //            }

    //            auto indices = ut::size_t_range(pi_basis_d_v2.size());
    //            std::sort(indices.begin(), indices.end(), [&pi_basis_d_v2](size_t i, size_t j) { return pi_basis_d_v2[i] < pi_basis_d_v2[j]; });
    //            nodes_pi_basis.back()[deg].nodes_pi_basis.clear();
    //            nodes_pi_basis.back()[deg].Einf.clear();
    //            for (auto i : indices) {
    //                nodes_pi_basis.back().at(deg).nodes_pi_basis.push_back(pi_basis_d_v2[i]);
    //                nodes_pi_basis.back().at(deg).Einf.push_back(std::move(pi_basis_d_Einf[i]));
    //            }
    //        }
    //        else {
    //            /* Add new relations in homotopy */
    //            algZ::Poly1d pi_rels;
    //            for (auto& m : pi_basis_d)
    //                pi_rels.push_back(m + algZ::Mon::O(deg.s + 1));
    //            pi_gb.AddRels(pi_rels, t_max, rings_.basis_ss_possEinf);
    //        }
    //    }
    //}
    // if (count_ss_old != count_ss)
    //    UpdatePossEinf(rings_.nodes_ss, rings_.basis_ss_possEinf);
}

void Diagram::SyncCofHomotopy(int iCof, AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
{
    // auto& ssCof = modules_[iCof];
    // int t_max = ssCof.t_max;
    // auto& gb = ssCof.gb;
    // auto& basis = ssCof.basis;
    // auto& nodes_ss = ssCof.nodes_ss;
    // auto& pi_gb = ssCof.pi_gb;
    // auto& nodes_pi_basis = ssCof.nodes_pi_basis;
    // auto& pi_gen_Einf = ssCof.pi_gen_Einf;
    // auto& pi_qt = ssCof.nodes_pi_qt;

    // int count_ss_old = count_ss;

    // Mod tmp;
    // for (int t = 0; t <= t_max; ++t) {
    //     for (int s = 0; s <= t; ++s) {
    //         AdamsDeg deg(s, t);
    //         if (deg.stem() < deg_min.stem() || deg.s < deg_min.s)
    //             continue;
    //         algZ::MMod1d pi_basis_d = GenBasis(pi_gb, deg, rings_.nodes_pi_basis);
    //         if (basis.find(deg) != basis.end()) {
    //             /* Add new boundaries to ss */
    //             int2d boundaries = GetModuleGbEinf(iCof, deg);
    //             const int r = deg.s + 1;
    //             const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
    //             for (auto& boundary : boundaries) {
    //                 if (!boundaries.empty()) {
    //                     const int count1 = SetModuleDiffGlobal(iCof, deg_src, {}, boundary, r);
    //                     if (count1 > 0) {
    //                         Logger::LogDiffBoun(depth, enumReason::htpy2ss, ssCof.name, deg, boundary);
    //                         count_ss += count1;
    //                     }
    //                 }
    //             }

    //            /* Construct the projection map */
    //            auto& basis_d = basis.at(deg);
    //            int2d projs;
    //            {
    //                auto& sc = ut::GetRecentValue(nodes_ss, deg);
    //                size_t first_PC = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
    //                for (auto& m : pi_basis_d) {
    //                    int1d proj = Mod2Indices(gb.Reduce(Proj(m, rings_.pi_gen_Einf, pi_gen_Einf)), basis_d);
    //                    proj = lina::Residue(sc.basis.begin(), sc.basis.begin() + first_PC, proj);
    //                    projs.push_back(proj);
    //                }
    //            }
    //            int2d image, kernel, g;
    //            lina::SetLinearMap(projs, image, kernel, g);

    //            /* Compute pi_basis */
    //            int1d non_leads = lina::add(ut::int_range((int)pi_basis_d.size()), lina::GetLeads(kernel));
    //            algZ::MMod1d pi_basis_d_v2;
    //            int2d pi_basis_d_Einf;
    //            for (int i : non_leads) {
    //                pi_basis_d_v2.push_back(pi_basis_d[i]);
    //                pi_basis_d_Einf.push_back(projs[i]);
    //            }

    //            /* Add new relations in homotopy */
    //            algZ::Mod1d pi_rels;
    //            for (auto& k : kernel) {
    //                pi_rels.push_back(algZ::Indices2Mod(k, pi_basis_d) + algZ::MMod::O(deg.s + 1));
    //                Logger::LogHtpyRel(depth, enumReason::ss2htpy, ssCof.name, deg, pi_rels.back());
    //            }
    //            pi_gb.AddRels(pi_rels, t_max, ssCof.basis_ss_possEinf);
    //            count_homotopy += (int)pi_rels.size();

    //            /* Add new generators in homotopy */
    //            int2d Einf;
    //            {
    //                auto& sc = ut::GetRecentValue(nodes_ss, deg);
    //                size_t first_PC = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
    //                size_t last_PC = GetFirstIndexOnLevel(sc, LEVEL_PERM + 1);
    //                for (size_t i = first_PC; i < last_PC; ++i)
    //                    Einf.push_back(sc.basis[i]);
    //            }
    //            int2d new_generators = QuotientSpace(Einf, image); /* image is supposed to be a subspace of Einf */
    //            size_t gen_size_old = pi_gb.v_degs().size();
    //            for (size_t i = 0; i < new_generators.size(); ++i) {
    //                pi_gb.AddGen(deg);
    //                pi_gen_Einf.push_back(Indices2Mod(new_generators[i], basis_d));
    //                pi_basis_d_v2.push_back(algZ::MMod({}, (uint32_t)(gen_size_old + i), deg.s));
    //                pi_basis_d_Einf.push_back(new_generators[i]);
    //                Logger::LogHtpyGen(depth, enumReason::ss2htpy, ssCof.name, deg, pi_gb.v_degs().size() - 1, pi_gen_Einf.back());

    //                Mod x = Indices2Mod(new_generators[i], basis.at(deg));
    //                Poly fx = rings_.gb.Reduce(subs(x, ssCof.qt));
    //                AdamsDeg deg_S0 = deg - ssCof.deg_qt;
    //                if (fx) {
    //                    int1d ifx = Poly2Indices(fx, rings_.basis.at(deg_S0));
    //                    auto& sc = ut::GetRecentValue(rings_.nodes_ss, deg_S0);
    //                    size_t first = GetFirstIndexOnLevel(sc, LEVEL_MAX / 2);
    //                    ifx = lina::Residue(sc.basis.begin(), sc.basis.begin() + first, std::move(ifx));

    //                    int2d image, g;
    //                    auto& pi_basis_d = *GetRecentPiBasis(rings_.nodes_pi_basis, deg_S0);
    //                    lina::GetInvMap(pi_basis_d.Einf, image, g);
    //                    int1d pi_ifx = lina::GetImage(image, g, ifx);
    //                    pi_qt.back().push_back(algZ::Indices2Poly(pi_ifx, pi_basis_d.nodes_pi_basis) + algZ::Mon::O(deg_S0.s + 1));
    //                }
    //                else
    //                    pi_qt.back().push_back(algZ::Mon::O(deg_S0.s + 1));
    //            }
    //            count_homotopy += (int)new_generators.size();

    //            auto indices = ut::size_t_range(pi_basis_d_v2.size());
    //            std::sort(indices.begin(), indices.end(), [&pi_basis_d_v2](size_t i, size_t j) { return pi_basis_d_v2[i] < pi_basis_d_v2[j]; });
    //            nodes_pi_basis.back()[deg].nodes_pi_basis.clear();
    //            nodes_pi_basis.back()[deg].Einf.clear();
    //            for (auto i : indices) {
    //                nodes_pi_basis.back().at(deg).nodes_pi_basis.push_back(pi_basis_d_v2[i]);
    //                nodes_pi_basis.back().at(deg).Einf.push_back(std::move(pi_basis_d_Einf[i]));
    //            }
    //        }
    //        else {
    //            /* Add new relations in homotopy */
    //            algZ::Mod1d pi_rels;
    //            for (auto& m : pi_basis_d)
    //                pi_rels.push_back(m + algZ::MMod::O(deg.s + 1));
    //            pi_gb.AddRels(pi_rels, t_max, ssCof.basis_ss_possEinf);
    //        }
    //    }
    //}
    // if (count_ss_old != count_ss)
    //    UpdatePossEinf(ssCof.nodes_ss, ssCof.basis_ss_possEinf);
}

int Diagram::DeduceTrivialExtensions(int depth)
{
    int count_homotopy = 0;

    ///* top cell maps */
    //{
    //    algZ::Poly1d new_rels_S0; /* By exactness qi * h = 0 */
    //    for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //        auto& ssCof = modules_[iCof];
    //        auto& q = ssCof.nodes_pi_qt.back();
    //        for (size_t i = 0; i < q.size(); ++i) {
    //            auto& qi = q[i];
    //            AdamsDeg deg = ssCof.pi_gb.v_degs()[i] - ssCof.deg_qt;
    //            algZ::Poly qi_extended;
    //            if (qi && ExtendRelRing(deg.stem(), qi, qi_extended)) {
    //                Logger::LogHtpyMap(depth, enumReason::degree, ssCof.name, ssCof.pi_gb.v_degs()[i], "q", i, qi, qi_extended);
    //                qi = std::move(qi_extended);
    //                ++count_homotopy;
    //                algZ::Poly h = rings_.pi_gb.Gen((uint32_t)iCof);
    //                algZ::Poly relS0 = rings_.pi_gb.ReduceForGbRel(qi * h);
    //                if (algZ::IsValidRel(relS0)) {
    //                    new_rels_S0.push_back(std::move(relS0));
    //                    Logger::LogHtpyRel(depth, enumReason::nat, "S0", deg, new_rels_S0.back());
    //                    ++count_homotopy;
    //                }
    //            }
    //        }
    //    }
    //    AddPiRelsRing(std::move(new_rels_S0));
    //}

    ///* multiplicative structures */
    // int old_count_homotopy = count_homotopy;
    // std::vector<size_t> indices_start(all_basis_ss_.size(), 0);
    // while (true) {
    //     { /* sphere */
    //         auto& pi_gb = rings_.pi_gb;
    //         size_t iSS = 0;

    //        algZ::Poly1d new_rels;
    //        size_t indices_end = pi_gb.data().size();
    //        for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
    //            auto& rel = pi_gb.data()[i];
    //            AdamsDeg deg = GetDeg(rel.GetLead(), pi_gb.gen_degs());
    //            algZ::Poly rel_extended;
    //            if (ExtendRelRing(deg.stem(), rel, rel_extended)) {
    //                rel_extended = pi_gb.ReduceForGbRel(std::move(rel_extended));
    //                if (algZ::IsValidRel(rel_extended)) {
    //                    Logger::LogHtpyRel(depth, enumReason::degree, "S0", deg, rel, rel_extended);
    //                    new_rels.push_back(std::move(rel_extended));
    //                    ++count_homotopy;
    //                }
    //            }
    //        }
    //        indices_start[iSS] = indices_end;
    //        AddPiRelsRing(std::move(new_rels));
    //    }
    //    for (size_t iCof = 0; iCof < modules_.size(); ++iCof) { /* module */
    //        auto& ssCof = modules_[iCof];
    //        auto& pi_gb = ssCof.pi_gb;
    //        auto& nodes_ss = ssCof.nodes_ss;
    //        size_t iSS = iCof + 1;

    //        algZ::Mod1d new_rels;
    //        size_t indices_end = pi_gb.data().size();
    //        for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
    //            auto& rel = pi_gb.data()[i];
    //            AdamsDeg deg = GetDeg(rel.GetLead(), rings_.pi_gb.gen_degs(), pi_gb.v_degs());
    //            algZ::Mod rel_extended;
    //            if (ExtendRelMod(iCof, deg.stem(), rel, rel_extended)) {
    //                rel_extended = pi_gb.ReduceForGbRel(std::move(rel_extended));
    //                if (algZ::IsValidRel(rel_extended)) {
    //                    Logger::LogHtpyRel(depth, enumReason::degree, ssCof.name, deg, rel, rel_extended);
    //                    new_rels.push_back(std::move(rel_extended));
    //                    ++count_homotopy;
    //                }
    //            }
    //        }
    //        indices_start[iSS] = indices_end;
    //        AddPiRelsCof(iCof, std::move(new_rels));
    //    }

    //    if (old_count_homotopy != count_homotopy) {
    //        old_count_homotopy = count_homotopy;
    //        continue;
    //    }

    //    {
    //        auto& nodes_pi_basis = rings_.nodes_pi_basis;
    //        auto& pi_gb = rings_.pi_gb;
    //        const int t_max = rings_.t_max;
    //        algZ::Poly1d new_rels_S0;
    //        algZ::Mod2d new_rels_Cofs(modules_.size());
    //        algZ::Poly tmp;
    //        algZ::Mod tmpm;
    //        for (auto& [deg, _] : nodes_pi_basis.front()) {
    //            for (auto& m : GetRecentPiBasis(nodes_pi_basis, deg)->nodes_pi_basis) {
    //                for (uint32_t gen_id = 0; gen_id < 4; ++gen_id) { /* Warning: in general gen_id might be incorrect */
    //                    AdamsDeg deg_prod = deg + rings_.pi_gb.gen_degs()[gen_id];
    //                    if (deg_prod.t <= t_max) {
    //                        algZ::Poly prod;
    //                        auto h = pi_gb.Gen(gen_id);
    //                        prod.iaddmulP(h, m, tmp, pi_gb.gen_2tor_degs());
    //                        auto prod_reduced = pi_gb.ReduceV2(prod);
    //                        algZ::Poly prod_extended;
    //                        if (prod_reduced && deg_prod.t + prod_reduced.GetLead().fil() - deg_prod.s <= t_max && ExtendRelRing(deg_prod.stem(), prod_reduced, prod_extended)) {
    //                            Logger::LogHtpyProd(depth, enumReason::degree, "S0", deg, h, m, prod_reduced, prod_extended);
    //                            prod.isubP(prod_extended, tmp, pi_gb.gen_2tor_degs());
    //                            new_rels_S0.push_back(std::move(prod));
    //                            ++count_homotopy;
    //                        }
    //                    }
    //                }
    //                for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //                    if (deg.t <= modules_[iCof].t_max) {
    //                        algZ::Mod x(m, 0, 0);
    //                        auto x_reduced = modules_[iCof].pi_gb.ReduceV2(x);
    //                        algZ::Mod x_extended;
    //                        if (x_reduced && deg.t + x_reduced.GetLead().fil() - deg.s <= t_max && ExtendRelMod(iCof, deg.stem(), x_reduced, x_extended)) {
    //                            Logger::LogHtpyProd(depth, enumReason::degree, modules_[iCof].name, deg, m, MOD_V0, x_reduced, x_extended);
    //                            x.isubP(x_extended, tmpm, pi_gb.gen_2tor_degs());
    //                            new_rels_Cofs[iCof].push_back(std::move(x));
    //                            ++count_homotopy;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        AddPiRelsRing(std::move(new_rels_S0));
    //        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
    //            AddPiRelsCof(iCof, std::move(new_rels_Cofs[iCof]));
    //    }
    //    for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //        auto& ssCof = modules_[iCof];
    //        auto& nodes_pi_basis = ssCof.nodes_pi_basis;
    //        auto& pi_gb = ssCof.pi_gb;
    //        const int t_max = ssCof.t_max;
    //        algZ::Mod1d new_rels_cof;
    //        algZ::Mod tmp;
    //        for (auto& [deg, _] : nodes_pi_basis.front()) {
    //            for (auto& m : GetRecentPiBasis(nodes_pi_basis, deg)->nodes_pi_basis) {
    //                for (uint32_t gen_id = 0; gen_id < 4; ++gen_id) {
    //                    AdamsDeg deg_prod = deg + rings_.pi_gb.gen_degs()[gen_id];
    //                    if (deg_prod.t <= t_max) {
    //                        auto h = rings_.pi_gb.Gen(gen_id);
    //                        auto prod = h * m;
    //                        auto prod_reduced = pi_gb.ReduceV2(prod);
    //                        algZ::Mod prod_extended;
    //                        if (prod_reduced && deg_prod.t + prod_reduced.GetLead().fil() - deg_prod.s <= t_max && ExtendRelMod(iCof, deg_prod.stem(), prod_reduced, prod_extended)) {
    //                            Logger::LogHtpyProd(depth, enumReason::degree, modules_[iCof].name, deg, h, m, prod_reduced, prod_extended);
    //                            prod.isubP(prod_extended, tmp, rings_.pi_gb.gen_2tor_degs());
    //                            new_rels_cof.push_back(std::move(prod));
    //                            ++count_homotopy;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        AddPiRelsCof(iCof, std::move(new_rels_cof));
    //    }

    //    if (old_count_homotopy != count_homotopy)
    //        old_count_homotopy = count_homotopy;
    //    else
    //        break;
    //}
    return count_homotopy;
}

int Diagram::DeduceExtensions2tor()
{
    int count_homotopy = 0;
    //{
    //    ut::map_seq2d<int, 0> torS0;
    //    for (size_t stem = 0; stem < rings_.basis_ss_possEinf.data.size(); ++stem) {
    //        for (size_t s = rings_.basis_ss_possEinf.data[stem].size(); s-- > 0;) {
    //            if (rings_.basis_ss_possEinf.data[stem][s] > 0)
    //                torS0(stem, s) = torS0.at(stem, s + 1) + 1;
    //            else
    //                torS0(stem, s) = torS0.at(stem, s + 1);
    //        }
    //    }

    //    algZ::Poly1d new_rels;
    //    int stem_max = (rings_.t_max * 2 - 2) / 3;  ////
    //    for (size_t i = 1; i < rings_.pi_gb.gen_degs().size(); ++i) {
    //        int stem = rings_.pi_gb.gen_degs()[i].stem();
    //        if (stem > stem_max)
    //            continue;
    //        algZ::Poly g = rings_.pi_gb.ReduceV2(rings_.pi_gb.Gen((uint32_t)i));
    //        if (!g)
    //            continue;
    //        algZ::Poly two = algZ::Poly::twoTo(1);
    //        uint32_t e = 0;
    //        while (true) {
    //            g = rings_.pi_gb.ReduceV2(two * g);
    //            ++e;
    //            if (!g)
    //                break;
    //            if (g.GetLead().IsUnKnown()) {
    //                int s = g.GetLead().fil();
    //                int tor = e + torS0.at(stem, s);
    //                algZ::Poly rel = algZ::Poly::twoTo(tor) * rings_.pi_gb.Gen((uint32_t)i);
    //                rel = rings_.pi_gb.Reduce(std::move(rel));
    //                if (algZ::IsValidRel(rel)) {
    //                    AdamsDeg deg = GetDeg(rel.GetLead(), rings_.pi_gb.gen_degs());
    //                    new_rels.push_back(std::move(rel));
    //                    Logger::LogHtpyRel(0, enumReason::degree, "S0", deg, new_rels.back());
    //                    ++count_homotopy;
    //                }
    //                break;
    //            }
    //        }
    //    }
    //    AddPiRelsRing(std::move(new_rels));
    //}
    // for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //    auto& ssCof = modules_[iCof];
    //    ut::map_seq2d<int, 0> torCof;
    //    for (size_t stem = 0; stem < ssCof.basis_ss_possEinf.data.size(); ++stem) {
    //        for (size_t s = ssCof.basis_ss_possEinf.data[stem].size(); s-- > 0;) {
    //            if (ssCof.basis_ss_possEinf.data[stem][s] > 0)
    //                torCof(stem, s) = torCof.at(stem, s + 1) + 1;
    //            else
    //                torCof(stem, s) = torCof.at(stem, s + 1);
    //        }
    //    }

    //    algZ::Mod1d new_rels;
    //    int stem_max = (ssCof.t_max * 2 - 2) / 3;  ////
    //    for (size_t i = 1; i < ssCof.pi_gb.v_degs().size(); ++i) {
    //        int stem = ssCof.pi_gb.v_degs()[i].stem();
    //        if (stem > stem_max)
    //            continue;
    //        algZ::Mod g = ssCof.pi_gb.ReduceV2(ssCof.pi_gb.Gen((uint32_t)i));
    //        if (!g)
    //            continue;
    //        algZ::Poly two = algZ::Poly::twoTo(1);
    //        uint32_t e = 0;
    //        while (true) {
    //            g = ssCof.pi_gb.ReduceV2(two * g);
    //            ++e;
    //            if (!g)
    //                break;
    //            if (g.GetLead().IsUnKnown()) {
    //                int s = g.GetLead().fil();
    //                int tor = e + torCof.at(stem, s);
    //                algZ::Mod rel = algZ::Poly::twoTo(tor) * ssCof.pi_gb.Gen((uint32_t)i);
    //                rel = ssCof.pi_gb.Reduce(std::move(rel));
    //                if (algZ::IsValidRel(rel)) {
    //                    AdamsDeg deg = GetDeg(rel.GetLead(), rings_.pi_gb.gen_degs(), ssCof.pi_gb.v_degs());
    //                    new_rels.push_back(std::move(rel));
    //                    Logger::LogHtpyRel(0, enumReason::degree, ssCof.name, deg, new_rels.back());
    //                    ++count_homotopy;
    //                }
    //                break;
    //            }
    //        }
    //    }
    //    AddPiRelsCof(iCof, std::move(new_rels));
    //}
    return 0;
}

int Diagram::DeduceExtensionsByExactness(int stem_min_para, int stem_max_para, int depth)
{
    int count_htpy = 0;

    /* Populate `basis_S0` and `basis_Cofs` */
    // algZ::Mon2d basis_S0(size_t(rings_.t_max + 1));
    // std::vector<algZ::MMod2d> basis_Cofs(modules_.size());
    // for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
    //     basis_Cofs[iCof].resize(size_t(modules_[iCof].t_max + 1));
    // if (depth == 0) {
    //     for (auto& [deg, pi_basis_d] : rings_.nodes_pi_basis.front())
    //         for (auto& m : pi_basis_d.nodes_pi_basis)
    //             basis_S0[deg.stem()].push_back(m);  //
    //     for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
    //         for (auto& [deg, pi_basis_d] : modules_[iCof].nodes_pi_basis.front())
    //             for (auto& m : pi_basis_d.nodes_pi_basis)
    //                 basis_Cofs[iCof][deg.stem()].push_back(m);
    // }
    // else {
    //     std::set<AdamsDeg> degs_S0;
    //     for (auto& pi_basis_i : rings_.nodes_pi_basis)
    //         for (auto& [deg, _] : pi_basis_i)
    //             degs_S0.insert(deg);
    //     for (AdamsDeg deg : degs_S0)
    //         for (auto& m : GetRecentPiBasis(rings_.nodes_pi_basis, deg)->nodes_pi_basis)
    //             basis_S0[deg.stem()].push_back(m);
    //     std::set<AdamsDeg> degs_Cof;
    //     for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //         degs_Cof.clear();
    //         for (auto& pi_basis_i : modules_[iCof].nodes_pi_basis)
    //             for (auto& [deg, _] : pi_basis_i)
    //                 degs_Cof.insert(deg);
    //         for (AdamsDeg deg : degs_Cof)
    //             for (auto& m : GetRecentPiBasis(modules_[iCof].nodes_pi_basis, deg)->nodes_pi_basis)
    //                 basis_Cofs[iCof][deg.stem()].push_back(m);
    //     }
    // }
    // int1d O1s_S0, O2s_S0, isSingle_S0;
    // PossMoreEinfFirstS_Ring(O1s_S0, O2s_S0, isSingle_S0);
    // algZ::Poly tmp;
    // algZ::Mod tmpm;

    // algZ::Poly1d new_rels_S0;
    // algZ::Mod2d new_rels_Cofs(modules_.size());
    // int1d f_changed(modules_.size(), 0);
    // AdamsDeg2d deg_perms_Cof(modules_.size());
    // for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //     auto& ssCof = modules_[iCof];
    //     auto& pi_f = ssCof.nodes_pi_qt;
    //     int d_f = ssCof.deg_qt.t;
    //     auto& h = rings_.pi_gb.Gen((uint32_t)iCof);

    //    int1d O1s_Cof, O2s_Cof, isSingle_Cof;
    //    PossMoreEinfFirstS_Mod(iCof, O1s_Cof, O2s_Cof, isSingle_Cof);

    //    /* exactness at h * q */
    //    int stem_min = std::max(0, stem_min_para - d_f);
    //    int stem_max = std::min({rings_.t_max, ssCof.t_max - d_f, stem_max_para});
    //    for (int stem = stem_min; stem <= stem_max; ++stem) {
    //        /* kernel of h */
    //        auto& basis_stem = basis_S0[stem];
    //        algZ::Poly1d x_h, hx;
    //        for (size_t i = 0; i < basis_stem.size(); ++i) {
    //            x_h.push_back(basis_stem[i]);
    //            hx.push_back(rings_.pi_gb.ReduceV2(h * basis_stem[i]));
    //        }
    //        ut::map_seq<int, 0> num_leads;
    //        ut::map_seq<int, algZ::FIL_MAX + 1> src_Os;
    //        GetImageLeads(x_h, hx, rings_.pi_gb, rings_.pi_gb, num_leads, src_Os);
    //        for (size_t i = 0; i < basis_stem.size(); ++i) { /* lift maps but yield unknown kernel */
    //            algZ::Poly hxi = hx[i];
    //            bool use_hxi = false;
    //            if (hxi.UnknownFil() < algZ::FIL_MAX) {
    //                ExtendRelRingV2(stem + d_f - 1, hxi, num_leads);
    //                if (hxi.UnknownFil() > algZ::FIL_MAX) {
    //                    int src_O = -1;
    //                    for (size_t j = hx[i].UnknownFil(); j < num_leads.data.size(); ++j) {
    //                        if (num_leads[j] > 0) {
    //                            src_O = src_Os.data[j];
    //                            break;
    //                        }
    //                    }
    //                    if (src_O != -1) {
    //                        if (src_O > basis_stem[i].fil()) {
    //                            x_h[i] += algZ::Poly::O(src_O);
    //                            hx[i] = std::move(hxi);
    //                            use_hxi = true;
    //                        }
    //                    }
    //                    else {
    //                        hx[i] = std::move(hxi);
    //                        use_hxi = true;
    //                    }
    //                }
    //            }
    //            if (!use_hxi)
    //                hx[i] = rings_.pi_gb.Reduce(h * basis_stem[i]);
    //        }
    //        algZ::Poly1d kernel_h;
    //        GetKernel(x_h, hx, rings_.pi_gb, rings_.pi_gb, kernel_h);
    //        ut::map_seq<int, 0> num_leads_kernel;
    //        CountByS(kernel_h, num_leads_kernel);

    //        /* image of q */
    //        algZ::Mod1d x_f;
    //        algZ::Poly1d fx;
    //        int stem_src = stem + d_f;
    //        auto& basis_cof_stem1 = basis_Cofs[iCof][stem_src];
    //        for (size_t i = 0; i < basis_cof_stem1.size(); ++i) {
    //            x_f.push_back(basis_cof_stem1[i]);
    //            auto fxi = algZ::subs(basis_cof_stem1[i], modules_[iCof].nodes_pi_qt.back(), modules_[iCof].pi_gb.v_degs());
    //            fxi = rings_.pi_gb.ReduceV2(std::move(fxi));
    //            fx.push_back(std::move(fxi));
    //        }
    //        algZ::Poly1d image1_f;
    //        int O1 = O1s_Cof[stem_src];
    //        int O2 = O2s_Cof[stem_src];
    //        algZ::Mod gO = algZ::Mod::O(isSingle_Cof[stem_src]);
    //        GetImage(x_f, fx, ssCof.pi_gb, rings_.pi_gb, image1_f, O1, O2, gO);

    //        /* Check exactness */
    //        for (size_t i = 0; i < kernel_h.size(); ++i) {
    //            auto k = Residue(image1_f, kernel_h[i], rings_.pi_gb);
    //            if (k && !k.GetLead().IsUnKnown()) {
    //                int s = k.GetLead().fil();
    //                AdamsDeg deg(s, stem + s);
    //                if (!IsPossTgt(rings_.nodes_ss, deg)) {
    //                    if (s < O1) {
    //                        Logger::LogException(depth, 0xdb0fa447U, "q cannot hit the kernel of h. {} {} k={} O1={}\n", ssCof.name, deg, k, O1);
    //                        throw SSPiException(0xdb0fa447U, "q cannot hit the kernel of h. ");  ////TODO: use reason as message
    //                    }
    //                    if (s < O2) {
    //                        if (!PossMoreEinf(rings_.nodes_ss, deg) && num_leads_kernel[s] == 1) {
    //                            if (gO.data.size() == 1 && !gO.GetLead().IsUnKnown() && !gO.GetLead().m) {
    //                                size_t v = (size_t)gO.GetLead().v;
    //                                algZ::Poly f_original = ssCof.nodes_pi_qt.back()[v];  ////
    //                                auto& f = ssCof.nodes_pi_qt.back()[v];
    //                                f.data.pop_back();
    //                                f.iaddP(k.LF(), tmp);
    //                                f.data.push_back(algZ::Mon::O(s + 1));
    //                                ExtendRelRing(stem, f);
    //                                Logger::LogHtpyMap(depth, enumReason::exact_hq, ssCof.name, deg, "q", v, f_original, f);
    //                                ++count_htpy;
    //                                f_changed[iCof] = 1;
    //                            }
    //                            else if (gO && gO.GetLead().IsUnKnown() && gO.GetLead().fil() == 1) {
    //                                AdamsDeg deg_perm(O1, stem_src + O1);
    //                                if (depth == 0) {
    //                                    deg_perms_Cof[iCof].push_back(deg_perm);
    //                                    Logger::LogDiffPerm(depth, enumReason::exact_hq, ssCof.name, deg_perm);  ////
    //                                    auto hx0 = rings_.pi_gb.ReduceV2(h * basis_stem[0]);
    //                                    ExtendRelRingV2(stem + d_f - 1, hx0, num_leads);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //                else if (s < O1) {
    //                    if (!IsPossTgt(rings_.nodes_ss, AdamsDeg(s - 1, stem + s - 1))) {
    //                        new_rels_S0.push_back(std::move(k));
    //                        Logger::LogHtpyRel(depth, enumReason::exact_hq, "S0", deg, new_rels_S0.back());
    //                        ++count_htpy;
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    /* exactness at i * h */
    //    stem_min = std::max(d_f - 1, stem_min_para);
    //    stem_max = std::min({rings_.t_max, ssCof.t_max, stem_max_para});
    //    for (int stem = stem_min; stem <= stem_max; ++stem) {
    //        /* kernel of i */
    //        auto& basis_stem = basis_S0[stem];
    //        algZ::Poly1d x_i;
    //        algZ::Mod1d ix;
    //        for (size_t i = 0; i < basis_stem.size(); ++i) {
    //            x_i.push_back(basis_stem[i]);
    //            ix.push_back(ssCof.pi_gb.ReduceV2(algZ::Mod(basis_stem[i], 0, 0)));
    //        }
    //        ut::map_seq<int, 0> num_leads;
    //        ut::map_seq<int, algZ::FIL_MAX + 1> src_Os;
    //        GetImageLeads(x_i, ix, rings_.pi_gb, ssCof.pi_gb, num_leads, src_Os);
    //        for (size_t i = 0; i < basis_stem.size(); ++i) {
    //            auto ixi = ssCof.pi_gb.ReduceV2(algZ::Mod(basis_stem[i], 0, 0));
    //            bool use_ixi = false;
    //            if (ixi.UnknownFil() < algZ::FIL_MAX) {
    //                int old_fil = ixi.UnknownFil();
    //                ExtendRelCofV2(iCof, stem, ixi, num_leads);
    //                if (ixi.UnknownFil() > algZ::FIL_MAX) {
    //                    int src_O = -1;
    //                    for (size_t j = old_fil; j < num_leads.data.size(); ++j) {
    //                        if (num_leads[j] > 0) {
    //                            src_O = src_Os.data[j];
    //                            break;
    //                        }
    //                    }
    //                    if (src_O != -1) {
    //                        if (src_O > basis_stem[i].fil()) {
    //                            x_i[i] += algZ::Poly::O(src_O);
    //                            ix[i] = std::move(ixi);
    //                            use_ixi = true;
    //                        }
    //                    }
    //                    else {
    //                        ix[i] = std::move(ixi);
    //                        use_ixi = true;
    //                    }
    //                }
    //            }
    //            if (!use_ixi)
    //                ix[i] = ssCof.pi_gb.Reduce(algZ::Mod(basis_stem[i], 0, 0));
    //        }
    //        algZ::Poly1d kernel_i;
    //        GetKernel(x_i, ix, rings_.pi_gb, modules_[iCof].pi_gb, kernel_i);
    //        ut::map_seq<int, 0> num_leads_kernel;
    //        CountByS(kernel_i, num_leads_kernel);

    //        /* image of h */
    //        algZ::Poly1d x_h;
    //        algZ::Poly1d image_h;
    //        int stem_src = stem - d_f + 1;
    //        auto& basis_S0_stem1 = basis_S0[stem_src];
    //        for (size_t i = 0; i < basis_S0_stem1.size(); ++i) {
    //            x_h.push_back(basis_S0_stem1[i]);
    //            auto fxi = rings_.pi_gb.ReduceV2(h * basis_S0_stem1[i]);
    //            image_h.push_back(std::move(fxi));
    //        }
    //        algZ::Poly1d image1_h;

    //        int O1 = O1s_S0[stem_src] + 1;
    //        int O2 = O2s_S0[stem_src] + 1;
    //        algZ::Poly gO = algZ::Poly::O(isSingle_S0[stem_src]);

    //        GetImage(x_h, image_h, rings_.pi_gb, rings_.pi_gb, image1_h, O1, O2, gO);

    //        /* Check exactness */
    //        for (size_t i = 0; i < kernel_i.size(); ++i) {
    //            auto k = Residue(image1_h, std::move(kernel_i[i]), rings_.pi_gb);
    //            if (k && !k.GetLead().IsUnKnown()) {
    //                int s = k.GetLead().fil();
    //                AdamsDeg deg(s, stem + s);
    //                if (!IsPossTgt(rings_.nodes_ss, deg)) {
    //                    if (s < O1) {
    //                        Logger::LogException(depth, 0xc927323U, "h cannot hit the kernel of i. {} {} k={} O1={}\n", ssCof.name, deg, k, O1);
    //                        throw SSPiException(0xc927323U, "h cannot hit the kernel of i");
    //                    }
    //                    if (s < O2) {
    //                        if (!PossMoreEinf(rings_.nodes_ss, deg) && num_leads_kernel[s] == 1) {
    //                            if (gO && !gO.GetLead().IsUnKnown()) {
    //                                auto& rel = algZ::Poly(h) * gO + k.LF() + algZ::Poly::O(s + 1);
    //                                new_rels_S0.push_back(std::move(rel));
    //                                Logger::LogHtpyRel(depth, enumReason::exact_ih, ssCof.name, deg, new_rels_S0.back());
    //                                ++count_htpy;
    //                            }
    //                            else if (gO && gO.GetLead().IsUnKnown() && gO.GetLead().fil() == 1) {
    //                                AdamsDeg deg_perm(O1 - 1, stem_src + O1 - 1);
    //                                if (depth == 0)
    //                                    fmt::print(fmt::fg(fmt::color::light_pink), "exact_ih - S0 {} needs one more permanent cycle to hit kernel of i(to {})\n", deg_perm, ssCof.name);
    //                            }
    //                        }
    //                    }
    //                }
    //                else if (s < O1) {
    //                    if (!IsPossTgt(rings_.nodes_ss, AdamsDeg(s - 1, stem + s - 1))) {
    //                        new_rels_S0.push_back(std::move(k));
    //                        Logger::LogHtpyRel(depth, enumReason::exact_ih, "S0", deg, new_rels_S0.back());
    //                        ++count_htpy;
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    /* exactness at q * i */
    //    stem_min = std::max(0, stem_min_para);
    //    stem_max = std::min({rings_.t_max, ssCof.t_max, stem_max_para});
    //    for (int stem = stem_min; stem <= stem_max; ++stem) {
    //        /* kernel of q */
    //        auto& basis_stem = basis_Cofs[iCof][stem];
    //        algZ::Mod1d x_f;
    //        algZ::Poly1d fx;
    //        for (size_t i = 0; i < basis_stem.size(); ++i) {
    //            x_f.push_back(basis_stem[i]);
    //            auto fxi = algZ::subs(basis_stem[i], modules_[iCof].nodes_pi_qt.back(), modules_[iCof].pi_gb.v_degs());
    //            fxi = rings_.pi_gb.ReduceV2(std::move(fxi));
    //            fx.push_back(fxi);
    //        }
    //        ut::map_seq<int, 0> num_leads;
    //        ut::map_seq<int, algZ::FIL_MAX + 1> src_Os;
    //        GetImageLeads(x_f, fx, ssCof.pi_gb, rings_.pi_gb, num_leads, src_Os);
    //        x_f.clear();
    //        fx.clear();
    //        for (size_t i = 0; i < basis_stem.size(); ++i) {
    //            auto fxi = algZ::subs(basis_stem[i], modules_[iCof].nodes_pi_qt.back(), modules_[iCof].pi_gb.v_degs());
    //            auto fxi_v2 = rings_.pi_gb.ReduceV2(fxi);
    //            bool use_fxi_v2 = false;
    //            if (fxi_v2.UnknownFil() < algZ::FIL_MAX) {
    //                int old_fil = fxi_v2.UnknownFil();
    //                ExtendRelRingV2(stem - d_f, fxi_v2, num_leads);
    //                if (fxi_v2.UnknownFil() > algZ::FIL_MAX) {
    //                    int src_O = -1;
    //                    for (size_t j = old_fil; j < num_leads.data.size(); ++j) {
    //                        if (num_leads[j] > 0) {
    //                            src_O = src_Os.data[j];
    //                            break;
    //                        }
    //                    }
    //                    if (src_O != -1) {
    //                        if (src_O > basis_stem[i].fil()) {
    //                            x_f.push_back(basis_stem[i]);
    //                            x_f.back().iaddP(algZ::Mod::O(src_O), tmpm);
    //                            fx.push_back(std::move(fxi_v2));
    //                            use_fxi_v2 = true;
    //                        }
    //                    }
    //                    else {
    //                        x_f.push_back(basis_stem[i]);
    //                        fx.push_back(std::move(fxi_v2));
    //                        use_fxi_v2 = true;
    //                    }
    //                }
    //            }
    //            if (!use_fxi_v2) {
    //                fxi = rings_.pi_gb.Reduce(std::move(fxi));
    //                ExtendRelRingV2(stem - d_f, fxi, num_leads);
    //                if (fxi.UnknownFil() > algZ::FIL_MAX) {
    //                    x_f.push_back(basis_stem[i]);
    //                    fx.push_back(std::move(fxi));
    //                }
    //            }
    //        }

    //        algZ::Mod1d kernel_f;
    //        GetKernel(x_f, fx, modules_[iCof].pi_gb, rings_.pi_gb, kernel_f);
    //        ut::map_seq<int, 0> num_leads_kernel;
    //        CountByS(kernel_f, num_leads_kernel);
    //        ut::RemoveIf(kernel_f, [](const algZ::Mod& k) { return k.GetLead().v == 0; });

    //        /* image of i */
    //        algZ::Poly1d x_i;
    //        algZ::Mod1d image_i;
    //        int stem_src = stem;
    //        auto& basis_S0_stem1 = basis_S0[stem_src];
    //        for (size_t i = 0; i < basis_S0_stem1.size(); ++i) {
    //            x_i.push_back(basis_S0_stem1[i]);
    //            auto fxi = ssCof.pi_gb.ReduceV2(algZ::Mod(basis_S0_stem1[i], 0, 0));
    //            image_i.push_back(std::move(fxi));
    //        }
    //        auto indices = ut::size_t_range(image_i.size());
    //        std::stable_sort(indices.begin(), indices.end(), [&image_i](size_t i, size_t j) { return image_i[i].UnknownFil() > image_i[j].UnknownFil(); });
    //        algZ::Poly1d x_i_sorted;
    //        algZ::Mod1d image_i_sorted;
    //        for (size_t i : indices) {
    //            x_i_sorted.push_back(std::move(x_i[i]));
    //            image_i_sorted.push_back(std::move(image_i[i]));
    //        }

    //        algZ::Mod1d image1_i;
    //        int O1 = O1s_S0[stem_src];
    //        int O2 = O2s_S0[stem_src];
    //        algZ::Poly gO = algZ::Poly::O(isSingle_S0[stem_src]);
    //        GetImage(x_i_sorted, image_i_sorted, rings_.pi_gb, ssCof.pi_gb, image1_i, O1, O2, gO);

    //        /* Check exactness */
    //        for (size_t i = 0; i < kernel_f.size(); ++i) {
    //            auto k = Residue(image1_i, std::move(kernel_f[i]), ssCof.pi_gb);
    //            if (k && !k.GetLead().IsUnKnown()) {
    //                int s = k.GetLead().fil();
    //                AdamsDeg deg(s, stem + s);
    //                if (!IsPossTgt(ssCof.nodes_ss, deg)) {
    //                    if (s < O1) {
    //                        Logger::LogException(depth, 0xb977ae09U, "i cannot hit the kernel of q. {} {} k={} O1={}\n", ssCof.name, deg, k, O1);
    //                        throw SSPiException(0xb977ae09U, "i can not hit the kernel of q");
    //                    }
    //                    if (s < O2) {
    //                        if (!PossMoreEinf(ssCof.nodes_ss, deg) && num_leads_kernel[s] == 1) {
    //                            if (gO && !gO.GetLead().IsUnKnown()) {
    //                                auto& rel = algZ::Mod(gO, 0, 0) + k.LF() + algZ::Mod::O(s + 1);
    //                                new_rels_Cofs[iCof].push_back(std::move(rel));
    //                                Logger::LogHtpyRel(depth, enumReason::exact_qi, ssCof.name, deg, new_rels_Cofs[iCof].back());
    //                                ++count_htpy;
    //                            }
    //                            else if (gO && gO.GetLead().IsUnKnown() && gO.GetLead().fil() == 1) {
    //                                AdamsDeg deg_perm(O1, stem_src + O1);
    //                                if (depth == 0)
    //                                    fmt::print(fmt::fg(fmt::color::light_pink), "exact_qi - S0 {} needs one more permanent cycle to hit {}\n", deg_perm, ssCof.name);
    //                            }
    //                        }
    //                    }
    //                }
    //                else if (s < O1) {
    //                    if (!IsPossTgt(ssCof.nodes_ss, AdamsDeg(s - 1, stem + s - 1))) {
    //                        new_rels_Cofs[iCof].push_back(std::move(k));
    //                        Logger::LogHtpyRel(depth, enumReason::exact_qi, ssCof.name, deg, new_rels_Cofs[iCof].back());
    //                        ++count_htpy;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}
    // AddPiRelsRing(std::move(new_rels_S0));
    // for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //    AddPiRelsCof(iCof, std::move(new_rels_Cofs[iCof]));
    //    if (f_changed[iCof])
    //        AddPiRelsByNat(iCof);
    //}
    // int sum_perms = 0;
    // for (size_t iCof = 0; iCof < deg_perms_Cof.size(); ++iCof) {
    //    for (auto& deg_x : deg_perms_Cof[iCof]) {
    //        SetPermanentCycle(depth, iCof, deg_x);
    //        ++sum_perms;
    //    }
    //}
    // if (count_htpy + sum_perms) {
    //    int count_ss = 0;
    //    SyncHomotopy(AdamsDeg(0, 0), count_ss, count_htpy, depth);
    //    if (sum_perms)
    //        count_htpy += DeduceTrivialExtensions(depth);
    //}

    return count_htpy;
}

unsigned Diagram::TryExtS0(algZ::Poly rel, AdamsDeg deg_change, int depth, DeduceFlag flag)
{
    AddNode(flag);
    unsigned error = 0;
    /*try {
        Logger::LogHtpyRel(depth + 1, enumReason::try_, "S0", deg_change, rel);
        AddPiRelsRing({std::move(rel)});
        int count_ss1 = 0, count_homotopy1 = 0;
        SyncHomotopy(deg_change, count_ss1, count_homotopy1, depth + 1);
        if (flag & DeduceFlag::pi_exact) {
            if (count_ss1)
                DeduceTrivialExtensions(depth + 1);
            DeduceExtensionsByExactness(deg_change.stem(), stem_max_exactness_, depth + 1);
        }
    }
    catch (SSException& e) {
        error = e.id();
    }*/
    PopNode(flag);

    return error;
}

unsigned Diagram::TryExtCof(size_t iCof, algZ::Mod rel, AdamsDeg deg_change, int depth, DeduceFlag flag)
{
    AddNode(flag);
    unsigned error = 0;
    /*try {
        Logger::LogHtpyRel(depth + 1, enumReason::try_, modules_[iCof].name, deg_change, rel);
        AddPiRelsCof(iCof, {std::move(rel)});
        int count_ss1 = 0, count_homotopy1 = 0;
        auto deg_min = deg_change - modules_[iCof].deg_qt;
        SyncHomotopy(deg_min, count_ss1, count_homotopy1, depth + 1);
        if (flag & DeduceFlag::pi_exact) {
            if (count_ss1)
                DeduceTrivialExtensions(depth + 1);
            DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth + 1);
        }
    }
    catch (SSException& e) {
        error = e.id();
    }
    PopNode(flag);*/

    return error;
}

unsigned Diagram::TryExtQ(size_t iCof, size_t gen_id, algZ::Poly q, AdamsDeg deg_change, int depth, DeduceFlag flag)
{
    AddNode(flag);
    unsigned error = 0;
    /*try {
        Logger::LogHtpyMap(depth + 1, enumReason::try_, modules_[iCof].name, deg_change, "q", gen_id, q);
        modules_[iCof].nodes_pi_qt.back()[gen_id] = std::move(q);
        AddPiRelsByNat(iCof);
        int count_ss1 = 0, count_homotopy1 = 0;
        auto deg_min = deg_change;
        SyncHomotopy(deg_min, count_ss1, count_homotopy1, depth + 1);
        if (flag & DeduceFlag::pi_exact) {
            DeduceTrivialExtensions(depth + 1);
            DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth + 1);
        }
    }
    catch (SSException& e) {
        error = e.id();
    }
    PopNode(flag);*/

    return error;
}

void Diagram::DeduceExtensions(int stem_min, int stem_max, int& count_ss, int& count_htpy, int depth, DeduceFlag flag)
{
    // SyncHomotopy(AdamsDeg(0, 0), count_ss, count_htpy, depth);
    // count_htpy += DeduceTrivialExtensions(depth);
    // if (flag & DeduceFlag::pi_exact)
    //     count_htpy += DeduceExtensionsByExactness(0, stem_max_exactness_, depth);

    // algZ::Poly tmp;
    // algZ::Mod tmpm;
    ///* top cell maps */
    //{
    //    algZ::Poly1d new_rels_S0; /* By exactness h * qi = 0 */
    //    for (size_t iCof = 0; iCof < modules_.size(); ++iCof) {
    //        auto& ssCof = modules_[iCof];
    //        auto& name = ssCof.name;
    //        size_t f_size = modules_[iCof].nodes_pi_qt.back().size();

    //        for (size_t i = 0; i < f_size; ++i) {
    //            std::cout << name << "  qi=" << i << '/' << f_size << "                     \r";
    //            int s = modules_[iCof].nodes_pi_qt.back()[i].UnknownFil();
    //            if (s > algZ::FIL_MAX)
    //                continue;
    //            int stem = modules_[iCof].pi_gb.v_degs()[i].stem() - modules_[iCof].deg_qt.stem();
    //            if (stem < stem_min || stem > stem_max)
    //                continue;
    //            AdamsDeg deg(s, stem + s);
    //            if (deg.t > rings_.t_max || PossMoreEinf(rings_.nodes_ss, deg) > 0)
    //                continue;

    //            algZ::Poly O1 = algZ::Poly::O(s + 1);
    //            ExtendRelRing(stem, O1);
    //            auto& pi_basis_d = GetRecentPiBasis(rings_.nodes_pi_basis, deg)->nodes_pi_basis;
    //            int ne_cout = (int)pi_basis_d.size();
    //            int count_pass = 0;
    //            unsigned i_max = 1 << ne_cout;
    //            algZ::Poly q, q_pass;
    //            /*std::vector<unsigned> error_list;*/
    //            bool bNewExt = false;
    //            for (unsigned b = 1; b < i_max; ++b) {
    //                q = modules_[iCof].nodes_pi_qt.back()[i];
    //                q.data.pop_back();
    //                q.iaddP(O1, tmp);
    //                for (int j : two_expansion(b))
    //                    q.iaddP(pi_basis_d[j], tmp);

    //                if (unsigned error = TryExtQ(iCof, i, q, deg, depth, flag)) {
    //                    /*error_list.push_back(error);*/
    //                }
    //                else {
    //                    ++count_pass;
    //                    if (count_pass > 1)
    //                        break;
    //                    q_pass = std::move(q);
    //                }
    //            }

    //            if (count_pass == 0) {
    //                q_pass = modules_[iCof].nodes_pi_qt.back()[i];
    //                q_pass.data.pop_back();
    //                q_pass.iaddP(O1, tmp);
    //                bNewExt = true;
    //            }
    //            else if (count_pass == 1) {
    //                q = modules_[iCof].nodes_pi_qt.back()[i];
    //                q.data.pop_back();
    //                q.iaddP(O1, tmp);
    //                if (unsigned error = TryExtQ(iCof, i, q, deg, depth, flag)) {
    //                    /*error_list.push_back(error);*/
    //                    bNewExt = true;
    //                }
    //            }
    //            if (bNewExt) {
    //                ++count_htpy;
    //                Logger::PrintDepth();

    //                auto& q1 = modules_[iCof].nodes_pi_qt.back()[i];
    //                Logger::LogHtpyMap(depth, enumReason::deduce, name, deg, "q", i, modules_[iCof].nodes_pi_qt.back()[i], q_pass);
    //                q1 = std::move(q_pass);

    //                AddPiRelsByNat(iCof);
    //                AdamsDeg deg_min = deg;
    //                SyncHomotopy(deg_min, count_ss, count_htpy, depth);
    //                count_htpy += DeduceTrivialExtensions(depth);
    //                if (flag & DeduceFlag::pi_exact)
    //                    count_htpy += DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth);

    //                algZ::Poly h = rings_.pi_gb.Gen((uint32_t)iCof);
    //                algZ::Poly relS0 = rings_.pi_gb.ReduceForGbRel(h * q1);
    //                if (algZ::IsValidRel(relS0)) {
    //                    new_rels_S0.push_back(std::move(relS0));
    //                    Logger::LogHtpyRel(depth, enumReason::nat, "S0", deg, new_rels_S0.back());
    //                    ++count_htpy;
    //                }
    //            }
    //            else
    //                Logger::ClearDepth();
    //        }
    //    }
    //    AddPiRelsRing(std::move(new_rels_S0));

    //    if (depth == 0)
    //        SimplifyPiRels();
    //}

    ///* multiplicative structures */
    // std::vector<size_t> i_start(all_basis_ss_.size(), 0);
    // int t_max = rings_.t_max;
    // for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
    //     if (t_max < modules_[iCof].t_max)
    //         t_max = modules_[iCof].t_max;
    // int t = 0;
    // while (t <= t_max) {
    //     if (ut::has(rings_.pi_gb.leads_group_by_t(), t)) { /* sphere */
    //         auto& pi_gb = rings_.pi_gb;
    //         size_t iSS = 0;
    //         auto& indices = pi_gb.leads_group_by_t().at(t);
    //         size_t indices_size = indices.size();

    //        for (size_t i = i_start[iSS]; i < indices_size; ++i) {
    //            std::cout << "S0  t=" << t << '/' << t_max << "  i=" << i << '/' << indices.size() << "                \r";
    //            algZ::Poly rel = pi_gb.data()[indices[i]];
    //            int s = rel.UnknownFil();
    //            if (s > algZ::FIL_MAX)
    //                continue;
    //            int stem = GetDeg(rel.GetLead(), pi_gb.gen_degs()).stem();
    //            if (stem < stem_min || stem > stem_max)
    //                continue;
    //            AdamsDeg deg(s, stem + s);
    //            if (deg.t > rings_.t_max || PossMoreEinf(rings_.nodes_ss, deg) > 0 || GetRecentPiBasis(rings_.nodes_pi_basis, deg) == nullptr)
    //                continue;
    //            algZ::Poly O1 = algZ::Poly::O(s + 1);

    //            ExtendRelRing(stem, O1);
    //            auto& pi_basis_d = GetRecentPiBasis(rings_.nodes_pi_basis, deg)->nodes_pi_basis;
    //            int ne_cout = (int)pi_basis_d.size();
    //            int count_pass = 0;
    //            unsigned i_max = 1 << ne_cout;
    //            algZ::Poly rel1, rel_pass;
    //            bool bNewExt = false;
    //            for (unsigned b = 1; b < i_max; ++b) {
    //                rel1 = rel;
    //                rel1.data.pop_back();
    //                rel1.iaddP(O1, tmp);
    //                for (int j : two_expansion(b))
    //                    rel1.iaddP(pi_basis_d[j], tmp);

    //                if (!TryExtS0(rel1, deg, depth, flag)) {
    //                    ++count_pass;
    //                    if (count_pass > 1)
    //                        break;
    //                    rel_pass = std::move(rel1);
    //                }
    //            }

    //            if (count_pass == 0) {
    //                rel_pass = rel;
    //                rel_pass.data.pop_back();
    //                rel_pass.iaddP(O1, tmp);
    //                bNewExt = true;
    //            }
    //            else if (count_pass == 1) {
    //                rel1 = rel;
    //                rel1.data.pop_back();
    //                rel1.iaddP(O1, tmp);
    //                if (TryExtS0(rel1, deg, depth, flag))
    //                    bNewExt = true;
    //            }
    //            if (bNewExt) {
    //                rel_pass = rings_.pi_gb.Reduce(std::move(rel_pass));
    //                if (algZ::IsValidRel(rel_pass)) {
    //                    ++count_htpy;
    //                    Logger::PrintDepth();

    //                    Logger::LogHtpyRel(depth, enumReason::deduce, "S0", deg, rel, rel_pass);
    //                    AddPiRelsRing({std::move(rel_pass)});
    //                    SyncHomotopy(deg, count_ss, count_htpy, depth);
    //                    count_htpy += DeduceTrivialExtensions(depth);
    //                    if (flag & DeduceFlag::pi_exact)
    //                        count_htpy += DeduceExtensionsByExactness(deg.stem(), stem_max_exactness_, depth);
    //                }
    //                else
    //                    Logger::ClearDepth();
    //            }
    //            else
    //                Logger::ClearDepth();
    //        }
    //        i_start[iSS] = indices_size;
    //    }

    //    for (size_t iCof = 0; iCof < modules_.size(); ++iCof) { /* module */
    //        auto& ssCof = modules_[iCof];
    //        if (ut::has(ssCof.pi_gb.leads_group_by_t(), t)) {
    //            auto& name = ssCof.name;
    //            auto& nodes_ss = ssCof.nodes_ss;
    //            auto& pi_gb = ssCof.pi_gb;
    //            auto& nodes_pi_basis = ssCof.nodes_pi_basis;
    //            size_t iSS = iCof + 1;
    //            auto& indices = pi_gb.leads_group_by_t().at(t);
    //            size_t indices_size = indices.size();

    //            for (size_t i = i_start[iSS]; i < indices_size; ++i) {
    //                std::cout << name << "  t=" << t << '/' << t_max << "  i=" << i << '/' << indices.size() << "                \r";
    //                algZ::Mod rel = pi_gb.data()[indices[i]];
    //                int s = rel.UnknownFil();
    //                if (s > algZ::FIL_MAX)
    //                    continue;
    //                int stem = GetDeg(rel.GetLead(), rings_.pi_gb.gen_degs(), pi_gb.v_degs()).stem();
    //                if (stem < stem_min || stem > stem_max)
    //                    continue;
    //                AdamsDeg deg(s, stem + s);
    //                if (deg.t > ssCof.t_max || PossMoreEinf(nodes_ss, deg) > 0 || GetRecentPiBasis(ssCof.nodes_pi_basis, deg) == nullptr)
    //                    continue;

    //                algZ::Mod O1 = algZ::Mod::O(s + 1);
    //                ExtendRelMod(iCof, stem, O1);
    //                auto& pi_basis_d = GetRecentPiBasis(ssCof.nodes_pi_basis, deg)->nodes_pi_basis;
    //                int ne_cout = (int)pi_basis_d.size();
    //                int count_pass = 0;
    //                unsigned i_max = 1 << ne_cout;
    //                algZ::Mod rel1, rel_pass;
    //                bool bNewExt = false;
    //                for (unsigned b = 1; b < i_max; ++b) {
    //                    rel1 = rel;
    //                    rel1.data.pop_back();
    //                    rel1.iaddP(O1, tmpm);
    //                    for (int j : two_expansion(b))
    //                        rel1.iaddP(pi_basis_d[j], tmpm);

    //                    if (!TryExtCof(iCof, rel1, deg, depth, flag)) {
    //                        ++count_pass;
    //                        if (count_pass > 1)
    //                            break;
    //                        rel_pass = std::move(rel1);
    //                    }
    //                }

    //                if (count_pass == 0) {
    //                    rel_pass = rel;
    //                    rel_pass.data.pop_back();
    //                    rel_pass.iaddP(O1, tmpm);
    //                    bNewExt = true;
    //                }
    //                else if (count_pass == 1) {
    //                    rel1 = rel;
    //                    rel1.data.pop_back();
    //                    rel1.iaddP(O1, tmpm);
    //                    if (TryExtCof(iCof, rel1, deg, depth, flag))
    //                        bNewExt = true;
    //                }
    //                if (bNewExt) {
    //                    rel_pass = ssCof.pi_gb.Reduce(std::move(rel_pass));
    //                    if (algZ::IsValidRel(rel_pass)) {
    //                        ++count_htpy;
    //                        Logger::PrintDepth();

    //                        Logger::LogHtpyRel(depth, enumReason::deduce, name, deg, rel, rel_pass);
    //                        AddPiRelsCof(iCof, {std::move(rel_pass)});
    //                        auto deg_min = deg - ssCof.deg_qt;
    //                        SyncHomotopy(deg_min, count_ss, count_htpy, depth);
    //                        count_htpy += DeduceTrivialExtensions(depth);
    //                        if (flag & DeduceFlag::pi_exact)
    //                            count_htpy += DeduceExtensionsByExactness(deg_min.stem(), stem_max_exactness_, depth);
    //                    }
    //                    else
    //                        Logger::ClearDepth();
    //                }
    //                else
    //                    Logger::ClearDepth();
    //            }
    //            i_start[iSS] = indices_size;
    //        }
    //    }

    //    bool all_at_end = true;
    //    for (size_t iSS = 0; iSS < i_start.size(); ++iSS) {
    //        if (iSS == 0) {
    //            if (ut::has(rings_.pi_gb.leads_group_by_t(), t) && i_start[iSS] != rings_.pi_gb.leads_group_by_t().at(t).size()) {
    //                all_at_end = false;
    //                break;
    //            }
    //        }
    //        else {
    //            size_t iCof = iSS - 1;
    //            auto& ssCof = modules_[iCof];
    //            if (ut::has(ssCof.pi_gb.leads_group_by_t(), t) && i_start[iSS] != ssCof.pi_gb.leads_group_by_t().at(t).size()) {
    //                all_at_end = false;
    //                break;
    //            }
    //        }
    //    }
    //    if (all_at_end) {
    //        ++t;
    //        for (auto& i : i_start)
    //            i = 0;
    //    }
    //}
}

int main_deduce_ext(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";
    int stem_min = 0, stem_max = 30;
    std::vector<std::string> strFlags;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"stem_min", &stem_min}, {"stem_max", &stem_max}, {"diagram", &diagram_name}, {"flags...", &strFlags}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    DeduceFlag flag = DeduceFlag::pi;
    for (auto& f : strFlags) {
        // if (f == "exact")
        //     flag = flag | DeduceFlag::pi_exact;
    }

    Diagram diagram(diagram_name, flag);

    try {
        int count_ss = 0, count_homotopy = 0;
        try {
            diagram.DeduceExtensions(stem_min, stem_max, count_ss, count_homotopy, 0, flag);
        }
        catch (TerminationException&) {
        }
        diagram.SimplifyPiRels();
        /*std::cout << "S0 pi_gb.size()=" << diagram.GetRings().pi_gb.data().size() << '\n';
        for (size_t iMod = 0; iMod < diagram.GetModules().size(); ++iMod)
            std::cout << diagram.GetModules()[iMod].name << " pi_gb.size()=" << diagram.GetModules()[iMod].pi_gb.data().size() << '\n';*/

        diagram.save(diagram_name, flag);
        Logger::LogSummary("Changed differentials", count_ss);
        Logger::LogSummary("Changed homotopy", count_homotopy);
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    return 0;
}

int main_deduce_ext_2tor(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    DeduceFlag flag = DeduceFlag::pi;
    Diagram diagram(diagram_name, flag);

    try {
        diagram.DeduceExtensions2tor();
        diagram.SimplifyPiRels();

        diagram.save(diagram_name, flag);
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    // bench::Counter::print();
    return 0;
}
