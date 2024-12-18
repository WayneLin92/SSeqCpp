/* This module deals with operations that might affect all spectra */

#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

/* Add a node */
void Category::AddNode(SSFlag flag)
{
    ++depth_;
    for (auto& ring : rings_)
        ring.nodes_ss.AddNode();
    for (auto& mod : modules_)
        mod.nodes_ss.AddNode();
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iTri = 0; iTri < 3; ++iTri)
                cofseq.nodes_cofseq[iTri].AddNode();
    }

    if (flag & SSFlag::pi) {
        for (auto& ring : rings_) {
            ring.pi_gb.AddNode();
            ring.nodes_pi_basis.push_back({});
        }
        for (auto& mod : modules_) {
            mod.pi_gb.AddNode();
            mod.nodes_pi_basis.push_back({});
        }

        //// TODO: Add pi_maps
    }
}

/* Pop the lastest node */
void Category::PopNode(SSFlag flag)
{
    --depth_;
    for (auto& ring : rings_)
        ring.nodes_ss.pop_back();
    for (auto& mod : modules_)
        mod.nodes_ss.pop_back();
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iTri = 0; iTri < 3; ++iTri)
                cofseq.nodes_cofseq[iTri].pop_back();
    }

    if (flag & SSFlag::pi) {
        for (auto& ring : rings_) {
            ring.pi_gb.PopNode();
            ring.nodes_pi_basis.pop_back();
        }
        for (auto& mod : modules_) {
            mod.pi_gb.PopNode();
            mod.nodes_pi_basis.pop_back();
        }
        //// TODO: Add pi_maps
        // UpdateAllPossEinf();
    }
}

int1d MapRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = category.GetRings();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg)
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_x));
    }
    return result;
}

int1d MapMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = category.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Mod::Verify(const Category& category, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = category.GetModules();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (AdamsDeg deg_x : mods[from.index].degs_ss) {
        auto& basis_d = mods[from.index].basis.at(deg_x);
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, category);
                int1d fx = map(int1d{(int)i}, deg_x, category);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_gfx = mods[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Mod2Indices(alg_gfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw ErrorIdMsg(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

auto subs(const Mod& x, const Poly1d& map1, const Mod1d& map2)
{
    Mod result{}, tmp{};
    Poly tmp_prod_p, tmf_p;
    for (const MMod& m : x.data) {
        Poly fm = Poly::Unit();
        for (auto p = m.m.begin(); p != m.m.end(); ++p) {
            powP(map1[p->g()], p->e_masked(), tmp_prod_p, tmf_p);
            fm.imulP(tmp_prod_p, tmf_p);
        }
        result.iaddP(fm * map2[m.v], tmp);
    }
    return result;
}

int1d MapMod2ModV2::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = category.GetModules();
        auto& maps = category.GetMaps();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2ModV2::Verify(const Category& category, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = category.GetModules();
    auto& maps = category.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (AdamsDeg deg_x : mods[from.index].degs_ss) {
        auto& basis_d = mods[from.index].basis.at(deg_x);
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, category);
                int1d fx = map(int1d{(int)i}, deg_x, category);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_fgfx = mods[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Mod2Indices(alg_fgfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw ErrorIdMsg(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMod2Ring::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = category.GetRings();
        auto& mods = category.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Ring::Verify(const Category& category, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = category.GetModules();
    auto& rings = category.GetRings();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (AdamsDeg deg_x : mods[from.index].degs_ss) {
        auto& basis_d = mods[from.index].basis.at(deg_x);
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, category);
                int1d fx = map(int1d{(int)i}, deg_x, category);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_gfx = rings[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Poly2Indices(alg_gfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw ErrorIdMsg(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

auto subs(const Mod& x, const Poly1d& map1, const Poly1d& map2)
{
    Poly result, tmp_prod, tmp;
    for (const MMod& m : x.data) {
        Poly fm = Poly::Unit();
        for (auto p = m.m.begin(); p != m.m.end(); ++p) {
            powP(map1[p->g()], p->e_masked(), tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        result.iaddP(fm * map2[m.v], tmp);
    }
    return result;
}

int1d MapMod2RingV2::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = category.GetRings();
        auto& mods = category.GetModules();
        auto& maps = category.GetMaps();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2RingV2::Verify(const Category& category, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = category.GetModules();
    auto& rings = category.GetRings();
    auto& maps = category.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (AdamsDeg deg_x : mods[from.index].degs_ss) {
        auto& basis_d = mods[from.index].basis.at(deg_x);
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, category);
                int1d fx = map(int1d{(int)i}, deg_x, category);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_fgfx = rings[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Poly2Indices(alg_fgfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw ErrorIdMsg(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMulRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = category.GetRings();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulRing2Mod::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = category.GetRings();
        auto& mods = category.GetModules();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Category& category) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = category.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(factor * x_alg);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

SSRet Category::SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    SSRet rt;
    auto& ring = rings_[iRing];
    auto& nodes_ss = ring.nodes_ss;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        int depth = depth_;
        AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
        if (x.empty()) {
            if (rt += SetRingBoundaryLeibniz(iRing, deg_dx, dx, r - 1, flag))
                return rt;
        }
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, ring.t_max, deg_x, r + 1) - 1;
            if (rt += SetRingDiffLeibniz(iRing, deg_x, x, dx, r, r_min, flag))
                return rt;
        }

        for (size_t iMap : ring.ind_maps) {
            auto map = (MapRing2Ring*)maps_[iMap].get();
            if (deg_x.t <= map->t_max) {
                int1d fdx;
                if (deg_dx.t <= map->t_max)
                    fdx = map->map(dx, deg_dx, *this);
                else if (!dx.empty())
                    continue;
                auto fx = map->map(x, deg_x, *this);
                if (IsNewDiff(rings_[map->to.index].nodes_ss, deg_x, fx, fdx, r)) {
                    Logger::LogDiff(depth, EnumReason::nat, rings_[map->to.index].name, deg_x, r, fx, fdx, map->name, flag);
                    if (rt += SetCwDiffGlobal(map->to, deg_x, fx, fdx, r, true, flag | SSFlag::stacked)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Consider the map `{}`.\n{}", rings_[iRing].name, deg_x, r, myio::Serialize(x), myio::Serialize(dx), map->name, rt.err_msg);
                        return rt;
                    }
                }
            }
        }

        /* x^2 are cycles in E_r */
        auto deg_xx = deg_x * 2;
        if (deg_xx.t <= ring.t_max && !x.empty() && dx.empty() && r < R_PERM - 1) {
            Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
            Poly poly_xx;
            poly_x.frobP(poly_xx);
            poly_xx = ring.gb.Reduce(std::move(poly_xx));
            if (poly_xx) {
                int1d xx = Poly2Indices(poly_xx, ring.basis.at(deg_xx));
                if (IsNewDiff(nodes_ss, deg_xx, xx, {}, r + 1)) {
                    Logger::LogDiff(depth, EnumReason::deduce_xx, ring.name, deg_xx, r + 1, xx, {}, "", flag);
                    if (rt += SetRingDiffGlobal(iRing, deg_xx, xx, {}, r + 1, true, flag | SSFlag::stacked)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Consider the square of `{} {} [{}]`", rings_[iRing].name, deg_x, myio::Serialize(x));
                        return rt;
                    }
                }
            }
        }
    }
    if ((flag & SSFlag::cofseq) && !(flag & SSFlag::stacked)) {
        if (rt += SyncToCofseq(flag))
            return rt;
    }
    return rt;
}

SSRet Category::SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    SSRet rt;
    auto& mod = modules_[iMod];
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (x.empty()) {
            if (rt += SetModuleBoundaryLeibniz(iMod, deg_dx, dx, r - 1, flag))
                return rt;
        }
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, t_max, deg_x, r + 1) - 1;
            if (rt += SetModuleDiffLeibniz(iMod, deg_x, x, dx, r, r_min, flag))
                return rt;
        }

        for (size_t iMap : mod.ind_maps) {
            auto& map = maps_[iMap];
            if (deg_x.t <= map->t_max) {
                int1d fdx;
                if (deg_dx.t <= map->t_max)
                    fdx = map->map(dx, deg_dx, *this);
                else if (!dx.empty())
                    continue;
                auto fx = map->map(x, deg_x, *this);
                AdamsDeg deg_fx = deg_x + map->deg;
                if (IsNewDiff(GetNodesSS(map->to), deg_fx, fx, fdx, r)) {
                    Logger::LogDiff(depth_, EnumReason::nat, GetCwName(map->to), deg_fx, r, fx, fdx, map->name, flag);
                    if (rt += SetCwDiffGlobal(map->to, deg_fx, fx, fdx, r, true, flag | SSFlag::stacked)) {
                        if (flag & SSFlag::log_proof)
                            rt.err_msg = fmt::format("Get `{} {} d_{}[{}]=[{}]`. Consider the map `{}`.\n{}", mod.name, deg_x, r, myio::Serialize(x), myio::Serialize(dx), map->name, rt.err_msg);
                        return rt;
                    }
                }
            }
        }
    }
    if ((flag & SSFlag::cofseq) && !(flag & SSFlag::stacked)) {
        if (rt += SyncToCofseq(flag))
            return rt;
    }
    return rt;
}

SSRet Category::SetCwDiffGlobal(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)  //// TODO: write two functions into one
{
    SSRet rt;
    if (iCw.isRing())
        rt += SetRingDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);
    else
        rt += SetModuleDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);
    if (rt)
        return rt;
    if (depth_ <= 1 && (flag & SSFlag::synthetic) && r <= 6 && x.size() && dx != NULL_DIFF && dx.size()) {
        bool has_cross = GetCrossR(GetNodesSS(iCw), deg_x, GetTMax(iCw), r) <= r;
        if (rt += SetCwDiffSynthetic(iCw, deg_x, x, dx, r, has_cross, flag))
            return rt;
    }
    return rt;
}

/* */
SSRet Category::SetCwDiffSynthetic(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool hasCross, SSFlag flag)
{
    SSRet rt;
    if (dx.empty())
        return rt;
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    bool printed_x = false;
    auto& name = GetCwName(iCw);
    int r_fx_tmp = -1;
    int1d dfx_tmp;
    std::string proof;
    for (IndexUniv iCof : GetIndexCofs(iCw)) {
        auto& cof = cofseqs_[iCof.index];
        size_t iTri_fx = NextiTri(iCof.iTri);
        auto& nodes_ss_fx = *cof.nodes_ss[iTri_fx];
        int t_max_fx = cof.t_max[iTri_fx];
        auto iCw_fx = cof.indexCw[iTri_fx];
        auto& f = maps_[cof.indexMap[iCof.iTri]];
        AdamsDeg deg_fx;
        int1d fx;
        if (GetSynImage(iCof, deg_x, x, LEVEL_MAX - r, deg_fx, fx, -1, hasCross ? 0 : R_PERM) == 0) {
            GetRAndDiff(nodes_ss_fx, deg_fx, fx, r_fx_tmp, dfx_tmp);
            int cross_dfx_tmp = GetCrossR(nodes_ss_fx, deg_fx, t_max_fx, deg_fx.s - deg_x.s - f->deg.s + 2);
            if (r_fx_tmp >= cross_dfx_tmp) {
                r_fx_tmp = cross_dfx_tmp - 1;
                dfx_tmp.clear();
            }
            int cross_fdx_min = deg_fx.s + r_fx_tmp + (dfx_tmp.empty() ? 1 : 0);
            /*if (GetCwName(iUniv) == "S0" && deg_x == AdamsDeg(9, 124 + 9) && x == int1d{1, 2} && iCof.iTri == 2) {
                fmt::print("cross_dfx_tmp={}, cross_fdx_min={}, r_fx_tmp={}, dfx_tmp={}\n", cross_dfx_tmp, cross_fdx_min, r_fx_tmp, myio::Serialize(dfx_tmp));
                 fmt::print("debug\n");
            }*/

            AdamsDeg deg_dfx;
            int1d dfx;
            if (!dx.empty()) {
                if (GetSynImage(iCof, deg_dx, dx, r, deg_dfx, dfx, deg_fx.s, cross_fdx_min) == 0) {
                    if (deg_dfx.s < cross_fdx_min)
                        continue;
                    int r_fx = deg_dfx.s - deg_fx.s;
                    if (IsNewDiff(GetNodesSS(iCw_fx), deg_fx, fx, dfx, r_fx)) {
                        if (depth_ == 0 && !printed_x) {
                            printed_x = true;
                            Logger::LogGreyDiff(depth_, name, deg_x, r, x, dx, flag);
                        }
                        auto map_name = maps_[cofseqs_[iCof.index].indexMap[iCof.iTri]]->name;
                        if (depth_ == 0 && flag & SSFlag::log_proof) {
                            proof.clear();
                            proof = fmt::format("Apply General Leibniz with map `{}` and k={} where the Adams differential has{} crossing", map_name, r_fx_tmp + (dfx_tmp.empty() ? 1 : 0), hasCross ? "" : " no");
                        }
                        Logger::LogDiff(depth_, EnumReason::syn, fmt::format("{}", GetCwName(iCw_fx)), deg_fx, r_fx, fx, dfx, proof, flag);
                        if (rt += SetCwDiffGlobal(iCw_fx, deg_fx, fx, dfx, r_fx, true, flag)) {
                            if (flag & SSFlag::log_proof)
                                rt.err_msg = fmt::format("Apply General Leibniz with map `{}` and k={} where the Adams differential has{} crossing.\n{}", map_name, r_fx_tmp + (dfx_tmp.empty() ? 1 : 0), hasCross ? "" : " no", rt.err_msg);
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

int Category::GetSynImage(IndexUniv iCof, AdamsDeg deg_x, const int1d& x, int level_x, AdamsDeg& deg_fx, int1d& fx, int s_f_dinv_x, int cross_min)
{
    auto& cof = cofseqs_[iCof.index];
    auto iTri = size_t(iCof.iTri);
    const size_t iTri_next = NextiTri(iTri);
    const size_t iTri_prev = PreviTri(iTri);
    auto& f = maps_[cof.indexMap[iTri]];
    auto& f_next = maps_[cof.indexMap[iTri_next]];
    auto& f_prev = maps_[cof.indexMap[iTri_prev]];

    if (deg_x.t > f->t_max)
        return 1;
    fx = f->map(x, deg_x, *this);

    AdamsDeg deg_xtop, deg_dxtop;
    int1d xtop, dxtop;
    int level_xtop = -1, r_xtop = -1;
    auto& nodes_ss_prev = *cof.nodes_ss[iTri_prev];
    if (fx.size()) {
        deg_fx = deg_x + f->deg;
        auto [dinv_fx, level_fx] = GetDiffAndLevel(*cof.nodes_ss[iTri_next], deg_fx, fx);
        int r_new = deg_fx.s - s_f_dinv_x;
        if (level_fx >= r_new)
            return 0; /* no hidden extension */

        /* compute xtop
         * dinv_fx --f_next--> xtop
         */
        if (dinv_fx == NULL_DIFF)
            return 2;
        if (deg_fx.t > f_next->t_max)
            return 3;
        auto deg_dinv_fx = deg_fx - AdamsDeg(level_fx, level_fx - 1);
        xtop = f_next->map(dinv_fx, deg_dinv_fx, *this);
        if (xtop.empty())
            return 4;
        deg_xtop = deg_dinv_fx + f_next->deg;
        std::tie(dxtop, level_xtop) = GetDiffAndLevel(nodes_ss_prev, deg_xtop, xtop);
        if (dxtop == NULL_DIFF)
            return 5;
        if (level_xtop <= level_x - (deg_x.s - deg_xtop.s))
            return 6;

        AdamsDeg deg_x1;
        int1d x1;
        if (GetSynImage(iCof.prev(), deg_xtop, xtop, level_xtop, deg_x1, x1, -1, 0) != 0)
            return 7; /* failed to confirm that xtop extends to x */
        if (deg_x1 != deg_x)
            return 8; /* failed to confirm that xtop extends to x */
        auto& nodes_ss = *cof.nodes_ss[iTri];
        auto& sc_x = nodes_ss.GetRecentSc(deg_x);
        size_t iFirst_level_x = GetFirstIndexOnLevel(sc_x, level_x);
        if (lina::Residue(sc_x.basis.begin(), sc_x.basis.begin() + iFirst_level_x, lina::add(x, x1)).size())
            return 9;
        r_xtop = LEVEL_MAX - level_xtop;
    }

    if (xtop.empty()) {
        /* compute xtop
         * xtop --f_prev--> x
         */
        deg_xtop = deg_x - f_prev->deg;
        if (deg_xtop.t > f_prev->t_max)
            return 10;
        auto& sc_xtop = nodes_ss_prev.GetRecentSc(deg_xtop);
        int2d l_x, l_fx, l_domain, l_f, image, g, kernel;
        for (size_t i = 0; i < sc_xtop.basis.size(); ++i) { /* Make sure that level(xtop) is as small as possible */
            l_x.push_back(sc_xtop.basis[i]);
            l_fx.push_back(f_prev->map(sc_xtop.basis[i], deg_xtop, *this));
        }
        lina::SetLinearMapV3(l_x, l_fx, l_domain, l_f, image, g, kernel);
        xtop = lina::GetImage(image, g, x);

        /* compute dxtop
         * xtop --d--> dxtop
         */
        std::tie(dxtop, level_xtop) = GetDiffAndLevel(nodes_ss_prev, deg_xtop, xtop);
        r_xtop = LEVEL_MAX - level_xtop;
        if (dxtop == NULL_DIFF)
            return 11;
        if (level_xtop <= level_x)
            return 12;
    }
    deg_dxtop = deg_xtop + AdamsDeg(r_xtop, r_xtop - 1);
    if (deg_dxtop.t > f_prev->t_max)
        return 13;

    /* compute fx
     * fx --f_next--> dxtop
     */
    int2d l_fx, image, kernel, g;
    deg_fx = deg_dxtop - f_next->deg;
    if (deg_fx.t > f_next->t_max)
        return 14;
    auto& nodes_ss_next = *cof.nodes_ss[iTri_next];
    auto& sc_fx = nodes_ss_next.front().has(deg_fx) ? nodes_ss_next.GetRecentSc(deg_fx) : EmptyStaircase;
    auto& sc_dxtop = dxtop.size() ? nodes_ss_prev.GetRecentSc(deg_dxtop) : EmptyStaircase;
    size_t iFirst_level_dxtop = GetFirstIndexOnLevel(sc_dxtop, r_xtop);
    if (level_x <= LEVEL_PERM) {
        int2d l_x, l_domain, l_f;
        size_t iLast = sc_fx.basis.size() ? GetFirstIndexOfFixedLevels(nodes_ss_next, deg_fx, LEVEL_PERM + 1) : 0;
        for (size_t i = 0; i < iLast; ++i) {
            l_x.push_back(sc_fx.basis[i]);
            auto l_fxi = f_next->map(sc_fx.basis[i], deg_fx, *this);
            l_fxi = lina::Residue(sc_dxtop.basis.begin(), sc_dxtop.basis.begin() + iFirst_level_dxtop, l_fxi);
            l_fx.push_back(std::move(l_fxi));
        }
        lina::SetLinearMapV3(l_x, l_fx, l_domain, l_f, image, g, kernel);
        int r_new = deg_fx.s - s_f_dinv_x;
        size_t iFirst = GetFirstIndexOnLevel(sc_fx, r_new);
        for (auto& k : kernel) {
            if (lina::Residue(sc_fx.basis.begin(), sc_fx.basis.begin() + iFirst, k).size())
                return 15; /* For permanent or boundary elements the hidden extension must be determined uniquely modulo d_1,...,d_{r_new} */
        }
        if (lina::Residue(image, dxtop).size())  //// cound let fx = {} or this case never happens
            return 16;
    }
    else {
        for (size_t i = 0; i < sc_fx.basis.size(); ++i) {
            auto l_fxi = f_next->map({(int)i}, deg_fx, *this);
            l_fxi = lina::Residue(sc_dxtop.basis.begin(), sc_dxtop.basis.begin() + iFirst_level_dxtop, l_fxi);
            l_fx.push_back(std::move(l_fxi));
        }
        lina::SetLinearMap(l_fx, image, kernel, g);

        if (deg_fx.s >= cross_min && kernel.size())
            return 17;
    }

    if (lina::Residue(image, dxtop).size())
        return 18; /* dxtop does not come from bottom */

    if (deg_fx.s >= cross_min) { /* Make sure that there is no crossing extensions */
        auto& nodes_ss = *cof.nodes_ss[iTri];
        AdamsDeg deg_fx1;
        int1d fx1;
        for (int s1 = deg_x.s + 1; s1 < deg_fx.s - f->deg.s; ++s1) {
            auto deg_x1 = AdamsDeg(s1, deg_x.stem() + s1);
            if (!nodes_ss.front().has(deg_x1))
                continue;
            auto& sc_x1 = nodes_ss.GetRecentSc(deg_x1);
            size_t iFirst_x1 = 0, iLast_x1 = sc_x1.basis.size();
            /*if (level_x <= LEVEL_PERM && sc_x1.basis.size() && deg_fx.s == cross_min) {
                iFirst_x1 = GetFirstIndexOnLevel(sc_x1, s1 - deg_x.s + 2);
                iLast_x1 = GetFirstIndexOfFixedLevels(nodes_ss, deg_x1, LEVEL_PERM + 1);
            }*/

            for (size_t i = iFirst_x1; i < iLast_x1; ++i) {
                if (sc_x1.levels[i] < level_x)
                    continue;
                if (int error = GetSynImage(iCof, deg_x1, sc_x1.basis[i], sc_x1.levels[i], deg_fx1, fx1, s_f_dinv_x, R_PERM); error == 0) {
                    if (deg_fx1.s <= deg_fx.s && deg_fx1.s >= cross_min) {
                        if (level_x <= LEVEL_PERM) {
                            auto& sc_fx1 = nodes_ss_next.GetRecentSc(deg_fx1);
                            size_t iLast_fx1 = GetFirstIndexOfFixedLevels(nodes_ss_next, deg_fx1, LEVEL_PERM + 1);
                            if (fx1.size() && lina::Residue(sc_fx1.basis.begin(), sc_fx1.basis.begin() + iLast_fx1, fx1).empty()) {
                                return 19;
                            }
                        }
                        else
                            return 20;
                    }
                }
                else if (error != 12)  ////
                    return 100 + error;
            }
        }
    }
    fx = lina::GetImage(image, g, dxtop);
    return 0;
}
