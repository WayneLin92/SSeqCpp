/* This module deals with operations that might affect all spectra */

#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

/* Add a node */
void Diagram::AddNode(SSFlag flag)
{
    ++depth_;
    for (auto& ring : rings_)
        ring.nodes_ss.push_back({});
    for (auto& mod : modules_)
        mod.nodes_ss.push_back({});
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iCs = 0; iCs < 3; ++iCs)
                cofseq.nodes_cofseq[iCs].push_back({});
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
void Diagram::PopNode(SSFlag flag)
{
    --depth_;
    for (auto& ring : rings_)
        ring.nodes_ss.pop_back();
    for (auto& mod : modules_)
        mod.nodes_ss.pop_back();
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iCs = 0; iCs < 3; ++iCs)
                cofseq.nodes_cofseq[iCs].pop_back();
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

int1d MapRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg)
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_x));
    }
    return result;
}

int1d MapMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Mod::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_gfx = mods[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Mod2Indices(alg_gfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
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

int1d MapMod2ModV2::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto& maps = diagram.GetMaps();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2ModV2::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& maps = diagram.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_fgfx = mods[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Mod2Indices(alg_fgfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMod2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Ring::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& rings = diagram.GetRings();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_gfx = rings[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Poly2Indices(alg_gfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
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

int1d MapMod2RingV2::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto& maps = diagram.GetMaps();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2RingV2::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& rings = diagram.GetRings();
    auto& maps = diagram.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_fgfx = rings[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Poly2Indices(alg_fgfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMulRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulRing2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(factor * x_alg);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int Diagram::SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    auto& nodes_ss = ring.nodes_ss;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        int depth = depth_;
        AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
        if (x.empty())
            count += SetRingBoundaryLeibniz(iRing, deg_dx, dx, r - 1, flag);
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, ring.t_max, deg_x, r + 1) - 1;
            count += SetRingDiffLeibniz(iRing, deg_x, x, dx, r, r_min, flag);
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
                    Logger::LogDiff(depth, EnumReason::nat, fmt::format("({}) {}", map->display, rings_[map->to.index].name), deg_x, fx, fdx, r);
                    count += SetCwDiffGlobal(map->to, deg_x, fx, fdx, r, true, flag);
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
                    Logger::LogDiff(depth, EnumReason::deduce_xx, ring.name, deg_xx, xx, {}, r + 1);
                    count += SetRingDiffGlobal(iRing, deg_xx, xx, {}, r + 1, true, flag);
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    int count = 0;
    auto& mod = modules_[iMod];
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (x.empty())
            count += SetModuleBoundaryLeibniz(iMod, deg_dx, dx, r - 1, flag);
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, t_max, deg_x, r + 1) - 1;
            count += SetModuleDiffLeibniz(iMod, deg_x, x, dx, r, r_min, flag);
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
                if (IsNewDiff(GetSS(map->to), deg_fx, fx, fdx, r)) {
                    Logger::LogDiff(depth_, EnumReason::nat, fmt::format("({}) {}", map->display, GetCwName(map->to)), deg_fx, fx, fdx, r);
                    count += SetCwDiffGlobal(map->to, deg_fx, fx, fdx, r, true, flag);
                }
            }
        }
    }
    return count;
}

int Diagram::SetCwDiffGlobal(IndexCw iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)  //// TODO: write two functions into one
{
    int count = 0;
    if (iCw.isRing)
        count += SetRingDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);
    else
        count += SetModuleDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);

    if (depth_ <= 1 && (flag & SSFlag::synthetic) && r <= 6 && x.size() && dx != NULL_DIFF && dx.size()) {
        bool has_cross = GetCrossR(GetSS(iCw), deg_x, GetTMax(iCw), r) <= r;
        SetCwDiffSynthetic(iCw, deg_x, x, dx, r, has_cross, flag);
    }
    return count;
}

int Diagram::SetCwDiffSynthetic(IndexCw iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool hasCross, SSFlag flag)
{
    int count = 0;
    if (dx.empty())
        return count;
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    bool printed_x = false;
    auto& name = GetCwName(iCw);
    int r_fx_tmp = -1;
    int1d dfx_tmp;
    for (IndexCof iCof : GetIndexCof(iCw)) {
        auto& cof = cofseqs_[iCof.iCof];
        size_t iCs_fx = (iCof.iCs + 1) % 3;
        auto& nodes_ss_fx = *cof.nodes_ss[iCs_fx];
        int t_max_fx = cof.t_max[iCs_fx];
        auto iCw_fx = cof.indexCw[iCs_fx];
        auto& f = maps_[cof.indexMap[iCof.iCs]];
        AdamsDeg deg_fx;
        int1d fx;
        if (GetSynImage(iCof, deg_x, x, LEVEL_MAX - r, deg_fx, fx, -1, hasCross ? 0 : R_PERM) == 0) {
            GetRAndDiff(nodes_ss_fx, deg_fx, fx, r_fx_tmp, dfx_tmp);
            int cross_dfx_tmp = GetCrossR(nodes_ss_fx, deg_fx, t_max_fx, deg_fx.s - deg_x.s - f->deg.s + 2);
            if (r_fx_tmp >= cross_dfx_tmp) {
                r_fx_tmp = cross_dfx_tmp - 1;
                dfx_tmp.clear();
            }
            int cross_fd_min = deg_fx.s + r_fx_tmp + (dfx_tmp.empty() ? 1 : 0);
            /*if (GetCwName(iCw) == "S0" && deg_x == AdamsDeg(9, 124 + 9) && x == int1d{1, 2} && iCof.iCs == 2) {
                fmt::print("cross_dfx_tmp={}, cross_fd_min={}, r_fx_tmp={}, dfx_tmp={}\n", cross_dfx_tmp, cross_fd_min, r_fx_tmp, myio::Serialize(dfx_tmp));
                 fmt::print("debug\n");
            }*/

            AdamsDeg deg_dfx;
            int1d dfx;
            if (!dx.empty()) {
                if (GetSynImage(iCof, deg_dx, dx, r, deg_dfx, dfx, deg_fx.s, cross_fd_min) == 0) {
                    if (deg_dfx.s < cross_fd_min)
                        continue;
                    int r_fx = deg_dfx.s - deg_fx.s;
                    if (IsNewDiff(GetSS(iCw_fx), deg_fx, fx, dfx, r_fx)) {
                        if (!printed_x) {
                            printed_x = true;
                            Logger::LogDiff(depth_, name, deg_x, x, dx, r);
                        }
                        auto map_name = maps_[cofseqs_[iCof.iCof].indexMap[iCof.iCs]]->display;
                        Logger::LogDiff(depth_, EnumReason::synnat, fmt::format("({}) {}", map_name, GetCwName(iCw_fx)), deg_fx, fx, dfx, r_fx);
                        SetCwDiffGlobal(iCw_fx, deg_fx, fx, dfx, r_fx, true, flag);
                        ++count;
                    }
                }
            }
        }
    }
    return count;
}

int Diagram::GetSynImage(IndexCof iCof, AdamsDeg deg_x, const int1d& x, int level_x, AdamsDeg& deg_fx, int1d& fx, int s_f_dinv_x, int cross_min)
{
    auto& cof = cofseqs_[iCof.iCof];
    auto iCs = size_t(iCof.iCs);
    auto iCs_next = (iCs + 1) % 3;
    auto iCs_prev = (iCs + 2) % 3;
    auto& f = maps_[cof.indexMap[iCs]];
    auto& f_next = maps_[cof.indexMap[iCs_next]];
    auto& f_prev = maps_[cof.indexMap[iCs_prev]];

    if (deg_x.t > f->t_max)
        return 1;
    fx = f->map(x, deg_x, *this);

    AdamsDeg deg_xtop, deg_dxtop;
    int1d xtop, dxtop;
    int level_xtop = -1, r_xtop = -1;
    auto& nodes_ss_prev = *cof.nodes_ss[iCs_prev];
    if (fx.size()) {
        deg_fx = deg_x + f->deg;
        int level_fx;
        auto dinv_fx = GetLevelAndDiff(*cof.nodes_ss[iCs_next], deg_fx, fx, level_fx);
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
        dxtop = GetLevelAndDiff(nodes_ss_prev, deg_xtop, xtop, level_xtop);
        if (dxtop == NULL_DIFF)
            return 5;
        if (level_xtop <= level_x - (deg_x.s - deg_xtop.s))
            return 6;

        AdamsDeg deg_x1;
        int1d x1;
        if (GetSynImage(IndexCof{iCof.iCof, iCs_prev}, deg_xtop, xtop, level_xtop, deg_x1, x1, -1, 0) != 0)
            return 7; /* failed to confirm that xtop extends to x */
        if (deg_x1 != deg_x)
            return 8; /* failed to confirm that xtop extends to x */
        auto& nodes_ss = *cof.nodes_ss[iCs];
        auto& sc_x = ut::GetRecentValue(nodes_ss, deg_x);
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
        auto& sc_xtop = ut::GetRecentValue(nodes_ss_prev, deg_xtop);
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
        dxtop = GetLevelAndDiff(nodes_ss_prev, deg_xtop, xtop, level_xtop);
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
    Staircase sc_empty;
    auto& nodes_ss_next = *cof.nodes_ss[iCs_next];
    auto& sc_fx = ut::has(nodes_ss_next.front(), deg_fx) ? ut::GetRecentValue(nodes_ss_next, deg_fx) : sc_empty;
    auto& sc_dxtop = dxtop.size() ? ut::GetRecentValue(nodes_ss_prev, deg_dxtop) : sc_empty;
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
        auto& nodes_ss = *cof.nodes_ss[iCs];
        AdamsDeg deg_fx1;
        int1d fx1;
        for (int s1 = deg_x.s + 1; s1 < deg_fx.s - f->deg.s; ++s1) {
            auto deg_x1 = AdamsDeg(s1, deg_x.stem() + s1);
            if (!ut::has(nodes_ss.front(), deg_x1))
                continue;
            auto& sc_x1 = ut::GetRecentValue(nodes_ss, deg_x1);
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
                            auto& sc_fx1 = ut::GetRecentValue(nodes_ss_next, deg_fx1);
                            size_t iLast_fx1 = GetFirstIndexOfFixedLevels(nodes_ss_next, deg_fx1, LEVEL_PERM + 1);
                            if (fx1.size() && lina::Residue(sc_fx1.basis.begin(), sc_fx1.basis.begin() + iLast_fx1, fx1).empty()) {
                                return 19;
                            }
                        }
                        else
                            return 20;
                    }
                }
                else if (error != 12) ////
                    return 100 + error;
            }
        }
    }
    fx = lina::GetImage(image, g, dxtop);
    return 0;
}
