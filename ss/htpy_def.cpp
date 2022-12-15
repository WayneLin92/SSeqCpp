#include "algebras/linalg.h"
#include "main.h"

bool EndWithIndAndO(const algZ::Poly& p, uint32_t& gen_id)
{
    if (p.data.size() >= 2 && p.data.back().IsUnKnown()) {
        auto& m = p.data[p.data.size() - 2];
        if (m.c() == 0 && m.m0().size() + m.m1().size() == 1) {
            if ((m.m0().size() == 1 && m.m0().begin()->e() == 1) || m.m1().begin()->e() == 1) {
                gen_id = m.m0().size() ? m.m0().begin()->g() : m.m1().begin()->g();
                return true;
            }
        }
    }
    return false;
}

bool EndWithIndAndO(const algZ::Mod& p, uint32_t& v_id)
{
    if (p.data.size() >= 2 && p.data.back().IsUnKnown()) {
        auto& m = p.data[p.data.size() - 2];
        if (m.c() == 0 && !m.m) {
            v_id = m.v;
            return true;
        }
    }
    return false;
}

algZ::Poly DefRel(algZ::Poly p, const algZ::Groebner& gb)
{
    algZ::Poly tmp;
    p.data.pop_back();
    /*algZ::Poly gen = p.data.back();
    p.data.pop_back();
    p.isubP(gen, tmp, gb.gen_2tor_degs());
    p = gb.Reduce(std::move(p));*/
    return p;
}

algZ::Mod DefRel(algZ::Mod p, const algZ::GroebnerMod& gb)
{
    algZ::Mod tmp;
    p.data.pop_back();
    /*algZ::Mod gen = p.data.back();
    p.data.pop_back();
    p.isubP(gen, tmp, gb.gen_2tor_degs());
    p = gb.Reduce(std::move(p));*/
    return p;
}

int Diagram::DefineDependenceInExtensions(int depth)
{
    int count_homotopy = 0;

    /* multiplicative structures */
    { /* sphere */
        auto& pi_gb = ssS0_.pi_gb;
        if (ssS0_.pi_gen_defs.empty())
            ssS0_.pi_gen_defs.resize(ssS0_.pi_gb.gen_degs().size(), DefFlag::no_def);
        if (ssS0_.pi_gen_def_mons.empty())
            ssS0_.pi_gen_def_mons.resize(ssS0_.pi_gb.gen_degs().size());

        algZ::Poly1d new_rels;
        for (size_t i = 0; i < pi_gb.data().size(); ++i) {
            auto& rel = pi_gb.data()[i];
            uint32_t gen_id = -1;
            if (EndWithIndAndO(rel, gen_id)) {
                AdamsDeg deg = GetDeg(rel.GetLead(), pi_gb.gen_degs());
                if (ssS0_.pi_gen_defs[gen_id] == DefFlag::no_def) {
                    ssS0_.pi_gen_defs[gen_id] = DefFlag::dec;
                    algZ::Poly rel1 = DefRel(rel, pi_gb);

                    if (depth == 0)
                        myio::Logger::cout_fout2() << "By definition:  S0  stem=" << deg.stem() << "  " << rel << " --> " << rel1 << '\n';
                    new_rels.push_back(std::move(rel1));
                    ++count_homotopy;
                }
            }
        }
        AddPiRelsS0(std::move(new_rels));
    }
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
        auto& ssCof = ssCofs_[iCof];
        auto& pi_gb = ssCof.pi_gb;
        auto& basis_ss = ssCof.basis_ss;
        size_t iSS = iCof + 1;

        if (ssCof.pi_gen_defs.empty())
            ssCof.pi_gen_defs.resize(ssCof.pi_gb.v_degs().size(), DefFlag::no_def);
        if (ssCof.pi_gen_def_mons.empty())
            ssCof.pi_gen_def_mons.resize(ssCof.pi_gb.v_degs().size());

        algZ::Mod1d new_rels;
        for (size_t i = 0; i < pi_gb.data().size(); ++i) {
            auto& rel = pi_gb.data()[i];
            uint32_t v_id = -1;
            if (EndWithIndAndO(rel, v_id)) {
                AdamsDeg deg = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs(), pi_gb.v_degs());
                if (ssCof.pi_gen_defs[v_id] == DefFlag::no_def) {
                    ssCof.pi_gen_defs[v_id] = DefFlag::dec;
                    algZ::Mod rel1 = DefRel(rel, pi_gb);

                    if (depth == 0)
                        myio::Logger::cout_fout2() << "By definition:  " << ssCof.name << "  stem=" << deg.stem() << "  " << rel << " --> " << rel1 << '\n';
                    new_rels.push_back(std::move(rel1));
                    ++count_homotopy;
                }
            }
        }
        AddPiRelsCof(iCof, std::move(new_rels));
    }

    /* top cell maps */
    {
        algZ::Poly1d new_rels_S0; /* By exactness qi * h = 0 */
        int old_count_homotopy = count_homotopy;
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            auto& ssCof = ssCofs_[iCof];
            auto& q = ssCof.pi_qt.back();
            for (size_t i = 0; i < q.size(); ++i) {
                auto& qi = q[i];
                qi = ssS0_.pi_gb.Reduce(std::move(qi));
                uint32_t gen_id = -1;
                if (EndWithIndAndO(qi, gen_id)) {
                    AdamsDeg deg = GetDeg(qi.GetLead(), ssS0_.pi_gb.gen_degs());
                    if (ssS0_.pi_gen_defs[gen_id] == DefFlag::no_def) {
                        ssS0_.pi_gen_defs[gen_id] = DefFlag::dec;
                        algZ::Poly qi1 = DefRel(qi, ssS0_.pi_gb);

                        if (depth == 0)
                            myio::Logger::cout_fout2() << "By definition:  " << ssCof.name << "  stem=" << deg.stem() << "  q" << std::to_string(i) << "=" << qi << " --> " << qi1 << '\n';
                        qi = std::move(qi1);
                        ++count_homotopy;

                        algZ::Poly h = ssS0_.pi_gb.Gen((uint32_t)iCof);
                        algZ::Poly relS0 = ssS0_.pi_gb.Reduce(qi * h);
                        if (algZ::IsValidRel(relS0)) {
                            if (depth == 0)
                                myio::Logger::cout_fout2() << "    -> S0 rel: " << relS0 << '\n';
                            new_rels_S0.push_back(std::move(relS0));
                            ++count_homotopy;
                        }
                    }
                }
            }
            if (count_homotopy > old_count_homotopy) {
                AddPiRelsCof2S0(iCof);
                old_count_homotopy = count_homotopy;
            }
        }
        AddPiRelsS0(std::move(new_rels_S0));
    }

    return count_homotopy;
}

int main_deduce_extdef(int argc, char** argv, int index)
{
    std::string selector = "default";
    std::vector<std::string> strFlags;
    DeduceFlag flag = DeduceFlag::no_op;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce extdef [selector] [flags...]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    if (myio::load_args(argc, argv, ++index, "flags", strFlags))
        return index;
    auto dbnames = GetDbNames(selector);
    for (auto& f : strFlags) {
        if (f == "exact")
            flag = flag | DeduceFlag::check_exactness;
    }

    Diagram diagram(dbnames);

    try {
        int count_ss = 0, count_homotopy = 0;
        try {
            diagram.DefineDependenceInExtensions(0);
        }
        catch (TerminationException&) {
        }
        diagram.SimplifyPiRels();

        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            auto pi_table = GetComplexName(dbnames[k]);
            db.begin_transaction();

            db.drop_and_create_pi_relations(pi_table);
            db.drop_and_create_pi_basis(pi_table);
            db.drop_and_create_pi_definitions(pi_table);

            if (k == 0) {
                db.drop_and_create_pi_generators(pi_table);
                db.save_pi_generators(pi_table, diagram.GetS0().pi_gb.gen_degs(), diagram.GetS0().pi_gen_Einf);
                db.save_pi_gb(pi_table, diagram.GetS0().pi_gb.OutputForDatabase(), diagram.GetS0GbEinf());
                db.save_pi_basis(pi_table, diagram.GetS0().pi_basis.front());
                db.save_pi_def(pi_table, diagram.GetS0().pi_gen_defs, diagram.GetS0().pi_gen_def_mons);
            }
            else {
                db.drop_and_create_pi_generators_mod(pi_table);
                auto& Cof = diagram.GetCofs()[k - 1];
                if (Cof.pi_qt.size() != 1)
                    throw MyException(0x925afecU, "Not on the initial node");
                db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_qt.front());
                db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), diagram.GetCofGbEinf(int(k - 1)));
                db.save_pi_basis_mod(pi_table, Cof.pi_basis.front());
                db.save_pi_def_mod(pi_table, Cof.pi_gen_defs, Cof.pi_gen_def_mons);
            }

            db.end_transaction();
        }
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