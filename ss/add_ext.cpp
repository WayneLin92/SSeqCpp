#include "main.h"

int main_add_ext(int argc, char** argv, int index)
{
    std::string diagram_name = "debug";
    std::string cw;
    std::string strRel, strExt;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Add extension\n";
        std::cout << "Usage:\n  ss add_ext <cw> <rel> <ext> [diagram]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  diagram = " << diagram_name << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return index;
    if (myio::load_arg(argc, argv, ++index, "rel", strRel))
        return index;
    if (myio::load_arg(argc, argv, ++index, "ext", strExt))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "diagram", diagram_name))
        return index;

    Diagram diagram(diagram_name, DeduceFlag::homotopy);

    /*int1d arr_rel = myio::Deserialize<int1d>(strRel);
    algZ::Poly tmp;
    algZ::Mod tmpm;
    if (iSS == 0) {
        if (arr_rel.empty() || arr_rel.size() % 2 == 1) {
            std::cout << "Invalid input rel\n";
            return 1;
        }
        int stem = 0;
        algZ::Poly rel = algZ::Poly::Unit();
        for (size_t i = 0; i < arr_rel.size(); i += 2) {
            auto power = diagram.GetRings().pi_gb.Gen((uint32_t)arr_rel[i], (uint32_t)arr_rel[i + 1]);
            rel = rel * power;
            stem += (diagram.GetRings().pi_gb.gen_degs()[arr_rel[i]] * arr_rel[i + 1]).stem();
        }
        algZ::Poly ext;
        if (!strExt.empty()) {
            int1d arr_ext = myio::Deserialize<int1d>(strExt);
            ext = algZ::Poly::Unit();
            int stem_ext = 0;
            for (size_t i = 0; i < arr_ext.size(); i += 2) {
                auto power = diagram.GetRings().pi_gb.Gen((uint32_t)arr_ext[i], (uint32_t)arr_ext[i + 1]);
                ext = ext * power;
                stem_ext += (diagram.GetRings().pi_gb.gen_degs()[arr_ext[i]] * arr_ext[i + 1]).stem();
            }
            if (stem_ext != stem) {
                std::cout << "rel and ext not in the same stem!\n";
                return 2;
            }
        }
        auto rel_reduced = diagram.GetRings().pi_gb.ReduceV2(rel);
        if (rel_reduced && rel_reduced.data.back().IsUnKnown()) {
            int s = rel_reduced.UnknownFil();
            if (ext && s != ext.GetLead().fil()) {
                std::cout << "ext is in incorrect filtration: " << ext.GetLead().fil() << "!=" << s << '\n';
                return 3;
            }
            algZ::Poly rel_reduced1 = rel_reduced;
            rel_reduced1.data.pop_back();
            rel_reduced1 += ext;
            rel_reduced1 += algZ::Poly::O(s + 1);
            algZ::Poly rel1 = rel;
            rel1.isubP(rel_reduced1, tmp, diagram.GetRings().pi_gb.gen_2tor_degs());
            
            diagram.ExtendRelRing(stem, rel1);
            std::cout << "S0  stem=" << stem << "  " << rel << '=' << rel_reduced << " --> " << rel_reduced1 << '\n';
            std::cout << rel1 << "=0\n";
            if (myio::UserConfirm()) {
                diagram.AddPiRelsRing({std::move(rel1)});
            }
            else {
                std::cout << "Aborted.\n";
                return 4;
            }
        }
        else {
            std::cout << "There is no extension to solve\n";
            return 5;
        }
    }
    else {
        auto& ssCof = diagram.GetModules()[iSS - 1];
        if (arr_rel.empty() || arr_rel.size() % 2 == 0) {
            std::cout << "Invalid input rel\n";
            return 6;
        }
        int stem = 0;
        algZ::Poly p = algZ::Poly::Unit();
        for (size_t i = 0; i + 2 < arr_rel.size(); i += 2) {
            auto power = diagram.GetRings().pi_gb.Gen((uint32_t)arr_rel[i], (uint32_t)arr_rel[i + 1]);
            stem += (diagram.GetRings().pi_gb.gen_degs()[arr_rel[i]] * arr_rel[i + 1]).stem();
            p = p * power;
        }
        algZ::Mod rel = p * ssCof.pi_gb.Gen(arr_rel.back());
        stem += ssCof.pi_gb.v_degs()[arr_rel.back()].stem();

        algZ::Mod ext;
        if (!strExt.empty()) {
            int1d arr_ext = myio::Deserialize<int1d>(strExt);
            algZ::Poly pext;
            pext = algZ::Poly::Unit();
            int stem_ext = 0;
            for (size_t i = 0; i + 2 < arr_ext.size(); i += 2) {
                auto power = diagram.GetRings().pi_gb.Gen((uint32_t)arr_ext[i], (uint32_t)arr_ext[i + 1]);
                pext = pext * power;
                stem_ext += (diagram.GetRings().pi_gb.gen_degs()[arr_ext[i]] * arr_ext[i + 1]).stem();
            }
            ext = pext * ssCof.pi_gb.Gen(arr_ext.back());
            stem_ext += ssCof.pi_gb.v_degs()[arr_ext.back()].stem();
            if (stem_ext != stem) {
                std::cout << "rel and ext not in the same stem!\n";
                return 7;
            }
        }
        auto rel_reduced = ssCof.pi_gb.ReduceV2(rel);
        if (rel_reduced && rel_reduced.data.back().IsUnKnown()) {
            int s = rel_reduced.UnknownFil();
            if (ext && s != ext.GetLead().fil()) {
                std::cout << "ext is in incorrect filtration: " << ext.GetLead().fil() << "!=" << s << '\n';
                return 8;
            }
            algZ::Mod rel_reduced1 = rel_reduced;
            rel_reduced1.data.pop_back();
            rel_reduced1.iaddP(ext, tmpm);
            rel_reduced1.iaddP(algZ::Mod::O(s + 1), tmpm);
            algZ::Mod rel1 = rel;
            rel1.isubP(rel_reduced1, tmpm, diagram.GetRings().pi_gb.gen_2tor_degs());

            diagram.ExtendRelMod(iSS - 1, stem, rel1);
            std::cout << ssCof.name << "  stem=" << stem << "  " << rel << '=' << rel_reduced << " --> " << rel_reduced1 << '\n';
            std::cout << rel1 << "=0\n";
            if (myio::UserConfirm()) {
                diagram.AddPiRelsCof(iSS - 1, {std::move(rel1)});
            }
            else {
                std::cout << "Aborted.\n";
                return 9;
            }
        }
        else {
            std::cout << "There is no extension to solve\n";
            return 10;
        }
    }

    try {
        diagram.SimplifyPiRels();
        int count_ss = 0, count_homotopy = 0;
        diagram.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_homotopy, 0);

        diagram.save(dbnames, DeduceFlag::homotopy);
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
    }*/

    // bench::Counter::print();
    return 0;
}