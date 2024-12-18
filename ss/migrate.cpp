#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <fstream>

/* Add the differentials from category1 to category2 */
SSRet Migrate(const Category& category1, Category& category2, SSFlag flag)
{
    int1d diff;
    auto rings_size = category1.GetRings().size();
    const size_t num_cw = rings_size + category1.GetModules().size();
    for (size_t i_Cw = 0; i_Cw < num_cw; ++i_Cw) {
        auto iCw1 = i_Cw < rings_size ? IndexRing(i_Cw) : IndexMod(i_Cw - rings_size);

        auto& name = category1.GetCwName(iCw1);
        auto iCw2 = category2.GetIndexCwByName(name);
        if (!iCw2)
            continue;

        auto& nodes_ss1 = category1.GetNodesSS(iCw1);
        auto& degs_ss1 = category1.GetSSDegs(iCw1);
        auto& nodes_ss2 = category2.GetNodesSS(iCw1);

        int t_max2 = nodes_ss2.front().t_max();
        for (AdamsDeg deg : degs_ss1) {
            if (deg.t > t_max2)
                break;
            const auto& sc1 = nodes_ss1.GetRecentSc(deg);
            for (size_t i = 0; i < sc1.levels.size(); ++i) {
                if (sc1.levels[i] > LEVEL_MAX / 2) {
                    int r = LEVEL_MAX - sc1.levels[i];
                    if (deg.t + r - 1 > t_max2 || sc1.diffs[i] == NULL_DIFF) {
                        --r;
                        diff = int1d{};
                    }
                    else
                        diff = sc1.diffs[i];
                    if (category2.IsNewDiff(nodes_ss2, deg, sc1.basis[i], diff, r)) {
                        Logger::LogDiff(0, EnumReason::migrate, name, deg, r, sc1.basis[i], diff, "", flag);
                        if (auto rt = category2.SetCwDiffGlobal(iCw2, deg, sc1.basis[i], diff, r, true, flag))
                            return rt;
                    }
                }
                else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                    int r = sc1.levels[i] + 1;
                    const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                    if (category2.IsNewDiff(nodes_ss2, deg_src, {}, sc1.basis[i], r)) {
                        Logger::LogDiffInv(0, EnumReason::migrate, name, deg, r, {}, sc1.basis[i], "", flag);
                        if (auto rt = category2.SetCwDiffGlobal(iCw2, deg_src, {}, sc1.basis[i], r, true, flag))
                            return rt;
                    }
                }
            }
        }
    }
    auto cofseqs_size = category1.GetCofSeqs().size();
    for (auto i_Cof1 = 0; i_Cof1 < cofseqs_size; ++i_Cof1) {
        auto& cofseq1 = category1.GetCofSeqs()[i_Cof1];
        auto& name = cofseq1.name;
        auto i_Cof2 = category2.GetCofSeqIndexByName(name);
        if (i_Cof2 == -1)
            continue;
        auto& cofseq2 = category2.GetCofSeqs()[i_Cof2];
        for (size_t iTri = 0; iTri < 3; ++iTri) {
            int stem_map = cofseq1.degMap[iTri].stem();
            int stem_map_prev = cofseq1.degMap[PreviTri(iTri)].stem();
            int t_max2 = category2.GetCofSeqs()[i_Cof2].t_max[iTri];
            int t_max2_next = category2.GetCofSeqs()[i_Cof2].t_max[NextiTri(iTri)];
            auto& nodes_cofseq = cofseq1.nodes_cofseq[iTri];
            for (AdamsDeg deg : nodes_cofseq.front().degs()) {
                if (deg.t > t_max2)
                    continue;
                const auto& sc1 = nodes_cofseq.GetRecentSc(deg);
                for (size_t i = 0; i < sc1.levels.size(); ++i) {
                    if (sc1.levels[i] > LEVEL_MAX / 2) {
                        int r = LEVEL_MAX - sc1.levels[i];
                        if (deg.t + stem_map + r > t_max2_next || sc1.diffs[i] == NULL_DIFF) {
                            --r;
                            diff = int1d{};
                        }
                        else
                            diff = sc1.diffs[i];
                        if (auto rt = category2.SetCwDiffGlobal(cofseq2.indexCw[iTri], deg, sc1.basis[i], {}, R_PERM - 1, false, flag))  //// log
                            throw RunTimeError("Failed to category2.SetCwDiffGlobal()");
                        if (diff.size()) {
                            if (auto rt = category2.SetCwDiffGlobal(cofseq2.indexCw[NextiTri(iTri)], deg + AdamsDeg(r, stem_map + r), diff, {}, R_PERM - 1, false, flag))
                                throw RunTimeError("Failed to category2.SetCwDiffGlobal()");
                        }
                        if (category2.IsNewDiffCofseq(cofseq2, iTri, deg, sc1.basis[i], diff, r)) {
                            Logger::LogDiff(0, EnumReason::migrate, fmt::format("{}:{}", name, iTri), deg, r, sc1.basis[i], diff, "", flag);
                            if (auto rt = category2.SetDiffGlobalCofseq(cofseq2, iTri, deg, sc1.basis[i], diff, r, true, flag))
                                return rt;
                        }
                    }
                    else if (sc1.levels[i] < LEVEL_MAX / 2 && sc1.diffs[i] == NULL_DIFF) {
                        int r = sc1.levels[i] + 1;
                        const AdamsDeg deg_src = deg - AdamsDeg(r, stem_map_prev + r);
                        if (auto rt = category2.SetCwDiffGlobal(cofseq2.indexCw[iTri], deg, sc1.basis[i], {}, R_PERM - 1, false, flag))  //// log
                            throw RunTimeError("Failed to category2.SetCwDiffGlobal()");
                        if (category2.IsNewDiffCofseq(cofseq2, PreviTri(iTri), deg_src, {}, sc1.basis[i], r)) {
                            Logger::LogDiffInv(0, EnumReason::migrate, name, deg, r, {}, sc1.basis[i], "", flag);
                            if (auto rt = category2.SetDiffGlobalCofseq(cofseq2, PreviTri(iTri), deg_src, {}, sc1.basis[i], r, true, flag))
                                return rt;
                        }
                    }
                }
            }
        }
    }
    // Logger::LogSummary("Changed differentials", count);
    return SSRet::NUL();
}

/* Add the homotopy from category1 to category2
 * Currently we only support migration from small range to bigger range
 */
// void Migrate_htpy(const Category& category1, Category& category2)
//{
//  for (size_t iSS = 0; iSS < category1.GetAllBasisSs().size(); ++iSS) {
//      if (iSS == 0) {
//          auto& pi_gb1 = category1.GetRings().pi_gb;
//          auto& pi_gb2 = category2.GetRings().pi_gb;
//          int t_max2 = category2.GetRings().nodes_ss.front().rbegin()->first.t;
//          if (pi_gb2.gen_degs().size() > 1) {  //// TODO: Support merge of homotopy information
//              Logger::LogException(0, 0x43dbbfcaU, "Diagram2 should be reset first.\n");
//              throw MyException(0x43dbbfcaU, "Diagram2 should be reset first.");
//          }
//          algZ::Poly1d rels = pi_gb1.data();
//          uint32_t i = 0;
//          for (auto& deg : pi_gb1.gen_degs()) {
//              if (deg.t <= t_max2) {
//                  pi_gb2.AddGen(deg);
//                  if (deg.stem() % 2 == 1)
//                      rels.push_back(algZ::Mon::two_x_square(i, deg.s));
//              }
//              ++i;
//          }
//          category2.GetRings().pi_gen_Einf = category1.GetRings().pi_gen_Einf;
//          category2.GetRings().pi_gen_defs = category1.GetRings().pi_gen_defs;
//          category2.GetRings().pi_gen_def_mons = category1.GetRings().pi_gen_def_mons;
//          category2.AddPiRelsRing(std::move(rels));
//      }
//      else {
//          auto& ssCof1 = category1.GetModules()[iSS - 1];
//          auto& ssCof2 = category2.GetModules()[iSS - 1];
//          auto& pi_gb1 = ssCof1.pi_gb;
//          auto& pi_gb2 = ssCof2.pi_gb;
//          int t_max2 = ssCof2.nodes_ss.front().rbegin()->first.t;
//          if (pi_gb2.v_degs().size() > 1) {  //// TODO: Support merge of homotopy information
//              Logger::LogException(0, 0x92cb2691, "Diagram2 should be reset first.\n");
//              throw MyException(0x92cb2691, "Diagram2 should be reset first.");
//          }
//          for (auto& deg : pi_gb1.v_degs())
//              if (deg.t <= t_max2)
//                  pi_gb2.AddGen(deg);
//          ssCof2.pi_gen_Einf = ssCof1.pi_gen_Einf;
//          ssCof2.pi_gen_defs = ssCof1.pi_gen_defs;
//          ssCof2.pi_gen_def_mons = ssCof1.pi_gen_def_mons;
//          ssCof2.nodes_pi_qt = ssCof1.nodes_pi_qt;
//          category2.AddPiRelsCof(iSS - 1, pi_gb1.data());
//      }
//  }
//  int count_ss = 0, count_htpy = 0;
//  category2.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_htpy, 0);
//}

int main_migrate(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name1, cat_name2;

    myio::CmdArg1d args = {{"category1", &cat_name1}, {"category2", &cat_name2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    auto [cat_root, ckpt] = ParseCatName(cat_name1);
    MyException::Assert(cat_name1 != cat_name2, "Need cat_name1 != cat_name2");

    DbLog log2(fmt::format("{}/log.db", cat_name2));

    auto flag = SSFlag::cofseq;
    Category category1(cat_root, ckpt, flag, false);
    Category category2(cat_name2, "", flag, true);  // log to category2 only when we do not want to migrate log
    try {
        if (auto rt = Migrate(category1, category2, flag))
            throw ErrorMsg("Migrate failed");
        category2.SaveNodes(cat_name2, "", true, flag);
    }
    /*catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}

void migrate_log(const std::string& log1, const std::string& log2)
{
    DbLog dblog2(log2);
    int id_max = dblog2.get_int("SELECT COALESCE(MAX(id), -1) FROM log");
    dblog2.execute_cmd(fmt::format("ATTACH '{}' as orig;"
                                   " INSERT INTO log (depth, reason, name, s, t, r, x, dx, info) SELECT depth, reason, name, s, t, r, x, dx, info FROM orig.log;"
                                   " INSERT INTO exclusions (id) SELECT id+{} FROM orig.exclusions;",
                                   log1, id_max, id_max, id_max));
}

int main_migrate_from_logs(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    myio::string1d cats_of_logs;

    myio::CmdArg1d args = {{"cat_name", &cat_name}, {"cats_of_logs", &cats_of_logs}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::cofseq;
    Category category(cat_name, "", flag, false);

    for (auto& cat : cats_of_logs) {
        fmt::print("cat={}\n", cat);

        DbLog log(fmt::format("{}/log.db", cat));
        myio::Statement stmt(log, fmt::format("SELECT id, name, s, t, r, x, dx, name like '%' || ':' || '%', reason like '%' || 'I' FROM log WHERE depth=0 AND reason IS NOT NULL AND reason NOT IN ('G', 'GI', 'N', 'CsCm')"));  ////
        while (stmt.step() == MYSQLITE_ROW) {
            auto [id, name, deg, r, x, dx, isCs, isDInv] = columns_diff(stmt);
            fmt::print("id={}     \r", id);
            category.SetUnivDiffGlobal(name, deg, r, x, dx, isCs, isDInv, flag);
        }
        if (auto rt = category.DeduceTrivialDiffs(flag))
            throw RunTimeError("Failed to category.DeduceTrivialDiffs()");

        migrate_log(fmt::format("{}/log.db", cat), fmt::format("{}/log.db", cat_name));
    }
    category.SaveNodes(cat_name, "", true, flag);

    return 0;
}

int main_migrate_pi(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name1, cat_name2;

    myio::CmdArg1d args = {{"category1", &cat_name1}, {"category2", &cat_name2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    auto [cat_root, ckpt] = ParseCatName(cat_name1);

    SSFlag flag = SSFlag::pi | SSFlag::pi_def;
    Category category1(cat_root, ckpt, flag);
    Category category2(cat_name2, "", flag, false);

    try {
        // Migrate_htpy(category1, category2);
        category2.SaveNodes(cat_name2, "", true, flag);
    }
    /* catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}

int main_resetfrom(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    std::string cat_name_from;

    myio::CmdArg1d args = {{"category", &cat_name}, {"category_from", &cat_name_from}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    std::vector<std::string> names, paths;
    std::vector<int> isRing;
    GetAllDbNames(cat_name, names, paths, isRing);
    std::vector<std::string> names2, paths2;
    std::vector<int> isRing2;
    GetAllDbNames(cat_name_from, names2, paths2, isRing2);

    for (size_t k = 0; k < names.size(); ++k) {
        std::ifstream src(paths2[k], std::ios::binary);
        std::ofstream dst(paths[k], std::ios::binary);
        dst << src.rdbuf();
    }

    return 0;
}