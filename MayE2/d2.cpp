#include "MayE2.h"
#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/myexception.h"
#include <iostream>

alg::PolyRevlex Decompose(const std::vector<std::string>& gen_names, const alg::array& S)
{
    alg::PolyRevlex result = alg::PolyRevlex::Unit();
    int count_t_minus_s = 0;
    size_t prev_i = 0;
    for (size_t i = 1; i <= S.size(); ++i) {
        count_t_minus_s += 1;
        if (i < S.size())
            count_t_minus_s -= S[i] - S[i - 1] - 1;
        if (count_t_minus_s == 0 || i == S.size()) {
            alg::array S1(S.begin() + prev_i, S.begin() + i);
            prev_i = i;
            result = result * GetPolyByName<alg::CmpRevlex>(gen_names, get_name(S1));
        }
        else if (count_t_minus_s < 0)
            throw MyException(0x210be4aaU, "Incorrect S");
    }
    return result;
}

void GetD2(int n)
{
    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/HX9.db");
    std::string table_prefix = "HX" + std::to_string(n);
    std::vector<std::string> gen_names = db.load_gen_names(table_prefix);
    alg::PolyRevlex null = alg::PolyRevlex::Gen(-1);
    alg::PolyRevlex1d gen_diffs(gen_names.size(), null);

    /* dh_{S, T} = sum_i h_{i+1}h_{S-{i}+{i+1},T-{i+1}+{i}}*/
    for (int k = 0; k < n; ++k) {
        for (int m = k + 1; m <= n; ++m) {
            for (const auto& [S, T] : H(k, m + 1)) {
                //std::cout << S << " " << T << '\n';
                alg::PolyRevlex diff;
                for (int i = 0; i < m; ++i) {
                    alg::array S1 = S;
                    alg::array T1 = T;
                    auto p1 = std::find(S1.begin(), S1.end(), i);
                    auto p2 = std::find(T1.begin(), T1.end(), i + 1);
                    if (p1 != S1.end() && p2 != T1.end()) {
                        *p1 = i + 1;
                        *p2 = i;
                        // std::cout << "diff " << S1 << " " << T1 << '\n';
                        if (S1.front() < T1.front())
                            diff += GetPolyByName<alg::CmpRevlex>(gen_names, get_name({i + 1})) * Decompose(gen_names, S1);
                    }
                }
                auto p_gen_names = std::find(gen_names.begin(), gen_names.end(), get_name(S));
                if (p_gen_names == gen_names.end())
                    throw MyException(0xc34463fdU, "Gen name not found");
                gen_diffs[p_gen_names - gen_names.begin()] = std::move(diff);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 2; j < n; ++j) {
            alg::PolyRevlex h1 = GetPolyByName<alg::CmpRevlex>(gen_names, "h_" + std::to_string(i + 1));
            alg::PolyRevlex h2 = GetPolyByName<alg::CmpRevlex>(gen_names, "h_" + std::to_string(j));
            alg::PolyRevlex b1, b2;
            if (j - i > 2) {
                b1 = GetPolyByName<alg::CmpRevlex>(gen_names, "b_{" + std::to_string(i + 1) + std::to_string(j) + "}");
                b2 = GetPolyByName<alg::CmpRevlex>(gen_names, "b_{" + std::to_string(i) + std::to_string(j - 1) + "}");
            }
            else {
                b1 = pow(GetPolyByName<alg::CmpRevlex>(gen_names, "h_" + std::to_string(i + 1)), 2);
                b2 = pow(GetPolyByName<alg::CmpRevlex>(gen_names, "h_" + std::to_string(i)), 2);
            }
            alg::PolyRevlex diff = h1 * b1 + h2 * b2;
            std::string name = "b_{" + std::to_string(i) + std::to_string(j) + "}";
            //std::cout << "d" << name << "=" << StrPoly(diff.data, gen_names) << '\n';
            auto p_gen_names = std::find(gen_names.begin(), gen_names.end(), name);
            if (p_gen_names == gen_names.end())
                throw MyException(0x1a06fa82U, "Gen name not found");
            gen_diffs[p_gen_names - gen_names.begin()] = std::move(diff);
        }
    }

    for (size_t i = 0; i < gen_names.size(); ++i) {
        if (gen_diffs[i] != null)
        std::cout << "d" << gen_names[i] << "=" << StrPoly(gen_diffs[i].data, gen_names) << '\n';
    }

    db.begin_transaction();

    {
        myio::Statement stmt(db, "UPDATE " + table_prefix + "_generators SET gen_diff = ?1 WHERE gen_id= ?2;");
        for (size_t i = 0; i < gen_diffs.size(); ++i) {
            if (gen_diffs[i] != null) {
                stmt.bind_str(1, db.Serialize(gen_diffs[i].data));
                stmt.bind_int(2, (int)i);
                stmt.step_and_reset();
            }
        }
        std::cout << gen_diffs.size() << ' ' << "gen_diffs are inserted into " + table_prefix + "_generators!\n";
    }

    db.end_transaction();
}