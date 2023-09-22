#include "main.h"

/* generate the table of the spectral sequence */
void generate_db(const std::string& diagram, const std::string& name, int t_max, int r = 2)
{
    using namespace alg2;

    std::string path = fmt::format("{}/{}_AdamsSS.db", diagram, name);
    DBSS db(path);
    std::string table_prefix = fmt::format("{}_AdamsE2", name);
    auto gen_degs = db.load_gen_adamsdegs(table_prefix);
    auto rels = db.load_gb(table_prefix, 200);
    std::map<AdamsDeg, Mon1d> basis;
    std::map<AdamsDeg, Staircase> nodes_ss;

    int1d gen_degs_t;
    for (auto& deg : gen_degs)
        gen_degs_t.push_back(deg.t);
    Groebner gb(t_max, gen_degs_t);
    gb.AddRels(rels, t_max);
    auto leads = gb.GetLeadings();

    /* Generate basis */
    basis[AdamsDeg{0, 0}].push_back(Mon());
    for (int t = 0; t <= t_max; ++t) {
        std::map<AdamsDeg, Mon1d> basis_new;
        for (size_t gen_id = gen_degs.size(); gen_id-- > 0;) {
            int t1 = t - gen_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    auto deg_mon = p->first + gen_degs[gen_id];
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        Mon& m = p->second[i];
                        if (!m || (int)gen_id <= m[0].g()) {
                            Mon mon = m * Mon::Gen((uint32_t)gen_id);
                            if (gen_id >= leads.size() || std::none_of(leads[gen_id].begin(), leads[gen_id].end(), [&mon](const Mon& _m) { return divisible(_m, mon); })) {
                                basis_new[deg_mon].push_back(std::move(mon));
                            }
                        }
                    }
                }
            }
        }
        for (auto it = basis_new.begin(); it != basis_new.end(); ++it)
            basis.insert(std::move(*it));
    }

    /* Fill nodes_ss */
    int prev_t = 0;
    for (auto& [d, basis_d] : basis) {
        for (int i = 0; i < (int)basis_d.size(); ++i) {
            nodes_ss[d].basis.push_back({i});
            nodes_ss[d].diffs.push_back({-1});
            nodes_ss[d].levels.push_back(LEVEL_MAX - r);
        }
    }

    nodes_ss[AdamsDeg{0, 0}].diffs = {{}};
    nodes_ss[AdamsDeg{0, 0}].levels = {LEVEL_PERM};

    /* insert into the database */
    db.begin_transaction();
    db.create_generators(table_prefix);
    db.create_relations(table_prefix);
    db.create_basis(table_prefix);
    db.drop_and_create_ss(table_prefix);
    db.save_ss(table_prefix, nodes_ss);

    db.drop_and_create_pi_relations(name);
    db.drop_and_create_pi_basis(name);
    db.drop_and_create_pi_definitions(name);
    db.drop_and_create_pi_generators(name);
    db.end_transaction();
}

int main_add_basis(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram, name;
    int t_max = 0;

    myio::CmdArg1d args = {{"diagram", &diagram}, {"name", &name}, {"t_max", &t_max}};  // TODO: with diagram only
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    generate_db(diagram, name, t_max, 2);

    return 0;
}