#include "main.h"

/* generate the table of the spectral sequence */
void generate_db(const std::string& name, int t_max, int r = 2)
{
    using namespace alg2;

    std::string path = fmt::format("{}_AdamsSS.db", name);
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

    /* Generate basis */
    for (int t = 0; t <= t_max; ++t) {
        for (size_t gen_id = gen_degs.size(); gen_id-- > 0;) {
            
        }
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

int main_generate_basis(int argc, char** argv, int& index, const char* desc)
{
    std::string name;
    int t_max = 0;

    myio::CmdArg1d args = {{"name", &name}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    generate_db(name, t_max, 2);

    return 0;
}