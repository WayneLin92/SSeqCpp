#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/myexception.h"

int main()
{
    int t_max = 511;

    Timer timer;

    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/HX9.db");
    myio::DbAlg db1("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/E4.db");

    std::string table_prefix_E2 = "HX9";
    std::string table_prefix_E4 = "E4";

    db1.create_generators_and_delete(table_prefix_E4);
    db1.create_relations_and_delete(table_prefix_E4);

    auto gen_degs = db.load_gen_maydegs(table_prefix_E2);
    auto gen_names = db.load_gen_names(table_prefix_E2);
    auto gen_diffs = db.load_gen_diffs<alg::CmpRevlex>(table_prefix_E2);
    auto gb = db.load_gb<alg::CmpRevlex>(table_prefix_E2, t_max);

    alg::GroebnerRevlex gb_h;
    alg::MayDeg1d gen_degs_h;
    alg::PolyRevlex1d gen_repr_h;
    std::vector<std::string> gen_names_h;
    
    alg::Homology(gb, gen_degs, gen_diffs, alg::MayDeg{1, 0, -2}, gb_h, gen_degs_h, gen_names_h, gen_repr_h, 300);

    db1.begin_transaction();
    db1.save_gen_maydegs(table_prefix_E4, gen_degs_h);
    db1.save_gen_names(table_prefix_E4, gen_names_h);
    db1.save_gen_reprs(table_prefix_E4, gen_repr_h);
    db1.save_gb(table_prefix_E4, gb_h, gen_degs);
    db1.end_transaction();

    return 0;
}