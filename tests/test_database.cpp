#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "algebras/dbalg.h"

TEST_CASE( "Save and load from database", "[DbAlg]" ) {
    using FnCmp = alg::CmpRevlex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;

    const myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/test.db");
    std::vector<alg::MayDeg> gen_degs = {{1, 1, 1}, {1, 2, 3}};
    std::vector<std::string> gen_names = {"x_1", "x_2"};
    Poly1d gen_reprs = {GenExp(1, 1), GenExp(2, 2)};
    Gb gb(Poly1d{GenExp(0, 2) + GenExp(1, 2), GenExp(1, 4) + GenExp(2, 1)});

    db.execute_cmd("delete from generators;");
    db.execute_cmd("delete from relations;");
    db.save_gen_maydegs("generators", gen_degs);
    db.save_gen_names("generators", gen_names);
    db.save_gen_reprs("generators", gen_reprs);

    db.save_gb("relations", gb, gen_degs, 0);

    auto gen_degs1 = db.load_gen_maydegs("generators");
    auto gen_names1 = db.load_gen_names("generators");
    auto gen_reprs1 = db.load_gen_reprs<FnCmp>("generators");

    auto gb1 = db.load_gb<FnCmp>("relations", -1);

    REQUIRE(gen_degs == gen_degs1);
    REQUIRE(gen_names == gen_names1);
    REQUIRE(gen_reprs == gen_reprs1);
    REQUIRE(gb.size() == gb1.size());
}
