#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "algebras/dbalg.h"

TEST_CASE( "Save and load from database", "[DbAlg]" )
{
    using FnCmp = alg2::CmpRevlex;
    using Poly = alg2::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg2::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;

    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/TestMilnorProduct.db");
    std::vector<alg2::MayDeg> gen_degs = {{1, 1, 1}, {1, 2, 3}, {3, 4, 5}};
    std::vector<std::string> gen_names = {"x_1", "x_2", "x_3"};
    Poly1d gen_reprs = {GenExp(1, 1), GenExp(2, 2), GenExp(0, 1)};
    Gb gb(100, Poly1d{GenExp(0, 2) + GenExp(1, 2), GenExp(1, 4) + GenExp(2, 1), GenExp(2, 1)});

    db.create_generators_and_delete("table_bcce9d67");
    db.create_relations_and_delete("table_bcce9d67");
    
    db.save_gen_maydegs("table_bcce9d67", gen_degs);
    db.save_gen_names("table_bcce9d67", gen_names);
    db.save_gen_reprs("table_bcce9d67", gen_reprs);

    db.save_gb("table_bcce9d67", gb, gen_degs, 0);

    auto gen_degs1 = db.load_gen_maydegs("table_bcce9d67");
    auto gen_names1 = db.load_gen_names("table_bcce9d67");
    auto gen_reprs1 = db.load_gen_reprs<FnCmp>("table_bcce9d67");

    auto gb1 = db.load_gb<FnCmp>("table_bcce9d67", 100);

    REQUIRE(gen_degs == gen_degs1);
    REQUIRE(gen_names == gen_names1);
    REQUIRE(gen_reprs == gen_reprs1);
    REQUIRE(gb.size() == gb1.size());
}

TEST_CASE("Save and load blobs from database", "[DbAlg]")
{
    using T = std::array<int, 10>;

    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/TestMilnorProduct.db");
    db.create_generators_and_delete("table_57fc659c");
    std::vector<alg2::MayDeg> gen_degs = {{1, 1, 1}, {1, 2, 3}};
    
    T a = {1, 2, 3, 5, 9, 1, 3, 4, 1, 2};
    T b = {1, 2, 3, 5, 9, 1, 3, 4, 1, 2};
    std::vector<T> reprs = {a, b};

    db.save_gen_maydegs("table_57fc659c", gen_degs);
    db.update_blob_column(
        "table_57fc659c_generators", "repr", "gen_id", reprs, [](const T& a) { return std::make_pair((const void*)a.data(), a.size() * sizeof(int)); }, 0);
    std::vector<T> reprs1 = db.get_column_from_blob<T>("table_57fc659c_generators", "repr", "", [](const void* data, size_t n) {
        T result;
        memcpy(result.data(), data, n);
        return result;
    });
    REQUIRE(reprs == reprs1);
}
