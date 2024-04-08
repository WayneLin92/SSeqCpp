#include "main.h"
#include "mylog.h"

int main_mul(int argc, char** argv, int& index, const char* desc)
{
    std::string cw;
    int stem1, s1, stem2, s2;
    std::string str_x1, str_x2;
    int1d x1, x2;

    std::string diagram_name = "default";

    myio::CmdArg1d args = {{"cw", &cw}, {"stem1", &stem1}, {"s1", &s1}, {"x1", &str_x1}, {"stem2", &stem2}, {"s2", &s2}, {"x2", &str_x2}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    AdamsDeg d1(s1, stem1 + s1), d2(s2, stem2 + s2), d3(d1 + d2);
    x1 = myio::Deserialize<int1d>(str_x1);
    x2 = myio::Deserialize<int1d>(str_x2);

    SSFlag flag = SSFlag::no_op;
    Diagram diagram(diagram_name, flag);

    if (auto iCw = diagram.GetIndexCwByName(cw); iCw.isRing) {
        auto& ring = diagram.GetRingByName(cw);
        if (d3.t > ring.t_max) {
            fmt::print("degree out of range");
            return 0;
        }
        Poly poly_x1 = Indices2Poly(x1, ring.basis.at(d1));
        Poly poly_x2 = Indices2Poly(x2, ring.basis.at(d2));
        Poly poly_x3 = ring.gb.Reduce(poly_x1 * poly_x2);
        int1d x3 = Poly2Indices(poly_x3, ring.basis.at(d3));
        fmt::print("{} [{}]*[{}]=[{}]\n", d3, myio::Serialize(x1), myio::Serialize(x2), myio::Serialize(x3));
    }
    else {
        auto& mod = diagram.GetModuleByName(cw);
        auto& ring = diagram.GetRings()[mod.iRing];
        if (d3.t > ring.t_max) {
            fmt::print("degree out of range");
            return 0;
        }
        Poly poly_x1 = Indices2Poly(x1, ring.basis.at(d1));
        Mod mod_x2 = Indices2Mod(x2, mod.basis.at(d2));
        Mod mod_x3 = mod.gb.Reduce(poly_x1 * mod_x2);
        int1d x3 = Mod2Indices(mod_x3, mod.basis.at(d3));
        fmt::print("{} [{}]*[{}]=[{}]\n", d3, myio::Serialize(x1), myio::Serialize(x2), myio::Serialize(x3));
    }

    return 0;
}