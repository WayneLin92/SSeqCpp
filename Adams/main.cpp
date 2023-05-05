#include "main.h"
#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main_res(int argc, char** argv, int index);
int main_res_RP(int argc, char** argv, int index);
int main_prod(int argc, char** argv, int index);
int main_prod_mod(int argc, char** argv, int index);
int main_prod_hi(int argc, char** argv, int index);
int main_export(int argc, char** argv, int index);
int main_export_mod(int argc, char** argv, int index);
int main_2cell(int argc, char** argv, int index);
int main_2cell_export(int argc, char** argv, int index);
int main_generators_to_csv(const std::string& db_filename, const std::string& tablename);

int main(int argc, char** argv)
{
    std::string cmd;

    int index = 0;
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Usage:\n  Adams <cmd> [-h] ...\n\n");

        fmt::print("<cmd> can be one of the following:\n");

        fmt::print("  res: Compute an A-resolution\n");
        fmt::print("  prod: Compute the multiplications\n");
        fmt::print("  prod_mod: Compute the multiplications of a module over a ring spectra\n");
        fmt::print("  2cell: Generate the basis of Ext for a cell complex\n");
        fmt::print("  export: Export AdamsE2 of a ring spectrum\n");
        fmt::print("  export_mod: Export AdamsE2 of a module over a ring spectrum\n");
        fmt::print("  plot: Generate an html file (feature not supported yet)\n");

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    bench::Timer timer;

    if (cmd == "res")
        return main_res(argc, argv, index);
    else if (cmd == "res_RP")
        return main_res_RP(argc, argv, index);
    else if (cmd == "prod")
        return main_prod(argc, argv, index);
    else if (cmd == "prod_mod")
        return main_prod_mod(argc, argv, index);
    else if (cmd == "2cell")
        return main_2cell(argc, argv, index);
    else if (cmd == "prod_hi")
        return main_prod_hi(argc, argv, index);
    else if (cmd == "export")
        return main_export(argc, argv, index);
    else if (cmd == "export_mod")
        return main_export_mod(argc, argv, index);
    else if (cmd == "res_csv")
        return main_generators_to_csv("X:/Nan/AdamsDB/res/S0_Adams_res_d261.db", "S0_Adams_res");
    else
        std::cerr << "Invalid cmd: " << cmd << '\n';

    return 0;
}