#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main(int argc, char** argv)
{
    myio::Logger::Init("ss.log");
    for (int i = 0; i < argc; ++i)
        myio::Logger::fout() << "cmd: " << argv[i] << ' ';
    myio::Logger::fout() << '\n';

    bench::Timer timer;

    std::string cmd;
    int index = 0;
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  reset: Initialize the ss table\n";
        std::cout << "  resetpi: Initialize the homotopy data\n";
        std::cout << "  basis_prod: Generate the basis_prod table for S0\n";
        std::cout << "  plot: Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting\n";
        std::cout << "  add_diff: Manually input a differential into the ss table\n";
        std::cout << "  try_add_diff: Try to input a differential into the ss table and detect contradictions without changing the database\n";
        std::cout << "  deduce: Deduce differentials and extensions\n\n";
        std::cout << "  mod: Deal with a complex where the its AdamsSS is a module over S0_AdamsSS\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "reset")
        return main_reset(argc, argv, index);
    else if (cmd == "resetpi")
        return main_resetpi(argc, argv, index);
    else if (cmd == "truncate")
        return main_truncate(argc, argv, index);
    else if (cmd == "basis_prod")
        return main_basis_prod(argc, argv, index);
    else if (cmd == "plot")
        return main_plot(argc, argv, index);
    else if (cmd == "plotpi")
        return main_plotpi(argc, argv, index);
    else if (cmd == "add_diff")
        return main_add_diff(argc, argv, index);
    else if (cmd == "add_ext")
        return main_add_ext(argc, argv, index);
    else if (cmd == "deduce")
        return main_deduce(argc, argv, index);
    else if (cmd == "mod")
        return main_mod(argc, argv, index);
    else {
        std::cerr << "Invalid cmd: " << cmd << '\n';
        return 0;
    }
}