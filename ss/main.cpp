#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main(int argc, char** argv)
{
    std::string cmd;

    int index = 0;
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  init: Initialize the ss table\n";
        std::cout << "  basis_prod: Generate the basis_prod table for S0\n";
        std::cout << "  plot: Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting\n";
        std::cout << "  add_diff: Manually input a differential into the ss table\n";
        std::cout << "  try_add_diff: Try to input a differential into the ss table and detect contradictions without changing the database\n";
        std::cout << "  deduce: Deduce differentials by Leibniz rules\n\n";
        std::cout << "  mod: Deal with a complex where the its AdamsSS is a module over S0_AdamsSS\n\n";

        std::cout << "Version:\n  1.0 (2022-7-29)" << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "init")
        return main_generate_ss(argc, argv, index);
    else if (cmd == "basis_prod")
        return main_basis_prod(argc, argv, index);
    else if (cmd == "plot")
        return main_plot(argc, argv, index);
    else if (cmd == "add_diff")
        return main_add_diff(argc, argv, index);
    else if (cmd == "try_add_diff")
        return main_try_add_diff(argc, argv, index);
    else if (cmd == "deduce")
        return main_deduce(argc, argv, index);
    else if (cmd == "mod")
        return main_mod(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << cmd << '\n';
}