#include "main.h"
#include "algebras/utility.h"
#include "algebras/myio.h"
#include <iostream>

int main(int argc, char** argv)
{
    std::string cmd;

    int index = 0;
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  basis_prod: Generate the basis_prod table\n";
        std::cout << "  ss_prod: Generate the ss_prod table and ss_diff table\n";
        std::cout << "  ss: Generate the ss table\n";
        std::cout << "  add_diff: Manually input a differential into the ss table\n";
        std::cout << "  deduce: Deduce differentials by Leibniz rules\n\n";

        std::cout << "Version:\n  1.0 (2022-7-10)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "ss")
        return main_generate_ss(argc, argv, index);
    else if (cmd == "basis_prod")
        return main_basis_prod(argc, argv, index);
    else if (cmd == "ss_prod")
        return main_ss_prod(argc, argv, index);
    else if (cmd == "add_diff")
        return main_add_diff(argc, argv, index);
    else if (cmd == "deduce")
        return main_deduce(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << argv[1] << '\n';
}