#include "main.h"
#include "algebras/utility.h"
#include "algebras/myio.h"
#include <iostream>

int main(int argc, char** argv)
{
    std::string cmd;

    int index = 0;
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  Adams <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  res: Compute an A-resolution\n";
        std::cout << "  cell: Generate the basis of Ext for a cell complex\n";
        std::cout << "  prod: Compute the multiplications\n";
        std::cout << "  export: Extract info of Ext from the resolution data\n";
        std::cout << "  plot: Generate an html file (feature not supported yet)\n";

        std::cout << "Version:\n  2.0 (2022-08-07)" << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "res")
        return main_res(argc, argv, index);
    else if (cmd == "prod")
        return main_prod(argc, argv, index);
    else if (cmd == "export")
        return main_export(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << argv[1] << '\n';
}