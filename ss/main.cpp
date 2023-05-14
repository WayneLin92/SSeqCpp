#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "mylog.h"
#include <iostream>

int main_basis_prod(int argc, char** argv, int index);
int main_plot(int argc, char** argv, int index);
int main_plotpi(int argc, char** argv, int index);
int main_add_diff(int argc, char** argv, int index);
int main_add_ext(int argc, char** argv, int index);

int main_reset(int argc, char** argv, int index);
int main_resetpi(int argc, char** argv, int index);
int main_resetfrom(int argc, char** argv, int index);
int main_migrate_ss(int argc, char** argv, int index);
int main_migrate_htpy(int argc, char** argv, int index);

int main(int argc, char** argv)
{
    Logger::SetOutMain("ss.log");
    Logger::LogCmd(argc, argv);

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

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    int rt = 0;
    if (cmd == "reset")
        rt = main_reset(argc, argv, index);
    else if (cmd == "resetpi")
        rt = main_resetpi(argc, argv, index);
    else if (cmd == "resetfrom")
        rt = main_resetfrom(argc, argv, index);
    else if (cmd == "migrate_ss")
        rt = main_migrate_ss(argc, argv, index);
    else if (cmd == "migrate_htpy")
        rt = main_migrate_htpy(argc, argv, index);
    else if (cmd == "plot")
        rt = main_plot(argc, argv, index);
    else if (cmd == "plotpi")
        rt = main_plotpi(argc, argv, index);
    else if (cmd == "add_diff")
        rt = main_add_diff(argc, argv, index);
    else if (cmd == "add_ext")
        rt = main_add_ext(argc, argv, index);
    else if (cmd == "deduce")
        rt = main_deduce(argc, argv, index);
    else {
        std::cerr << "Invalid cmd: " << cmd << '\n';
        return -1;
    }

    Logger::LogTime(timer.print2str());
    return rt;
}