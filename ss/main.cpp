#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "mylog.h"
#include <iostream>


int main_reset(int, char**, int&, const char*);
int main_resetpi(int, char**, int&, const char*);
int main_resetfrom(int, char**, int&, const char*);
int main_migrate_ss(int, char**, int&, const char*);
int main_migrate_htpy(int, char**, int&, const char*);

int main_plot(int, char**, int&, const char*);
int main_plotpi(int, char**, int&, const char*);
int main_add_diff(int, char**, int&, const char*);
int main_add_ext(int, char**, int&, const char*);

int main(int argc, char** argv)
{
    Logger::SetOutMain("ss.log");
    Logger::LogCmd(argc, argv);

    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"reset_ss", "Initialize the ss tables", main_reset},
        {"reset_pi", "Initialize the homotopy tables", main_resetpi},
        {"resetfrom", "Initialize from another diagram", main_resetfrom},
        {"migrate_ss", "Migrate ss between diagrams", main_migrate_ss},
        {"migrate_htpy", "Migrate homotopy between diagrams", main_migrate_htpy},
        {"plot_ss", "Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting", main_plot},
        {"plot_pi", "Compute the multiplications by hi", main_plotpi},
        {"add_diff", "Manually input a differential into the ss table", main_add_diff},
        {"add_ext", "Export the Adams E2 page of a module", main_add_ext},
        {"deduce", "Deduce differentials and extensions", main_deduce},
    };
    int index = 1;
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, "Manage spectral sequences and homotopy groups", VERSION, subcmds))
        return error;

    Logger::LogTime(timer.print2str());
    return 0;
}