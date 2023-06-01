#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "mylog.h"
#include <iostream>


int main_reset_ss(int, char**, int&, const char*);
int main_reset_pi(int, char**, int&, const char*);
int main_resetfrom(int, char**, int&, const char*);
int main_migrate_ss(int, char**, int&, const char*);
int main_migrate_pi(int, char**, int&, const char*);

int main_plot(int, char**, int&, const char*);
int main_plot_htpy(int, char**, int&, const char*);
int main_add_diff(int, char**, int&, const char*);
int main_add_pi(int, char**, int&, const char*);

int main(int argc, char** argv)
{
    Logger::SetOutMain("ss.log");
    Logger::LogCmd(argc, argv);

    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"reset_ss", "Initialize the ss tables", main_reset_ss},
        {"reset_pi", "Initialize the homotopy tables", main_reset_pi},
        {"resetfrom", "Initialize from another diagram", main_resetfrom},
        {"migrate_ss", "Migrate ss between diagrams", main_migrate_ss},
        {"migrate_pi", "Migrate homotopy between diagrams", main_migrate_pi},
        {"plot_ss", "Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting", main_plot},
        {"plot_htpy", "Compute the multiplications by hi", main_plot_htpy},
        {"add_diff", "Manually input a differential into the ss table", main_add_diff},
        {"add_pi", "Manually input a homotopy into the pi table", main_add_pi},
        {"deduce", "Deduce differentials and extensions", main_deduce},
    };
    int index = 1;
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, "Manage spectral sequences and homotopy groups", VERSION, subcmds))
        return error;

    Logger::LogTime(timer.print2str());
    return 0;
}
