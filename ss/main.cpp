#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "mylog.h"
#include <iostream>


int main_reset_ss(int, char**, int&, const char*);
int main_reset_pi(int, char**, int&, const char*);
int main_resetfrom(int, char**, int&, const char*);
int main_add_diff(int, char**, int&, const char*);
int main_add_pi(int, char**, int&, const char*);

int main_deduce(int, char**, int&, const char*);

int main_migrate_ss(int, char**, int&, const char*);
int main_migrate_pi(int, char**, int&, const char*);
int main_import_chua_d2(int, char**, int&, const char*);

int main_plot_ss(int, char**, int&, const char*);
int main_plot_pi(int, char**, int&, const char*);
int main_rename_gen(int, char**, int&, const char*);


int main(int argc, char** argv)
{
    Logger::SetOutMain("ss.log");
    Logger::LogCmd(argc, argv);

    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"reset_ss", "Initialize the ss tables", main_reset_ss},
        {"reset_pi", "Initialize the homotopy tables", main_reset_pi},
        {"resetfrom", "Initialize from another diagram", main_resetfrom},
        {"add_diff", "Manually input a differential into the ss table", main_add_diff},
        {"add_pi", "Manually input a homotopy into the pi table", main_add_pi},
        {"deduce", "Deduce differentials and extensions", main_deduce},
        {"migrate_ss", "Migrate ss between diagrams", main_migrate_ss},
        {"migrate_pi", "Migrate homotopy between diagrams", main_migrate_pi},
        {"import_chua_d2", "Import the d2 data of Chua", main_import_chua_d2},
        {"plot_ss", "Generate the json data for plotting", main_plot_ss},
        {"plot_pi", "Compute the multiplications by hi", main_plot_pi},
        {"rename_gen", "Manage generator names", main_rename_gen},
    };
    int index = 1;
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, "Manage spectral sequences and homotopy groups", VERSION, subcmds))
        return error;

    Logger::LogTime(timer.print2str());
    return 0;
}
