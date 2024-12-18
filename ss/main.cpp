#include "main.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "mylog.h"
#include <iostream>

int main_reset(int, char**, int&, const char*);
int main_reset_cofseq(int, char**, int&, const char*);
int main_reset_pi(int, char**, int&, const char*);
int main_checkpoint(int, char**, int&, const char*);
int main_rollback(int, char**, int&, const char*);
int main_resetfrom(int, char**, int&, const char*);
int main_add_diff(int, char**, int&, const char*);
int main_add_diff_from_file(int, char**, int&, const char*);
int main_add_diff_from_log(int, char**, int&, const char*);
int main_add_cofseq_diff(int, char**, int&, const char*);

int main_deduce(int, char**, int&, const char*);

int main_migrate(int, char**, int&, const char*);
int main_migrate_from_logs(int, char**, int&, const char*);
int main_migrate_pi(int, char**, int&, const char*);

int main_plot_ss(int, char**, int&, const char*);
int main_plot_cofseq(int, char**, int&, const char*);
int main_plot_pi(int, char**, int&, const char*);
int main_name(int, char**, int&, const char*);

int main_ut(int argc, char** argv, int& index, const char* desc);
int main_mul(int argc, char** argv, int& index, const char* desc);
int main_copy(int argc, char** argv, int& index, const char* desc);

int main_test(int argc, char** argv, int& index, const char* desc);

int main(int argc, char** argv)
{
    Logger::SetOutMain("ss.log");
    Logger::LogCmd(argc, argv);

    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"reset", "Initialize the category", main_reset},
        {"reset_cofseq", "Initialize the cofseq tables", main_reset_cofseq},
        //{"reset_pi", "Initialize the homotopy tables", main_reset_pi},
        {"checkpoint", "Make a checkpoint", main_checkpoint},
        {"rollback", "Roll back to a checkpoint", main_rollback},
        {"resetfrom", "Initialize from another category", main_resetfrom},
        {"add_diff", "Manually input a differential into the spectral sequences", main_add_diff},
        {"add_diff_from_file", "Input differentials from a file", main_add_diff_from_file},
        {"add_diff_from_log", "Input differentials from a log database", main_add_diff_from_log},
        {"add_cofseq_diff", "Manually input a differential into the cofiber sequence spectral sequence", main_add_cofseq_diff},
        {"deduce", "Deduce differentials and extensions", main_deduce},
        {"migrate", "Migrate data between categories", main_migrate},
        {"migrate_from_logs", "Migrate data from logs", main_migrate_from_logs},
        //{"migrate_pi", "Migrate homotopy between categories", main_migrate_pi},
        {"plot_ss", "Generate the json data for plotting ss", main_plot_ss},
        {"plot_cofseq", "Generate the json data for plotting cofiber sequences", main_plot_cofseq},
        {"plot_pi", "Compute the multiplications by hi", main_plot_pi},
        {"name", "Manage generator names", main_name},
        {"ut", "Utilities", main_ut},
        {"mul", "Display the product", main_mul},
        {"copy", "Copy a category", main_copy},
        {"test", "Test", main_test},
    };
    int index = 1;
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, "Manage spectral sequences and homotopy groups", VERSION, subcmds))
        return error;

    Logger::LogTime(timer.print2str());
    return 0;
}
