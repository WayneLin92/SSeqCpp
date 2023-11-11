#include "main.h"
#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main_res(int, char**, int&, const char*);
int main_map_coh(int, char**, int&, const char*);
int main_map_res(int, char**, int&, const char*);
int main_verify_map(int, char**, int&, const char*);

int main_prod(int, char**, int&, const char*);
int main_prod_mod(int, char**, int&, const char*);
int main_prod_hi(int, char**, int&, const char*);

int main_export(int, char**, int&, const char*);
int main_export_mod(int, char**, int&, const char*);
int main_export_map(int, char**, int&, const char*);

int main_2cell(int, char**, int&, const char*);
int main_res_csv(int, char**, int&, const char*);

int main_status(int, char**, int&, const char*);
int main_verify_status(int, char**, int&, const char*);
int main_ut(int, char**, int&, const char*);

int main(int argc, char** argv)
{
    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"res", "Compute a minimal A-resolution", main_res},
        {"map_coh", "Set up the cohomology map in file", main_map_coh},
        {"map_res", "Compute a chain map between resolutions", main_map_res},
        {"verify_map", "Verify the correctness of a chain map", main_verify_map},
        {"prod", "Compute the multiplications for a ring", main_prod},
        {"prod_mod", "Compute the multiplications for a module", main_prod_mod},
        {"prod_hi", "Compute the multiplications by hi", main_prod_hi},
        {"export", "Export the Adams E2 page of a ring", main_export},
        {"export_mod", "Export the Adams E2 page of a module", main_export_mod},
        {"export_map", "Export the map between Adams E2 pages", main_export_map},
        {"2cell", "Functions for Cofibers of Hopf elements", main_2cell},
        {"res_csv", "Export the resolution data to a csv file", main_res_csv},
        {"status", "Display the computation status in the current directory", main_status},
        {"verify_status", "Display the verification status in the current directory", main_verify_status},
        {"ut", "Utilities", main_ut},
    };
    int index = 1;
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, "Build A-resolutions and chain maps.", VERSION, subcmds))
        return error;
    return 0;
}