#include "main.h"
#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main_test(int, char**, int&, const char*);

int main(int argc, char** argv)
{
    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"test", "test", main_test}
    };
    int index = 1;
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, "Build R-motivic A-resolutions and chain maps.", VERSION, subcmds))
        return error;
    return 0;
}