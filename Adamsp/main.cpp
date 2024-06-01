#include "main.h"
#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <iostream>

int main_test_associativity(int, char**, int&, const char*);
int main_mul(int, char**, int&, const char*);
int main_mul_v2(int, char**, int&, const char*);


int main(int argc, char** argv)
{
    bench::Timer timer;
    myio::SubCmdArg1d subcmds = {
        {"test_assoc", "Test associativity", main_test_associativity},
        {"mul", "Multiply two Milnor basis", main_mul},
        {"mul_v2", "Multiply two Milnor basis using (c: E: R:) formatted input", main_mul_v2}
    };
    int index = 1;
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, "Build R-motivic A-resolutions and chain maps.", VERSION, subcmds))
        return error;
    return 0;
}