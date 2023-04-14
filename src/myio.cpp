#include "myio.h"
#include "utility.h"
#include <iostream>

/*********** FUNCTIONS **********/

namespace myio {

/*
** Consume and ignore string `pattern` from istream.
** Set badbit error if pattern is not matched.
*/
void consume(std::istream& sin, const char* pattern)
{
    size_t i;
    for (i = 0; pattern[i] != '\0' && sin.peek() == int(pattern[i]); ++i)
        sin.ignore();
    if (pattern[i] != '\0')
        sin.setstate(std::ios_base::badbit);
}

bool UserConfirm()
{
    std::string input;
    while (true) {
        fmt::print("Input to confirm: [Y/N]\n");
        std::cin >> input;
        if (input == "Y" || input == "y")
            return true;
        else if (input == "N" || input == "n")
            return false;
        else
            fmt::print("Invalid input!\n");
    }
}

}  // namespace myio