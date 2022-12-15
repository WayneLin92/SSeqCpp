#include "myio.h"
#include "utility.h"

/*********** FUNCTIONS **********/

namespace myio {

std::ofstream Logger::fout1_;
std::ofstream Logger::fout2_;
bool Logger::bInitialized = false;


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

void Logger::Init(const char* filename, const char* filename1)
{
    bInitialized = true;
    fout1_.open(filename, std::ofstream::app);
    fout2_.open(filename1, std::ofstream::app);
    fout1_ << "\nLogging start at " << ut::get_time() << std::endl;
    fout2_ << "\nLogging start at " << ut::get_time() << std::endl;
}

std::ostream& Logger::smart_stream()
{
    return bInitialized ? fout1_ : std::cout;
}

bool UserConfirm()
{
    std::string input;
    while (true) {
        std::cout << "Input to confirm: [Y/N]" << std::endl;
        std::cin >> input;
        if (input == "Y")
            return true;
        else if (input == "N")
            return false;
        else
            std::cout << "Invalid input!\n";
    }
}

}  // namespace myio