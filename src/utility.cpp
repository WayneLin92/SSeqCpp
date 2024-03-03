#include "utility.h"
#include <fmt/format.h>
#include "myexception.h"
#include <fmt/core.h>

void MyException::Assert(bool statement, const char* message)
{
    if (!statement) {
        fmt::print("Assert failed: {}\n", message);
        throw MyException(0xd0dec985, "Assert failed");
    }
}

void MyException::Assert(bool statement, const std::string& message)
{
    if (!statement) {
        fmt::print("Assert failed: {}\n", message);
        throw MyException(0xd0dec985, "Assert failed");
    }
}
