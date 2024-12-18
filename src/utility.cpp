#include "utility.h"
#include <fmt/format.h>
#include "myexception.h"
#include <fmt/core.h>

static_assert(std::ranges::range<ut::Iter2d<std::vector<std::vector<int>>>>);

constexpr std::string_view light_yellow = "";

RunTimeError::RunTimeError(const char* message, const std::source_location location)
	: ErrorMsg(fmt::format("File:\n{}({}:{})\nFunction:\n`{}`:\nError:\n{}", location.file_name(), location.line(), location.column(), location.function_name(), message))
{
    fmt::print("\n\033[38;2;200;255;255mFile: \033[0m\n{}({}:{})\n\033[38;2;200;255;255mFunction: \033[0m\n`{}`:\n\033[38;2;200;255;255mError: \033[0m\n{}\n", location.file_name(), location.line(), location.column(), location.function_name(), message);
    fmt::print("\a");
}

RunTimeError::RunTimeError(const std::string& message, const std::source_location location)
	: ErrorMsg(fmt::format("File:\n{}({}:{})\nFunction:\n`{}`:\nError:\n{}", location.file_name(), location.line(), location.column(), location.function_name(), message))
{
    fmt::print("\n\033[38;2;200;255;255mFile: \033[0m\n{}({}:{})\n\033[38;2;200;255;255mFunction: \033[0m\n`{}`:\n\033[38;2;200;255;255mError: \033[0m\n{}\n", location.file_name(), location.line(), location.column(), location.function_name(), message);
    fmt::print("\a");
}

void MyException::Assert(bool statement, const char* message, const std::source_location location)
{
    if (!statement) {
        throw RunTimeError(fmt::format("Assert failed: {}", message), location);
    }
}

void MyException::Assert(bool statement, const std::string& message, const std::source_location location)
{
    if (!statement) {
        throw RunTimeError(fmt::format("Assert failed: {}", message), location);
    }
}

void utility_tmp()
{
    
}