#include "myio.h"
#include "myexception.h"
#include "utility.h"
#include <fmt/ranges.h>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <sys/stat.h>

/*********** FUNCTIONS **********/

namespace myio {

std::string join(const std::string& sep, const string1d& strs)
{
    return fmt::format("{}", fmt::join(strs.begin(), strs.end(), sep));
}

/* split comma-delimited string */
std::vector<std::string> split(const std::string& str, char delim)
{
    std::vector<std::string> result;
    if (str.empty())
        return result;
    std::stringstream ss(str);
    while (ss.good()) {
        std::string substr;
        std::getline(ss, substr, delim);
        result.push_back(substr);
    }
    return result;
}

/*
 * Consume and ignore string `pattern` from istream.
 * Set badbit error if pattern is not matched.
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
    for (int i = 0; i < 5; ++i) {
        fmt::print("Input to confirm: [Y/N]\n");
        std::cin >> input;
        if (input == "Y" || input == "y")
            return true;
        else if (input == "N" || input == "n")
            return false;
        else
            fmt::print("Invalid input!\n");
    }
    fmt::print("Use default answer No.\n");
    return false;
}

bool FileExists(const std::string& filename)
{
    std::ifstream f(filename);
    return f.good();
}

void AssertFileExists(const std::string& filename)
{
    if (!FileExists(filename)) {
        fmt::print("Error: File {} does not exist.\n", filename);
        throw MyException(0x821119ec, "File does not exists.");
    }
}

void AssertFolderExists(const std::string& foldername)
{
    struct stat sb;
    if (stat(foldername.c_str(), &sb) != 0) {
        fmt::print("Error: Folder {} does not exist.\n", foldername);
        throw MyException(0xcfb27b11, "Folder does not exists.");
    }
}

template <typename T>
int TplLoadArg(int argc, char** argv, int& index, const char* name, T* x, bool optional)
{
    if (index < argc) {
        if constexpr (std::is_same<T, std::string>::value) {
            *x = argv[index];
            ++index;
            return 0;
        }
        std::istringstream ss(argv[index]);
        if (!(ss >> *x)) {
            fmt::print("Invalid: {}={}\n", name, argv[index]);
            return -index;
        }
        ++index;
        return 0;
    }
    else if (!optional) {
        fmt::print("Mising argument <{}>\n", name);
        return -index;
    }
    return 0;
}

std::string myio::CmdArg::StrValue()
{
    if (std::holds_alternative<int*>(value))
        return std::to_string(*std::get<int*>(value));
    else if (std::holds_alternative<double*>(value))
        return std::to_string(*std::get<double*>(value));
    else if (std::holds_alternative<std::string*>(value))
        return *std::get<std::string*>(value);
    else if (std::holds_alternative<std::vector<std::string>*>(value)) {
        auto& strs = *std::get<std::vector<std::string>*>(value);
        return join(",", strs);
    }
    return "?";
}

void PrintHelp(const std::string& cmd, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args)
{
    fmt::print("{}\nUsage:\n  {}", description, cmd);
    for (auto& x : args)
        fmt::print(" <{}>", x.name);
    for (auto& x : op_args)
        fmt::print(" [{}]", x.name);
    fmt::print("\n\nDefault values:\n");
    for (auto& x : op_args)
        fmt::print("  {} = {}\n", x.name, x.StrValue());
    fmt::print("\n{}\n", version);
}

int LoadCmdArgs_(int argc, char** argv, int& index, CmdArg1d& args, bool optional)
{
    for (size_t i = 0; i < args.size(); ++i) {
        if (std::holds_alternative<int*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<int*>(args[i].value), optional))
                return error;
        }
        else if (std::holds_alternative<double*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<double*>(args[i].value), optional))
                return error;
        }
        else if (std::holds_alternative<std::string*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<std::string*>(args[i].value), optional))
                return error;
        }
        else if (std::holds_alternative<std::vector<std::string>*>(args[i].value)) {
            if (i == args.size() - 1) {
                auto& strs = *std::get<std::vector<std::string>*>(args[i].value);
                while (index < argc)
                    strs.push_back(argv[index++]);
            }
            else
                throw MyException(0x58f9c91a, "Multiple-string command-line argument should be the last");
        }
        else if (std::holds_alternative<std::map<std::string, std::vector<std::string>>*>(args[i].value)) {
            if (i == args.size() - 1) {
                auto& strs = *std::get<std::map<std::string, std::vector<std::string>>*>(args[i].value);
                std::regex key_value("(\\w+)=((?:\\w|,|/|\\)+)");
                std::smatch match;
                while (index < argc) {
                    std::string arg(argv[index++]);
                    if (std::regex_search(arg, match, key_value); match[0].matched)
                        strs[match[1].str()] = split(match[2].str());
                    else {
                        fmt::print("Need key=value pairs as arguments\n");
                        throw MyException(0x238f3e02, "Need key=value pairs as arguments");
                    }
                }
            }
            else {
                fmt::print("Multiple-string command-line argument should be the last\n");
                throw MyException(0x58f9c91a, "Multiple-string command-line argument should be the last");
            }
        }
        else
            return -1023674697;
    }
    return 0;
}

int LoadCmdArgs(int argc, char** argv, int& index, const char* program, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args)
{
    if (index < argc && strcmp(argv[index], "-h") == 0) {
        auto cmd = fmt::format("{} {}", program, fmt::join(argv + 1, argv + index, " "));
        PrintHelp(cmd, description, version, args, op_args);
        return 1;
    }
    if (int error = LoadCmdArgs_(argc, argv, index, args, false))
        return error;
    if (int error = LoadCmdArgs_(argc, argv, index, op_args, true))
        return error;

    return 0;
}

void PrintHelp(const std::string& cmd, const char* description, const char* version, SubCmdArg1d& cmds)
{
    fmt::print("{}\nUsage:\n  {} <cmd> ...\n\ncmd can be one of the following:\n", description, cmd);
    for (auto& cmd : cmds)
        fmt::print("  {}: {}\n", cmd.name, cmd.description);
    fmt::print("\n{}\n", version);
}

int LoadSubCmd(int argc, char** argv, int& index, const char* program, const char* description, const char* version, SubCmdArg1d& cmds)
{
    if (index < argc && strcmp(argv[index], "-h") == 0) {
        std::string cmd = program;
        if (index > 1)
            cmd += fmt::format(" {}", fmt::join(argv + 1, argv + index, " "));
        PrintHelp(cmd, description, version, cmds);
        return 1;
    }
    std::string arg_cmd;
    if (int error = TplLoadArg(argc, argv, index, "cmd", &arg_cmd, false))
        return error;

    for (auto& cmd : cmds) {
        if (arg_cmd == cmd.name) {
            return cmd.f(argc, argv, index, cmd.description);
        }
    }
    fmt::print("Invalid cmd={}\n", arg_cmd);
    return -index;
}

nlohmann::json load_json(const std::string& file_name)
{
    AssertFileExists(file_name);
    std::ifstream ifs(file_name);
    try {
        return nlohmann::json::parse(ifs, nullptr, true, true);
    }
    catch (nlohmann::detail::parse_error& e) {
        fmt::print("json error: {}: {:#x} - {}\n", file_name, e.id, e.what());
        throw e;
    }
}

}  // namespace myio