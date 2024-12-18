#include "myio.h"
#include "myexception.h"
#include "utility.h"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <sys/stat.h>

inline const char* COMPILE_DATE = __DATE__;  //// TODO: move to caller
inline const char* COMPILE_TIME = __TIME__;

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

void AssertFileExists(const std::string& filename, const std::source_location location)
{
    if (!FileExists(filename))
        throw RunTimeError(fmt::format("File {} does not exist.\n", filename), location);
}

void AssertFolderExists(const std::string& foldername, const std::source_location location)
{
    struct stat sb;
    if (stat(foldername.c_str(), &sb) != 0)
        throw RunTimeError(fmt::format("Folder {} does not exist.\n", foldername), location);
}

/* x = argv[index] */
template <typename T>
[[nodiscard]] int TplLoadArg(int argc, char** argv, int& index, const char* name, T* x)
{
    if (index < argc) {
        if constexpr (std::is_same<T, std::string>::value) {
            *x = argv[index];
            ++index;
            return 0;
        }
        else {
            std::istringstream ss(argv[index]);
            if (!(ss >> *x)) {
                fmt::print("Invalid: {}={}\n", name, argv[index]);
                return -index;
            }
            ++index;
            return 0;
        }
    }
    else {
        fmt::print("Mising argument <{}>\n", name);
        return -index;
    }
}

/* x = argv[index + i] for some i >= 0 */
template <typename T>
[[nodiscard]] int TplLoadOpArg(int argc, char** argv, int index, const char* name, T* x)
{
    auto len_name = std::strlen(name);
    for (; index < argc; ++index) {
        if (std::strlen(argv[index]) >= len_name + 1 && strncmp(argv[index], name, len_name) == 0 && argv[index][len_name] == '=') {
            auto arg = argv[index] + len_name + 1;
            if constexpr (std::is_same<T, std::string>::value) {
                *x = arg;
                return 1;
            }
            else {
                std::istringstream ss(arg);
                if (!(ss >> *x)) {
                    fmt::print("Invalid: {}={}\n", name, arg);
                    return -index;
                }
                return 1;
            }
        }
    }
    return 0;
}

std::string myio::CmdArg::StrValue()
{
    if (std::holds_alternative<int*>(value))
        return std::to_string(*std::get<int*>(value));
    if (std::holds_alternative<bool*>(value))
        return std::to_string(*std::get<bool*>(value));
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

static void PrintHelp(const std::string& cmd, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args)
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

int1d Deserialize(const std::string& str)  ////
{
    int1d result;
    if (str.empty())
        return result;
    std::stringstream ss(str);
    while (ss.good()) {
        int i;
        ss >> i;
        result.push_back(i);
        if (ss.peek() == ',')
            ss.ignore();
    }
    return result;
}

static int LoadCmdArgs_(int argc, char** argv, int& index, CmdArg1d& args)
{
    for (size_t i = 0; i < args.size(); ++i) {
        if (std::holds_alternative<int*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<int*>(args[i].value)))
                return error;
        }
        else if (std::holds_alternative<double*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<double*>(args[i].value)))
                return error;
        }
        else if (std::holds_alternative<std::string*>(args[i].value)) {
            if (int error = TplLoadArg(argc, argv, index, args[i].name, std::get<std::string*>(args[i].value)))
                return error;
        }
        else if (std::holds_alternative<int1d*>(args[i].value)) {
            std::string nums;
            if (int error = TplLoadArg(argc, argv, index, args[i].name, &nums))
                return error;
            *std::get<int1d*>(args[i].value) = Deserialize(nums);
        }
        else if (std::holds_alternative<string1d*>(args[i].value)) {
            std::string strs;
            if (int error = TplLoadArg(argc, argv, index, args[i].name, &strs))
                return error;
            *std::get<string1d*>(args[i].value) = split(strs);
        }
        else
            return -43;
    }
    return 0;
}

static int LoadOpCmdArgs_(int argc, char** argv, int& index, CmdArg1d& op_args)
{
    int op_count = 0;
    for (size_t i = 0; i < op_args.size(); ++i) {
        if (std::holds_alternative<int*>(op_args[i].value)) {
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, std::get<int*>(op_args[i].value)); rt < 0)
                return rt;
            else
                op_count += rt;
        }
        else if (std::holds_alternative<bool*>(op_args[i].value)) {
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, std::get<bool*>(op_args[i].value)); rt < 0)
                return rt;
            else
                op_count += rt;
        }
        else if (std::holds_alternative<double*>(op_args[i].value)) {
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, std::get<double*>(op_args[i].value)); rt < 0)
                return rt;
            else
                op_count += rt;
        }
        else if (std::holds_alternative<std::string*>(op_args[i].value)) {
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, std::get<std::string*>(op_args[i].value)); rt < 0)
                return rt;
            else
                op_count += rt;
        }
        else if (std::holds_alternative<int1d*>(op_args[i].value)) {
            std::string nums;
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, &nums); rt < 0)
                return rt;
            else
                op_count += rt;
            *std::get<int1d*>(op_args[i].value) = Deserialize(nums);
        }
        else if (std::holds_alternative<string1d*>(op_args[i].value)) {
            std::string strs;
            if (int rt = TplLoadOpArg(argc, argv, index, op_args[i].name, &strs); rt < 0)
                return rt;
            else
                op_count += rt;
            *std::get<string1d*>(op_args[i].value) = split(strs);
        }
        else
            return -143;
    }
    if (index + op_count != argc) {
        fmt::print("unknown option key\n");
        std::fflush(stdout);
        return -144;
    }
    return 0;
}

int ParseArguments(int argc, char** argv, int& index, const char* program, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args)
{
    if (index < argc && strcmp(argv[index], "-h") == 0) {
        auto cmd = fmt::format("{} {}", program, fmt::join(argv + 1, argv + index, " "));
        PrintHelp(cmd, description, version, args, op_args);
        return 1;
    }
    if (int error = LoadCmdArgs_(argc, argv, index, args))
        return error;
    if (int error = LoadOpCmdArgs_(argc, argv, index, op_args))
        return error;

    return 0;
}

static void PrintHelp(const std::string& cmd, const char* description, const char* version, SubCmdArg1d& cmds)
{
    fmt::print("{}\nUsage:\n  {} <cmd> ...\n\n<cmd> can be one of the following:\n", description, cmd);
    for (auto& c : cmds)
        fmt::print("  {}: {}\n", c.name, c.description);
    fmt::print("\nVersion: {} ({} {}) {}\n", version, COMPILE_DATE, COMPILE_TIME, COMPILE_MODE);
}

int ParseSubCmd(int argc, char** argv, int& index, const char* program, const char* description, const char* version, SubCmdArg1d& cmds)
{
    if (index < argc && strcmp(argv[index], "-h") == 0) {
        std::string cmd = program;
        if (index > 1)
            cmd += fmt::format(" {}", fmt::join(argv + 1, argv + index, " "));
        PrintHelp(cmd, description, version, cmds);
        return 1;
    }
    std::string arg_cmd;
    if (int error = TplLoadArg(argc, argv, index, "cmd", &arg_cmd))
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

std::string load_text(const std::string& file_name)
{
    AssertFileExists(file_name);
    std::ifstream ifs(file_name);
    return std::string((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
}

void save_text(const std::string& file_name, const std::string& text)
{
    std::ofstream ofs(file_name);
    ofs << text;
}

}  // namespace myio