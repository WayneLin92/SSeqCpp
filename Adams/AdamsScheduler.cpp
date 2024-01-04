#include "algebras/myio.h"
#include "algebras/utility.h"
#include "main.h"
#include <array>
#include <filesystem>
#include <fmt/core.h>
#include <map>
#include <regex>

void SchedulerStatus()
{
    for (const auto& entry : std::filesystem::directory_iterator("/proc")) {
        std::string filename = entry.path().filename().string();
        std::string filepath = entry.path().string();

        fmt::print("{}\n", filename);
    }
}

int main_scheduler_status(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";
    myio::string1d options;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SchedulerStatus();
    return 0;
}


int main_scheduler(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"status", "Print tasks", main_scheduler_status},
    };
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;
    return 0;
}