#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "main.h"
#include <array>
#include <chrono>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <map>
#include <regex>
#include <thread>
#include <unordered_set>

using json = nlohmann::json;

/* Read CPU stat via /proc/stat */
void GetProcStat(long long& total, long long& work)
{
    std::ifstream file("/proc/stat");
    std::string line;
    while (std::getline(file, line)) {
        std::string cpu_stat = line.substr(5);
        std::vector<std::string> cpu_stat_split = myio::split(cpu_stat, ' ');
        std::vector<long long> cpu_stat_int;
        for (auto& stat : cpu_stat_split)
            cpu_stat_int.push_back(std::stoll(stat));

        auto user = cpu_stat_int[0];
        auto nice = cpu_stat_int[1];
        auto system = cpu_stat_int[2];
        auto idle = cpu_stat_int[3];
        auto iowait = cpu_stat_int[4];
        auto irq = cpu_stat_int[5];
        auto softirq = cpu_stat_int[6];
        auto steal = cpu_stat_int[7];
        auto guest = cpu_stat_int[8];
        auto guest_nice = cpu_stat_int[9];

        total = user + nice + system + idle + iowait + irq + softirq + steal + guest + guest_nice;
        work = user + nice + system;

        break;
    }
}

double GetCpuUsage()
{
    long long total1, work1, total2, work2;
    GetProcStat(total1, work1);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    GetProcStat(total2, work2);
    return (work2 - work1) * 1.0 / (total2 - total1);
}

double GetMemUsage()
{
    std::ifstream file("/proc/meminfo");
    std::string line;
    std::map<std::string, long long> meminfo;
    int i = 0;
    while (std::getline(file, line) && i++ < 10) {
        std::vector<std::string> meminfo_split = myio::split(line, ':');
        meminfo_split[1].pop_back();
        meminfo_split[1].pop_back();
        meminfo[meminfo_split[0]] = std::stoll(meminfo_split[1]);
    }
    return (meminfo["MemTotal"] - meminfo["MemFree"] - meminfo["Buffers"] - meminfo["Cached"]) * 1.0 / meminfo["MemTotal"];
}

uint64_t GetDiskAvailable()
{
    const std::filesystem::space_info si = std::filesystem::space(".");
    return si.available;
}

std::string GetThisPid()
{
    /* Get the pid of this program */
    std::string pid;
    std::ifstream file("/proc/self/stat");
    std::string line;
    std::getline(file, line);
    std::vector<std::string> line_split = myio::split(line, ' ');
    pid = line_split[0];
    return pid;
}

/* return {pid: cmd} */
std::map<std::string, std::string> GetRunningAdams()
{
    std::map<std::string, std::string> running_adams;
    std::regex is_Adams("^./Adams");
    std::regex is_num("^[0-9]+$");
    std::smatch match;
    auto this_pid = GetThisPid();
    for (const auto& entry : std::filesystem::directory_iterator("/proc")) {
        if (entry.is_directory()) {
            std::string filename_pid = entry.path().filename().string();
            if (std::regex_search(filename_pid, match, is_num); match[0].matched) {
                std::string filepath = entry.path().string();
                std::string filenameCmd = fmt::format("{}/cmdline", filepath);
                std::ifstream fileCmd(filenameCmd);
                std::string line;
                if (std::getline(fileCmd, line)) {
                    if (std::regex_search(line, match, is_Adams); match[0].matched) {
                        auto cmd = myio::join(" ", myio::split(line, '\0'));
                        cmd = cmd.substr(0, cmd.size() - 1);
                        if (filename_pid != this_pid)
                            running_adams[cmd] = filename_pid;
                    }
                }
            }
        }
    }
    return running_adams;
}

/* return if some Adams instance is already running  */
bool IsAdamsRunning(const std::string& cmd_prefix)
{
    std::regex is_Adams("^./Adams");
    std::regex is_num("^[0-9]+$");
    std::smatch match;
    auto this_pid = GetThisPid();
    for (const auto& entry : std::filesystem::directory_iterator("/proc")) {
        if (entry.is_directory()) {
            std::string filename_pid = entry.path().filename().string();
            if (std::regex_search(filename_pid, match, is_num); match[0].matched) {
                std::string filepath = entry.path().string();
                std::string filenameCmd = fmt::format("{}/cmdline", filepath);
                std::ifstream fileCmd(filenameCmd);
                std::string line;
                if (std::getline(fileCmd, line)) {
                    if (std::regex_search(line, match, is_Adams); match[0].matched) {
                        auto cmd = myio::join(" ", myio::split(line, '\0'));
                        cmd = cmd.substr(0, cmd.size() - 1);
                        if (filename_pid != this_pid && myio::starts_with(cmd, cmd_prefix)) {
                            /* Get the working dir of pid in linux */
                            std::string filenameCwd = fmt::format("{}/cwd", filepath);
                            if (std::filesystem::read_symlink(filenameCwd) == std::filesystem::current_path())
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

std::map<std::string, json::json_pointer> GetTasks(const json& js)
{
    std::map<std::string, json::json_pointer> tasks;
    auto tasks_flat = js.at("tasks").flatten();
    for (auto it : tasks_flat.items()) {
        if (myio::ends_with(it.key(), "cmd")) {
            auto ptr = json::json_pointer(it.key()).parent_pointer();
            tasks[ptr.back()] = ptr;
        }
    }
    return tasks;
}

void replace_all(std::string& str, const std::string& from, const std::string& to)
{
    std::string::size_type n = 0;
    while ((n = str.find(from, n)) != std::string::npos) {
        str.replace(n, from.size(), to);
        n += to.size();
    }
}

void tpl_expand(nlohmann::json& js)
{
    /* expand js.at("templates") */
    auto& js_templates = js.at("templates");
    auto js_templates_flat = js_templates.flatten();
    std::unordered_set<std::string> tpl_keys;
    for (auto it : js_templates_flat.items()) {
        auto& key = it.key();
        if (key.find("tpl123") != std::string::npos)
            tpl_keys.insert(json::json_pointer(key).parent_pointer().to_string());
    }
    std::string template_str;
    for (auto& str_ptr : tpl_keys) {
        auto ptr = json::json_pointer(str_ptr);
        auto& tpl = js_templates.at(ptr);
        std::vector<std::string> params;
        for (auto& k : tpl)
            params.push_back(k.get<std::string>());
        MyException::Assert(!params.empty(), "params should be nonempty");
        if (js_templates.contains(params[0])) {
            auto& template_ = js_templates.at(params[0]);
            template_str = template_.at("json").dump();
            for (size_t i = 1; i < params.size(); i++) {
                auto k = template_.at("params")[i - 1].get<std::string>();
                replace_all(template_str, fmt::format("{{{}}}", k), params[i]);
            }
            js_templates.at(ptr.parent_pointer()) = json::parse(template_str);
        }
        else {
            fmt::print("Template {} not found\n", params[0]);
            throw MyException(0xb1583f9e, "template not found");
        }
    }

    std::vector<std::string> tpl_names;
    for (auto it : js_templates.items())
        tpl_names.push_back(it.key());
    for (auto& k : tpl_names) {
        js_templates.at(k)["parse"] = js_templates.at(k).at("json").dump();
    }

    /* expand all */
    auto js_flat = js.flatten();
    tpl_keys.clear();
    for (auto it : js_flat.items()) {
        auto& key = it.key();
        if (key.find("tpl123") != std::string::npos)
            tpl_keys.insert(json::json_pointer(key).parent_pointer().to_string());
    }
    for (auto& str_ptr : tpl_keys) {
        auto ptr = json::json_pointer(str_ptr);
        auto& tpl = js.at(ptr);
        std::vector<std::string> params;
        for (auto& k : tpl)
            params.push_back(k.get<std::string>());
        MyException::Assert(!params.empty(), "params should be nonempty");
        if (js_templates.contains(params[0])) {
            auto& template_ = js_templates.at(params[0]);
            template_str = template_.at("parse").get<std::string>();
            for (size_t i = 1; i < params.size(); i++) {
                auto k = template_.at("params")[i - 1].get<std::string>();
                replace_all(template_str, fmt::format("{{{}}}", k), params[i]);
            }
            js.at(ptr.parent_pointer()) = json::parse(template_str);
        }
        else {
            fmt::print("Template {} not found\n", params[0]);
            throw MyException(0xb1583f9e, "template not found");
        }
    }
}

struct FinishedTasks
{
    std::unordered_set<std::string> set;
    std::vector<std::string> list;
    void push(const std::string& t)
    {
        set.insert(t);
        list.push_back(t);
    }
};

FinishedTasks GetFinishedTasksFromJson()
{
    FinishedTasks finished_tasks;
    if (myio::FileExists("tasks_status.json")) {
        auto js = myio::load_json("tasks_status.json");
        for (auto& t : js.at("finished_tasks"))
            finished_tasks.push(t.get<std::string>());
    }
    return finished_tasks;
}

/* Return tasks that are not finished but all their prerequisites are finished */
std::vector<std::string> GetPreparedTasks(const json& js, const std::map<std::string, json::json_pointer>& tasks, const std::unordered_set<std::string>& finished_tasks)
{
    std::vector<std::string> prepared_tasks;
    for (auto it : tasks) {
        if (!ut::has(finished_tasks, it.first)) {
            bool prepared = true;
            if (js.at("tasks").at(it.second).contains("pre")) {
                for (auto& pre : js.at("tasks").at(it.second).at("pre")) {
                    if (!ut::has(finished_tasks, pre.get<std::string>())) {
                        prepared = false;
                        break;
                    }
                }
            }
            if (prepared)
                prepared_tasks.push_back(it.first);
        }
    }
    return prepared_tasks;
}

int get_db_t_max(const myio::Database& db);
int get_db_metadata_int(const myio::Database& db, std::string_view key);
std::string get_db_metadata_str(const myio::Database& db, std::string_view key);

bool IsFinished(const std::smatch& match, const char* postfix)
{
    auto cw = match[1].str();
    auto t_max = std::stoi(match[2].str());
    auto filename = fmt::format("{}_{}.db", cw, postfix);
    if (!myio::FileExists(filename))
        return false;
    myio::Database db(filename);
    return get_db_t_max(db) >= t_max;
}

bool IsFinished(const std::string& task)
{
    std::regex is_Adams_res_regex("^res (\\w+) ([0-9]+)$"); /* match example: ./Adams res S0 200 */
    std::smatch match;
    if (std::regex_search(task, match, is_Adams_res_regex); match[0].matched)
        return IsFinished(match, "Adams_res");

    std::regex is_Adams_prod_regex("^prod (\\w+) ([0-9]+)$"); /* match example: ./Adams prod S0 200 */
    if (std::regex_search(task, match, is_Adams_prod_regex); match[0].matched)
        return IsFinished(match, "Adams_res_prod");

    std::regex is_Adams_export_regex("^export (\\w+) ([0-9]+)$"); /* match example: ./Adams export S0 200 */
    if (std::regex_search(task, match, is_Adams_export_regex); match[0].matched)
        return IsFinished(match, "AdamsSS");

    std::regex is_Adams_map_res_regex("^map_res (\\w+) (\\w+) ([0-9]+)$"); /* match example: ./Adams map_res C2 S0 200 */
    if (std::regex_search(task, match, is_Adams_map_res_regex); match[0].matched) {
        auto cw1 = match[1].str();
        auto cw2 = match[2].str();
        std::string from, to;
        auto t_max = std::stoi(match[3].str());
        int t_max_cw1, t_max_cw2, t_max_map, sus, fil;
        {
            auto filename = fmt::format("map_Adams_res_{}__{}.db", cw1, cw2);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_map = get_db_t_max(db);
            if (t_max_map == -3)
                return false;
            sus = get_db_metadata_int(db, "suspension");
            fil = get_db_metadata_int(db, "filtration");
            from = get_db_metadata_str(db, "from");
            to = get_db_metadata_str(db, "to");
            if (myio::starts_with(from, "Error:") || myio::starts_with(to, "Error:"))
                return false;
        }
        {
            auto filename = fmt::format("{}_Adams_res.db", from);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_cw1 = get_db_t_max(db);
        }
        {
            auto filename = fmt::format("{}_Adams_res.db", to);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_cw2 = get_db_t_max(db);
        }
        return t_max_map >= std::min({t_max, t_max_cw1, t_max_cw2 + sus - fil});
    }

    std::regex is_Adams_export_map_regex("^export_map (\\w+) (\\w+) ([0-9]+)$"); /* match example: ./Adams map_res C2 S0 200 */
    if (std::regex_search(task, match, is_Adams_export_map_regex); match[0].matched) {
        auto cw1 = match[1].str();
        auto cw2 = match[2].str();
        std::string from, to;
        auto t_max = std::stoi(match[3].str());
        int t_max_cw1, t_max_cw2, t_max_map, sus, fil;
        {
            auto filename = fmt::format("map_AdamsSS_{}__{}.db", cw1, cw2);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_map = get_db_t_max(db);
            if (t_max_map == -3)
                return false;
            sus = get_db_metadata_int(db, "suspension");
            fil = get_db_metadata_int(db, "filtration");
            from = get_db_metadata_str(db, "from");
            to = get_db_metadata_str(db, "to");
            if (myio::starts_with(from, "Error:") || myio::starts_with(to, "Error:"))
                return false;
        }
        {
            auto filename = fmt::format("{}_AdamsSS.db", from);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_cw1 = get_db_t_max(db);
        }
        {
            auto filename = fmt::format("{}_AdamsSS.db", to);
            if (!myio::FileExists(filename))
                return false;
            myio::Database db(filename);
            t_max_cw2 = get_db_t_max(db);
        }
        return t_max_map >= std::min({t_max, t_max_cw1, t_max_cw2 + sus - fil});
    }

    return false;
}

void AddNewFinishedTasks(std::vector<std::string>& prepared_tasks, FinishedTasks& finished_tasks)
{
    for (auto& t : prepared_tasks) {
        if (IsFinished(t)) {
            fmt::print("Finished: {}\n", t);
            finished_tasks.push(t);
        }
    }
    ut::RemoveIf(prepared_tasks, [&](const std::string& t) { return ut::has(finished_tasks.set, t); });
}

void SchedulerStatus()
{
    /* Get CPU and memory usage */
    double cpu_usage = GetCpuUsage();
    fmt::print("CPU usage: {:.2f}%\n", cpu_usage * 100);
    double mem_usage = GetMemUsage();
    fmt::print("Memory usage: {:.2f}%\n", mem_usage * 100);
    uint64_t disk_available = GetDiskAvailable();
    fmt::print("Disk available: {:.0f} GB\n", disk_available / 1.074e9);

    /* Get running Adams */
    auto running_adams = GetRunningAdams();
    fmt::print("\nRunning tasks:\n");
    for (auto& [cmd, pid] : running_adams)
        fmt::print("  {}\n", cmd);
}

void SchedulerKillAll()
{
    /* Get running Adams */
    auto running_adams = GetRunningAdams();
    fmt::print("Killing the following tasks:\n");
    for (auto& [cmd, pid] : running_adams) {
        if (int error = system(fmt::format("kill {}", pid).c_str()))
            fmt::print("  Error ({}): failed to kill {}\n", error, pid);
        else
            fmt::print("  {} ({})\n", pid, cmd);
    }
}

void SchedulerKillScheduler()
{
    /* Get running Adams */
    auto running_adams = GetRunningAdams();
    fmt::print("Killing the following tasks:\n");
    for (auto& [cmd, pid] : running_adams) {
        if (myio::starts_with(cmd, "./Adams scheduler loop")) {
            if (int error = system(fmt::format("kill {}", pid).c_str()))
                fmt::print("  Error ({}): failed to kill {}\n", error, pid);
            else
                fmt::print("  {} ({})\n", pid, cmd);
        }
    }
}

int GetProperty(const json& js, json::json_pointer ptr, std::string_view key)
{
    do {
        if (js.at("tasks").at(ptr).contains(key))
            return js.at("tasks").at(ptr).at(key).get<int>();
        ptr = ptr.parent_pointer();
    } while (!ptr.empty());
    return 0;
}

int GetPriority(const json& js, json::json_pointer ptr)
{
    return GetProperty(js, ptr, "priority");
}

void update_finished_tasks(const std::vector<std::string>& finished_tasks)
{
    /* Save finished task list */
    json js_status;
    if (myio::FileExists("tasks_status.json"))
        js_status = myio::load_json("tasks_status.json");
    std::ofstream outfile("tasks_status.json");
    js_status["finished_tasks"] = finished_tasks;
    outfile << js_status.dump(4) << std::endl;
}

std::unordered_set<std::string> ConvertToTaskName(const std::map<std::string, std::string>& running_adams)
{
    std::unordered_set<std::string> result;
    std::regex is_Adams_res_regex("^./Adams res (\\w+) ([0-9]+)$"); /* match example: ./Adams res S0 200 */
    std::regex is_Adams_prod_regex("^./Adams prod(?:_mod|) (\\w+)(?:\\s\\w+|) ([0-9]+)$"); /* match example: ./Adams prod S0 200 */
    std::regex is_Adams_export_regex("^./Adams export(?:_mod|) (\\w+)(?:\\s\\w+|) ([0-9]+)$"); /* match example: ./Adams export S0 200 */
    std::regex is_Adams_map_res_regex("^./Adams map_res (\\w+) (\\w+) ([0-9]+)$"); /* match example: ./Adams map_res C2 S0 200 */
    std::regex is_Adams_export_map_regex("^./Adams export_map (\\w+) (\\w+) ([0-9]+)$"); /* match example: ./Adams map_res C2 S0 200 */
    std::smatch match;
    for (auto& [cmd, pid] : running_adams) {
        if (std::regex_search(cmd, match, is_Adams_res_regex); match[0].matched)
            result.insert(fmt::format("res {} {}", match[1].str(), match[2].str()));
        else if (std::regex_search(cmd, match, is_Adams_prod_regex); match[0].matched)
            result.insert(fmt::format("prod {} {}", match[1].str(), match[2].str()));
        else if (std::regex_search(cmd, match, is_Adams_export_regex); match[0].matched)
            result.insert(fmt::format("export {} {}", match[1].str(), match[2].str()));
        else if (std::regex_search(cmd, match, is_Adams_map_res_regex); match[0].matched)
            result.insert(fmt::format("map_res {} {} {}", match[1].str(), match[2].str(), match[3].str()));
        else if (std::regex_search(cmd, match, is_Adams_export_map_regex); match[0].matched)
            result.insert(fmt::format("export_map {} {} {}", match[1].str(), match[2].str(), match[3].str()));
    }
    return result;
}

int SchedulerRunOnce(const json& js)
{
    /* Get CPU and memory usage */
    double cpu_usage = GetCpuUsage();
    fmt::print("CPU usage: {:.2f}%\n", cpu_usage * 100);
    double mem_usage = GetMemUsage();
    fmt::print("Memory usage: {:.2f}%\n", mem_usage * 100);
    uint64_t disk_available = GetDiskAvailable();
    fmt::print("Disk available: {:.0f} GB\n\n", disk_available / 1.074e9);


    auto running_adams = GetRunningAdams();
    auto num_tasks_limit = myio::get(js, "num_tasks_limit", 64);
    if (running_adams.size() > num_tasks_limit) {
		fmt::print("Too many tasks running. kill {} ({})\n", running_adams.begin()->second, running_adams.begin()->first);
		if (int error = system(fmt::format("kill {}", running_adams.begin()->second).c_str()))
			return error;
		return 0;
    }
    else if (running_adams.size() == num_tasks_limit) {
		fmt::print("Too many tasks running. No new task added.\n");
		return 0;
	}
    if (mem_usage > js.at("MEM_limit")) {
        if (running_adams.size() > 0) {
            fmt::print("Memory usage exceeds limit. kill {} ({})\n", running_adams.begin()->second, running_adams.begin()->first);
            if (int error = system(fmt::format("kill {}", running_adams.begin()->second).c_str()))
                return error;
            return 0;
        }
        else {
            fmt::print("Memory usage exceeds limit. No running tasks.\n");
            return 0;
        }
    }

    if (mem_usage > js.at("MEM_threshold")) {
        fmt::print("Memory usage is too high to add another task.\n");
        return 0;
    }

    if (cpu_usage > js.at("CPU_threshold")) {
        fmt::print("cpu usage is too high to add another task.\n");
        return 0;
    }

    if (disk_available / 1.074e9 < js.at("DISK_threshold")) {
        fmt::print("Disk available is too low to add another task.\n");
        return 0;
    }

    auto tasks = GetTasks(js);
    auto finished_tasks = GetFinishedTasksFromJson();
    auto prepared_tasks = GetPreparedTasks(js, tasks, finished_tasks.set);
    AddNewFinishedTasks(prepared_tasks, finished_tasks);
    update_finished_tasks(finished_tasks.list);
    auto running_tasks = ConvertToTaskName(running_adams);
    ut::RemoveIf(prepared_tasks, [&](const std::string& t) { return ut::has(running_tasks, t); });
    ut::RemoveIf(prepared_tasks, [&](const std::string& t) { return GetProperty(js, tasks.at(t), "disabled") == 1; });
    std::sort(prepared_tasks.begin(), prepared_tasks.end(), [&](const std::string& t1, const std::string& t2) { return GetPriority(js, tasks.at(t1)) > GetPriority(js, tasks.at(t2)); });

    /* Run all blocking tasks and one nohup task */
    bool nohup_task = false;
    for (auto& task : prepared_tasks) {
        auto cmd = js.at("tasks").at(tasks.at(task)).at("cmd").get<std::string>();
        if (!myio::starts_with(cmd, "nohup")) {
            fmt::print("Run: {}\n", cmd);
            std::fflush(stdout);
            if (int error = system(cmd.c_str()))
                return error;
            finished_tasks.push(task);
            update_finished_tasks(finished_tasks.list);
        }
        else if (!nohup_task) {
            fmt::print("Run: {}\n", cmd);
            std::fflush(stdout);
            if (int error = system(cmd.c_str()))
                return error;
            nohup_task = true;
        }
    }
    if (!nohup_task)
        fmt::print("No available nohup task to run.\n");

    fflush(stdout);
    return 0;
}

int main_scheduler_status(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SchedulerStatus();
    return 0;
}

int main_scheduler_killall(int argc, char** argv, int& index, const char* desc)
{
    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SchedulerKillAll();
    return 0;
}

int main_scheduler_killsch(int argc, char** argv, int& index, const char* desc)
{
    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    SchedulerKillScheduler();
    return 0;
}

int main_scheduler_run_once(int argc, char** argv, int& index, const char* desc)
{
    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    if (!myio::FileExists("tasks.json")) {
        fmt::print("tasks.json not found\n");
        return -1;
    }
    auto js = myio::load_json("tasks.json");
    tpl_expand(js);
    if (int error = SchedulerRunOnce(js)) {
        fmt::print("SchedulerRunOnce failed\n");
        return error;
    }
    return 0;
}

int main_scheduler_loop(int argc, char** argv, int& index, const char* desc)
{
    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
        
/* Prevent double run on linux */
#ifdef __linux__
    if (IsAdamsRunning(fmt::format("./Adams scheduler loop"))) {
        fmt::print("Error: ./Adams scheduler loop is already running.\n");
        return -1;
    }
#endif

    if (!myio::FileExists("tasks.json")) {
        fmt::print("tasks.json not found\n");
        return -2;
    }
    while (true) {
        bench::Timer timer;
        auto js = myio::load_json("tasks.json");
        tpl_expand(js);
        if (int error = SchedulerRunOnce(js)) {
            fmt::print("SchedulerRunOnce failed\n");
            return error;
        }
        timer.print("Scheduler");
        fmt::print("\n----------------------------\n\n");
        fflush(stdout);
        std::this_thread::sleep_for(std::chrono::seconds(js.at("interval").get<int>()));
    }
    return 0;
}

int main_scheduler(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"status", "Print tasks status", main_scheduler_status},
        {"run_once", "Run scheduler once", main_scheduler_run_once},
        {"loop", "Run scheduler infinitely", main_scheduler_loop},
        {"killall", "Kill all running Adams tasks", main_scheduler_killall},
        {"killsch", "Kill scheduler loop", main_scheduler_killsch}
    };
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;
    return 0;
}

int main_test(int argc, char** argv, int& index, const char* desc)
{
    std::string cmd_prefix;
    myio::CmdArg1d args = {{"cmd_prefix", &cmd_prefix}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    fmt::print("IsAdamsRunning({}) = {}\n", cmd_prefix, IsAdamsRunning(cmd_prefix));
    return 0;
}