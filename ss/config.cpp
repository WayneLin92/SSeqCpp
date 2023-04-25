/* @file config.cpp
 * Load settings from json
 */

#include "json.h"
#include "main.h"
#include "mylog.h"
#include <fstream>

std::vector<std::string> GetDbNames(const std::string& selector, bool log)
{
    using json = nlohmann::json;
    std::vector<std::string> result;

    json js;
    {
        std::ifstream ifs("ss.json");
        if (ifs.is_open())
            ifs >> js;
        else {
            Logger::LogException(0, 0xb8525e9bU, "File ss.json found\n");
            throw MyException(0xb8525e9bU, "File ss.json found");
        }
    }
    try {
        json& diagram = js["diagrams"][selector];
        if (log) {
            auto path_log = fmt::format("{}/{}", diagram["dir"].get<std::string>(), diagram["log"].get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }
        std::vector<std::string> names = {"S0", "C2", "Ceta", "Cnu", "Csigma"};
        for (auto& name : names) {
            auto path = fmt::format("{}/{}", diagram["dir"].get<std::string>(), diagram["spectra"][name].get<std::string>());
            result.push_back(path);
        }
    }
    catch (nlohmann::detail::exception e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
    return result;
}
