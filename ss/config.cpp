/* @file config.cpp
 * Load settings from json
 */

#include "json.h"
#include "main.h"
#include <fstream>

std::vector<std::string> GetDbNames(const std::string& selector)
{
    using json = nlohmann::json;
    std::vector<std::string> result;

    json js;
    {
        std::ifstream ifs("ss.json");
        if (ifs.is_open())
            ifs >> js;
        else
            throw MyException(0xb8525e9bU, "File ss.json found");
    }
    try {
        json& databases = js["databases"][selector];
        for (auto& db : databases)
            result.push_back(db["path"].get<std::string>());
    }
    catch (nlohmann::detail::exception e) {
        std::cout << "Error:" << e.what() << '\n';
        throw e;
    }
    return result;
}
