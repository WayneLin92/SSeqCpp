#ifndef MYIO_H
#define MYIO_H

#include "json.h"
#include <cstring>
#include <fmt/format.h>
#include <fmt/core.h>
#include <sstream>
#include <variant>
#include <vector>

namespace myio {

using int1d = std::vector<int>;
using int2d = std::vector<int1d>;
using int3d = std::vector<int2d>;
using string1d = std::vector<std::string>;

std::vector<std::string> split(const std::string& str, char delim);
/* split comma-delimited string */
inline std::vector<std::string> split(const std::string& str)
{
    return split(str, ',');
}

inline bool starts_with(const std::string& str, std::string_view start)
{
    return str.size() >= start.size() && str.compare(0, start.size(), start) == 0;
}

inline bool ends_with(const std::string& str, std::string_view end)
{
    return str.size() >= end.size() && str.compare(str.size() - end.size(), end.size(), end) == 0;
}

template <typename FwdIt, typename FnStr>                                                                                            //// Deprecate
std::string TplStrCont(const char* left, const char* sep, const char* right, const char* empty, FwdIt first, FwdIt last, FnStr str)  // TODO: make a performant version
{
    std::string result;
    if (first == last)
        result += empty;
    else {
        result += left;
        result += str(*first);
        while (++first != last) {
            result += sep;
            result += str(*first);
        }
        result += right;
    }
    return result;
}

std::string join(const std::string& sep, const string1d& strs);

template <typename FwdIt, typename FnStr>
std::string Tpljoin(const std::string& sep, FwdIt first, FwdIt last, FnStr str)
{
    string1d strs;
    for (auto it = first; it != last; ++it)
        strs.push_back(str(*it));
    return join(sep, strs);
}

template <typename Container, typename FnStr>
std::string StrCont(const char* left, const char* sep, const char* right, const char* empty, const Container& cont, FnStr str)
{
    return TplStrCont(left, sep, right, empty, cont.begin(), cont.end(), str);
}

/**
 * Consume and ignore string `pattern` from istream.
 * Set badbit error if pattern is not matched.
 */
void consume(std::istream& sin, const char* pattern);

inline std::istream& operator>>(std::istream& sin, const char* pattern)
{
    consume(sin, pattern);
    return sin;
}

/* Load container from an istream */
template <typename Container>
void load_vector(std::istream& sin, Container& cont, const char* left, const char* sep, const char* right)
{
    cont.clear();
    sin >> std::ws;
    consume(sin, left);
    typename Container::value_type ele;
    while (sin.good()) {
        sin >> ele;
        cont.push_back(ele);
        sin >> std::ws;
        consume(sin, sep);
    }
    if (sin.bad()) {
        sin.clear();
        consume(sin, right);
    }
}

/* Load the container from an istream */
template <typename Container, typename _Fn_load>
void load_vector(std::istream& sin, Container& cont, const char* left, const char* sep, const char* right, _Fn_load load)
{
    cont.clear();
    sin >> std::ws;
    consume(sin, left);
    typename Container::value_type ele;
    while (sin.good()) {
        load(sin, ele);
        cont.push_back(std::move(ele));
        sin >> std::ws;
        consume(sin, sep);
    }
    if (sin.bad()) {
        sin.clear();
        consume(sin, right);
    }
}

inline void load_array(std::istream& sin, int1d& a)
{
    load_vector(sin, a, "(", ",", ")");
}
inline void load_array2d(std::istream& sin, int2d& a)
{
    load_vector(sin, a, "[", ",", "]", load_array);
}
inline void load_array3d(std::istream& sin, int3d& a)
{
    load_vector(sin, a, "{", ",", "}", load_array2d);
}
inline std::istream& operator>>(std::istream& sin, int1d& a)
{
    load_array(sin, a);
    return sin;
}
inline std::istream& operator>>(std::istream& sin, int2d& a)
{
    load_array2d(sin, a);
    return sin;
}
inline std::istream& operator>>(std::istream& sin, int3d& a)
{
    load_array3d(sin, a);
    return sin;
}

bool FileExists(const std::string& filename);
void AssertFileExists(const std::string& filename);
void AssertFolderExists(const std::string& foldername);

/*********************************************************
                 Command line utilities
 *********************************************************/

struct CmdArg
{
    const char* name;
    std::variant<int*, double*, std::string*, std::vector<std::string>*, std::map<std::string, std::vector<std::string>>*> value;
    std::string StrValue();
};
using CmdArg1d = std::vector<CmdArg>;

int LoadCmdArgs(int argc, char** argv, int& index, const char* program, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args);

/* Return true if user inputs Y; Return false if user inputs N */
bool UserConfirm();

using MainFnType = int (*)(int, char**, int&, const char*);

struct SubCmdArg
{
    const char* name;
    const char* description;
    MainFnType f;
};
using SubCmdArg1d = std::vector<SubCmdArg>;

int LoadSubCmd(int argc, char** argv, int& index, const char* program, const char* description, const char* version, SubCmdArg1d& cmds);

/*********************************************************
                     json
 *********************************************************/
nlohmann::json load_json(const std::string& file_name);

}  // namespace myio

/*********************************************************
                    Formatters
 *********************************************************/

template <typename T>
struct fmt::formatter<T, char, std::enable_if_t<std::is_same_v<decltype(T().Str()), std::string>>>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const T& x, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};

#endif /* MYIO_H */