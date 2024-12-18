#ifndef MYIO_H
#define MYIO_H

#include "json.h"
#include <cstring>
#include <fmt/core.h>
#include <fmt/format.h>
#include <source_location>
#include <sstream>
#include <variant>
#include <vector>

namespace myio {

using int1d = std::vector<int>;
using int2d = std::vector<int1d>;
using int3d = std::vector<int2d>;
using string1d = std::vector<std::string>;

string1d split(const std::string& str, char delim);
/* split comma-delimited string */
inline string1d split(const std::string& str)
{
    return split(str, ',');
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

/*********************************************************
                        Files
 *********************************************************/

bool FileExists(const std::string& filename);
void AssertFileExists(const std::string& filename, const std::source_location location = std::source_location::current());
void AssertFolderExists(const std::string& foldername, const std::source_location location = std::source_location::current());
std::string load_text(const std::string& file_name);
void save_text(const std::string& file_name, const std::string& text);

/*********************************************************
                 Command line utilities
 *********************************************************/

struct CmdArg
{
    const char* name;
    std::variant<int*, bool*, double*, std::string*, std::vector<int>*, string1d*> value;
    std::string StrValue();
};
using CmdArg1d = std::vector<CmdArg>;

int ParseArguments(int argc, char** argv, int& index, const char* program, const char* description, const char* version, CmdArg1d& args, CmdArg1d& op_args);

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

int ParseSubCmd(int argc, char** argv, int& index, const char* program, const char* description, const char* version, SubCmdArg1d& cmds);

/*********************************************************
                     json
 *********************************************************/
nlohmann::json load_json(const std::string& file_name);

inline int get(const nlohmann::json& js, std::string key, int default_)
{
    return js.contains(key) ? js.at(key).get<int>() : default_;
}

inline std::string get(const nlohmann::json& js, std::string key, const std::string& default_)
{
    return js.contains(key) ? js.at(key).get<std::string>() : default_;
}

}  // namespace myio

/*********************************************************
                 Custom format string
----------------------------------------------------------
 * Example:
 \begin{equation}
 <arg>
 \end{equation}
 *********************************************************/
constexpr size_t Count(const std::string_view str, const std::string_view arg)
{
    size_t occurrences = 0;
    std::string::size_type pos = 0;
    while ((pos = str.find(arg, pos)) != std::string::npos) {
        ++occurrences;
        pos += arg.length();
    }
    return occurrences;
}

constexpr size_t FmtLen(const std::string_view str, const std::string_view arg)
{
    auto num_leftbrace = Count(str, "{");
    auto num_rightbrace = Count(str, "}");
    if (num_leftbrace != num_rightbrace)
        return 1 / size_t(num_leftbrace == num_rightbrace);
    auto num_args = Count(str, arg);
    return str.size() + num_leftbrace * 2 + num_args * 2 + 1 - num_args * arg.length();
}

template <size_t N>
constexpr std::array<char, N> ToFmtStr(const std::string_view str, const std::string_view arg)
{
    std::array<char, N> fmt_str;
    for (size_t i = 0, j = 0; i < str.size(); ++i) {
        if (str[i] == '{') {
            fmt_str[j++] = '{';
            fmt_str[j++] = '{';
        }
        else if (str[i] == '}') {
            fmt_str[j++] = '}';
            fmt_str[j++] = '}';
        }
        else if (str.substr(i, arg.size()) == arg) {
            fmt_str[j++] = '{';
            fmt_str[j++] = '}';
            i += arg.size() - 1;
        }
        else
            fmt_str[j++] = str[i];
    }
    fmt_str[N - 1] = '\0';
    return fmt_str;
}
template <size_t N>
constexpr std::string_view FmtView(const std::array<char, N>& arr)
{
    return std::string_view(arr.data(), N - 1);
}
#define MYFMT(Fmt, MyFmt)                                                         \
    constexpr auto MyFmt##Arr = ToFmtStr<FmtLen(MyFmt, "<arg>")>(MyFmt, "<arg>"); \
    constexpr auto Fmt = FmtView(MyFmt##Arr);


/*********************************************************
                    Formatters
 *********************************************************/

/* This enables formatting for classes with Str() method */
template <typename T>
struct fmt::formatter<T, char, std::enable_if_t<std::is_same_v<decltype(T().Str()), std::string>>>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) const
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const T& x, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};

#endif /* MYIO_H */