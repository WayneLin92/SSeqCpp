#ifndef MYEXCEPTION_H
#define MYEXCEPTION_H

#ifndef NDEBUG
#define MYDEBUG
#define COMPILE_MODE "debug"
#else
//#define MYDEBUG
#define COMPILE_MODE "release"
#endif

#include <string>
#include <source_location>

/* All exception class should be derived from this class */
class MyException : public std::exception
{
public:
    virtual ~MyException() noexcept {}

    static void Assert(bool statement, const char* message, const std::source_location location = std::source_location::current());
    static void Assert(bool statement, const std::string& message, const std::source_location location = std::source_location::current());
};

class ErrorMsg : public MyException
{
protected:
    std::string msg_;

public:
    ErrorMsg(const char* message) : msg_(message) {}
    ErrorMsg(std::string message) : msg_(std::move(message)) {}
    virtual ~ErrorMsg() noexcept {}
    virtual const char* what() const noexcept
    {
        return msg_.c_str();
    }
};

class ErrorIdMsg : public ErrorMsg
{
protected:
    unsigned int id_;

public:
    ErrorIdMsg(unsigned int error_id, const char* message) : id_(error_id), ErrorMsg(message) {}
    ErrorIdMsg(unsigned int error_id, std::string message) : id_(error_id), ErrorMsg(std::move(message)) {}
    virtual ~ErrorIdMsg() noexcept {}
    unsigned int id() const noexcept
    {
        return id_;
    }

};

class RunTimeError : public ErrorMsg
{
public:
    RunTimeError(const char* message, const std::source_location location = std::source_location::current());
    RunTimeError(const std::string& message, const std::source_location location = std::source_location::current());
	virtual ~RunTimeError() noexcept {}
};

#endif /* MYEXCEPTION_H */