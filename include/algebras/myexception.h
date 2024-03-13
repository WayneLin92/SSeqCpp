#ifndef MYEXCEPTION_H
#define MYEXCEPTION_H

#ifndef NDEBUG
#define MYDEBUG
#else
//#define MYDEBUG 1
#endif

#include <string>

/* A custom Exception class */
class MyException : public std::exception
{
public:
    MyException(unsigned int error_id, const char* message) : id_(error_id), msg_(message) {}
    MyException(unsigned int error_id, const std::string& message) : id_(error_id), msg_(message) {}
    virtual ~MyException() noexcept {}
    virtual const char* what() const noexcept
    {
        return msg_.c_str();
    }
    unsigned int id() const noexcept
    {
        return id_;
    }

    static void Assert(bool statement, const char* message);
    static void Assert(bool statement, const std::string& message);

protected:
    unsigned int id_;
    std::string msg_;
};

#endif /* MYEXCEPTION_H */