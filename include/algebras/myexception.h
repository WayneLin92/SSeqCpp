#ifndef MYEXCEPTION_H
#define MYEXCEPTION_H

#include <string>


/* A custom Exception class */
class MyException : public std::exception
{
public:
	MyException(unsigned int error_id, const char* message) : id_(error_id), msg_(message) {}
	MyException(unsigned int error_id, const std::string& message) : id_(error_id), msg_(message) {}
	virtual ~MyException() noexcept {}
	virtual const char* what() const noexcept { return msg_.c_str(); }
	unsigned int id() const noexcept { return id_; }
protected:
	unsigned int id_;
	std::string msg_;
};


#endif /* MYEXCEPTION_H */