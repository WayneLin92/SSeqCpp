#include "myparser.h"

/*
** Consume and ignore string `pattern` from istream. 
** Set badbit error if pattern is not matched.
*/
void consume(std::istream& sin, const char* pattern)
{
	size_t i;
	for (i = 0; pattern[i] != '\0' && sin.peek() == int(pattern[i]); ++i)
		sin.ignore();
	if (pattern[i] != '\0')
		sin.setstate(std::ios_base::badbit);
}

void load_py_str(std::string& s, std::istream& sin) // partial feature
{
	s.clear();
	sin >> std::ws; consume(sin, "'");
	sin >> std::noskipws;
	char c;
	while (sin.good() && sin.peek() != int('\'')) {
		sin >> c;
		s.push_back(c);
	}
	sin >> std::skipws;
	if (sin.good())
		consume(sin, "'");
}