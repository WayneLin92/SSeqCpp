#include "myio.h"

/*********** FUNCTIONS **********/

namespace myio {
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
}  // namespace myio