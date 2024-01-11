#include "algebras/steenrod.h"

namespace steenrod {

/* Secondary steenrod algebra in the Milnor basis.
 * When 0<=k<l, this is $Y_{k,l}$.
 * When k=-1, this is Sq(m).
 * When k=-2, this is 2Sq(m).
 */
class MMilnorSec
{
public:
    MMilnor m;
    int k, l;

public:
    MMilnorSec() : m(), k(-1), l(0) {}
};

}  // namespace steenrod
