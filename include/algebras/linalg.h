/*****************************************************************
** Linear Algebra Mod 2
**
** All row vectors are compressed
** Spaces are upper triangular matrices with some column ordering
*****************************************************************/

#ifndef LINALG_H
#define LINALG_H

#include <algorithm>
#include <iterator>
#include <vector>

namespace lina {

using int1d = std::vector<int>;
using int2d = std::vector<int1d>;
using int2dIt = int2d::const_iterator;

/* Add two compressed vectors */
inline int1d AddVectors(const int1d& v1, const int1d& v2)
{
    int1d result;
    std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
    return result;
};

/* Return the space spanned by `vectors` */
int2d GetSpace(const int2d& vectors);
int1d GetLeads(const int2d& spaceV);

/* Reduce the space to rref form */
int2d& SimplifySpace(int2d& spaceV);

/* Return a newer v such that (spaceV\\ v) is triangular */
int1d Residue(int2dIt spaceV_first, int2dIt spaceV_last, int1d v);
inline int1d Residue(const int2d& spaceV, int1d v)
{
    return Residue(spaceV.begin(), spaceV.end(), std::move(v));
}

/* g[i] = f^{-1}image[i] */
void GetInvMap(const int2d& fx, int2d& image, int2d& g);
/* Setup a linear map */
void SetLinearMap(const int2d& fx, int2d& image, int2d& kernel, int2d& g);
void SetLinearMapV2(int2dIt x_first, int2dIt x_last, int2dIt fx_first, int2dIt fx_last, int2d& image, int2d& kernel, int2d& g);
void SetLinearMapV2(const int1d& x, const int2d& fx, int2d& image, int2d& kernel, int2d& g);
void SetLinearMapV3(const int2d& x, const int2d& fx, int2d& domain, int2d& f, int2d& image, int2d& g, int2d& kernel);

/* Return f(v) for v\\in V. fi=f(vi) */
int1d GetImage(int2dIt spaceV_first, int2dIt spaceV_last, int2dIt f_first, int2dIt f_last, int1d v);
inline int1d GetImage(const int2d& spaceV, const int2d& f, int1d v)
{
    return GetImage(spaceV.begin(), spaceV.end(), f.begin(), f.end(), std::move(v));
}
/* Return f^{-1}(v) where fi=vi */
int1d GetInvImage(const int2d& spaceV, int1d v);

/* Compute the quotient of linear spaces V/W assuming that W is a subspace of V */
int2d QuotientSpace(const int2d& spaceV, const int2d& spaceW);

}  // namespace lina

#endif /* LINALG_H */