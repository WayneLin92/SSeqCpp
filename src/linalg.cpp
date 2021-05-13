#include "linalg.h"
#include "myexception.h"

#ifdef _DEBUG
#include <iostream>
#endif

namespace lina {

array2d GetSpace(const array2d& vectors)
{
	array2d result;
	for (array v : vectors) { /* Create copy on purpose */
		for (const auto& vi : result)
			if (std::binary_search(v.begin(), v.end(), vi[0]))
				v = AddVectors(v, vi);
		if (!v.empty())
			result.push_back(std::move(v));
	}
	return result;
}

array2d& SimplifySpace(array2d& spaceV)
{
	for (size_t i = spaceV.size() - 1; i != -1; i--)
		for (size_t j = 0; j < i; j++)
			if (std::binary_search(spaceV[j].begin(), spaceV[j].end(), spaceV[i][0]))
				spaceV[j] = AddVectors(spaceV[j], spaceV[i]);
	return spaceV;
}

array Residue(array2dIt spaceV_first, array2dIt spaceV_last, array v)
{
	for (auto p_vi = spaceV_first; p_vi != spaceV_last; ++p_vi)
		if (std::binary_search(v.begin(), v.end(), p_vi->front()))
			v = AddVectors(v, *p_vi);
	return v;
}

inline void AddToSpace(array2d& spaceV, const array& v)
{
	array v1 = Residue(spaceV, v);
	if (!v1.empty())
		spaceV.push_back(std::move(v1));
}

void GetInvMap(const array2d& fx, array2d& image, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { int(i) };
		array tgt = fx[i];
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = AddVectors(tgt, image[j]);
				src = AddVectors(src, g[j]);
			}
		}
		if (!tgt.empty()) {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void SetLinearMap(const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { int(i) };
		array tgt = fx[i];
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = AddVectors(tgt, image[j]);
				src = AddVectors(src, g[j]);
			}
		}
		if (tgt.empty())
			AddToSpace(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void SetLinearMapV2(array2dIt x_first, array2dIt x_last, array2dIt fx_first, array2dIt fx_last, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (auto p_x = x_first, p_fx = fx_first; p_x != x_last; ++p_x, ++p_fx) {
		array src = *p_x;
		array tgt = *p_fx;
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = AddVectors(tgt, image[j]);
				src = AddVectors(src, g[j]);
			}
		}
		if (tgt.empty())
			AddToSpace(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void SetLinearMapV2(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { x[i] };
		array tgt = fx[i];
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = AddVectors(tgt, image[j]);
				src = AddVectors(src, g[j]);
			}
		}
		if (tgt.empty())
			AddToSpace(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void SetLinearMapV3(const array2d& x, const array2d& fx, array2d& domain, array2d& f, array2d& image, array2d& g, array2d& kernel)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = x[i];
		array tgt = fx[i];
		for (size_t k = 0; k < domain.size(); ++k) {
			if (std::binary_search(src.begin(), src.end(), domain[k][0])) {
				src = AddVectors(src, domain[k]);
				tgt = AddVectors(tgt, f[k]);
			}
		}
		if (src.empty()) {
			if (!tgt.empty())
				throw MyException(0x67b8b67dU, "conflicting linear map definition");
		}
		else {
			domain.push_back(src);
			f.push_back(tgt);
			for (size_t j = 0; j < image.size(); j++) {
				if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
					tgt = AddVectors(tgt, image[j]);
					src = AddVectors(src, g[j]);
				}
			}
			if (tgt.empty())
				AddToSpace(kernel, src);
			else {
				image.push_back(std::move(tgt));
				g.push_back(std::move(src));
			}
		}
	}
}

array GetImage(array2dIt spaceV_first, array2dIt spaceV_last, array2dIt f_first, array2dIt f_last, array v)
{
	array result;
	for (auto p_Vi = spaceV_first, p_fi = f_first; p_Vi != spaceV_last && !v.empty(); ++p_Vi, ++p_fi)
		if (std::binary_search(v.begin(), v.end(), p_Vi->front())) {
			v = AddVectors(v, *p_Vi);
			result = AddVectors(result, *p_fi);
		}
#ifdef _DEBUG
	if (!v.empty()) {
		throw MyException(0x6a4fe8a1U, "v is not in spaceV");
	}
#endif
	return result;
}

array2d QuotientSpace(const array2d& spaceV, const array2d& spaceW)
{
	array2d quotient;
	size_t dimQuo = spaceV.size() - spaceW.size();
#ifdef _DEBUG
	for (size_t i = 0; i < spaceV.size(); i++)
#else
	for (size_t i = 0; i < spaceV.size() && quotient.size() < dimQuo; i++)
#endif
	{
		auto v1 = Residue(quotient, Residue(spaceW, spaceV[i]));
		if (!v1.empty())
			quotient.push_back(std::move(v1));
	}
#ifdef _DEBUG
	if (quotient.size() != dimQuo) {
		std::cerr << "W is not a subspace of V!\n";
		throw "cec7f701";
	}
#endif
	return quotient;
}

} /* namespace lina */