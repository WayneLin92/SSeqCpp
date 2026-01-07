#ifndef MILNOR_BASIS_HPP
#define MILNOR_BASIS_HPP

#include<vector>
#include<array>
#include<iostream>
#include<cassert>
#include<unordered_set>
#include <stdexcept>

const size_t CACHE_SIZE = 128; //The actual cache size <= pow(2, MATRIX_SIZE-1), set to pow(2, MATRIX_SIZE-1) for safety
const size_t MATRIX_SIZE = 8; //the rows of matrice are stored using uint8_t to exploit every bit of memory. If bigger matrices are needed, uint8_t should be replaced by a larger type
const size_t MAX_DEGREE = 276;

/**
The class defined the hash function required for using unordered set to store summations in product
*/
template<std::size_t N>
struct std::hash<std::array<unsigned int, N>> 
{
	std::size_t operator()(const std::array<unsigned int, N>& a) const 
	{
		std::size_t h = 0;
		std::size_t j = 0;

		for (int i : a) {
			h += (static_cast<std::size_t>(i) << (8 * j));
			j = (j + 1) % 8;
		}
		return h;
	}
};

/**
The class used to calculate product structure with milnor basis 
*/
class prod_calculator 
{
public:	
	prod_calculator() : cache(get_cache()), prod({}), no_summand_exceeding_degree(true) {}

	/**
	The function that calculates the product given sequences R and S
	@param R: R sequence
	@param S: S sequence
	@returns a set containing all sequences T in the summation of product
	*/
	std::unordered_set<std::array<unsigned int, MATRIX_SIZE>> get_prod(const std::array<int, MATRIX_SIZE>& R, const std::array<int, MATRIX_SIZE>& S);

	//std::unordered_set<std::array<unsigned int, MATRIX_SIZE>> get_prod() { return prod; }

private:
	std::unordered_set<std::array<unsigned int, MATRIX_SIZE>> prod; //want to use unordered set, need a hash function
	std::array<std::vector<std::array<uint8_t, MATRIX_SIZE - 1>>, CACHE_SIZE> cache;
	bool no_summand_exceeding_degree;

	void search_product(std::array<int, MATRIX_SIZE>& R, std::array<int, MATRIX_SIZE>& S, std::array<unsigned int, MATRIX_SIZE>& TX, int depth);

	std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> generate_matrix(const  std::array<int, MATRIX_SIZE - 1>& seq);
	
	/**
	The procedure that generates cache for computing Milnor basis
	@param cache: stores the generated cache, should be initialized to an array of empty vectors
	@param seq: auxlliary sequence for generating matrix
	@param depth: depth for search, should input 1
	*/
	void generate_cache(std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> >, CACHE_SIZE>& cache, std::array<int, MATRIX_SIZE - 1>& seq, int depth);

	std::array<std::vector<std::array<uint8_t, MATRIX_SIZE - 1>>, CACHE_SIZE> get_cache();

};

#endif // !MILNOR_BASIS_HPP