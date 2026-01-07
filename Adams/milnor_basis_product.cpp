#include "milnor_basis_product.hpp"

std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> prod_calculator::generate_matrix(const  std::array<int, MATRIX_SIZE-1>& seq)
{
	std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> matrix = { 0 };
	for (int i = 0; i < MATRIX_SIZE - 1; ++i)
	{
		if (seq[i] > 0)
		{
			matrix[i + 1 - seq[i]][seq[i] - 1] = 1;
		}
	}
	return matrix;
}

void prod_calculator::generate_cache(std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> >, CACHE_SIZE>& cache, std::array<int, MATRIX_SIZE-1>& seq, int depth = 1)
{
	if (depth > MATRIX_SIZE - 1)
	{
		int index = 0;
		int degree = 0;
		for (int i = 0; i < MATRIX_SIZE - 1; ++i)
		{
			if (seq[i] > 0)
			{
				index += pow(2, i);
				degree += pow(2, i + 1)-1;
			}
		}
		if (degree <= MAX_DEGREE)
			assert(index < CACHE_SIZE);
			cache[index].push_back(generate_matrix(seq));
		return;
	}
	for (int i = 0; i <= depth; ++i)
	{
		seq[depth-1] = i;
		generate_cache(cache, seq, depth + 1);
	}
	return;
}

std::array<std::vector<std::array<uint8_t, MATRIX_SIZE - 1>>,CACHE_SIZE> prod_calculator::get_cache()
{
	std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE - 1>, MATRIX_SIZE - 1> >, CACHE_SIZE> cache = { {} };
	std::array<int, MATRIX_SIZE - 1> temp_seq;
	generate_cache(cache, temp_seq, 1);
	std::array<std::vector<std::array<uint8_t, MATRIX_SIZE-1>>, CACHE_SIZE> result = { {} };

	for (int i = 0; i < CACHE_SIZE; ++i)
	{
		if (cache[i].size() > 0)
		{
			for (auto matrix : cache[i])
			{
				std::array<uint8_t, MATRIX_SIZE-1> new_matrix = { 0 };
				for (int j = 0; j < MATRIX_SIZE-1; ++j)
					for (int k = 0; k < MATRIX_SIZE-1; ++k)
						new_matrix[j] = new_matrix[j] | (matrix[j][k] << k); //the jk spot is stored at the k-th bit (2^k) of j-th int8 (the 8-th bit is always 0)
				result[i].push_back(new_matrix);
			}
		}
	}
	return result;
}

void prod_calculator::search_product(std::array<int, MATRIX_SIZE>& R, std::array<int, MATRIX_SIZE>& S, std::array<unsigned int, MATRIX_SIZE>& TX, int depth = 0)
{
	int r_sum = 0, s_sum = 0;
	for (int i : R)
		r_sum += i;
	for (int i : S)
		s_sum += i;
	if (r_sum == 0 && s_sum == 0)
	{//Find a term in summation
		if (prod.count(TX))
			prod.erase(TX);
		else
			prod.insert(TX); 
		return;
	}

	unsigned int r_index = 0; //r_index is colunm 0
	for (int i = 0; i < MATRIX_SIZE; ++i)
		r_index += (R[i] % 2) << i;
	
	std::array<int, MATRIX_SIZE> R_copy = R;
	std::array<int, MATRIX_SIZE> S_copy = S;
	std::array<unsigned int, MATRIX_SIZE> TX_copy = TX;

	for (size_t i = 0; i < CACHE_SIZE; ++i)
		if ((i | (r_index >> 1)) == (i + (r_index >> 1))) //the last 7 bits of the index of r and index of matrix has no repeated bits
			for (auto matrix : cache[i])
			{				
				uint8_t row_0 = 0;
				int degree_count = 0;

				for (int j = 0; j < MATRIX_SIZE - 1; ++j)
				{
					R[j] -= (matrix[j] << 1);
					if (R[j] < 0) //maybe use try assert here, but it is not more efficeint, and we may need to turn off assertions
						goto next_matrix; //if we sort the matrices within an index in some way (increasningly by row), we may skip the whole index 
					R[j] /= 2; //the least significant bit of R is recored in r_index and discarded here
				}

				for (int k = 0; k < MATRIX_SIZE; ++k)
				{
					for (int j = 0; j < MATRIX_SIZE - 1; ++j)
					{
						S[k] -= (matrix[j] >> k) % 2;
					}
					if (S[k] < 0) 
						goto next_matrix;
					else if (S[k] % 2 != 0)
					{
						if (((r_index | (i << 1)) >> k) % 2 )
							goto next_matrix;
						row_0 += (1 << k);
						S[k] -= 1;
					}
					S[k] /= 2;
				}

				TX[0] += (r_index % 2 + row_0 % 2) << depth;
				for (int j = 1; j < MATRIX_SIZE; ++j)
				{
					TX[j] += ((i >> (j - 1)) % 2 + (r_index >> j)% 2 + (row_0 >> j) % 2) << depth; //i is the matrix contribution
				}
				
				for (int j = 0; j < MATRIX_SIZE; ++j)
					degree_count += TX[j] << j;
				
				if (degree_count > MAX_DEGREE)
				{
					no_summand_exceeding_degree = false;
					goto next_matrix;
				}

				search_product(R, S, TX, depth + 1);

				next_matrix:;
				R = R_copy;
				S = S_copy;
				TX = TX_copy;
			}

}

std::unordered_set<std::array<unsigned int, MATRIX_SIZE>> prod_calculator::get_prod(const std::array<int, MATRIX_SIZE>& R, const std::array<int, MATRIX_SIZE>& S)
{
	prod = {};
	no_summand_exceeding_degree = true;
	std::array<unsigned int, MATRIX_SIZE> TX = { 0 };
	std::array<int, MATRIX_SIZE> R_copy = R;
	std::array<int, MATRIX_SIZE> S_copy = S;
	search_product(R_copy, S_copy, TX, 0);
	
	assert(no_summand_exceeding_degree); //The assertion fails if there is a possible summand whose dgree is too high.

	return prod;
}