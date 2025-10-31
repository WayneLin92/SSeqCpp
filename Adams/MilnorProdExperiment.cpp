#include<vector>
#include<array>
#include<iostream>
#include<fstream>

const int CACHE_SIZE = 193;
const int MATRIX_SIZE = 7;
const int MAX_DEGREE = 276;

std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> generate_matrix(const  std::array<int, MATRIX_SIZE>& seq)
{
	std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> matrix = { 0 };
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		if (seq[i] > 0)
		{
			matrix[i + 1 - seq[i]][seq[i] - 1] = 1;
		}
	}
	return matrix;
}

/**
The procedure that generates cache for computing Milnor basis
@param cache : stores the generated cache, should be initialized to an array of empty vectors
@param seq: auxlliary sequence for generating matrix
@param depth: depth for search, should input 1
*/
void generate_cache(std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> >, CACHE_SIZE>& cache, std::array<int, MATRIX_SIZE>& seq, int depth = 1)
{
	if (depth > MATRIX_SIZE)
	{
		int index = 0;
		int degree = 0;
		for (int i = 0; i < MATRIX_SIZE; ++i)
		{
			if (seq[i] > 0)
			{
				index += pow(2, i);
				degree += pow(2, i + 1)-1;
			}
		}
		if(degree<=MAX_DEGREE)
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

void print_matrix(const std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE>& matrix)
{
	for (int j = 0; j < MATRIX_SIZE; ++j)
	{
	
		for (int k = 0; k < MATRIX_SIZE; ++k)
			std::cout <<matrix[j][k] << ' ';
		std::cout << std::endl;
	}
}


// /**
//  * Column0[1] x_{11} x_{12} x_{13} x_{14} x_{15} x_{16} x_{17}
//  * Column0[2] x_{21} x_{22} x_{23} x_{24} x_{25} x_{26}
//  * Column0[3] x_{31} x_{32} x_{33} x_{34} x_{35}
//  * Column0[4] x_{41} x_{42} x_{43} x_{44}
//  * Column0[5] x_{51} x_{52} x_{53}
//  * Column0[6] x_{61} x_{62}
//  * Column0[7] x_{71}
//  * Column0[8]
//  */
// void MulMilnorV4(const std::array<uint32_t, 8>& R, const std::array<uint32_t, 8>& S, std::vector<const std::array<uint32_t, 8>>& result_app)
// {
//     std::array<std::array<uint32_t, 9>, 9> jDiag; /* Each X represents the bits of a row in columns 1-7 */
//     std::array<std::array<uint8_t, 9>, 9> Column0;
//     const auto deg_max = 276;
//     for (uint8_t iBit = 0; iBit <= 8; ++iBit) { /* We consider the iBit-th bit of Milnor's matrix X */
//         for (size_t i = 1; i <= 8; ++i)
//             Column0[iBit][i] = R[i - 1] % 2;
//         for (jDiag[iBit][3] = 0; jDiag[iBit][3] <= 2; ++jDiag[iBit][3]) {
//             for (jDiag[iBit][4] = 0; jDiag[iBit][4] <= 3; ++jDiag[iBit][4]) {
//                 for (jDiag[iBit][5] = 0; jDiag[iBit][5] <= 4; ++jDiag[iBit][5]) {
//                     for (jDiag[iBit][6] = 0; jDiag[iBit][6] <= 5; ++jDiag[iBit][6]) {
//                         for (jDiag[iBit][7] = 0; jDiag[iBit][7] <= 6; ++jDiag[iBit][7]) {
//                             for (jDiag[iBit][8] = 0; jDiag[iBit][8] <= 7; ++jDiag[iBit][8]) {

//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

std::array<uint8_t,MATRIX_SIZE> convert_to_packed(std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> milnor_matrix) {

    std::array<uint8_t,MATRIX_SIZE> result = {0};

    for (int i = 0; i < MATRIX_SIZE; ++i) {
        
        for (int j = 0; j < MATRIX_SIZE; ++j) {

            uint8_t val = (uint8_t) milnor_matrix[i][j];
            result[i] = result[i] | (val << j);

            // if (milnor_matrix[i][j]) {
            //     std::cout << "(" << i << "," << j << ") true" << std::endl;
            // }

            // std::cout << (result[i] == 0) << std::endl;
            
        }
    }


    return result;
}

void write_cache(std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> >, CACHE_SIZE> cache, std::string filename) {
    std::ofstream outfile(filename, std::ios::binary);

    if (!outfile) {
        std::cout << "Error opening file" << std::endl;
    }

    for (int i = 0; i < CACHE_SIZE; ++i) {

        auto chunk = cache[i];
        int size = chunk.size();
        
        outfile.write(reinterpret_cast<const char*>(&i), sizeof(i));
        outfile.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (int j = 0; j < chunk.size(); ++j) {
        
            std::array<uint8_t, MATRIX_SIZE> triangle = 
                convert_to_packed(chunk[j]);
            outfile.write(reinterpret_cast<const char*>(triangle.data()), 
                    triangle.size());
        }
    }

    outfile.close();
}

int main()
{

	std::array<std::vector<std::array<std::array<bool, MATRIX_SIZE>, MATRIX_SIZE> >, CACHE_SIZE> cache = { {} };
	std::array<int, MATRIX_SIZE> temp_seq;
	generate_cache(cache, temp_seq, 1);
	int sum = 0;
	for (int i = 0; i < CACHE_SIZE; ++i)
	{
		
		if (cache[i].size() > 0)
		{
			std::cout << "index" << i << ": " << cache[i].size() << " matrices" << std::endl;

			sum += cache[i].size();
			//for (std::array<std::array<int, MATRIX_SIZE>, MATRIX_SIZE> matrix : cache[i])
			//	print_matrix(matrix);
			std::cout << std::endl;
		}
	}
	std::cout << sum << std::endl << std::endl;

    std::array<uint8_t, MATRIX_SIZE> sample = convert_to_packed(cache[18][2]);

    for (int i = 0; i < MATRIX_SIZE; ++i) {
        std::cout << std::to_string(static_cast<int>(sample[i])) << std::endl;
        // std::cout << (sample[i] == 0) << std::endl;
    }

    write_cache(cache, "cachedata/milnor.cache");
}
