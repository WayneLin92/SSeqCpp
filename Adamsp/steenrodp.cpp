
#include <iostream>
#include <vector>
#include <array>
#include <random>  
#include <ctime>
#include <algorithm>

namespace steenrodp
{
    // Constants for the project
    constexpr uint P = 3;  // A prime number
    constexpr uint N = 5;  // A fixed size

    // Divide a number a by P i times
    constexpr uint DivPtoi(uint a, uint i)
    {
        while (i--)
            a /= P;
        return a;
    }

    // Multiply a number a by P i times
    constexpr uint MulPtoi(uint a, uint i)
    {
        while (i--)
            a *= P;
        return a;
    }

    // Type aliases for arrays

    // Type R sequence R = (r_1, r_2, ..., r_N)
    using RSeq = std::array<uint, N>;
    // Type E sequence E = (e_0, e_1, ..., e_N)
    using ESeq = std::array<bool, N + 1>;
    // Type Y sequence Y = (y_0, y_1, ..., y_N)
    using YSeq = std::array<int, N + 1>;
    // Type X matrix X = (X_{ij}), i,j = 0, 1, ..., N
    using XMatrix = std::array<uint, (N + 1) * (N + 1)>;

    namespace detail
    {
        // For 0 < n < P, CacheInv[n] computes the mod p inverse of n. I.e. CacheInv[n] = m if n*m = 1 mod p.  This result will be stored in the array INV_P.
        constexpr std::array<uint, P> CacheInv()
        {
            std::array<uint, P> result = {0};
            for (size_t i = 1; i < P; ++i)
            {
                for (uint j = 1; j < P; ++j)
                {
                    if ((i * j) % P == 1)
                        result[i] = j;
                }
            }
            return result;
        }

        // for n < P, Factorials[n] computes the factorial of n mod P. This result will be stored in the array FACTORIALS.     
        constexpr std::array<uint, P> Factorials()
        {
            std::array<uint, P> result = {1};
            for (uint i = 1; i < P; ++i)
            {
                result[i] = (result[i - 1] * i) % P;
            }
            return result;
        }

        // InnerDegreeOfRSeq[i] gives the inner degree of xi_{i+1}, which is 2p^{i+1} - 2
        constexpr std::array<uint, N> InnerDegreeOfRSeq()
        {
            std::array<uint, N> result = {};
            for (size_t i = 0; i < N; ++i)
                result[i] = 2 * MulPtoi(1, i + 1) - 2;
            return result;
        }

        // InnerDegreeOfESeq[i] gives the inner degree of tau_{i}, which is 2p^{i} - 1
        constexpr std::array <uint, N + 1> InnerDegreeOfESeq()
        {
            std::array<uint, N + 1> result = {};
            for (size_t i = 0; i < N + 1; ++i)
                result[i] = 2 * MulPtoi(1, i) - 1;
            return result;
        }

        // Maximum inner degree of the Milnor basis element we will consider in this project, which is |xi_{N+1}| - 1 =2p^{N+1}-3
        constexpr uint UpperBoundOfInnerDegree()
        {
            return 2 * MulPtoi(1, N + 1) - 3;
        }
    }

    // Memoize the results from detail functions
    constexpr std::array<uint, P> INV_P = detail::CacheInv();
    constexpr std::array<uint, P> FACTORIALS = detail::Factorials();
    constexpr std::array<uint, N> RInnerDegree = detail::InnerDegreeOfRSeq();
    constexpr std::array<uint, N + 1> EInnerDegree = detail::InnerDegreeOfESeq();
    constexpr uint MaxInnerDegree = detail::UpperBoundOfInnerDegree();


    // Struct to represent a monomial in the Milnor basis
    struct MMilnor
    {
        uint c;    // A coefficient value from 0 to P-1
        ESeq E;    // An array of size N+1 for the E part
        RSeq R;    // Another array of size N for the R part
    };

    // Struct to represent selected information of a Y sequence for computational convienence
    struct YInfo
    {
        uint sigma; // sigma(Y), the number of pairs (i,j) such that i < j and Y[i] > Y[j] >= 0
        ESeq T;     // T(Y) = (t_0, t_1, ..., t_N)
        RSeq R;     // R(Y) = (r_1, r_2, ..., r_N)
    };

    // Struct to represent selected information of a X matrix for computational convienence
    struct XInfo
    {
        RSeq T;       // T(X) = (t_1, t_2, ..., t_N)
        uint b;       // b(X), a coefficient value from 0 to P-1
    };

    using VectorOfMMilnor = std::vector<MMilnor>;
    using VectorOfYInfo = std::vector<YInfo>;
    using VectorOfXInfo = std::vector<XInfo>;

    //define an ordering of the monomials in the Milnor basis
    bool MMilnorLessThan(const MMilnor &lhs, const MMilnor &rhs)
    {
        for (size_t i = 0; i < N + 1; ++i)
        {
            if (lhs.E[i] < rhs.E[i])
                return true;
            if (lhs.E[i] > rhs.E[i])
                return false;
        }
        for (size_t i = 0; i < N; ++i)
        {
            if (lhs.R[i] < rhs.R[i])
                return true;
            if (lhs.R[i] > rhs.R[i])
                return false;
        }
        return false;
    }

    // Sort the monomials in a vector of MMilnor using the ordering defined by MMilnorLessThan
    void SortMMilnor(VectorOfMMilnor &v)
    {
        std::sort(v.begin(), v.end(), MMilnorLessThan);
    }

    // Sort and combine the monomials in a vector of MMilnor using the ordering defined by MMilnorLessThan
    void CombineMMilnor(VectorOfMMilnor &v)
    {
        if (v.size() == 0)
            return;
        SortMMilnor(v);
        VectorOfMMilnor result;
        MMilnor current = v[0];
        for (size_t i = 1; i < v.size(); ++i)
        {
            if (MMilnorLessThan(current, v[i]))
            {
                if (current.c)
                    result.push_back(current);
                current = v[i];
            }
            else
            {
                current.c = (current.c + v[i].c) % P;
            }
        }
        if (current.c)
            result.push_back(current);
        v = result;
    }

    // Generate a random Milnor basis element such that its inner degree is less than or equal to max
    MMilnor RandomMMilnor(const uint &max)
    {
        MMilnor m;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint> dist(0, max);
        //make the coefficient nonzero
        m.c = dist(gen) % (P-1) + 1;
        uint CurrentInnerDegree = 0;
        for (size_t i = 0; i < N + 1; ++i)
        {
            m.E[i] = 0;
            int maxvalue = (max - CurrentInnerDegree)/EInnerDegree[i];
            if (maxvalue > 0)
            {
                m.E[i] = dist(gen) % 2;
                CurrentInnerDegree += m.E[i] * EInnerDegree[i];
            }
        }
        for (size_t i = 0; i < N; ++i)
        {
            m.R[i] = 0;
            int maxvalue = (max - CurrentInnerDegree)/RInnerDegree[i];
            if (maxvalue > 0)
            {
                m.R[i] = dist(gen) % (maxvalue + 1);
                CurrentInnerDegree += m.R[i] * RInnerDegree[i];
            }
        }
        return m;
    }

    // Compute the inner degree of a Milnor basis element
    uint InnerDegreeOfMMilnor(const MMilnor &m)
    {
        uint result = 0;
        if (m.c == 0)
            return 0;
        for (size_t i = 0; i < N + 1; ++i)
            result += m.E[i] * EInnerDegree[i];
        for (size_t i = 0; i < N; ++i)
            result += m.R[i] * RInnerDegree[i];
        return result;
    }

    // Given a matrix X without carriers, compute b(X)
    uint b(XMatrix X)
    {
        int result = 1;
        for (int k = 1; k <= N; ++k)
        {
            /* Compute the multinomial (X[0][k], X[1][k-1], ..., X[k][0]) */
            while (true)
            {
                size_t j = 0;
                while (j <= k && X[j * (N + 1) + k - j] == 0)
                    ++j; //j records the first nonzero entry along the diagonal
                if (j == k + 1)
                    break; //If all entries have become zero (all digits have been computed), then move to the next k

                //Otherwise, compute the multinomial formed by the last digits of the entries along the diagonal    
                uint prod_fact = 1;
                for (int i = j; i <= k; ++i)
                    prod_fact = (prod_fact * FACTORIALS[X[i * (N + 1) + k - i] % P]) % P;
                uint sum = 0;
                for (int i = j; i <= k; ++i)
                    sum += X[i * (N + 1) + k - i];
                result = (result * FACTORIALS[sum % P] * INV_P[prod_fact]) % P;

                //Delete the last digit of each entry. In other words, move to the next digit.
                for (size_t i = j; i <= k; ++i)      
                    X[i * (N + 1) + k - i] /= P;  
            }                   
        }
        return result;
    }

    // Check if there's a carrier in the addition of two numbers in P-ary system
    bool HasCarrier(uint a, uint b)
    {
        while (a && b)
        {
            if ((a % P) + (b % P) >= P)
                return true;
            a /= P;
            b /= P;
        }
        return false;
    }

    // Given that Xij > 0, Xij + t has no carrier, find the largest number n such that n < Xij and n + t has no carrier
    uint NextXij(uint Xij, uint t)
    {
        uint pow = 1;
        uint x = Xij;
        //find the lowest nonzero digit of Xij, we will decrease that digit by 1.
        while (x % P == 0)
        {
            x /= P;
            pow *= P;
        }
        return Xij - 1 - (t % pow); //all higher digits are unchanged, hence won't have a carrier
    }
    
    // find the largest number n such that n <= upper_bound and n has no carrier with mask
    uint32_t max_mask(uint32_t upper_bound, uint32_t mask)
    {
        if (HasCarrier(upper_bound, mask))
        {
            uint32_t a = upper_bound; 
            uint32_t b = mask;
            uint32_t d = 0;
            uint32_t j = 0;
            while (a || b)
            {
                if ((a % P) + (b % P) >= P)
                    d = j; // d records the highest position of the digit that causes the carrier
                j++; // j is the current position of the digit    
                a /= P;
                b /= P;
            }
            return upper_bound - (upper_bound % (MulPtoi(1, d+1))) + (MulPtoi(1, d+1) - 1 - (mask % (MulPtoi(1, d+1)))); //all higher digits are unchanged, hence won't have a carrier, and the lower digits are set to the maximum possible value
        }
        else
            return upper_bound;        
    }

    // Find all type X matrix such that R(X) = R1 and S(X) = R2, record the corresponding b(X), T(X) in result
    VectorOfXInfo RPartMulMilnor(const RSeq &R1, const RSeq &R2)
    {
        //result records the result of the computation.  Each element of result is a struct XInfo, which contains the b(X) and T(X) of a type X matrix X.
        VectorOfXInfo result = {};

        //R1 corresponds to R in the product formula, R2 corresponds to S.
        if (HasCarrier(R1[N - 1], R2[N - 1]))
            return result;
        //R1[N-1] = r_N = X_{N,0}, R2[N-1] = s_N = X_{0,N}.  If X_{N,0} + X_{0,N} has a carrier, then the product is forced to be zero.    

        XMatrix X, XR, XS, XT;
        //XR, XS, XT stands for the remaining unused weighted row sum, unused column sum and current diagonal sum at a certain matrix position respectively.

        //initialize X to 0
        for (size_t i = 0; i < (N + 1) * (N + 1); ++i)
            X[i] = 0;

        for (size_t i = 1; i <= N; ++i)
            XS[(N - i) * (N + 1) + i] = R2[i - 1];
        for (size_t i = 1; i <= N - 1; ++i)
            XR[i * (N + 1) + (N - i)] = R1[i - 1];
        for (size_t i = 1; i <= N - 1; ++i)
            XT[i * (N + 1)] = 0;
        X[N] = R2[N - 1];
        X[N*(N + 1)] = R1[N - 1];
        //in the odd prime version, we need to compute b(X).  Hence we not only need to have T(X), we also need to have all relavent entries of X filled in. 

        XT[(N - 1) * (N + 1) + 1] = R1[N - 1] + R2[N - 1];

        size_t i = N - 1, j = 1;
        bool decrease = false;
        while (true)
        {                      
            bool move_right = false;
            if (j)
            {
                const size_t index = i * (N + 1) + j;
                const size_t index_up = (i - 1) * (N + 1) + j;
                const size_t index_up_right = (i - 1) * (N + 1) + j + 1;
                const size_t index_left = i * (N + 1) + j - 1;
                if (i == 1)
                {
                    if (decrease)
                    {
                        if (j == N - 1)
                            X[index] = NextXij(X[index], XT[index]);
                        else
                            X[index] = NextXij(X[index], XT[index] + X[index_up_right]);
                        decrease = false;
                    }
                    else
                    {
                        if (j == N-1)
                            X[index] = max_mask(std::min(DivPtoi(XR[index], j), XS[index]), XT[index]); // X_{0,N} has already been determined and added to the diagonal sum from the very beginning
                        else
                            X[index] = max_mask(std::min(DivPtoi(XR[index], j), XS[index]), XT[index] + X[index_up_right]);
                    }
                    X[index_up] = XS[index] - X[index]; //the entry X_{i,0} is determined as the remaining column sum
                    if (HasCarrier(X[index_up], XT[index_left]))
                    {
                        if (X[index])
                            decrease = true;
                        else
                            move_right = true;
                    }
                    else
                    {
                        XR[index_left] = XR[index] - MulPtoi(X[index], j);
                        //record the diagonal sum in XT[0,j+1]
                        if (j == N - 1)
                            XT[index_up_right] = XT[index] + X[index];
                        else
                            XT[index_up_right] = XT[index] + X[index] + X[index_up_right];
                        --j;  
                    }
                }
                else
                {
                    if (decrease)
                    {
                        X[index] = NextXij(X[index], XT[index]);
                        decrease = false;
                    }
                    else
                        X[index] = max_mask(std::min(DivPtoi(XR[index], j), XS[index]), XT[index]);
                    XR[index_left] = XR[index] - MulPtoi(X[index],j) ;
                    XS[index_up] = XS[index] - X[index];
                    XT[index_up_right] = XT[index] + X[index];
                    --j;
                }
            }
            else
            {
                if (i == 1)
                {
                    if (HasCarrier(XR[N + 1], X[1]) == false)
                    {   
                        X[N + 1] = XR[N + 1];                        
                        XT[1] = X[N + 1] + X[1];
                        //form the summand of result corresponding to X
                        XInfo m;
                        m.b = b(X);
                        for (size_t i = 0; i < N; ++i)
                        {                            
                            m.T[i] = XT[i + 1];
                        }
                        /* Add to result. */
                        result.push_back(m);
                    }
                    move_right = true;
                }
                else
                {
                    const size_t index = i * (N + 1);
                    const size_t index_up_right = (i - 1) * (N + 1) + 1;
                    X[index] = XR[index];//the entry X_{i,0} is determined as the remaining row sum
                    XT[index_up_right] = XR[index];
                    j = N - (--i);
                }
            }
            if (move_right)
            {
                size_t index = i * (N + 1) + j;
                do
                {
                    if (i + j < N)
                    {
                        ++j;
                        ++index;
                    }
                    else
                    {
                        ++i;
                        index += (N + 2) - j;
                        j = 1;
                    }
                } while (i < N && X[index] == 0);
                if (i >= N)
                    break;
                decrease = true;
            }
        }
        return result;
    }

    // Find the number of pairs (i,j) in a Y sequence such that i < j and Y[i] > Y[j] >= 0
    uint PermOfYseq (const YSeq &Y)
    {
        uint result = 0;
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = i + 1; j < N + 1; ++j)
            {
                if (Y[i] > Y[j] && Y[j] >= 0)
                    ++result;
            }
        }
        return result;
    }

    // Find the number of pairs (i,j) in two E sequences E1, E2 such that i > j and E1[i] = E2[j] = 1
    uint PermOfTwoEseq (const ESeq &E1, const ESeq &E2)
    {
        uint result = 0;
        for (size_t i = 0; i < N + 1; ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                if (E1[i] && E2[j])
                    ++result;
            }
        }
        return result;
    }

    // Find all Y sequences satisfying given conditions
    VectorOfYInfo DetermineY (const ESeq &E1, const RSeq &R1, const ESeq &E2)
    {
        //result records all Y sequences satisfying the conditions. 
        VectorOfYInfo result = {};
        
        YSeq Y;             //Y = (y_0, y_1, ..., y_N)
        uint sigma = 0;     //sigma(Y), initialized to 0
        ESeq T = {0};       //T(Y) = (t_0, t_1, ..., t_N)
        RSeq R = {0};       //R(Y) = (r_1, r_2, ..., r_N)

        //initialize Y
        for (size_t i = 0; i < N + 1; ++i)
            Y[i] = -1;

        //find all i such that S(Y)[i] = E2[i] = 1, record such i in the sequence Values in decreasing order
        std::vector<int> Values;
        for (int i = N; i >= 0; --i)
        {
            if (E2[i])
                Values.push_back(i);
        }

        //If Values is empty, then this forces Y to be -1 everywhere.
        if (Values.size() == 0)
        {
            YInfo m;
            m.sigma = 0;
            m.T = T;
            m.R = R;
            result.push_back(m);
            return result;
        }

        //If Values is not empty, then each of the nonnegative numbers in Values need to be assigned to a position i in Y such that E1[i] = 0.  Note, if E1[i] = 1, then T[i] = 0, Y[i] = -1.

        //count the number of 0's in E1
        size_t count = 0;
        for (size_t i = 0; i < N + 1; ++i)
        {
            if (E1[i] == 0)
                ++count;
        }

        //If the number of values is greater than the number of positions, then there's no such Y.
        if (Values.size() > count)
            return result;

        //Position[i] records the position of the i-th value in Values in Y, Y[Position[i]] = Values[i]. Initialize all positions to N+1.
        std::vector<uint> Position(Values.size(), N + 1);
        //index records the current index of the value in Values awaiting assignment
        int index = 0;
        //whether we need to decrease the position assigned to the current value
        bool decrease = false;     
        //whether we need to move to the previous value in Values
        bool move_left = false;
        
        while(true)
        {
            if (move_left == true)
            {
                if (index == 0)
                    break;
                //if the current value has been assigned a position, unsign the current value, note R[a-1] = R_{a}
                if (Position[index] < N + 1)
                {
                    Y[Position[index]] = -1;
                    T[Position[index]] = 0;
                    if (Position[index] > Values[index])
                        R[Position[index] - Values[index] - 1] -= MulPtoi(1, Values[index]);
                    Position[index] = N + 1;
                }
                //index moves left    
                --index;
                decrease = true;
                move_left = false;
            }
            else if (decrease == true && Position[index] == 0)
            {
                //if we need to decrease the position assigned to the current value, and the position is already 0, move left
                move_left = true;
            }
            else if (decrease == true && Position[index] != 0)
            {
                //find the largest position k such that k < Position[index], k >= Values[index], E1[k] = 0, T[k] = 0, and R[k - Values[index] - 1] + MulPtoi(1, Values[index]) <= R1[k - Values[index] - 1] if k - Values[index] > 0. 
                if(Position[index] > N)
                    std::cout << "Error: Position[index] > N and decrease = true.\n";

                int k = Position[index] - 1;
                while (k >= Values[index] && (E1[k] || T[k] || (k > Values[index] && R[k - Values[index] - 1] + MulPtoi(1, Values[index]) > R1[k - Values[index] - 1])))
                    --k;

                //if such k doesn't exist, move left    
                if (k < Values[index])
                {
                    move_left = true;
                    decrease = true;
                }

                //if such k exists, unsign the previous position, and assign k to the current value     
                else
                {
                    Y[Position[index]] = -1;
                    T[Position[index]] = 0;
                    if (Position[index] > Values[index])
                        R[Position[index] - Values[index] - 1] -= MulPtoi(1, Values[index]);
                                        
                    Position[index] = k;
                    Y[k] = Values[index];
                    T[k] = 1;
                    if (k > Values[index])
                        R[k - Values[index] - 1] += MulPtoi(1, Values[index]);

                    //if this is not the last value, move to the next value
                    if (index + 1 < Values.size())
                    {
                        ++index;
                        decrease = false;
                    }

                    //if this is already the last value, record the info in result, and set decrease to true to find other possible positions for the current value.
                    else
                    {
                        YInfo m;
                        m.sigma = PermOfYseq(Y);
                        m.T = T;
                        m.R = R;
                        result.push_back(m);
                        decrease = true;

                    }
                }    
            }
            else
            {
                //find the largest position k such that k >= Values[index], E1[k] = 0, T[k] = 0, and R[k - Values[index] - 1] + MulPtoi(1, Values[index]) <= R1[k - Values[index] - 1] if k - Values[index] > 0. 
                int k = N;
                while (k >= Values[index] && (E1[k] || T[k] || (k > Values[index] && R[k - Values[index] - 1] + MulPtoi(1, Values[index]) > R1[k - Values[index] - 1])))
                    --k;
                    
                //if such k doesn't exist, move left    
                if (k < Values[index])
                {
                    move_left = true;
                    decrease = true;
                }

                //if such k exists, assign the position to the current value     
                else
                {
                    Position[index] = k;
                    Y[k] = Values[index];
                    T[k] = 1;
                    if (k > Values[index])
                        R[k - Values[index] - 1] += MulPtoi(1, Values[index]);

                    //if this is not the last value, move to the next value
                    if (index + 1 < Values.size())
                    {
                        ++index;
                        decrease = false;
                    }

                    //if this is already the last value, record the info in result, and set decrease to true to find other possible positions for the current value.
                    else
                    {
                        YInfo m;
                        m.sigma = PermOfYseq(Y);
                        m.T = T;
                        m.R = R;
                        result.push_back(m);
                        decrease = true;
                    }
                }    
            }
        }
        return result;
    }

    // Find the product of two Milnor basis elements, sort the terms in the product, and combine the terms with the same E and R parts
    VectorOfMMilnor MulMilnor (const MMilnor &M1, const MMilnor &M2)
    {
        if (InnerDegreeOfMMilnor(M1) + InnerDegreeOfMMilnor(M2) > MaxInnerDegree)
        {
            std::cout << "Error: The inner degree of the product exceeds the maximum inner degree.\n";
            return {};
        }
        
        ESeq E1 = M1.E;
        RSeq R1 = M1.R;
        ESeq E2 = M2.E;
        RSeq R2 = M2.R;
        uint CoefficientP = M1.c * M2.c % P;
        VectorOfMMilnor product = {};

        // sum1 = sum E1, sum2 = sum E2. In the product formula, the first coefficient (from E part) is -1 to the power of sum1 * sum2
        uint sum1 = 0, sum2 = 0;
        for (size_t i = 0; i < N + 1; ++i)
        {
            sum1 += E1[i];
            sum2 += E2[i];
        }
        uint PowerE = sum1 * sum2;

        VectorOfYInfo ListOfY = DetermineY(E1, R1, E2);

        //read all Y info in ListOfY
        for (const YInfo &CurrentYInfo : ListOfY)
        {
            // The second coefficient in the product formula (from Y part) is -1 to the power of PowerE
            uint PowerY = PermOfTwoEseq(E1, CurrentYInfo.T) + CurrentYInfo.sigma;
            RSeq R1MinusRY;
            for (size_t i = 0; i < N; ++i)
                R1MinusRY[i] = R1[i] - CurrentYInfo.R[i];
            VectorOfXInfo ListOfX = RPartMulMilnor(R1MinusRY, R2);

            //read all X info in ListOfX
            for (const XInfo &CurrentXInfo : ListOfX)
            {
                MMilnor m;
                //if there is a negative sign in this summand, the coefficient is  - CoefficientP * b(X) mod P
                if ((PowerE + PowerY) % 2)
                    m.c = (P - CoefficientP * CurrentXInfo.b % P) % P;
                else
                    m.c = CoefficientP * CurrentXInfo.b % P;

                if (m.c == 0)
                    continue;    

                ESeq E1plusTY;
                for (size_t i = 0; i < N + 1; ++i)
                {
                    E1plusTY[i] = E1[i] + CurrentYInfo.T[i];
                    if (E1[i] && CurrentYInfo.T[i])
                        std::cout << "Error: E1[i] and T[i] are both 1.\n";
                }
                m.E = E1plusTY;
                m.R = CurrentXInfo.T;
                // Add this summand to the product
                product.push_back(m);
            }
        }
        CombineMMilnor(product);
        return product;
    }

    // Find the product of two vectors of Milnor basis elements, sort the terms in the product, and combine the terms with the same E and R parts
    VectorOfMMilnor MulMilnorVectors(const VectorOfMMilnor &v1, const VectorOfMMilnor &v2)
    {
        VectorOfMMilnor result;
        for (const auto &m1 : v1)
        {
            for (const auto &m2 : v2)
            {
                VectorOfMMilnor product = MulMilnor(m1, m2);
                for (const auto &m : product)
                    result.push_back(m);
            }
        }
        CombineMMilnor(result);
        return result;
    }

    // Determine whether two vectors of MMilnor (already sorted and combined) are equal
    bool VectorOfMMilnorEqual (const VectorOfMMilnor &v1, const VectorOfMMilnor &v2)
    {
        if (v1.size() != v2.size())
            return false;
        for (size_t i = 0; i < v1.size(); ++i)
        {
            if (v1[i].c != v2[i].c)
                return false;
            if (MMilnorLessThan(v1[i], v2[i]) || MMilnorLessThan(v2[i], v1[i]))
                return false;
        }
        return true;
    }

}


// Main function
int main_test(int, char**, int&, const char*)
{
    using namespace steenrodp;

 /*
    // Test the MulMilnor function for chosen milnor basis elements M1, M2, and cout M1 * M2.
    MMilnor M1 = {1, {0, 0, 1, 0, 0, 0}, {19, 6, 3, 0, 0}};
    MMilnor M2 = {1, {0, 0, 1, 0, 0, 0}, {3, 0, 0, 0, 0}};
    VectorOfMMilnor M1M2 = MulMilnor(M1, M2);
    std::cout << "P = " << P << ", N = " << N << "\n";
    std::cout << "The number of MMilnor elements in M1 * M2 is " << M1M2.size() << std::endl;
    for (const auto &m : M1M2)
    {
        std::cout << "c = " << m.c << ", E = (";
        for (size_t i = 0; i < N + 1; ++i)
            std::cout << m.E[i] << (i < N ? ", " : "), R = (");
        for (size_t i = 0; i < N; ++i)
            std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
    }
 */

 /*
    // Test the MulMilnor function, generate random milnor basis elements M1, M2 and cout M1 * M2.  We require M1 * M2 to be nonzero, the inner degree of M1 * M2 to be less than or equal to MaxInnerDegree, and the number of summands in M1 * M2 to be less than 10 so that the output is easy to read.
    MMilnor M1, M2;
    VectorOfMMilnor M1M2;
    uint degree1, degree2;
    uint count, MaxCount = 100;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint> dist(0, MaxInnerDegree);

    for (count = 0; count < MaxCount; ++count)
    {
        degree1 = dist(gen) % MaxInnerDegree;
        degree2 = dist(gen) % (MaxInnerDegree - degree1 + 1);
        M1 = RandomMMilnor(degree1);
        M2 = RandomMMilnor(degree2);
        M1M2 = MulMilnor(M1, M2);
        if (M1M2.size() != 0 && M1M2.size() < 10)
        {
            std::cout << "The product of the following two MMilnor elements is nonzero.\n";
            std::cout << "P = " << P << ", N = " << N << "\n";
            std::cout << "M1: c = " << M1.c << ", E = (";
            for (size_t i = 0; i < N + 1; ++i)
                std::cout << M1.E[i] << (i < N ? ", " : "), R = (");
            for (size_t i = 0; i < N; ++i)
                std::cout << M1.R[i] << (i < N - 1 ? ", " : ")\n");
            std::cout << "M2: c = " << M2.c << ", E = (";
            for (size_t i = 0; i < N + 1; ++i)
                std::cout << M2.E[i] << (i < N ? ", " : "), R = (");
            for (size_t i = 0; i < N; ++i)
                std::cout << M2.R[i] << (i < N - 1 ? ", " : ")\n");
            std::cout << "The number of MMilnor elements in M1 * M2 is " << M1M2.size() << std::endl;
            for (const auto &m : M1M2)
            {
                std::cout << "c = " << m.c << ", E = (";
                for (size_t i = 0; i < N + 1; ++i)
                    std::cout << m.E[i] << (i < N ? ", " : "), R = (");
                for (size_t i = 0; i < N; ++i)
                    std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
            }
            break;
        }
    }
    if (count == MaxCount)
        std::cout << "For all random MMilnor elements generated in the last " << MaxCount << " times, the product is either zero or has more than 10 summands. Please try again.\n";
 */

// /*
    // Test whether the MulMilnor function satisfies associativity.  Generate three random Milnor basis elements M1, M2, M3, and compute the product of M1*M2*M3 in two different ways. Compare the two results.
    MMilnor M1, M2, M3;
    VectorOfMMilnor M1M2, M2M3, result1, result2;
    uint degree1, degree2, degree3;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint> dist(0, MaxInnerDegree);

    //Generate random MMilnor elements. Generate at most MaxCount times.
    size_t count = 0;
    size_t MaxCount = 100;
    size_t nonzerotimes = 0;
    while (VectorOfMMilnorEqual (result1, result2) && count < MaxCount)
    {
        degree1 = dist(gen);
        degree2 = dist(gen) % (MaxInnerDegree - degree1 + 1);
        degree3 = dist(gen) % (MaxInnerDegree - degree1 - degree2 + 1);
        M1 = RandomMMilnor(degree1);
        M2 = RandomMMilnor(degree2);
        M3 = RandomMMilnor(degree3);
        M1M2 = MulMilnor(M1, M2);
        M2M3 = MulMilnor(M2, M3);
        result1 = MulMilnorVectors(M1M2, {M3});
        result2 = MulMilnorVectors({M1}, M2M3);
        ++count;
        if (result1.size() != 0)
            ++nonzerotimes;
    }
    if (count == MaxCount)
    {
        std::cout << "For all random MMilnor elements generated in the last " << MaxCount << " times, the two results are equal.\n";
        std::cout << "Amongst these examples, the number of times the product is nonzero is " << nonzerotimes << ".\n";
        std::cout << "The last example is:\n";
    }
    else
        std::cout << "The two results are not equal for the following example.\n";

    // Display the three random MMilnor elements
    std::cout << "P = " << P << ", N = " << N << "\n";
    std::cout << "M1: c = " << M1.c << ", E = (";
    for (size_t i = 0; i < N + 1; ++i)
        std::cout << M1.E[i] << (i < N ? ", " : "), R = (");
    for (size_t i = 0; i < N; ++i)
        std::cout << M1.R[i] << (i < N - 1 ? ", " : ")\n");
    std::cout << "M2: c = " << M2.c << ", E = (";
    for (size_t i = 0; i < N + 1; ++i)
        std::cout << M2.E[i] << (i < N ? ", " : "), R = (");
    for (size_t i = 0; i < N; ++i)
        std::cout << M2.R[i] << (i < N - 1 ? ", " : ")\n");
    std::cout << "M3: c = " << M3.c << ", E = (";
    for (size_t i = 0; i < N + 1; ++i)
        std::cout << M3.E[i] << (i < N ? ", " : "), R = (");
    for (size_t i = 0; i < N; ++i)
        std::cout << M3.R[i] << (i < N - 1 ? ", " : ")\n");

    std::cout << "The number of MMilnor elements in M1 * M2 is " << M1M2.size() << std::endl;
    for (const auto &m : M1M2)
    {
        std::cout << "c = " << m.c << ", E = (";
        for (size_t i = 0; i < N + 1; ++i)
            std::cout << m.E[i] << (i < N ? ", " : "), R = (");
        for (size_t i = 0; i < N; ++i)
            std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
    }

    std::cout << "The number of MMilnor elements in result1 is " << result1.size() << std::endl;
    for (const auto &m : result1)
    {
        std::cout << "c = " << m.c << ", E = (";
        for (size_t i = 0; i < N + 1; ++i)
            std::cout << m.E[i] << (i < N ? ", " : "), R = (");
        for (size_t i = 0; i < N; ++i)
            std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
    }

    std::cout << "The number of MMilnor elements in M2 * M3 is " << MulMilnor(M2, M3).size() << std::endl;
    for (const auto &m : MulMilnor(M2, M3))
    {
        std::cout << "c = " << m.c << ", E = (";
        for (size_t i = 0; i < N + 1; ++i)
            std::cout << m.E[i] << (i < N ? ", " : "), R = (");
        for (size_t i = 0; i < N; ++i)
            std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
    }

    std::cout << "The number of MMilnor elements in result2 is " << result2.size() << std::endl;
    for (const auto &m : result2)
    {
        std::cout << "c = " << m.c << ", E = (";
        for (size_t i = 0; i < N + 1; ++i)
            std::cout << m.E[i] << (i < N ? ", " : "), R = (");
        for (size_t i = 0; i < N; ++i)
            std::cout << m.R[i] << (i < N - 1 ? ", " : ")\n");
    }
// */

    return 0;
}