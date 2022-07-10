#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ut {

inline void hash_combine(uint64_t& seed, uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);
}

template <typename>
struct is_sequence : std::false_type
{
};
template <typename T>
struct is_sequence<std::vector<T>> : std::true_type
{
    using value_type = T;
};
template <typename T>
struct is_sequence<std::unordered_set<T>> : std::true_type
{
    using value_type = T;
};

template <typename>
struct is_map : std::false_type
{
};
template <typename K, typename T>
struct is_map<std::map<K, T>> : std::true_type
{
    using key_type = K;
    using value_type = T;
};

template <typename It>
uint64_t hash(It begin, It end)
{
    uint64_t seed = 0;
    for (auto it = begin; it != end; ++it)
        ut::hash_combine(seed, hash(*it));
    return seed;
}

template <typename T>
uint64_t hash(const T& x)
{
    if constexpr (is_sequence<T>::value) {
        hash(x.begin(), x.end());
    }
    else if constexpr (is_map<T>::value) {
        uint64_t seed = 0;
        for (auto it = x.begin(); it != x.end(); ++it) {
            ut::hash_combine(seed, hash(it->first));
            ut::hash_combine(seed, hash(it->second));
        }
        return seed;
    }
    else
        return x.hash();
}

template <>
inline uint64_t hash<int>(const int& x)
{
    return uint64_t(x);
}

template <>
inline uint64_t hash<uint64_t>(const uint64_t& x)
{
    return uint64_t(x);
}

}  // namespace ut