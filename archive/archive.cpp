

/* Find the indices of leading lm from [begin, end) such that when ignoring the exponents
 * lm[-j-1] divides mon[i] and lm[:-j-1] divides mon[:i].
 */
int IndexOfDivisibleLeading_(const Mon& mon, size_t i, size_t j, IterIndex begin, IterIndex end) const
{
    IterIndex begin1 = std::lower_bound(begin, end, mon[i].gen, [j](const EleIndex& p, int gen) { return p.first[j] < gen; });
    IterIndex end1 = std::lower_bound(begin1, end, mon[i].gen + 1, [j](const EleIndex& p, int gen) { return p.first[j] < gen; });
    if (begin1 != end1) {
        if (begin1->first.size() == j + 1) {
            for (int i : begin1->second) {
                if (divisible(data[i].GetLead(), mon))
                    return i;
            }
            ++begin1;
        }
        for (size_t l = 0; l < i; ++l) {
            int ind = IndexOfDivisibleLeading_(mon, l, j + 1, begin1, end1);
            if (ind != -1)
                return ind;
        }
    }
    return -1;
}