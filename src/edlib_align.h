#include "edlib/include/edlib.h"

EdlibEqualityPair additionalEqualities[5] = {
    {'a', 'A'}, {'c', 'C'}, {'g', 'G'}, {'t', 'T'}, {'n', 'N'}};

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start,
                   int *end, int k);