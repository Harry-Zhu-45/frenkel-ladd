#include "randomDouble.h"
#include <cstdlib>

double randomDouble()
{
    return rand() / (RAND_MAX + 1.0);
}
