#include <iostream>
#include <stdexcept>

#include "pairingedlib.hpp"
#include "tests.hpp"

int main()
{
    levels_test();
    occupation_test();
    basis_test();
    generateBoundedOccupations_test();
    countBoundedOccupations_test();
    return 0;
}
