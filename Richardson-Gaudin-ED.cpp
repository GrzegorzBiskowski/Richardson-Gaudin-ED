#include <fstream>
#include <iostream>
#include <memory>
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
    blockedConfiguration_test();
    calculateSectorDegeneracy_test();
    modelGenerator_test();
    return 0;
}
