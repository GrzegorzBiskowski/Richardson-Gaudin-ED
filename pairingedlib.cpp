#include "pairingedlib.hpp"

#include <cmath>
#include <iostream>
#include <vector>

std::size_t pairingModelParameters::totalElectronCount() const
{
	if (oddParity)
	{
		return 2 * pairCount + 1;
	}
	else
	{
		return 2 * pairCount;
	}
}

std::size_t pairingModelParameters::totalPairCapacity() const
{
	std::size_t result = 0;
	for (std::size_t i = 0; i < levels.size(); i++)
	{
		result += levels[i].pairCapacity;
	}
	return result;
}