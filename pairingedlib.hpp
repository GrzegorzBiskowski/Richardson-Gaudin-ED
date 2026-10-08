#pragma once

#include<vector>

struct level
{
	double energy;
	std::size_t pairCapacity;
};

struct pairingModelParameters
{
	std::vector<level> levels;
	double pairingStrength;
	std::size_t pairCount;
	bool oddParity;
	std::size_t totalElectronCount() const;
	std::size_t totalPairCapacity() const;
};