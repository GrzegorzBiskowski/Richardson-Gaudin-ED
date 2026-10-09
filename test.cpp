#include "pairingedlib.hpp"
#include "tests.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void levels_test()
{
	Level something1 = {2.0, 1 };
	Level something2 = {2.0, 2 };
	Level something3 = {3.0, 1 };
	Level something4 = {1.0, 3 };
	std::vector<Level> testVec = { something1, something2, something3, something4 };
	PairingModelParameters systemPar = { testVec, 0.1, 4, false };
	std::cout << "systemPar.levels[i].energy = ";
	for (std::size_t i = 0; i < systemPar.levels.size(); i++)
	{
		std::cout << systemPar.levels[i].energy << ' ';
	}
	std::cout << '\n' << "systemPar.levels[i].pairCapacity = ";
	for (std::size_t i = 0; i < systemPar.levels.size(); i++)
	{
		std::cout << systemPar.levels[i].pairCapacity << ' ';
	}
	std::cout << '\n';
	std::cout << "totalElectronCount = " << systemPar.totalElectronCount() << '\n';
	std::cout << "totalPairCapacity = " << systemPar.totalPairCapacity() << '\n' << '\n';
}

void occupation_test()
{
	std::vector<std::uint32_t> occupations = { 1,2,3,4,0 };
	OccupationState myState(occupations);
	std::cout << "myState.size() = " << myState.size() << '\n';
	std::cout << "myState.totalOccupation() = " << myState.totalOccupation() << '\n';
	std::cout << "myState.occupations()[i] = ";
	for (std::size_t i = 0; i < myState.size(); i++)
	{
		std::cout << myState.occupations()[i] << ' ';
	}
	std::cout << '\n';
	std::cout << "myState[2] = " << myState[2] << '\n' << '\n';
}

void basis_test()
{
	std::vector<OccupationState> myStates = { OccupationState({ 1,1,0,0 }),OccupationState({ 1,0,1,0 }),
	OccupationState({ 1,0,0,1 }), OccupationState({ 0,1,1,0 }), OccupationState({ 0,1,0,1 }) , OccupationState({0,0,1,1})};
	Basis myBasis;
	for (std::size_t i = 0; i < myStates.size(); i++)
	{
		myBasis.addState(myStates[i]);
	}
	std::cout << "myBasis.size() = " << myBasis.size() << '\n';
	std::cout << "myBasis[1].totalOccupation() = " << myBasis[1].totalOccupation() << '\n';
	std::cout << "myBasis[2].size() = " << myBasis[2].size() << '\n';
	std::cout << "myBasis[i][j]:" << '\n';
	for (std::size_t i = 0; i < myBasis.size(); i++)
	{
		for (std::size_t j = 0; j < myBasis[i].size(); j++)
		{
			std::cout << myBasis[i][j] << ' ';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void generateBoundedOccupations_test()
{
	std::vector<OccupationState> myStates = generateBoundedOccupations({ 1,1,1,1 }, 2);
	Basis myBasis;
	std::cout << "Basis generated from generateBoundedOccupations({ 1,1,1,1 }, 2):" << '\n';
	for (std::size_t i = 0; i < myStates.size(); i++)
	{
		myBasis.addState(myStates[i]);
	}
	for (std::size_t i = 0; i < myBasis.size(); i++)
	{
		for (std::size_t j = 0; j < myBasis[i].size(); j++)
		{
			std::cout << myBasis[i][j] << ' ';
		}
		std::cout << '\n';
	}
	std::vector<OccupationState> myStates1 = generateBoundedOccupations({ 1,2,2,1 }, 3);
	Basis myBasis1;
	std::cout << "Basis generated from generateBoundedOccupations({ 1,2,2,1 }, 3):" << '\n';
	for (std::size_t i = 0; i < myStates1.size(); i++)
	{
		myBasis1.addState(myStates1[i]);
	}
	for (std::size_t i = 0; i < myBasis1.size(); i++)
	{
		for (std::size_t j = 0; j < myBasis1[i].size(); j++)
		{
			std::cout << myBasis1[i][j] << ' ';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void countBoundedOccupations_test()
{
	std::vector<std::uint32_t> upperBounds1 = { 1,1,1,1,1,1,1,1 };
	std::vector<std::uint32_t> upperBounds2 = { 1,2,2,1 };
	std::vector<std::uint32_t> upperBounds3 = { 1,2,3,4,5 };
	std::cout << "countBound... dla {1,1,1,1,1, 1,1,1} M=4: " << countBoundedOccupations(upperBounds1, 4) << '\n';
	std::cout << "countBound... dla {1,2,2,1} M=3: " << countBoundedOccupations(upperBounds2, 3) << '\n';
	std::cout << "countBound... dla {1,2,3,4,5} M=6: " << countBoundedOccupations(upperBounds3, 6) << '\n';
	std::cout << '\n';
}

void blockedConfiguration_test()
{
	BlockedConfiguration blocked({ 0, 1, 0, 1 });
	std::cout << "BlockedConfiguration blocked({0,1,0,1})" << '\n';

	std::cout << "Seniority: "
		<< blocked.seniority() << '\n';

	for (std::size_t j = 0;
		j < blocked.blockedElectrons().size();
		++j)
	{
		std::cout << "Level " << j
			<< ": "
			<< blocked.blockedElectrons()[j]
			<< " blocked electrons\n";
	}
	std::cout << '\n';
}

void calculateSectorDegeneracy_test()
{
	std::vector<Level> levels = {
	{ 1.0, 2 }
	};

	BlockedConfiguration blocked({ 1 });

	const auto d =
		calculateSectorDegeneracy(levels, blocked);
	std::cout << "calculateSectorDegeneracy({1.0,2}, {1}) = " << d << '\n';

	std::vector<Level> levels2 = {
	{ 1.0, 1 },
	{ 2.0, 2 },
	{ 3.0, 1 }
	};

	BlockedConfiguration blocked2({ 1, 1, 0 });
	const auto b = calculateSectorDegeneracy(levels2, blocked2);
	std::cout << "Second example = " << b << '\n' << '\n';
}

void modelGenerator_test()
{
	PairingModelParameters model{
	{
		{ 1.0, 1 },
		{ 2.0, 2 },
	},
	0.1,  // pairingStrength
	2,    // pairCount
	false // even parity
	};

	SectorGenerator generator(model);
	const auto sectors = generator.generate();
	std::cout << "Number of sectors: "
		<< sectors.size() << '\n';

	for (const auto& sector : sectors) {
		std::cout
			<< "ID: " << sector.id
			<< ", seniority: " << sector.seniority()
			<< ", pairs: " << sector.pairCount
			<< ", degeneracy: " << sector.degeneracy
			<< '\n';
	}
	std::cout << '\n';
}

void pairingHamiltonian_test()
{
	PairingModelParameters model{
		{
			{ 1.0, 1 },
			{ 2.0, 1 }
		},
		0.5,
		1,
		false
	};

	SenioritySector sector{
		0,
		BlockedConfiguration({ 0, 0 }),
		1,
		{ 1, 1 },
		0.0,
		1
	};

	Basis basis;

	basis.addState(OccupationState({ 1, 0 }));
	basis.addState(OccupationState({ 0, 1 }));

	PairingHamiltonian hamiltonian(
		model,
		sector,
		basis
	);

	DenseMatrix H = hamiltonian.buildMatrix();

	assert(H.rows() == 2);
	assert(H.cols() == 2);

	assert(std::abs(H(0, 0) - 1.5) < 1e-12);
	assert(std::abs(H(1, 1) - 3.5) < 1e-12);

	assert(std::abs(H(0, 1) + 0.5) < 1e-12);
	assert(std::abs(H(1, 0) + 0.5) < 1e-12);
}

void memoryEstimate_test()
{
	PairingModelParameters model;

	model.levels = {
		{ 1.0, 1 },
		{ 2.0, 1 }
	};

	model.pairingStrength = 0.5;
	model.pairCount = 1;
	model.oddParity = false;

	SectorGenerator generator(model);

	const auto sectors = generator.generate();

	const auto estimate =
		estimateMemory(model, sectors);

	assert(estimate.sectorCount == 2);

	assert(estimate.totalEigenvalues == 3);

	assert(estimate.largestBlockDimension == 2);

	assert(estimate.matrixBytes ==
		4 * sizeof(double));

	assert(estimate.spectrumBytes ==
		3 * sizeof(SpectrumEntry));

	assert(estimate.peakBytes >=
		estimate.matrixBytes +
		estimate.spectrumBytes);

	assert(estimate.compatibleWithLapackInt);
}
