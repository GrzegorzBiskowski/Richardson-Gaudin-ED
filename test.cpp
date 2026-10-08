#include "pairingedlib.hpp"
#include "tests.hpp"

#include <cmath>
#include <iostream>
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
	std::cout << "myBasis.size() = " << myBasis.size();
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