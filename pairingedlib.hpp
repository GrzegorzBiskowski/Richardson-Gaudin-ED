#pragma once

#include <vector>

struct Level
{
	double energy;
	std::size_t pairCapacity;
	void validate() const;
};

struct PairingModelParameters
{
	std::vector<Level> levels;
	double pairingStrength;
	std::size_t pairCount;
	bool oddParity;
	std::size_t totalElectronCount() const;
	std::size_t totalPairCapacity() const;
	void validate() const;
};

class OccupationState {
public:
	explicit OccupationState(std::vector<std::uint32_t> occupations)
		: occupations_(std::move(occupations)){}
	const std::vector<std::uint32_t>& occupations() const noexcept;

	std::size_t size() const noexcept;

	std::uint64_t totalOccupation() const noexcept;

	std::uint32_t operator[](std::size_t i) const;

private:
	std::vector<std::uint32_t> occupations_;
};

class Basis {
public:
	void addState(OccupationState state);

	const OccupationState& operator[](std::size_t i) const;

	std::size_t size() const noexcept;

	const std::vector<OccupationState>& states() const noexcept;

private:
	std::vector<OccupationState> states_;
};

std::vector<OccupationState> generateBoundedOccupations(const std::vector<std::uint32_t>& upperBounds,std::size_t totalOccupation);

std::size_t countBoundedOccupations(const std::vector<std::uint32_t>& upperBounds,std::size_t totalOccupation);