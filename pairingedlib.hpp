#pragma once

#include <iosfwd>
#include <vector>

struct Level
{
	double energy;
	std::size_t pairCapacity;
};

struct PairingModelParameters
{
	std::vector<Level> levels;
	double pairingStrength;
	std::size_t pairCount;
	bool oddParity;
	std::size_t totalElectronCount() const;
	std::size_t totalPairCapacity() const;
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

class Basis
{
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

class BlockedConfiguration
{
public:
	explicit BlockedConfiguration(std::vector<std::uint32_t> blockedElectrons);
	const std::vector<std::uint32_t>& blockedElectrons() const noexcept;
	std::size_t seniority() const noexcept;

private:
	std::vector<std::uint32_t> blockedElectrons_;
	std::size_t seniority_{ 0 };
};

std::uint64_t calculateSectorDegeneracy(const std::vector<Level>& levels, const BlockedConfiguration& blocked);

struct SenioritySector
{
	SenioritySector makeSenioritySector(std::size_t id, BlockedConfiguration blocked, const PairingModelParameters& model);
	std::size_t id;
	BlockedConfiguration blocked;
	std::size_t pairCount;
	std::vector<std::uint32_t> effectivePairCapacity;
	double blockedEnergy;
	std::uint64_t degeneracy;

	// Seniority wynika bezpośrednio z konfiguracji blokującej.
	std::size_t seniority() const noexcept 
	{
		return blocked.seniority();
	}
};

struct SpectrumEntry
{
	std::size_t sectorId;
	std::size_t seniority;
	double energy;
	std::uint64_t degeneracy;
};

class SpectrumWriter
{
public:
	explicit SpectrumWriter(std::ostream& output);
	void write(const SpectrumEntry& entry);
	void flush();
	std::size_t recordsWritten() const noexcept;
private:
	std::ostream& output_;
	std::size_t recordsWritten_{ 0 };
};

void validate(const PairingModelParameters& model);

std::vector<BlockedConfiguration> generateBlockedConfigurations(const std::vector<std::uint32_t>& pairCapacity,std::size_t totalSeniority);

class SectorGenerator {
public:
	explicit SectorGenerator(PairingModelParameters parameters);

	std::vector<SenioritySector> generate() const;

private:
	PairingModelParameters parameters_;
};