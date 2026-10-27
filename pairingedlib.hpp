#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
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

void validate(const PairingModelParameters& model);

std::vector<BlockedConfiguration> generateBlockedConfigurations(const std::vector<std::uint32_t>& pairCapacity,std::size_t totalSeniority);

class SectorGenerator
{
public:
	explicit SectorGenerator(PairingModelParameters parameters);
	std::vector<SenioritySector> generate() const;
private:
	PairingModelParameters parameters_;
};

class DenseMatrix {
public:
	DenseMatrix() noexcept = default;

	DenseMatrix(std::size_t rows, std::size_t cols);

	// Kopiowanie i przenoszenie.
	DenseMatrix(const DenseMatrix&) = default;
	DenseMatrix& operator=(const DenseMatrix&) = default;

	DenseMatrix(DenseMatrix&& other) noexcept;
	DenseMatrix& operator=(DenseMatrix&& other) noexcept;

	// Informacje o macierzy.
	std::size_t rows() const noexcept;
	std::size_t cols() const noexcept;
	std::size_t size() const noexcept;

	bool empty() const noexcept;
	bool isSquare() const noexcept;

	std::size_t memoryBytes() const;

	// Dostęp do elementów.
	double& operator()(
		std::size_t row,
		std::size_t col
		);

	const double& operator()(
		std::size_t row,
		std::size_t col
		) const;

	// Dostęp do ciągłego obszaru pamięci.
	double* data() noexcept;
	const double* data() const noexcept;

	// Operacje na macierzy.
	void fill(double value);

private:
	static std::size_t checkedElementCount(
		std::size_t rows,
		std::size_t cols
	);

	std::size_t offset(
		std::size_t row,
		std::size_t col
	) const;

	std::size_t rows_{ 0 };
	std::size_t cols_{ 0 };

	std::vector<double> data_;
};

class PairingHamiltonian {
public:
	PairingHamiltonian(
		const PairingModelParameters& model,
		const SenioritySector& sector,
		const Basis& basis
	);

	DenseMatrix buildMatrix() const;

private:
	const PairingModelParameters& model_;
	const SenioritySector& sector_;
	const Basis& basis_;
};

Basis generateSectorBasis(const SenioritySector& sector);

struct MemoryEstimateOptions {
	// Planowana liczba double w workspace LAPACK
	// przypadająca na jeden wymiar macierzy.
	std::size_t lapackWorkspaceDoublesPerDimension = 64;

	// Przybliżony narzut pojedynczej alokacji pamięci.
	std::size_t heapOverheadPerAllocationBytes = 32;

	// Rezerwa na pozostałe struktury, bibliotekę
	// standardową, dodatkowe bufory itd.
	std::size_t reserveBytes = 64ULL * 1024 * 1024;
};

struct MemoryEstimate {
	std::size_t sectorCount = 0;
	std::size_t totalEigenvalues = 0;

	std::size_t largestBlockDimension = 0;
	std::size_t largestBlockSectorId = 0;

	std::size_t matrixBytes = 0;
	std::size_t basisBytes = 0;
	std::size_t lookupBytes = 0;
	std::size_t lapackBytes = 0;
	std::size_t eigenvaluesBytes = 0;

	std::size_t spectrumBytes = 0;
	std::size_t sectorsBytes = 0;
	std::size_t parametersBytes = 0;

	std::size_t dpBytes = 0;
	std::size_t reserveBytes = 0;

	std::size_t peakBytes = 0;

	bool compatibleWithLapackInt = true;
};

MemoryEstimate estimateMemory(
	const PairingModelParameters& model,
	const std::vector<SenioritySector>& sectors,
	MemoryEstimateOptions options = {}
);

void printMemoryEstimate(
	const MemoryEstimate& estimate,
	std::ostream& output
);

class SymmetricDiagonalizer {
public:
	static constexpr std::size_t defaultWorkspaceMultiplier = 64;

	explicit SymmetricDiagonalizer(
		std::size_t workspaceMultiplier =
		defaultWorkspaceMultiplier
	);

	std::vector<double> diagonalize(
		DenseMatrix& matrix
	) const;

	std::size_t workspaceMultiplier() const noexcept;

private:
	std::size_t workspaceMultiplier_;
};

class Spectrum {
public:
	Spectrum() = default;

	// Zarządzanie pamięcią.
	void reserve(std::size_t count);
	void clear() noexcept;

	// Dodawanie nowych poziomów.
	void add(SpectrumEntry entry);

	// Informacje o przechowywanym widmie.
	std::size_t size() const noexcept;
	std::size_t capacity() const noexcept;
	bool empty() const noexcept;

	// Dostęp do wszystkich rekordów.
	const std::vector<SpectrumEntry>& entries() const noexcept;

	// Operacje na widmie.
	void sortByEnergy();
	double groundStateEnergy() const;

	// Liczba stanów z uwzględnieniem degeneracji.
	std::uint64_t physicalStateCount() const;

private:
	std::vector<SpectrumEntry> entries_;
};

class SpectrumWriter {
public:
	explicit SpectrumWriter(std::ostream& output);

	void write(const SpectrumEntry& entry);
	void flush();

	std::size_t recordsWritten() const noexcept;

private:
	std::ostream& output_;
	std::size_t recordsWritten_{ 0 };
};


class SpectrumReader {
public:
	explicit SpectrumReader(std::istream& input);

	// Odczytuje jeden rekord.
	// Zwraca false po osiągnięciu końca pliku.
	bool readNext(SpectrumEntry& entry);

	// Odczytuje wszystkie pozostałe rekordy.
	Spectrum readAll();

	std::size_t recordsRead() const noexcept;
	std::size_t lineNumber() const noexcept;

private:
	bool readLine(std::string& line);

	std::istream& input_;
	std::size_t lineNumber_{ 0 };
	std::size_t recordsRead_{ 0 };
};

struct ThermodynamicPoint {
	double temperature;

	long double logPartitionFunction;
	long double shiftedPartitionFunction;

	long double internalEnergy;
	long double heatCapacity;
	long double entropy;
	long double freeEnergy;

	// Z = exp(log Z). Może zwrócić 0 lub Inf,
	// jeśli Z nie mieści się w typie long double.
	long double partitionFunction() const;
};


class Thermodynamics {
public:
	explicit Thermodynamics(
		const Spectrum& spectrum,
		double boltzmannConstant = 1.0
	);

	ThermodynamicPoint evaluate(double temperature) const;

	double boltzmannConstant() const noexcept;

private:
	const Spectrum& spectrum_;

	long double groundEnergy_;
	double boltzmannConstant_;
};

struct HeatCapacityPeak
{
	double temperature;
	long double heatCapacity;
};

std::vector<HeatCapacityPeak> findHeatCapacityPeaks(
	const Thermodynamics& thermo,
	const std::vector<ThermodynamicPoint>& samples
);


// Ogólna funkcja numeryczna, niezależna od Thermodynamics.
std::vector<double> logarithmicGrid(
	double minimum,
	double maximum,
	std::size_t points
);

struct RunConfiguration {
	PairingModelParameters model;

	double temperatureMin;
	double temperatureMax;
	std::size_t temperaturePoints;

	double boltzmannConstant;
};

RunConfiguration readInput(std::istream& input);