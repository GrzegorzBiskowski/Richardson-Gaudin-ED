#include "pairingedlib.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <istream>
#include <iostream>
#include <cblas.h>
#include <limits>
#include <locale>
#include <numeric>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {

    void dsyev_(
        char* jobz,
        char* uplo,
        int* n,
        double* a,
        int* lda,
        double* w,
        double* work,
        int* lwork,
        int* info
    );

}

std::size_t PairingModelParameters::totalElectronCount() const
{
	if (oddParity)
	{
		return pairCount*2 + 1;
	}
	else
	{
		return pairCount*2;
	}
}

std::size_t PairingModelParameters::totalPairCapacity() const
{
	std::size_t result = 0;
	for (std::size_t i = 0; i < levels.size(); i++)
	{
		result += levels[i].pairCapacity;
	}
	return result;
}

const std::vector<std::uint32_t>& OccupationState::occupations() const noexcept
{
	return occupations_;
}

std::size_t OccupationState::size() const noexcept
{
	return occupations_.size();
}

std::uint64_t OccupationState::totalOccupation() const noexcept
{
	std::int32_t totalOccupation = 0;
	for (std::size_t i = 0; i < occupations_.size(); i++)
	{
		totalOccupation += occupations_[i];
	}
	return totalOccupation;
}

std::uint32_t OccupationState::operator[](std::size_t i) const
{
	return occupations_.at(i);
}

void Basis::addState(OccupationState state)
{
	states_.push_back(std::move(state));
}

const OccupationState& Basis::operator[](std::size_t i) const
{
		return states_.at(i);
}

std::size_t Basis::size() const noexcept
{
	return states_.size();
}

const std::vector<OccupationState>& Basis::states() const noexcept
{
	return states_;
}

namespace
{
    void generateRecursive(const std::vector<std::uint32_t>& upperBounds, const std::vector<std::size_t>& suffixCapacity, std::size_t level,
        std::size_t remaining, std::vector<std::uint32_t>& current, std::vector<OccupationState>& result)
    {
        const std::size_t numberOfLevels = upperBounds.size();
        // Doszliśmy do końca wszystkich poziomów.
        if (level == numberOfLevels) 
        {
            if (remaining == 0) 
            {
                result.emplace_back(current);
            }
            return;
        }

        // Ile maksymalnie możemy umieścić na obecnym poziomie?
        const std::size_t maxOccupation = std::min<std::size_t>(upperBounds[level], remaining);

        for (std::size_t occupation = 0; occupation <= maxOccupation; ++occupation)
        {
            // Po obsadzeniu obecnego poziomu musi być jeszcze
            // możliwe rozłożenie pozostałych par na kolejnych poziomach.
            const std::size_t newRemaining = remaining - occupation;
            if (newRemaining > suffixCapacity[level + 1])
            {
                continue;
            }
            current[level] = static_cast<std::uint32_t>(occupation);
            generateRecursive(upperBounds, suffixCapacity, level + 1, newRemaining, current, result );
        }
    }
} // namespace


std::vector<OccupationState> generateBoundedOccupations(const std::vector<std::uint32_t>& upperBounds, std::size_t totalOccupation)
{
    std::vector<OccupationState> result;
    // Nie ma poziomów.
    if (upperBounds.empty())
    {
        if (totalOccupation == 0)
        {
            result.emplace_back(std::vector<std::uint32_t>{});
        }
        return result;
    }
    // suffixCapacity[i] = maksymalna liczba par,
    // którą można jeszcze umieścić na poziomach i, i+1, ...
    std::vector<std::size_t> suffixCapacity(upperBounds.size() + 1, 0);
    for (std::size_t i = upperBounds.size(); i > 0; --i)
    {
        suffixCapacity[i - 1] = suffixCapacity[i] + static_cast<std::size_t>(upperBounds[i - 1]);
    }

    // Za dużo par na dostępną przestrzeń Hilberta.
    if (totalOccupation > suffixCapacity[0])
    {
        return result;
    }
    std::vector<std::uint32_t> current(upperBounds.size(), 0);
    generateRecursive(upperBounds, suffixCapacity, 0, totalOccupation, current, result);
    return result;
}

std::size_t countBoundedOccupations(const std::vector<std::uint32_t>& upperBounds, std::size_t totalOccupation)
{
    // Nie ma poziomów.
    if (upperBounds.empty())
    {
        return (totalOccupation == 0) ? 1 : 0;
    }

    // Obliczamy maksymalną możliwą liczbę par.
    std::size_t maximumOccupation = 0;

    for (std::uint32_t bound : upperBounds)
    {
        const std::size_t boundSize = static_cast<std::size_t>(bound);
        if (maximumOccupation > std::numeric_limits<std::size_t>::max() - boundSize)
        {
            throw std::overflow_error("Overflow while calculating total pair capacity.");
        }
        maximumOccupation += boundSize;
    }

    // Nie da się rozmieścić tylu par.
    if (totalOccupation > maximumOccupation)
    {
        return 0;
    }

    // dp[s] = liczba sposobów rozmieszczenia s par
    // na poziomach przetworzonych do tej pory.
    std::vector<std::size_t> dp(totalOccupation + 1, 0);
    std::vector<std::size_t> newDp(totalOccupation + 1, 0);

    // Na początku mamy zero poziomów:
    // istnieje dokładnie jeden sposób otrzymania 0 par.
    dp[0] = 1;

    for (std::uint32_t bound32 : upperBounds)
    {
        const std::size_t bound = static_cast<std::size_t>(bound32);
        std::fill(newDp.begin(), newDp.end(), 0);
        std::size_t window = 0;
        for (std::size_t s = 0; s <= totalOccupation; ++s)
        {
            // Dodajemy dp[s] do bieżącego okna.
            if (window > std::numeric_limits<std::size_t>::max() - dp[s])
            {
                throw std::overflow_error("Overflow while counting bounded occupations.");
            }
            window += dp[s];

            // Okno reprezentuje:
            //
            // dp[s] + dp[s-1] + ... + dp[s-bound]
            //
            // więc jeżeli przekroczyliśmy szerokość okna,
            // usuwamy dp[s-bound-1].
            if (s > bound)
            {
                window -= dp[s - bound - 1];
            }
            newDp[s] = window;
        }
        dp.swap(newDp);
    }
    return dp[totalOccupation];
}

namespace 
{
    std::size_t checkedAdd(std::size_t a, std::size_t b)
    {
        if (b > std::numeric_limits<std::size_t>::max() - a)
        {
            throw std::overflow_error("Integer overflow in checkedAdd().");
        }
        return a + b;
    }

    std::size_t checkedMultiply(std::size_t a, std::size_t b)
    {
        if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a)
        {
            throw std::overflow_error("Integer overflow in checkedMultiply().");
        }
        return a * b;
    }

} // namespace

BlockedConfiguration::BlockedConfiguration(std::vector<std::uint32_t> blockedElectrons)
    : blockedElectrons_(std::move(blockedElectrons)), seniority_(0)
{
    for (const std::uint32_t count : blockedElectrons_)
    {
        seniority_ = checkedAdd(seniority_, static_cast<std::size_t>(count));
    }
}

const std::vector<std::uint32_t>& BlockedConfiguration::blockedElectrons() const noexcept
{
    return blockedElectrons_;
}

std::size_t BlockedConfiguration::seniority() const noexcept
{
    return seniority_;
}

SenioritySector makeSenioritySector(std::size_t id, BlockedConfiguration blocked, const PairingModelParameters& model)
{
    const auto& levels = model.levels;
    const auto& v = blocked.blockedElectrons();

    if (v.size() != levels.size())
    {
        throw std::invalid_argument("Blocked configuration size does not match the number of energy levels.");
    }

    // Liczba elektronów w całym układzie.
    const std::size_t electronCount = checkedAdd(checkedMultiply(model.pairCount, 2), model.oddParity ? 1 : 0);
    const std::size_t nu = blocked.seniority();
    if (nu > electronCount)
    {
        throw std::invalid_argument("Seniority exceeds the total electron count.");
    }

    if ((electronCount - nu) % 2 != 0)
    {
        throw std::invalid_argument("Invalid seniority parity for this electron number.");
    }

    const std::size_t numberOfPairs = (electronCount - nu) / 2;

    std::vector<std::uint32_t> effectiveCapacity;
    effectiveCapacity.reserve(levels.size());

    std::size_t totalEffectiveCapacity = 0;
    long double blockedEnergySum = 0.0L;

    for (std::size_t j = 0; j < levels.size(); ++j)
    {
        const std::uint32_t omega = levels[j].pairCapacity;

        if (v[j] > omega)
        {
            throw std::invalid_argument("Blocked electron count exceeds pair capacity at level " + std::to_string(j) + ".");
        }
        const std::uint32_t effectiveOmega = omega - v[j];
        effectiveCapacity.push_back(effectiveOmega);
        totalEffectiveCapacity = checkedAdd(totalEffectiveCapacity, static_cast<std::size_t>(effectiveOmega));
        blockedEnergySum += static_cast<long double>(levels[j].energy) * static_cast<long double>(v[j]);
    }

    if (numberOfPairs > totalEffectiveCapacity)
    {
        throw std::invalid_argument("The sector has more pairs than its effective pair capacity allows.");
    }

    const double blockedEnergy = static_cast<double>(blockedEnergySum);

    if (!std::isfinite(blockedEnergySum) || !std::isfinite(blockedEnergy))
    {
        throw std::overflow_error("Blocked-electron energy is not finite.");
    }

    const std::uint64_t degeneracy = calculateSectorDegeneracy(levels, blocked);

    return SenioritySector{id, std::move(blocked), numberOfPairs, std::move(effectiveCapacity), blockedEnergy, degeneracy};
}

namespace
{
    std::uint64_t checkedMultiplyU64(std::uint64_t a,std::uint64_t b)
    {
        if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        {
            throw std::overflow_error("Overflow while calculating sector degeneracy.");
        }
        return a * b;
    }

    // Oblicza symbol Newtona n choose k.
    // Rzuca overflow_error, jeśli wynik nie mieści się
    // w uint64_t.
    std::uint64_t binomialCoefficientChecked(std::uint64_t n, std::uint64_t k)
    {
        if (k > n)
        {
            return 0;
        }
        // C(n,k) = C(n,n-k).
        k = (k < n - k) ? k : n - k;
        std::uint64_t result = 1;
        for (std::uint64_t i = 1; i <= k; ++i)
        {
            std::uint64_t numerator = n - k + i;
            std::uint64_t denominator = i;
            // Skracamy licznik i mianownik.
            const auto g1 = std::gcd(numerator, denominator);

            numerator /= g1;
            denominator /= g1;

            // Skracamy pozostały mianownik z dotychczasowym
            // wynikiem. Po tym kroku denominator powinien
            // wynosić 1, ponieważ współczynnik dwumianowy
            // jest liczbą całkowitą.
            const auto g2 = std::gcd(result, denominator);
            result /= g2;
            denominator /= g2;

            if (denominator != 1)
            {
                throw std::logic_error("Unexpected non-integer intermediate result in binomialCoefficientChecked().");
            }
            result = checkedMultiplyU64(result, numerator);
        }
        return result;
    }

} // namespace

std::uint64_t calculateSectorDegeneracy(const std::vector<Level>& levels, const BlockedConfiguration& blocked)
{
    const auto& v = blocked.blockedElectrons();

    if (levels.size() != v.size())
    {
        throw std::invalid_argument("Cannot calculate sector degeneracy: the number of levels does not match the blocked configuration size.");
    }

    std::uint64_t totalDegeneracy = 1;

    for (std::size_t j = 0; j < levels.size(); ++j)
    {
        const std::uint32_t omega = levels[j].pairCapacity;
        const std::uint32_t seniority = v[j];
        if (omega == 0)
        {
            throw std::invalid_argument("Cannot calculate sector degeneracy: pair capacity is zero at level " + std::to_string(j) + ".");
        }
        if (seniority > omega) {
            throw std::invalid_argument("Cannot calculate sector degeneracy: partial seniority exceeds pair capacity at level " + std::to_string(j) + ".");
        }
        // Poziom ma 2*Omega stanów jednocząstkowych.
        const std::uint64_t singleParticleStates = checkedMultiplyU64(2, static_cast<std::uint64_t>(omega));
        const std::uint64_t firstTerm = binomialCoefficientChecked(singleParticleStates,seniority);
        std::uint64_t secondTerm = 0;

        if (seniority >= 2) 
        {
            secondTerm = binomialCoefficientChecked(singleParticleStates, seniority - 2);
        }

        if (secondTerm > firstTerm)
        {
            throw std::logic_error("Invalid local degeneracy formula result at level " + std::to_string(j) + ".");
        }

        const std::uint64_t localDegeneracy = firstTerm - secondTerm;

        totalDegeneracy = checkedMultiplyU64(totalDegeneracy, localDegeneracy);
    }
    return totalDegeneracy;
}

void validate(const PairingModelParameters& model)
{
    if (model.levels.empty())
    {
        throw std::invalid_argument("Invalid model: at least one energy level is required.");
    }

    if (!std::isfinite(model.pairingStrength))
    {
        throw std::invalid_argument("Invalid model: pairingStrength must be finite.");
    }

    std::size_t totalPairCapacity = 0;

    for (std::size_t j = 0; j < model.levels.size(); ++j)
    {
        const auto& level = model.levels[j];

        if (!std::isfinite(level.energy))
        {
            throw std::invalid_argument("Invalid model: non-finite energy at level " + std::to_string(j) + ".");
        }

        if (level.pairCapacity == 0)
        {
            throw std::invalid_argument("Invalid model: pairCapacity must be positive at level " + std::to_string(j) + ".");
        }

        totalPairCapacity = checkedAdd(totalPairCapacity, static_cast<std::size_t>(level.pairCapacity));
    }

    // N_e = 2M + b
    const std::size_t electronCount = checkedAdd(checkedMultiply(model.pairCount, 2), model.oddParity ? 1 : 0);

    // Maksymalna liczba elektronów: 2 * sum_j Omega_j.
    const std::size_t maximumElectronCount = checkedMultiply(totalPairCapacity, 2);

    if (electronCount > maximumElectronCount)
    {
        throw std::invalid_argument("Invalid model: the number of electrons exceeds the total capacity of the energy levels.");
    }
}

std::vector<BlockedConfiguration> generateBlockedConfigurations(const std::vector<std::uint32_t>& pairCapacity, std::size_t totalSeniority)
{
    const auto states = generateBoundedOccupations(pairCapacity, totalSeniority);
    std::vector<BlockedConfiguration> result;
    result.reserve(states.size());
    for (const auto& state : states)
    {
        result.emplace_back(state.occupations());
    }
    return result;
}

SectorGenerator::SectorGenerator(PairingModelParameters parameters)
    : parameters_(std::move(parameters))
{
    validate(parameters_);
}

namespace
{
    SenioritySector buildSector(std::size_t id, BlockedConfiguration blocked, const PairingModelParameters& model, std::size_t electronCount)
    {
        const auto& levels = model.levels;
        const auto& v = blocked.blockedElectrons();

        if (v.size() != levels.size())
        {
            throw std::logic_error("Internal error: blocked configuration has incorrect size.");
        }

        const std::size_t nu = blocked.seniority();

        if (nu > electronCount || (electronCount - nu) % 2 != 0)
        {
            throw std::logic_error("Internal error: incompatible seniority.");
        }

        const std::size_t numberOfPairs =(electronCount - nu) / 2;

        std::vector<std::uint32_t> effectiveCapacity;
        effectiveCapacity.reserve(levels.size());

        std::size_t totalEffectiveCapacity = 0;
        long double blockedEnergySum = 0.0L;

        for (std::size_t j = 0; j < levels.size(); ++j)
        {
            const auto omega = levels[j].pairCapacity;
            if (v[j] > omega)
            {
                throw std::logic_error("Internal error: blocked electron count exceeds pair capacity.");
            }

            const std::uint32_t effectiveOmega = omega - v[j];
            effectiveCapacity.push_back(effectiveOmega);

            totalEffectiveCapacity = checkedAdd(totalEffectiveCapacity, static_cast<std::size_t>(effectiveOmega));

            blockedEnergySum += static_cast<long double>(levels[j].energy) * static_cast<long double>(v[j]);
        }

        if (numberOfPairs > totalEffectiveCapacity)
        {
            throw std::logic_error("Internal error: sector has more pairs than its effective capacity allows.");
        }

        const double blockedEnergy = static_cast<double>(blockedEnergySum);

        if (!std::isfinite(blockedEnergySum) || !std::isfinite(blockedEnergy))
        {
            throw std::overflow_error("Blocked-electron energy is not finite.");
        }

        // Oblicz degenerację, zanim przeniesiemy konfigurację.
        const std::uint64_t degeneracy = calculateSectorDegeneracy(levels, blocked);

        return SenioritySector{id, std::move(blocked), numberOfPairs, std::move(effectiveCapacity), blockedEnergy, degeneracy};
    }

} // namespace

std::vector<SenioritySector> SectorGenerator::generate() const
{
    const auto& levels = parameters_.levels;

    std::vector<std::uint32_t> pairCapacity;
    pairCapacity.reserve(levels.size());

    std::size_t totalPairCapacity = 0;

    for (const auto& level : levels)
    {
        pairCapacity.push_back(level.pairCapacity);

        totalPairCapacity = checkedAdd(totalPairCapacity, static_cast<std::size_t>(level.pairCapacity));
    }

    const std::size_t electronCount = checkedAdd(checkedMultiply(parameters_.pairCount, 2), parameters_.oddParity ? 1 : 0);
    const std::size_t maximumElectronCount = checkedMultiply(totalPairCapacity, 2);
    const std::size_t maximumSeniority = std::min(electronCount, maximumElectronCount - electronCount);
    const std::size_t firstSeniority = parameters_.oddParity ? 1 : 0;

    std::vector<SenioritySector> sectors;
    std::size_t nextId = 0;

    for (std::size_t nu = firstSeniority; nu <= maximumSeniority; )
    {
        auto configurations = generateBlockedConfigurations(pairCapacity, nu);

        for (auto& blocked : configurations)
        {
            if (nextId == std::numeric_limits<std::size_t>::max())
            {
                throw std::overflow_error("Too many seniority sectors.");
            }

            sectors.push_back(buildSector(nextId, std::move(blocked), parameters_, electronCount));
            ++nextId;
        }

        // Nie dodajemy 2, jeśli wyszlibyśmy poza zakres.
        if (maximumSeniority - nu < 2)
        {
            break;
        }
        nu = checkedAdd(nu, 2);
    }
    return sectors;
}

std::size_t DenseMatrix::checkedElementCount(
    std::size_t rows,
    std::size_t cols
) {
    if ((rows == 0) != (cols == 0)) {
        throw std::invalid_argument(
            "DenseMatrix: only a 0x0 matrix "
            "may have a zero dimension."
        );
    }

    return checkedMultiply(rows, cols);
}

DenseMatrix::DenseMatrix(
    std::size_t rows,
    std::size_t cols
)
    : rows_(rows),
    cols_(cols),
    data_(checkedElementCount(rows, cols), 0.0)
{}

DenseMatrix::DenseMatrix(DenseMatrix&& other) noexcept
    : rows_(std::exchange(other.rows_, 0)),
    cols_(std::exchange(other.cols_, 0)),
    data_(std::move(other.data_))
{
    other.data_.clear();
}

DenseMatrix& DenseMatrix::operator=(
    DenseMatrix&& other
    ) noexcept
{
    if (this != &other) {
        data_ = std::move(other.data_);

        rows_ = std::exchange(other.rows_, 0);
        cols_ = std::exchange(other.cols_, 0);

        other.data_.clear();
    }

    return *this;
}

std::size_t DenseMatrix::rows() const noexcept
{
    return rows_;
}

std::size_t DenseMatrix::cols() const noexcept
{
    return cols_;
}

std::size_t DenseMatrix::size() const noexcept
{
    return data_.size();
}

bool DenseMatrix::empty() const noexcept
{
    return data_.empty();
}

bool DenseMatrix::isSquare() const noexcept
{
    return rows_ == cols_;
}

std::size_t DenseMatrix::memoryBytes() const
{
    return checkedMultiply(
        data_.capacity(),
        sizeof(double)
    );
}

std::size_t DenseMatrix::offset(
    std::size_t row,
    std::size_t col
) const
{
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range(
            "DenseMatrix: index (" +
            std::to_string(row) + ", " +
            std::to_string(col) +
            ") out of range for " +
            std::to_string(rows_) + "x" +
            std::to_string(cols_) + " matrix."
        );
    }

    return col * rows_ + row;
}

double& DenseMatrix::operator()(
    std::size_t row,
    std::size_t col
    )
{
    return data_[offset(row, col)];
}

const double& DenseMatrix::operator()(
    std::size_t row,
    std::size_t col
    ) const
{
    return data_[offset(row, col)];
}

double* DenseMatrix::data() noexcept
{
    return data_.data();
}

const double* DenseMatrix::data() const noexcept
{
    return data_.data();
}

void DenseMatrix::fill(double value)
{
    std::fill(
        data_.begin(),
        data_.end(),
        value
    );
}

namespace {

    struct OccupationHash {

        std::size_t operator()(
            const std::vector<std::uint32_t>& occupations
            ) const noexcept
        {
            std::size_t seed = 0;

            for (std::uint32_t occupation : occupations) {
                seed ^= std::hash<std::uint32_t>{}(occupation)
                    +std::size_t{ 0x9e3779b9 }
                    + (seed << 6)
                    + (seed >> 2);
            }

            return seed;
        }
    };

} // namespace

PairingHamiltonian::PairingHamiltonian(
    const PairingModelParameters& model,
    const SenioritySector& sector,
    const Basis& basis
)
    : model_(model),
    sector_(sector),
    basis_(basis)
{}

DenseMatrix PairingHamiltonian::buildMatrix() const
{
    const std::size_t D = basis_.size();
    const std::size_t L = model_.levels.size();

    const auto& capacity = sector_.effectivePairCapacity;
    const double g = model_.pairingStrength;

    // Macierz D x D, początkowo wypełniona zerami.
    DenseMatrix H(D, D);

    // Odwzorowanie wektora obsadzeń na indeks stanu.
    std::unordered_map<
        std::vector<std::uint32_t>,
        std::size_t,
        OccupationHash
    > stateIndex;

    stateIndex.reserve(D);

    for (std::size_t alpha = 0; alpha < D; ++alpha) {
        stateIndex.emplace(
            basis_[alpha].occupations(),
            alpha
        );
    }

    // Przechodzimy po kolejnych stanach bazowych.
    for (std::size_t alpha = 0; alpha < D; ++alpha) {

        const auto& n = basis_[alpha].occupations();

        // ==========================================
        // 1. Element diagonalny
        // ==========================================

        double diagonal = 0.0;

        for (std::size_t j = 0; j < L; ++j) {

            const double nj =
                static_cast<double>(n[j]);

            const double omegaJ =
                static_cast<double>(capacity[j]);

            const double epsilonJ =
                model_.levels[j].energy;

            diagonal += 2.0 * epsilonJ * nj;

            diagonal -= g * nj * (omegaJ - nj + 1.0);
        }

        H(alpha, alpha) = diagonal;

        // ==========================================
        // 2. Elementy pozadiagonalne
        // ==========================================

        // Roboczy wektor, który będziemy modyfikować
        // podczas przenoszenia pojedynczej pary.
        auto moved = n;

        // Poziom, z którego zabieramy parę.
        for (std::size_t j = 0; j < L; ++j) {

            if (n[j] == 0) {
                continue;
            }

            // Poziom, na który przenosimy parę.
            for (std::size_t k = 0; k < L; ++k) {

                if (j == k || n[k] >= capacity[k]) {
                    continue;
                }

                // Tworzymy konfigurację po przeniesieniu.
                --moved[j];
                ++moved[k];

                // Indeks nowego stanu bazowego.
                const std::size_t beta =
                    stateIndex.at(moved);

                // Przywracamy pierwotny stan.
                ++moved[j];
                --moved[k];

                // Wypełniamy tylko jedną połowę macierzy,
                // a następnie kopiujemy wartość symetrycznie.
                if (beta <= alpha) {
                    continue;
                }

                const double nj =
                    static_cast<double>(n[j]);

                const double nk =
                    static_cast<double>(n[k]);

                const double omegaJ =
                    static_cast<double>(capacity[j]);

                const double omegaK =
                    static_cast<double>(capacity[k]);

                const double factor = std::sqrt(
                    nj *
                    (omegaJ - nj + 1.0) *
                    (nk + 1.0) *
                    (omegaK - nk)
                );

                const double element = -g * factor;

                H(alpha, beta) = element;
                H(beta, alpha) = element;
            }
        }
    }

    return H;
}

Basis generateSectorBasis(const SenioritySector& sector)
{
    auto states = generateBoundedOccupations(
        sector.effectivePairCapacity,
        sector.pairCount
    );

    Basis basis;

    for (auto& state : states) {
        basis.addState(std::move(state));
    }

    return basis;
}

MemoryEstimate estimateMemory(
    const PairingModelParameters& model,
    const std::vector<SenioritySector>& sectors,
    MemoryEstimateOptions options
) {
    MemoryEstimate estimate;

    estimate.sectorCount = sectors.size();
    estimate.reserveBytes = options.reserveBytes;

    const std::size_t L = model.levels.size();

    // --------------------------------------------------
    // 1. Koszt przechowywania pojedynczego stanu
    // --------------------------------------------------

    const std::size_t occupationBytes =
        checkedMultiply(L, sizeof(std::uint32_t));

    // Uwzględniamy OccupationState, jego wektor danych
    // oraz zapas na tymczasowe kontenery podczas
    // generowania Basis.
    const std::size_t basisBytesPerState = checkedAdd(
        checkedMultiply(4, sizeof(OccupationState)),
        checkedAdd(
            occupationBytes,
            options.heapOverheadPerAllocationBytes
        )
    );

    // --------------------------------------------------
    // 2. Koszt indeksu stanów w PairingHamiltonian
    // --------------------------------------------------

    using IndexRecord = std::pair<
        const std::vector<std::uint32_t>,
        std::size_t
    >;

    const std::size_t indexBytesPerState = checkedAdd(
        sizeof(IndexRecord),
        checkedAdd(
            checkedMultiply(4, sizeof(void*)),
            checkedAdd(
                occupationBytes,
                checkedMultiply(
                    2,
                    options.heapOverheadPerAllocationBytes
                )
            )
        )
    );

    // --------------------------------------------------
    // 3. Pamięć parametrów fizycznych
    // --------------------------------------------------

    const std::size_t modelBytes = checkedMultiply(
        model.levels.capacity(),
        sizeof(Level)
    );

    // Przewidujemy kopię w main oraz w SectorGenerator.
    estimate.parametersBytes =
        checkedMultiply(2, modelBytes);

    // --------------------------------------------------
    // 4. Pamięć wszystkich opisów sektorów
    // --------------------------------------------------

    estimate.sectorsBytes = checkedMultiply(
        sectors.capacity(),
        sizeof(SenioritySector)
    );

    for (const auto& sector : sectors) {
        const std::size_t effectiveCapacityBytes =
            checkedMultiply(
                sector.effectivePairCapacity.capacity(),
                sizeof(std::uint32_t)
            );

        const std::size_t blockedBytes =
            checkedMultiply(
                sector.blocked.blockedElectrons().capacity(),
                sizeof(std::uint32_t)
            );

        estimate.sectorsBytes = checkedAdd(
            estimate.sectorsBytes,
            effectiveCapacityBytes
        );

        estimate.sectorsBytes = checkedAdd(
            estimate.sectorsBytes,
            blockedBytes
        );

        estimate.sectorsBytes = checkedAdd(
            estimate.sectorsBytes,
            checkedMultiply(
                2,
                options.heapOverheadPerAllocationBytes
            )
        );
    }

    // --------------------------------------------------
    // 5. Analiza poszczególnych sektorów
    // --------------------------------------------------

    std::size_t largestWorkingBytes = 0;

    for (const auto& sector : sectors) {

        // Wymiar Hamiltonianu bez generowania Basis.
        const std::size_t D = countBoundedOccupations(
            sector.effectivePairCapacity,
            sector.pairCount
        );

        if (D == 0) {
            throw std::logic_error(
                "Memory estimation encountered "
                "a zero-dimensional sector."
            );
        }

        estimate.totalEigenvalues = checkedAdd(
            estimate.totalEigenvalues,
            D
        );

        // Gęsta macierz D x D typu double.
        const std::size_t matrixBytes =
            checkedMultiply(
                checkedMultiply(D, D),
                sizeof(double)
            );

        // Basis wraz z szacowanymi buforami.
        const std::size_t basisBytes =
            checkedMultiply(D, basisBytesPerState);

        // Tablica haszująca PairingHamiltonian.
        const std::size_t indexBytes =
            checkedMultiply(D, indexBytesPerState);

        // Workspace LAPACK.
        const std::size_t lapackDoubles =
            std::max<std::size_t>(
                1,
                checkedMultiply(
                    D,
                    options.lapackWorkspaceDoublesPerDimension
                )
            );

        const std::size_t lapackBytes =
            checkedMultiply(
                lapackDoubles,
                sizeof(double)
            );

        // Wektor wartości własnych sektora.
        const std::size_t eigenvaluesBytes =
            checkedMultiply(D, sizeof(double));

        // Konserwatywny koszt pracy z jednym sektorem.
        const std::size_t workingBytes = checkedAdd(
            checkedAdd(
                matrixBytes,
                basisBytes
            ),
            checkedAdd(
                checkedAdd(indexBytes, lapackBytes),
                eigenvaluesBytes
            )
        );

        // Pamięć wykorzystywana przez dwa wektory
        // programowania dynamicznego.
        const std::size_t dpBytes =
            checkedMultiply(
                checkedMultiply(
                    checkedAdd(sector.pairCount, 1),
                    sizeof(std::size_t)
                ),
                2
            );

        estimate.dpBytes = std::max(
            estimate.dpBytes,
            dpBytes
        );

        // Zachowujemy informacje o najdroższym bloku.
        if (workingBytes > largestWorkingBytes) {
            largestWorkingBytes = workingBytes;

            estimate.largestBlockSectorId = sector.id;
            estimate.largestBlockDimension = D;

            estimate.matrixBytes = matrixBytes;
            estimate.basisBytes = basisBytes;
            estimate.lookupBytes = indexBytes;
            estimate.lapackBytes = lapackBytes;
            estimate.eigenvaluesBytes = eigenvaluesBytes;
        }
    }

    // --------------------------------------------------
    // 6. Pamięć pełnego widma
    // --------------------------------------------------

    estimate.spectrumBytes = checkedMultiply(
        estimate.totalEigenvalues,
        sizeof(SpectrumEntry)
    );

    // --------------------------------------------------
    // 7. Zgodność z LP64 LAPACK
    // --------------------------------------------------

    const std::size_t lapackIntMax =
        static_cast<std::size_t>(
            std::numeric_limits<int>::max()
            );

    const std::size_t largestWorkspaceCount =
        checkedMultiply(
            estimate.largestBlockDimension,
            options.lapackWorkspaceDoublesPerDimension
        );

    estimate.compatibleWithLapackInt =
        estimate.largestBlockDimension <= lapackIntMax &&
        largestWorkspaceCount <= lapackIntMax;

    // --------------------------------------------------
    // 8. Całkowita estymacja maksymalnego zużycia RAM
    // --------------------------------------------------

    const std::size_t persistentBytes = checkedAdd(
        checkedAdd(
            estimate.sectorsBytes,
            estimate.parametersBytes
        ),
        estimate.reserveBytes
    );

    // Faza obliczeniowa: pełne widmo + najdroższy blok.
    const std::size_t calculationBytes = checkedAdd(
        estimate.spectrumBytes,
        largestWorkingBytes
    );

    // Faza przygotowawcza: pamięć DP.
    const std::size_t transientBytes = std::max(
        calculationBytes,
        estimate.dpBytes
    );

    estimate.peakBytes = checkedAdd(
        persistentBytes,
        transientBytes
    );

    return estimate;
}

void printMemoryEstimate(
    const MemoryEstimate& e,
    std::ostream& output
) {
    constexpr double MiB = 1024.0 * 1024.0;
    constexpr double GiB = MiB * 1024.0;

    output << std::fixed << std::setprecision(2);

    output << "\n=== MEMORY ESTIMATE ===\n\n";

    output << "Number of sectors: "
        << e.sectorCount << '\n';

    output << "Total eigenvalues: "
        << e.totalEigenvalues << '\n';

    output << "Largest Hamiltonian dimension: "
        << e.largestBlockDimension << '\n';

    output << "Largest block sector ID: "
        << e.largestBlockSectorId << "\n\n";

    output << "Largest Hamiltonian: "
        << e.matrixBytes / MiB << " MiB\n";

    output << "Basis and temporary storage: "
        << e.basisBytes / MiB << " MiB\n";

    output << "State lookup table: "
        << e.lookupBytes / MiB << " MiB\n";

    output << "LAPACK workspace: "
        << e.lapackBytes / MiB << " MiB\n";

    output << "Sector eigenvalues: "
        << e.eigenvaluesBytes / MiB << " MiB\n";

    output << "Complete spectrum: "
        << e.spectrumBytes / MiB << " MiB\n";

    output << "Sector descriptions: "
        << e.sectorsBytes / MiB << " MiB\n";

    output << "DP working memory: "
        << e.dpBytes / MiB << " MiB\n";

    output << "Additional reserve: "
        << e.reserveBytes / MiB << " MiB\n";

    output << "\nESTIMATED PEAK RAM: "
        << e.peakBytes / GiB << " GiB\n";

    output << "LP64 LAPACK compatible: "
        << (e.compatibleWithLapackInt ? "yes" : "no")
        << '\n';
}

SymmetricDiagonalizer::SymmetricDiagonalizer(
    std::size_t workspaceMultiplier
)
    : workspaceMultiplier_(workspaceMultiplier)
{
    if (workspaceMultiplier_ < 3) {
        throw std::invalid_argument(
            "SymmetricDiagonalizer: workspace multiplier "
            "must be at least 3."
        );
    }
}

std::vector<double> SymmetricDiagonalizer::diagonalize(
    DenseMatrix& matrix
) const
{
    // ==========================================
    // 1. Wymiary macierzy
    // ==========================================

    if (!matrix.isSquare() || matrix.empty()) {
        throw std::invalid_argument(
            "SymmetricDiagonalizer: matrix must be "
            "non-empty and square."
        );
    }

    const std::size_t D = matrix.rows();

    // ==========================================
    // 2. Kontrola symetrii i wartości skończonych
    // ==========================================

    constexpr double tolerance = 1e-12;

    for (std::size_t i = 0; i < D; ++i) {

        if (!std::isfinite(matrix(i, i))) {
            throw std::invalid_argument(
                "SymmetricDiagonalizer: matrix contains "
                "a non-finite diagonal element."
            );
        }

        for (std::size_t j = i + 1; j < D; ++j) {

            const double a = matrix(i, j);
            const double b = matrix(j, i);

            if (!std::isfinite(a) || !std::isfinite(b)) {
                throw std::invalid_argument(
                    "SymmetricDiagonalizer: matrix contains "
                    "NaN or Inf."
                );
            }

            const double scale = std::max({
                1.0,
                std::abs(a),
                std::abs(b)
                });

            if (std::abs(a - b) > tolerance * scale) {
                throw std::invalid_argument(
                    "SymmetricDiagonalizer: matrix is not symmetric."
                );
            }
        }
    }

    // ==========================================
    // 3. Zgodność z LP64 LAPACK
    // ==========================================

    const std::size_t maxLapackInt =
        static_cast<std::size_t>(
            std::numeric_limits<int>::max()
            );

    // Zapewnia, że zarówno N, jak i LWORK
    // mieszczą się w typie int.
    if (D > maxLapackInt / workspaceMultiplier_) {
        throw std::length_error(
            "SymmetricDiagonalizer: matrix dimension or "
            "workspace exceeds LP64 LAPACK integer range."
        );
    }

    int n = static_cast<int>(D);
    int lda = n;

    int lwork = static_cast<int>(
        D * workspaceMultiplier_
        );

    // ==========================================
    // 4. Alokacja pamięci
    // ==========================================

    std::vector<double> eigenvalues(D);

    std::vector<double> work(
        static_cast<std::size_t>(lwork)
    );

    // ==========================================
    // 5. Parametry LAPACK
    // ==========================================

    char jobz = 'N';  // Tylko wartości własne
    char uplo = 'U';  // Górny trójkąt macierzy

    int info = 0;

    // ==========================================
    // 6. Diagonalizacja
    // ==========================================

    dsyev_(
        &jobz,
        &uplo,
        &n,
        matrix.data(),
        &lda,
        eigenvalues.data(),
        work.data(),
        &lwork,
        &info
    );

    // ==========================================
    // 7. Weryfikacja wyniku
    // ==========================================

    if (info < 0) {
        throw std::runtime_error(
            "LAPACK DSYEV: invalid argument number " +
            std::to_string(-info) + "."
        );
    }

    if (info > 0) {
        throw std::runtime_error(
            "LAPACK DSYEV: diagonalization did not converge. "
            "INFO = " + std::to_string(info) + "."
        );
    }

    for (double E : eigenvalues) {
        if (!std::isfinite(E)) {
            throw std::runtime_error(
                "LAPACK DSYEV returned a non-finite eigenvalue."
            );
        }
    }

    return eigenvalues;
}

std::size_t
SymmetricDiagonalizer::workspaceMultiplier() const noexcept
{
    return workspaceMultiplier_;
}

void Spectrum::reserve(std::size_t count)
{
    entries_.reserve(count);
}

void Spectrum::clear() noexcept
{
    entries_.clear();
}

void Spectrum::add(SpectrumEntry entry)
{
    if (!std::isfinite(entry.energy)) {
        throw std::invalid_argument(
            "Spectrum::add(): energy must be finite."
        );
    }

    if (entry.degeneracy == 0) {
        throw std::invalid_argument(
            "Spectrum::add(): degeneracy must be positive."
        );
    }

    entries_.push_back(entry);
}

std::size_t Spectrum::size() const noexcept
{
    return entries_.size();
}

std::size_t Spectrum::capacity() const noexcept
{
    return entries_.capacity();
}

bool Spectrum::empty() const noexcept
{
    return entries_.empty();
}

const std::vector<SpectrumEntry>&
Spectrum::entries() const noexcept
{
    return entries_;
}

void Spectrum::sortByEnergy()
{
    std::sort(
        entries_.begin(),
        entries_.end(),

        [](const SpectrumEntry& a,
            const SpectrumEntry& b)
        {
            if (a.energy != b.energy) {
                return a.energy < b.energy;
            }

            if (a.sectorId != b.sectorId) {
                return a.sectorId < b.sectorId;
            }

            return a.seniority < b.seniority;
        }
    );
}

double Spectrum::groundStateEnergy() const
{
    if (entries_.empty()) {
        throw std::logic_error(
            "Spectrum::groundStateEnergy(): empty spectrum."
        );
    }

    const auto it = std::min_element(
        entries_.begin(),
        entries_.end(),

        [](const SpectrumEntry& a,
            const SpectrumEntry& b)
        {
            return a.energy < b.energy;
        }
    );

    return it->energy;
}

std::uint64_t Spectrum::physicalStateCount() const
{
    std::uint64_t total = 0;

    constexpr auto maxValue =
        std::numeric_limits<std::uint64_t>::max();

    for (const auto& entry : entries_) {
        if (entry.degeneracy > maxValue - total) {
            throw std::overflow_error(
                "Spectrum::physicalStateCount(): "
                "state count exceeds uint64_t."
            );
        }

        total += entry.degeneracy;
    }

    return total;
}

namespace {

    constexpr const char* spectrumMagic =
        "# pairingedlib spectrum v1";

    constexpr const char* spectrumColumns =
        "# columns: sector_id seniority energy degeneracy";

} // namespace

SpectrumWriter::SpectrumWriter(std::ostream& output)
    : output_(output)
{
    if (!output_) {
        throw std::runtime_error(
            "SpectrumWriter: output stream is not writable."
        );
    }

    output_ << spectrumMagic << '\n'
        << spectrumColumns << '\n';

    if (!output_) {
        throw std::runtime_error(
            "SpectrumWriter: cannot write file header."
        );
    }
}

void SpectrumWriter::write(const SpectrumEntry& entry)
{
    if (!std::isfinite(entry.energy)) {
        throw std::invalid_argument(
            "SpectrumWriter: energy must be finite."
        );
    }

    if (entry.degeneracy == 0) {
        throw std::invalid_argument(
            "SpectrumWriter: degeneracy must be positive."
        );
    }

    if (recordsWritten_ ==
        std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error(
            "SpectrumWriter: record counter overflow."
        );
    }

    // Przygotowanie pojedynczego rekordu.
    std::ostringstream line;

    line.imbue(std::locale::classic());

    line << entry.sectorId << ' '
        << entry.seniority << ' '
        << std::setprecision(
            std::numeric_limits<double>::max_digits10
        )
        << entry.energy << ' '
        << entry.degeneracy << '\n';

    if (!line) {
        throw std::runtime_error(
            "SpectrumWriter: record formatting failed."
        );
    }

    // Zapis gotowej linii.
    output_ << line.str();

    if (!output_) {
        throw std::runtime_error(
            "SpectrumWriter: record write failed."
        );
    }

    ++recordsWritten_;
}

void SpectrumWriter::flush()
{
    output_.flush();

    if (!output_) {
        throw std::runtime_error(
            "SpectrumWriter: flushing failed."
        );
    }
}

std::size_t SpectrumWriter::recordsWritten() const noexcept
{
    return recordsWritten_;
}

namespace {

    template <typename UInt>
    UInt parseUnsigned(
        const std::string& token,
        const char* field,
        std::size_t lineNumber
    ) {
        UInt value{};

        const char* first = token.data();
        const char* last = first + token.size();

        const auto result = std::from_chars(
            first,
            last,
            value,
            10
        );

        if (result.ec != std::errc{} ||
            result.ptr != last) {
            throw std::runtime_error(
                "SpectrumReader: line " +
                std::to_string(lineNumber) +
                ": invalid " + field +
                " ('" + token + "')."
            );
        }

        return value;
    }


    double parseEnergy(
        const std::string& token,
        std::size_t lineNumber
    ) {
        std::istringstream parser(token);

        parser.imbue(std::locale::classic());

        double energy{};
        char extra;

        if (!(parser >> energy) ||
            (parser >> extra) ||
            !std::isfinite(energy)) {
            throw std::runtime_error(
                "SpectrumReader: line " +
                std::to_string(lineNumber) +
                ": invalid energy ('" + token + "')."
            );
        }

        return energy;
    }

} // namespace

bool SpectrumReader::readLine(std::string& line)
{
    try {
        if (std::getline(input_, line)) {
            ++lineNumber_;

            // Obsługa zakończeń linii Windows (CRLF).
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }

            return true;
        }
    }
    catch (const std::ios_base::failure&) {
        if (input_.bad() || !input_.eof()) {
            throw std::runtime_error(
                "SpectrumReader: I/O error while reading line."
            );
        }
    }

    if (input_.bad() || !input_.eof()) {
        throw std::runtime_error(
            "SpectrumReader: I/O error while reading line."
        );
    }

    return false;
}

SpectrumReader::SpectrumReader(std::istream& input)
    : input_(input)
{
    std::string line;

    if (!readLine(line) || line != spectrumMagic) {
        throw std::runtime_error(
            "SpectrumReader: unsupported or missing "
            "spectrum header."
        );
    }

    if (!readLine(line) || line != spectrumColumns) {
        throw std::runtime_error(
            "SpectrumReader: unexpected spectrum columns."
        );
    }
}

bool SpectrumReader::readNext(SpectrumEntry& entry)
{
    std::string line;

    while (readLine(line)) {

        // Ignorujemy puste linie i komentarze.
        const auto first =
            line.find_first_not_of(" \t\r\n");

        if (first == std::string::npos ||
            line[first] == '#') {
            continue;
        }

        // Rozdzielamy rekord na cztery pola.
        std::istringstream fields(line);

        fields.imbue(std::locale::classic());

        std::string sector;
        std::string seniority;
        std::string energy;
        std::string degeneracy;
        std::string extra;

        if (!(fields >> sector
            >> seniority
            >> energy
            >> degeneracy) ||
            (fields >> extra)) {
            throw std::runtime_error(
                "SpectrumReader: line " +
                std::to_string(lineNumber_) +
                ": expected exactly four columns."
            );
        }

        // Tworzymy i walidujemy rekord.
        SpectrumEntry parsed{
            parseUnsigned<std::size_t>(
                sector,
                "sector_id",
                lineNumber_
            ),

            parseUnsigned<std::size_t>(
                seniority,
                "seniority",
                lineNumber_
            ),

            parseEnergy(
                energy,
                lineNumber_
            ),

            parseUnsigned<std::uint64_t>(
                degeneracy,
                "degeneracy",
                lineNumber_
            )
        };

        if (parsed.degeneracy == 0) {
            throw std::runtime_error(
                "SpectrumReader: line " +
                std::to_string(lineNumber_) +
                ": degeneracy must be positive."
            );
        }

        if (recordsRead_ ==
            std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error(
                "SpectrumReader: record counter overflow."
            );
        }

        entry = parsed;
        ++recordsRead_;

        return true;
    }

    return false;
}

Spectrum SpectrumReader::readAll()
{
    Spectrum spectrum;
    SpectrumEntry entry{};

    while (readNext(entry)) {
        spectrum.add(entry);
    }

    return spectrum;
}

std::size_t SpectrumReader::recordsRead() const noexcept
{
    return recordsRead_;
}

std::size_t SpectrumReader::lineNumber() const noexcept
{
    return lineNumber_;
}

Thermodynamics::Thermodynamics(
    const Spectrum& spectrum,
    double boltzmannConstant
)
    : spectrum_(spectrum),
    groundEnergy_(0.0L),
    boltzmannConstant_(boltzmannConstant)
{
    if (spectrum_.empty()) {
        throw std::invalid_argument(
            "Thermodynamics: spectrum cannot be empty."
        );
    }

    if (!std::isfinite(boltzmannConstant_) ||
        boltzmannConstant_ <= 0.0) {
        throw std::invalid_argument(
            "Thermodynamics: Boltzmann constant must be "
            "finite and positive."
        );
    }

    groundEnergy_ = static_cast<long double>(
        spectrum_.groundStateEnergy()
        );
}


double Thermodynamics::boltzmannConstant() const noexcept
{
    return boltzmannConstant_;
}

ThermodynamicPoint Thermodynamics::evaluate(
    double temperature
) const
{
    // ----------------------------------------------------
    // 1. Temperatura i skala termiczna
    // ----------------------------------------------------

    if (!std::isfinite(temperature) ||
        temperature <= 0.0) {
        throw std::invalid_argument(
            "Thermodynamics: temperature must be "
            "finite and positive."
        );
    }

    const long double kB =
        static_cast<long double>(boltzmannConstant_);

    const long double thermalEnergy =
        kB * static_cast<long double>(temperature);

    if (!std::isfinite(thermalEnergy) ||
        thermalEnergy <= 0.0L) {
        throw std::overflow_error(
            "Thermodynamics: invalid thermal energy kB*T."
        );
    }

    // ----------------------------------------------------
    // 2. Akumulatory statystyczne
    // ----------------------------------------------------

    long double weightSum = 0.0L;

    // Ważona średnia bezwymiarowej energii wzbudzeń.
    long double meanX = 0.0L;

    // Nienormalizowany drugi moment centralny.
    long double moment2 = 0.0L;

    // ----------------------------------------------------
    // 3. Jedno przejście przez pełne widmo
    // ----------------------------------------------------

    for (const auto& level : spectrum_.entries()) {

        const long double deltaEnergy =
            static_cast<long double>(level.energy)
            - groundEnergy_;

        if (deltaEnergy < 0.0L) {
            throw std::logic_error(
                "Thermodynamics: spectrum changed "
                "after initialization."
            );
        }

        // x = (E_i - E_0) / (kB*T)
        const long double x =
            deltaEnergy / thermalEnergy;

        // Bardzo duża dodatnia energia wzbudzenia
        // daje zerowy wkład termiczny.
        if (!std::isfinite(x)) {
            continue;
        }

        // w_i = d_i * exp(-x_i)
        const long double weight =
            static_cast<long double>(level.degeneracy)
            * std::exp(-x);

        // Podprzepływ oznacza zaniedbywalny wkład.
        if (weight == 0.0L) {
            continue;
        }

        // ------------------------------------------------
        // 4. Stabilna aktualizacja średniej i wariancji
        // ------------------------------------------------

        const long double newWeightSum =
            weightSum + weight;

        const long double difference =
            x - meanX;

        meanX +=
            (weight / newWeightSum) * difference;

        moment2 +=
            weight * difference * (x - meanX);

        weightSum = newWeightSum;
    }

    if (!std::isfinite(weightSum) ||
        weightSum <= 0.0L) {
        throw std::runtime_error(
            "Thermodynamics: invalid partition sum."
        );
    }

    // ----------------------------------------------------
    // 5. Funkcja rozdziału
    // ----------------------------------------------------

    const long double logShiftedZ =
        std::log(weightSum);

    const long double logZ =
        logShiftedZ - groundEnergy_ / thermalEnergy;

    // ----------------------------------------------------
    // 6. Energia wewnętrzna
    // ----------------------------------------------------

    const long double internalEnergy =
        groundEnergy_ + thermalEnergy * meanX;

    // ----------------------------------------------------
    // 7. Pojemność cieplna
    // ----------------------------------------------------

    // Var(x) = <(x - <x>)^2>
    const long double varianceX =
        std::max(0.0L, moment2 / weightSum);

    const long double heatCapacity =
        kB * varianceX;

    // ----------------------------------------------------
    // 8. Entropia
    // ----------------------------------------------------

    const long double entropy =
        kB * (logShiftedZ + meanX);

    // ----------------------------------------------------
    // 9. Energia swobodna Helmholtza
    // ----------------------------------------------------

    const long double freeEnergy =
        groundEnergy_ - thermalEnergy * logShiftedZ;

    // ----------------------------------------------------
    // 10. Wynik
    // ----------------------------------------------------

    return ThermodynamicPoint{
        temperature,
        logZ,
        weightSum,
        internalEnergy,
        heatCapacity,
        entropy,
        freeEnergy
    };
}

long double ThermodynamicPoint::partitionFunction() const
{
    return std::exp(logPartitionFunction);
}

std::vector<double> logarithmicGrid(
    double minimum,
    double maximum,
    std::size_t points
)
{
    if (!std::isfinite(minimum) ||
        !std::isfinite(maximum) ||
        minimum <= 0.0 ||
        maximum < minimum ||
        points == 0) {
        throw std::invalid_argument(
            "logarithmicGrid: invalid arguments."
        );
    }

    if (points == 1) {
        return { minimum };
    }

    std::vector<double> grid;
    grid.reserve(points);

    const double logMin = std::log(minimum);
    const double logMax = std::log(maximum);

    for (std::size_t i = 0; i < points; ++i) {

        const double fraction =
            static_cast<double>(i) /
            static_cast<double>(points - 1);

        const double logT =
            (1.0 - fraction) * logMin
            + fraction * logMax;

        grid.push_back(std::exp(logT));
    }

    // Dokładne końce zakresu.
    grid.front() = minimum;
    grid.back() = maximum;

    return grid;
}

namespace {

    // Usuwa białe znaki z początku i końca tekstu.
    std::string trim(const std::string& text)
    {
        const auto first =
            text.find_first_not_of(" \t\r\n");

        if (first == std::string::npos)
            return {};

        const auto last =
            text.find_last_not_of(" \t\r\n");

        return text.substr(first, last - first + 1);
    }


    // Tworzy komunikat z numerem linii.
    std::runtime_error inputError(
        std::size_t line,
        const std::string& message
    )
    {
        return std::runtime_error(
            "Input line " + std::to_string(line) +
            ": " + message
        );
    }


    // Konwersja nieujemnej liczby całkowitej.
    template <typename UInt>
    UInt parseUInt(
        const std::string& text,
        std::size_t line
    )
    {
        UInt value{};

        const char* first = text.data();
        const char* last = first + text.size();

        const auto result =
            std::from_chars(first, last, value);

        if (result.ec != std::errc{} ||
            result.ptr != last) {

            throw inputError(
                line,
                "Invalid unsigned integer: " + text
            );
        }

        return value;
    }


    // Konwersja liczby zmiennoprzecinkowej.
    double parseReal(
        const std::string& text,
        std::size_t line
    )
    {
        std::istringstream stream(text);

        stream.imbue(std::locale::classic());

        double value{};
        std::string extra;

        if (!(stream >> value) ||
            (stream >> extra) ||
            !std::isfinite(value)) {

            throw inputError(
                line,
                "Invalid floating-point value: " + text
            );
        }

        return value;
    }


    // Konwersja wartości logicznej.
    bool parseBool(
        const std::string& text,
        std::size_t line
    )
    {
        if (text == "true" || text == "1")
            return true;

        if (text == "false" || text == "0")
            return false;

        throw inputError(
            line,
            "Boolean value must be true, false, 1 or 0."
        );
    }


    // Zapobiega wielokrotnemu definiowaniu parametru.
    template <typename T>
    void assignOnce(
        std::optional<T>& target,
        T value,
        const std::string& key,
        std::size_t line
    )
    {
        if (target.has_value()) {
            throw inputError(
                line,
                "Duplicate parameter: " + key
            );
        }

        target = std::move(value);
    }

} // namespace

RunConfiguration readInput(std::istream& input)
{
    // Parametry, których obecność będziemy sprawdzać.
    std::optional<double> g;
    std::optional<std::size_t> pairs;
    std::optional<bool> odd;

    std::optional<double> Tmin;
    std::optional<double> Tmax;
    std::optional<std::size_t> NT;
    std::optional<double> kB;

    std::vector<Level> levels;

    std::string line;
    std::size_t lineNumber = 0;

    // ==================================================
    // 1. Odczyt kolejnych linii
    // ==================================================

    while (std::getline(input, line)) {

        ++lineNumber;

        // Usuwanie komentarza, również na końcu linii.
        const auto comment = line.find('#');

        if (comment != std::string::npos)
            line.erase(comment);

        line = trim(line);

        if (line.empty())
            continue;

        // Oczekujemy składni: key = value.
        const auto equal = line.find('=');

        if (equal == std::string::npos) {
            throw inputError(
                lineNumber,
                "Expected 'key = value'."
            );
        }

        const std::string key =
            trim(line.substr(0, equal));

        const std::string value =
            trim(line.substr(equal + 1));

        if (key.empty() || value.empty()) {
            throw inputError(
                lineNumber,
                "Empty parameter name or value."
            );
        }

        // ==================================================
        // 2. Poziomy energetyczne
        // ==================================================

        if (key == "level") {

            std::istringstream tokens(value);

            std::string energyToken;
            std::string capacityToken;
            std::string extra;

            if (!(tokens >> energyToken >> capacityToken) ||
                (tokens >> extra)) {

                throw inputError(
                    lineNumber,
                    "Level requires exactly two values: "
                    "energy pairCapacity."
                );
            }

            const double energy =
                parseReal(energyToken, lineNumber);

            const std::uint32_t capacity =
                parseUInt<std::uint32_t>(
                    capacityToken,
                    lineNumber
                );

            levels.push_back({ energy, capacity });
        }

        // ==================================================
        // 3. Parametry fizyczne
        // ==================================================

        else if (key == "pairing_strength") {
            assignOnce(
                g,
                parseReal(value, lineNumber),
                key,
                lineNumber
            );
        }

        else if (key == "pair_count") {
            assignOnce(
                pairs,
                parseUInt<std::size_t>(
                    value, lineNumber
                ),
                key,
                lineNumber
            );
        }

        else if (key == "odd_parity") {
            assignOnce(
                odd,
                parseBool(value, lineNumber),
                key,
                lineNumber
            );
        }

        // ==================================================
        // 4. Parametry termodynamiczne
        // ==================================================

        else if (key == "temperature_min") {
            assignOnce(
                Tmin,
                parseReal(value, lineNumber),
                key,
                lineNumber
            );
        }

        else if (key == "temperature_max") {
            assignOnce(
                Tmax,
                parseReal(value, lineNumber),
                key,
                lineNumber
            );
        }

        else if (key == "temperature_points") {
            assignOnce(
                NT,
                parseUInt<std::size_t>(
                    value, lineNumber
                ),
                key,
                lineNumber
            );
        }

        else if (key == "boltzmann_constant") {
            assignOnce(
                kB,
                parseReal(value, lineNumber),
                key,
                lineNumber
            );
        }

        // ==================================================
        // 5. Nieznany parametr
        // ==================================================

        else {
            throw inputError(
                lineNumber,
                "Unknown parameter: " + key
            );
        }
    }

    if (input.bad() || !input.eof()) {
        throw std::runtime_error(
            "I/O error while reading input configuration."
        );
    }

    // ==================================================
    // 6. Sprawdzenie wymaganych parametrów
    // ==================================================

    const auto require = [](
        bool present,
        const char* key
        ) {
            if (!present) {
                throw std::runtime_error(
                    std::string("Missing input parameter: ") +
                    key
                );
            }
        };

    require(g.has_value(), "pairing_strength");
    require(pairs.has_value(), "pair_count");
    require(odd.has_value(), "odd_parity");

    require(!levels.empty(), "level");

    require(Tmin.has_value(), "temperature_min");
    require(Tmax.has_value(), "temperature_max");
    require(NT.has_value(), "temperature_points");
    require(kB.has_value(), "boltzmann_constant");

    // ==================================================
    // 7. Utworzenie konfiguracji
    // ==================================================

    RunConfiguration config{
        PairingModelParameters{
            std::move(levels),
            *g,
            *pairs,
            *odd
        },
        *Tmin,
        *Tmax,
        *NT,
        *kB
    };

    // Walidacja istniejącej warstwy fizycznej.
    validate(config.model);

    // Walidacja termodynamiki.
    if (config.temperatureMin <= 0.0 ||
        config.temperatureMax < config.temperatureMin) {

        throw std::invalid_argument(
            "Invalid temperature range."
        );
    }

    if (config.temperaturePoints < 2) {
        throw std::invalid_argument(
            "temperature_points must be at least 2."
        );
    }

    if (config.boltzmannConstant <= 0.0) {
        throw std::invalid_argument(
            "boltzmann_constant must be positive."
        );
    }

    return config;
}