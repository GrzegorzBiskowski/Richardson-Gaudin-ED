#include "pairingedlib.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

namespace {
    void generateRecursive(
        const std::vector<std::uint32_t>& upperBounds,
        const std::vector<std::size_t>& suffixCapacity,
        std::size_t level,
        std::size_t remaining,
        std::vector<std::uint32_t>& current,
        std::vector<OccupationState>& result
    ) {
        const std::size_t numberOfLevels = upperBounds.size();

        // Doszliśmy do końca wszystkich poziomów.
        if (level == numberOfLevels) {
            if (remaining == 0) {
                result.emplace_back(current);
            }
            return;
        }

        // Ile maksymalnie możemy umieścić na obecnym poziomie?
        const std::size_t maxOccupation =
            std::min<std::size_t>(
                upperBounds[level],
                remaining
            );

        for (std::size_t occupation = 0;
            occupation <= maxOccupation;
            ++occupation)
        {
            // Po obsadzeniu obecnego poziomu musi być jeszcze
            // możliwe rozłożenie pozostałych par na kolejnych poziomach.
            const std::size_t newRemaining = remaining - occupation;

            if (newRemaining > suffixCapacity[level + 1]) {
                continue;
            }

            current[level] =
                static_cast<std::uint32_t>(occupation);

            generateRecursive(
                upperBounds,
                suffixCapacity,
                level + 1,
                newRemaining,
                current,
                result
            );
        }
    }

} // namespace


std::vector<OccupationState> generateBoundedOccupations(const std::vector<std::uint32_t>& upperBounds, std::size_t totalOccupation)
{
    std::vector<OccupationState> result;

    // Nie ma poziomów.
    if (upperBounds.empty()) {
        if (totalOccupation == 0) {
            result.emplace_back(
                std::vector<std::uint32_t>{}
            );
        }

        return result;
    }

    // suffixCapacity[i] = maksymalna liczba par,
    // którą można jeszcze umieścić na poziomach i, i+1, ...
    std::vector<std::size_t> suffixCapacity(
        upperBounds.size() + 1,
        0
    );

    for (std::size_t i = upperBounds.size(); i > 0; --i) {
        suffixCapacity[i - 1] =
            suffixCapacity[i] +
            static_cast<std::size_t>(upperBounds[i - 1]);
    }

    // Za dużo par na dostępną przestrzeń Hilberta.
    if (totalOccupation > suffixCapacity[0]) {
        return result;
    }

    std::vector<std::uint32_t> current(
        upperBounds.size(),
        0
    );

    generateRecursive(
        upperBounds,
        suffixCapacity,
        0,
        totalOccupation,
        current,
        result
    );

    return result;
}

std::size_t countBoundedOccupations(
    const std::vector<std::uint32_t>& upperBounds,
    std::size_t totalOccupation
) {
    // Nie ma poziomów.
    if (upperBounds.empty()) {
        return (totalOccupation == 0) ? 1 : 0;
    }

    // Obliczamy maksymalną możliwą liczbę par.
    std::size_t maximumOccupation = 0;

    for (std::uint32_t bound : upperBounds) {
        const std::size_t boundSize =
            static_cast<std::size_t>(bound);

        if (maximumOccupation >
            std::numeric_limits<std::size_t>::max() - boundSize)
        {
            throw std::overflow_error(
                "Overflow while calculating total pair capacity."
            );
        }

        maximumOccupation += boundSize;
    }

    // Nie da się rozmieścić tylu par.
    if (totalOccupation > maximumOccupation) {
        return 0;
    }

    // dp[s] = liczba sposobów rozmieszczenia s par
    // na poziomach przetworzonych do tej pory.
    std::vector<std::size_t> dp(totalOccupation + 1, 0);
    std::vector<std::size_t> newDp(totalOccupation + 1, 0);

    // Na początku mamy zero poziomów:
    // istnieje dokładnie jeden sposób otrzymania 0 par.
    dp[0] = 1;

    for (std::uint32_t bound32 : upperBounds) {

        const std::size_t bound =
            static_cast<std::size_t>(bound32);

        std::fill(newDp.begin(), newDp.end(), 0);

        std::size_t window = 0;

        for (std::size_t s = 0;
            s <= totalOccupation;
            ++s)
        {
            // Dodajemy dp[s] do bieżącego okna.
            if (window >
                std::numeric_limits<std::size_t>::max() - dp[s])
            {
                throw std::overflow_error(
                    "Overflow while counting bounded occupations."
                );
            }

            window += dp[s];

            // Okno reprezentuje:
            //
            // dp[s] + dp[s-1] + ... + dp[s-bound]
            //
            // więc jeżeli przekroczyliśmy szerokość okna,
            // usuwamy dp[s-bound-1].
            if (s > bound) {
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

const std::vector<std::uint32_t>&
BlockedConfiguration::blockedElectrons() const noexcept
{
    return blockedElectrons_;
}

std::size_t BlockedConfiguration::seniority() const noexcept
{
    return seniority_;
}

SenioritySector makeSenioritySector(
    std::size_t id,
    BlockedConfiguration blocked,
    const PairingModelParameters& model
) {
    const auto& levels = model.levels;
    const auto& v = blocked.blockedElectrons();

    if (v.size() != levels.size()) {
        throw std::invalid_argument(
            "Blocked configuration size does not match "
            "the number of energy levels."
        );
    }

    // Liczba elektronów w całym układzie.
    const std::size_t electronCount = checkedAdd(
        checkedMultiply(model.pairCount, 2),
        model.oddParity ? 1 : 0
    );

    const std::size_t nu = blocked.seniority();

    if (nu > electronCount) {
        throw std::invalid_argument(
            "Seniority exceeds the total electron count."
        );
    }

    if ((electronCount - nu) % 2 != 0) {
        throw std::invalid_argument(
            "Invalid seniority parity for this electron number."
        );
    }

    const std::size_t numberOfPairs =
        (electronCount - nu) / 2;

    std::vector<std::uint32_t> effectiveCapacity;
    effectiveCapacity.reserve(levels.size());

    std::size_t totalEffectiveCapacity = 0;
    long double blockedEnergySum = 0.0L;

    for (std::size_t j = 0; j < levels.size(); ++j) {
        const std::uint32_t omega = levels[j].pairCapacity;

        if (v[j] > omega) {
            throw std::invalid_argument(
                "Blocked electron count exceeds pair capacity "
                "at level " + std::to_string(j) + "."
            );
        }

        const std::uint32_t effectiveOmega = omega - v[j];

        effectiveCapacity.push_back(effectiveOmega);

        totalEffectiveCapacity = checkedAdd(
            totalEffectiveCapacity,
            static_cast<std::size_t>(effectiveOmega)
        );

        blockedEnergySum +=
            static_cast<long double>(levels[j].energy) *
            static_cast<long double>(v[j]);
    }

    if (numberOfPairs > totalEffectiveCapacity) {
        throw std::invalid_argument(
            "The sector has more pairs than its effective "
            "pair capacity allows."
        );
    }

    const double blockedEnergy =
        static_cast<double>(blockedEnergySum);

    if (!std::isfinite(blockedEnergySum) ||
        !std::isfinite(blockedEnergy)) {
        throw std::overflow_error(
            "Blocked-electron energy is not finite."
        );
    }

    const std::uint64_t degeneracy =
        calculateSectorDegeneracy(levels, blocked);

    return SenioritySector{
        id,
        std::move(blocked),
        numberOfPairs,
        std::move(effectiveCapacity),
        blockedEnergy,
        degeneracy
    };
}

namespace {

    std::uint64_t checkedMultiplyU64(
        std::uint64_t a,
        std::uint64_t b
    ) {
        if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a) {
            throw std::overflow_error(
                "Overflow while calculating sector degeneracy."
            );
        }

        return a * b;
    }


    // Oblicza symbol Newtona n choose k.
    // Rzuca overflow_error, jeśli wynik nie mieści się
    // w uint64_t.
    std::uint64_t binomialCoefficientChecked(
        std::uint64_t n,
        std::uint64_t k
    ) {
        if (k > n) {
            return 0;
        }

        // C(n,k) = C(n,n-k).
        k = (k < n - k) ? k : n - k;

        std::uint64_t result = 1;

        for (std::uint64_t i = 1; i <= k; ++i) {
            std::uint64_t numerator = n - k + i;
            std::uint64_t denominator = i;

            // Skracamy licznik i mianownik.
            const auto g1 = std::gcd(
                numerator,
                denominator
            );

            numerator /= g1;
            denominator /= g1;

            // Skracamy pozostały mianownik z dotychczasowym
            // wynikiem. Po tym kroku denominator powinien
            // wynosić 1, ponieważ współczynnik dwumianowy
            // jest liczbą całkowitą.
            const auto g2 = std::gcd(
                result,
                denominator
            );

            result /= g2;
            denominator /= g2;

            if (denominator != 1) {
                throw std::logic_error(
                    "Unexpected non-integer intermediate result "
                    "in binomialCoefficientChecked()."
                );
            }

            result = checkedMultiplyU64(
                result,
                numerator
            );
        }

        return result;
    }

} // namespace

std::uint64_t calculateSectorDegeneracy(
    const std::vector<Level>& levels,
    const BlockedConfiguration& blocked
) {
    const auto& v = blocked.blockedElectrons();

    if (levels.size() != v.size()) {
        throw std::invalid_argument(
            "Cannot calculate sector degeneracy: "
            "the number of levels does not match "
            "the blocked configuration size."
        );
    }

    std::uint64_t totalDegeneracy = 1;

    for (std::size_t j = 0; j < levels.size(); ++j) {
        const std::uint32_t omega =
            levels[j].pairCapacity;

        const std::uint32_t seniority =
            v[j];

        if (omega == 0) {
            throw std::invalid_argument(
                "Cannot calculate sector degeneracy: "
                "pair capacity is zero at level "
                + std::to_string(j) + "."
            );
        }

        if (seniority > omega) {
            throw std::invalid_argument(
                "Cannot calculate sector degeneracy: "
                "partial seniority exceeds pair capacity "
                "at level " + std::to_string(j) + "."
            );
        }

        // Poziom ma 2*Omega stanów jednocząstkowych.
        const std::uint64_t singleParticleStates =
            checkedMultiplyU64(
                2,
                static_cast<std::uint64_t>(omega)
            );

        const std::uint64_t firstTerm =
            binomialCoefficientChecked(
                singleParticleStates,
                seniority
            );

        std::uint64_t secondTerm = 0;

        if (seniority >= 2) {
            secondTerm = binomialCoefficientChecked(
                singleParticleStates,
                seniority - 2
            );
        }

        if (secondTerm > firstTerm) {
            throw std::logic_error(
                "Invalid local degeneracy formula result "
                "at level " + std::to_string(j) + "."
            );
        }

        const std::uint64_t localDegeneracy =
            firstTerm - secondTerm;

        totalDegeneracy = checkedMultiplyU64(
            totalDegeneracy,
            localDegeneracy
        );
    }

    return totalDegeneracy;
}

SpectrumWriter::SpectrumWriter(std::ostream& output)
    : output_(output)
{
    if (!output_) {
        throw std::runtime_error(
            "Cannot initialize SpectrumWriter: "
            "output stream is not writable."
        );
    }

    // Zachowanie wystarczającej precyzji do odtworzenia double.
    output_ << std::defaultfloat
        << std::setprecision(
            std::numeric_limits<double>::max_digits10
        );

    output_ << "# pairingedlib spectrum v1\n"
        << "# columns: sector_id seniority energy degeneracy\n";

    if (!output_) {
        throw std::runtime_error(
            "Failed to write spectrum file header."
        );
    }
}


void SpectrumWriter::write(const SpectrumEntry& entry)
{
    if (!std::isfinite(entry.energy)) {
        throw std::invalid_argument(
            "Cannot write spectrum entry: energy is not finite."
        );
    }

    if (entry.degeneracy == 0) {
        throw std::invalid_argument(
            "Cannot write spectrum entry: degeneracy must be positive."
        );
    }

    output_ << entry.sectorId << ' '
        << entry.seniority << ' '
        << entry.energy << ' '
        << entry.degeneracy << '\n';

    if (!output_) {
        throw std::runtime_error(
            "Failed while writing a spectrum entry."
        );
    }

    // Przepełnienie licznika nie powinno być ignorowane.
    if (recordsWritten_ ==
        std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error(
            "SpectrumWriter record counter overflow."
        );
    }

    ++recordsWritten_;
}


void SpectrumWriter::flush()
{
    output_.flush();

    if (!output_) {
        throw std::runtime_error(
            "Failed to flush spectrum output stream."
        );
    }
}


std::size_t SpectrumWriter::recordsWritten() const noexcept
{
    return recordsWritten_;
}

void validate(const PairingModelParameters& model)
{
    if (model.levels.empty()) {
        throw std::invalid_argument(
            "Invalid model: at least one energy level is required."
        );
    }

    if (!std::isfinite(model.pairingStrength)) {
        throw std::invalid_argument(
            "Invalid model: pairingStrength must be finite."
        );
    }

    std::size_t totalPairCapacity = 0;

    for (std::size_t j = 0; j < model.levels.size(); ++j) {
        const auto& level = model.levels[j];

        if (!std::isfinite(level.energy)) {
            throw std::invalid_argument(
                "Invalid model: non-finite energy at level "
                + std::to_string(j) + "."
            );
        }

        if (level.pairCapacity == 0) {
            throw std::invalid_argument(
                "Invalid model: pairCapacity must be positive "
                "at level " + std::to_string(j) + "."
            );
        }

        totalPairCapacity = checkedAdd(
            totalPairCapacity,
            static_cast<std::size_t>(level.pairCapacity)
        );
    }

    // N_e = 2M + b
    const std::size_t electronCount = checkedAdd(
        checkedMultiply(model.pairCount, 2),
        model.oddParity ? 1 : 0
    );

    // Maksymalna liczba elektronów: 2 * sum_j Omega_j.
    const std::size_t maximumElectronCount =
        checkedMultiply(totalPairCapacity, 2);

    if (electronCount > maximumElectronCount) {
        throw std::invalid_argument(
            "Invalid model: the number of electrons exceeds "
            "the total capacity of the energy levels."
        );
    }
}

std::vector<BlockedConfiguration> generateBlockedConfigurations(
    const std::vector<std::uint32_t>& pairCapacity,
    std::size_t totalSeniority
) {
    const auto states = generateBoundedOccupations(
        pairCapacity,
        totalSeniority
    );

    std::vector<BlockedConfiguration> result;
    result.reserve(states.size());

    for (const auto& state : states) {
        result.emplace_back(state.occupations());
    }

    return result;
}

SectorGenerator::SectorGenerator(
    PairingModelParameters parameters
)
    : parameters_(std::move(parameters))
{
    validate(parameters_);
}

namespace {

    SenioritySector buildSector(
        std::size_t id,
        BlockedConfiguration blocked,
        const PairingModelParameters& model,
        std::size_t electronCount
    ) {
        const auto& levels = model.levels;
        const auto& v = blocked.blockedElectrons();

        if (v.size() != levels.size()) {
            throw std::logic_error(
                "Internal error: blocked configuration has "
                "incorrect size."
            );
        }

        const std::size_t nu = blocked.seniority();

        if (nu > electronCount ||
            (electronCount - nu) % 2 != 0) {
            throw std::logic_error(
                "Internal error: incompatible seniority."
            );
        }

        const std::size_t numberOfPairs =
            (electronCount - nu) / 2;

        std::vector<std::uint32_t> effectiveCapacity;
        effectiveCapacity.reserve(levels.size());

        std::size_t totalEffectiveCapacity = 0;
        long double blockedEnergySum = 0.0L;

        for (std::size_t j = 0; j < levels.size(); ++j) {
            const auto omega = levels[j].pairCapacity;

            if (v[j] > omega) {
                throw std::logic_error(
                    "Internal error: blocked electron count "
                    "exceeds pair capacity."
                );
            }

            const std::uint32_t effectiveOmega = omega - v[j];

            effectiveCapacity.push_back(effectiveOmega);

            totalEffectiveCapacity = checkedAdd(
                totalEffectiveCapacity,
                static_cast<std::size_t>(effectiveOmega)
            );

            blockedEnergySum +=
                static_cast<long double>(levels[j].energy) *
                static_cast<long double>(v[j]);
        }

        if (numberOfPairs > totalEffectiveCapacity) {
            throw std::logic_error(
                "Internal error: sector has more pairs "
                "than its effective capacity allows."
            );
        }

        const double blockedEnergy =
            static_cast<double>(blockedEnergySum);

        if (!std::isfinite(blockedEnergySum) ||
            !std::isfinite(blockedEnergy)) {
            throw std::overflow_error(
                "Blocked-electron energy is not finite."
            );
        }

        // Oblicz degenerację, zanim przeniesiemy konfigurację.
        const std::uint64_t degeneracy =
            calculateSectorDegeneracy(levels, blocked);

        return SenioritySector{
            id,
            std::move(blocked),
            numberOfPairs,
            std::move(effectiveCapacity),
            blockedEnergy,
            degeneracy
        };
    }

} // namespace

std::vector<SenioritySector> SectorGenerator::generate() const
{
    const auto& levels = parameters_.levels;

    std::vector<std::uint32_t> pairCapacity;
    pairCapacity.reserve(levels.size());

    std::size_t totalPairCapacity = 0;

    for (const auto& level : levels) {
        pairCapacity.push_back(level.pairCapacity);

        totalPairCapacity = checkedAdd(
            totalPairCapacity,
            static_cast<std::size_t>(level.pairCapacity)
        );
    }

    const std::size_t electronCount = checkedAdd(
        checkedMultiply(parameters_.pairCount, 2),
        parameters_.oddParity ? 1 : 0
    );

    const std::size_t maximumElectronCount =
        checkedMultiply(totalPairCapacity, 2);

    const std::size_t maximumSeniority = std::min(
        electronCount,
        maximumElectronCount - electronCount
    );

    const std::size_t firstSeniority =
        parameters_.oddParity ? 1 : 0;

    std::vector<SenioritySector> sectors;
    std::size_t nextId = 0;

    for (std::size_t nu = firstSeniority;
        nu <= maximumSeniority; )
    {
        auto configurations = generateBlockedConfigurations(
            pairCapacity,
            nu
        );

        for (auto& blocked : configurations) {
            if (nextId ==
                std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error(
                    "Too many seniority sectors."
                );
            }

            sectors.push_back(
                buildSector(
                    nextId,
                    std::move(blocked),
                    parameters_,
                    electronCount
                )
            );

            ++nextId;
        }

        // Nie dodajemy 2, jeśli wyszlibyśmy poza zakres.
        if (maximumSeniority - nu < 2) {
            break;
        }

        nu = checkedAdd(nu, 2);
    }

    return sectors;
}