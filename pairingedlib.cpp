#include "pairingedlib.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
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