#include "pairingedlib.hpp"
#include "tests.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <openblas/cblas.h>
#include <stdexcept>
#include <string>
#include <vector>


int main(int argc, char* argv[])
{
    /*levels_test();
    occupation_test();
    basis_test();
    generateBoundedOccupations_test();
    countBoundedOccupations_test();
    blockedConfiguration_test();
    calculateSectorDegeneracy_test();
    modelGenerator_test();
    pairingHamiltonian_test();*/

    try {

        // ============================================================
        // 1. DEFINICJA UKŁADU FIZYCZNEGO
        // ============================================================

        const std::string inputFilename =
            (argc >= 2) ? argv[1] : "input.dat";

        std::ifstream inputFile(inputFilename);

        if (!inputFile) {
            throw std::runtime_error(
                "Cannot open input file: " + inputFilename
            );
        }

        // Odczyt wszystkich parametrów.
        RunConfiguration config = readInput(inputFile);

        // Istniejący kod może nadal korzystać
        // z model dokładnie tak samo jak wcześniej.
        const PairingModelParameters& model = config.model;

        std::cout << "Input file: "
            << inputFilename << '\n';

        std::cout << "Energy levels: "
            << model.levels.size() << '\n';

        std::cout << "Pairing strength: "
            << model.pairingStrength << '\n';

        std::cout << "Initial pair count: "
            << model.pairCount << '\n';

        const std::string outputFilename = "spectrum.dat";

        // Kontrole diagnostyczne.
        const bool verifyTrace = true;
        const bool verifyOutputFile = true;

        // ============================================================
        // 2. WALIDACJA PARAMETRÓW
        // ============================================================

        validate(model);

        std::cout << "====================================\n";
        std::cout << "       PAIRING MODEL SOLVER\n";
        std::cout << "====================================\n\n";

        std::cout << "Energy levels: "
            << model.levels.size() << '\n';

        std::cout << "Pairing strength: "
            << model.pairingStrength << '\n';

        std::cout << "Initial pair count: "
            << model.pairCount << '\n';

        std::cout << "Particle parity: "
            << (model.oddParity ? "odd" : "even")
            << "\n\n";

        // ============================================================
        // 3. GENEROWANIE SEKTORÓW SENIORITY
        // ============================================================

        std::cout << "Generating seniority sectors...\n";

        SectorGenerator generator(model);

        const std::vector<SenioritySector> sectors =
            generator.generate();

        std::cout << "Generated sectors: "
            << sectors.size() << "\n\n";

        // ============================================================
        // 4. ESTYMACJA ZAPOTRZEBOWANIA NA RAM
        // ============================================================

        MemoryEstimateOptions memoryOptions;

        const MemoryEstimate memory =
            estimateMemory(model, sectors, memoryOptions);

        printMemoryEstimate(memory, std::cout);

        if (!memory.compatibleWithLapackInt) {
            throw std::runtime_error(
                "The calculation exceeds the LP64 LAPACK "
                "integer limits."
            );
        }

        std::cout << "\nContinue calculations? [y/n]: ";

        char answer = 'n';
        std::cin >> answer;

        if (answer != 'y' && answer != 'Y') {
            std::cout << "Calculation cancelled.\n";
            return 0;
        }

        // ============================================================
        // 5. INICJALIZACJA WIDMA I DIAGONALIZATORA
        // ============================================================

        Spectrum spectrum;

        spectrum.reserve(memory.totalEigenvalues);

        SymmetricDiagonalizer diagonalizer(
            memoryOptions.lapackWorkspaceDoublesPerDimension
        );

        std::size_t completedSectors = 0;

        // ============================================================
        // 6. OTWARCIE PLIKU WYNIKOWEGO
        // ============================================================

        {
            std::ofstream spectrumFile(
                outputFilename,
                std::ios::out | std::ios::trunc
            );

            if (!spectrumFile) {
                throw std::runtime_error(
                    "Cannot open spectrum.dat for writing."
                );
            }

            SpectrumWriter writer(spectrumFile);

            std::cout << "\nStarting diagonalization...\n";

            // ========================================================
            // 7. OBLICZENIA DLA KOLEJNYCH SEKTORÓW
            // ========================================================

            for (const auto& sector : sectors) {

                // ----------------------------------------------------
                // 7.1. Generowanie bazy sektora
                // ----------------------------------------------------

                Basis basis = generateSectorBasis(sector);

                // ----------------------------------------------------
                // 7.2. Budowa Hamiltonianu
                // ----------------------------------------------------

                PairingHamiltonian hamiltonian(
                    model,
                    sector,
                    basis
                );

                DenseMatrix H = hamiltonian.buildMatrix();

                // ----------------------------------------------------
                // 7.3. Kontrola śladu przed diagonalizacją
                // ----------------------------------------------------

                long double originalTrace = 0.0L;
                long double traceScale = 0.0L;

                if (verifyTrace) {
                    for (std::size_t i = 0; i < H.rows(); ++i) {
                        const long double value = H(i, i);

                        originalTrace += value;
                        traceScale += std::abs(value);
                    }
                }

                // ----------------------------------------------------
                // 7.4. Diagonalizacja
                // ----------------------------------------------------

                std::vector<double> eigenvalues =
                    diagonalizer.diagonalize(H);

                if (eigenvalues.size() != basis.size()) {
                    throw std::runtime_error(
                        "Unexpected eigenvalue count in sector " +
                        std::to_string(sector.id)
                    );
                }

                // ----------------------------------------------------
                // 7.5. Test niezmienniczości śladu
                // ----------------------------------------------------

                if (verifyTrace) {

                    long double eigenvalueSum = 0.0L;

                    for (double E : eigenvalues) {
                        eigenvalueSum += E;
                    }

                    const long double difference =
                        std::abs(originalTrace - eigenvalueSum);

                    const long double tolerance =
                        1.0e-10L * (1.0L + traceScale);

                    if (difference > tolerance) {
                        throw std::runtime_error(
                            "Trace verification failed in sector " +
                            std::to_string(sector.id)
                        );
                    }
                }

                // ----------------------------------------------------
                // 7.6. Tworzenie pełnych rekordów widma
                // ----------------------------------------------------

                for (double E : eigenvalues) {

                    SpectrumEntry entry{
                        sector.id,
                        sector.seniority(),

                        // Pełna energia fizyczna:
                        E + sector.blockedEnergy,

                        sector.degeneracy
                    };

                    // Widmo w pamięci RAM.
                    spectrum.add(entry);

                    // Zapis do pliku.
                    writer.write(entry);
                }

                // ----------------------------------------------------
                // 7.7. Zakończenie sektora
                // ----------------------------------------------------

                writer.flush();

                ++completedSectors;

                // Raport postępu co 10 sektorów.
                if (completedSectors % 10 == 0 ||
                    completedSectors == sectors.size()) {

                    std::cout
                        << "Completed "
                        << completedSectors
                        << " / "
                        << sectors.size()
                        << " sectors; eigenvalues: "
                        << spectrum.size()
                        << '\n';
                }

                // Basis, H oraz eigenvalues zostaną
                // zwolnione po zakończeniu iteracji.
            }

            // ========================================================
            // 8. WERYFIKACJA LICZBY ZAPISANYCH REKORDÓW
            // ========================================================

            if (spectrum.size() != memory.totalEigenvalues) {
                throw std::runtime_error(
                    "Incomplete spectrum in memory."
                );
            }

            if (writer.recordsWritten() !=
                memory.totalEigenvalues) {
                throw std::runtime_error(
                    "Incorrect number of written records."
                );
            }

            writer.flush();

            spectrumFile.close();

            if (!spectrumFile) {
                throw std::runtime_error(
                    "Error while closing spectrum output file."
                );
            }
        }

        std::cout << "\nDiagonalization completed.\n";

        // ============================================================
        // 9. KONTROLA PLIKU PRZEZ SPECTRUMREADER
        // ============================================================

        if (verifyOutputFile) {

            std::cout << "Verifying spectrum file...\n";

            std::ifstream inputFile(outputFilename);

            if (!inputFile) {
                throw std::runtime_error(
                    "Cannot reopen spectrum.dat for verification."
                );
            }

            SpectrumReader reader(inputFile);

            SpectrumEntry loaded{};

            // Porównujemy plik z oryginalnym widmem
            // przed jego sortowaniem.
            for (const auto& original : spectrum.entries()) {

                if (!reader.readNext(loaded)) {
                    throw std::runtime_error(
                        "Spectrum file ended prematurely."
                    );
                }

                if (loaded.sectorId != original.sectorId ||
                    loaded.seniority != original.seniority ||
                    loaded.energy != original.energy ||
                    loaded.degeneracy != original.degeneracy) {

                    throw std::runtime_error(
                        "Spectrum file verification failed "
                        "at record " +
                        std::to_string(reader.recordsRead())
                    );
                }
            }

            // Sprawdzamy, czy plik nie zawiera nadmiarowych danych.
            if (reader.readNext(loaded)) {
                throw std::runtime_error(
                    "Spectrum file contains extra records."
                );
            }

            std::cout << "Spectrum file verified successfully.\n";
        }

        // ============================================================
        // 10. GLOBALNE SORTOWANIE WIDMA
        // ============================================================

        spectrum.sortByEnergy();

        // ============================================================
        // 11. RAPORT KOŃCOWY
        // ============================================================

        std::cout << "\n====================================\n";
        std::cout << "        CALCULATION SUMMARY\n";
        std::cout << "====================================\n\n";

        std::cout << "Completed sectors: "
            << completedSectors << '\n';

        std::cout << "Calculated eigenvalues: "
            << spectrum.size() << '\n';

        std::cout << "Physical state count: "
            << spectrum.physicalStateCount() << '\n';

        std::cout << std::defaultfloat
            << std::setprecision(17);

        std::cout << "Ground-state energy: "
            << spectrum.groundStateEnergy() << '\n';

        std::cout << "Output file: "
            << outputFilename << '\n';

        // ============================================================
        // 12. WYŚWIETLENIE PIERWSZYCH WARTOŚCI WŁASNYCH
        // ============================================================

        const std::size_t displayCount =
            std::min<std::size_t>(10, spectrum.size());

        std::cout << "\nLowest energy levels:\n";
        std::cout << "Energy\tSeniority\tDegeneracy\tSector\n";

        const auto& entries = spectrum.entries();

        for (std::size_t i = 0; i < displayCount; ++i) {

            const auto& entry = entries[i];

            std::cout
                << entry.energy << '\t'
                << entry.seniority << '\t'
                << entry.degeneracy << '\t'
                << entry.sectorId << '\n';
        }

        std::cout << "\nCalculation finished successfully.\n";

        std::cout << "\nCalculating thermodynamics...\n";

        // Stała Boltzmanna w jednostkach używanych przez model.

        Thermodynamics thermo(
            spectrum,
            config.boltzmannConstant
        );

        const auto temperatures = logarithmicGrid(
            config.temperatureMin,
            config.temperatureMax,
            config.temperaturePoints
        );

        // Plik wynikowy.
        std::ofstream thermoFile("thermodynamics.dat");

        if (!thermoFile) {
            throw std::runtime_error(
                "Cannot open thermodynamics.dat."
            );
        }

        thermoFile.imbue(std::locale::classic());

        thermoFile << std::scientific
            << std::setprecision(
                std::numeric_limits<long double>::max_digits10
            );

        thermoFile
            << "# Thermodynamics from pairingedlib\n"
            << "# T logZ Z Zshift U Cv S F\n";

        // Obliczenia.
        for (double T : temperatures) {

            const ThermodynamicPoint result =
                thermo.evaluate(T);

            thermoFile
                << result.temperature << '\t'
                << result.logPartitionFunction << '\t'
                << result.partitionFunction() << '\t'
                << result.shiftedPartitionFunction << '\t'
                << result.internalEnergy << '\t'
                << result.heatCapacity << '\t'
                << result.entropy << '\t'
                << result.freeEnergy << '\n';
        }

        thermoFile.close();

        if (!thermoFile) {
            throw std::runtime_error(
                "Failed while writing thermodynamics.dat."
            );
        }

        std::cout << "Thermodynamics calculated for "
            << temperatures.size()
            << " temperatures.\n";

        std::cout << "Output: thermodynamics.dat\n";

        return 0;
    }
    catch (const std::exception& e) {

        std::cerr
            << "\nERROR: "
            << e.what()
            << '\n';

        return 1;
    }
}
