#include "pairingedlib.hpp"
#include "tests.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <optional>
#include <cblas.h>
#include <stdexcept>
#include <string>
#include <vector>

Spectrum calculateSpectrum(
    const PairingModelParameters& model,
    const std::vector<SenioritySector>& sectors)
{
    Spectrum spectrum;
    SymmetricDiagonalizer diagonalizer;

    for (const auto& sector : sectors)
    {
        // Generowanie bazy danego sektora.
        Basis basis = generateSectorBasis(sector);

        // Konstrukcja Hamiltonianu.
        PairingHamiltonian hamiltonian(
            model, sector, basis
        );

        DenseMatrix matrix = hamiltonian.buildMatrix();

        // Diagonalizacja.
        const auto eigenvalues =
            diagonalizer.diagonalize(matrix);

        // Zapis wartości własnych w pamięci.
        for (const double eigenvalue : eigenvalues)
        {
            spectrum.add({
                sector.id,
                sector.seniority(),
                eigenvalue + sector.blockedEnergy,
                sector.degeneracy
            });
        }
    }

    spectrum.sortByEnergy();

    return spectrum;
}

void scanHeatCapacityPeaks(
    const RunConfiguration& config,
    const std::vector<SenioritySector>& sectors)
    {

    if (!config.couplingScan.has_value())
    {
        throw std::invalid_argument(
            "Coupling scan is not configured."
        );
    }

    const auto& scan = *config.couplingScan;

    const double gMin = config.model.pairingStrength;
    const double gMax = scan.endStrength;

    const auto temperatures = logarithmicGrid(
        config.temperatureMin,
        config.temperatureMax,
        config.temperaturePoints
    );

    const std::string filename =
        config.seniorityZeroOnly
        ? "peaks_vs_g_zero.dat"
        : "peaks_vs_g_all.dat";

    std::ofstream output(filename);

    if (!output)
    {
        throw std::runtime_error(
            "Cannot open peak scan output file."
        );
    }

    output.imbue(std::locale::classic());

    output << std::scientific
           << std::setprecision(
                  std::numeric_limits<long double>::max_digits10
              );

    output << "# g peak_index T_peak C_V_peak\n";

    struct PeakScanRow
    {
        double g;
        std::vector<HeatCapacityPeak> peaks;
    };

    std::vector<PeakScanRow> scanRows;
    scanRows.reserve(scan.points);

    std::size_t maxPeakCount = 0;

    for (std::size_t i = 0; i < scan.points; ++i)
    {
        const double alpha =
            static_cast<double>(i) /
            static_cast<double>(scan.points - 1);

        const double g = (i + 1 == scan.points)
            ? gMax
            : gMin + alpha * (gMax - gMin);

        // Kopia parametrow modelu dla nowego g.
        PairingModelParameters model = config.model;
        model.pairingStrength = g;

        // NOWA diagonalizacja dla kazdego g.
        Spectrum spectrum =
            calculateSpectrum(model, sectors);

        Thermodynamics thermo(
            spectrum,
            config.boltzmannConstant
        );

        // Obliczenie C_V na siatce temperatur.
        std::vector<ThermodynamicPoint> samples;
        samples.reserve(temperatures.size());

        for (const double T : temperatures)
        {
            samples.push_back(thermo.evaluate(T));
        }

        // Nasza istniejaca funkcja znajdowania maksimow.
        const auto peaks =
            findHeatCapacityPeaks(thermo, samples);
        scanRows.push_back({g, peaks});

        maxPeakCount = std::max(
            maxPeakCount,
            peaks.size()
        );

        // Zapis wszystkich wykrytych maksimow.
        for (std::size_t j = 0; j < peaks.size(); ++j)
        {
            output << g << ' '
                   << (j + 1) << ' '
                   << peaks[j].temperature << ' '
                   << peaks[j].heatCapacity
                   << '\n';
        }

        // Pusta linia oddzielajaca kolejne g.
        output << '\n';

        std::cout
            << "[" << (i + 1) << "/"
            << scan.points << "] "
            << "g = " << g
            << ", peaks = " << peaks.size()
            << '\n';
    }

    const auto writeWideFile =
    [&](const std::string& filename,
        bool writeTemperatures)
    {
        std::ofstream file(filename);

        if (!file)
        {
            throw std::runtime_error(
                "Cannot open file: " + filename
            );
        }

        file.imbue(std::locale::classic());

        file << std::scientific
             << std::setprecision(
                 std::numeric_limits<long double>::max_digits10
             );

        // Nagłówek.
        file << "# g";

        for (std::size_t j = 0; j < maxPeakCount; ++j)
        {
            file << ' ';

            if (writeTemperatures)
                file << "T_peak_" << j + 1;
            else
                file << "CV_peak_" << j + 1;
        }

        file << '\n';

        // Wiersze odpowiadające kolejnym wartościom g.
        for (const auto& row : scanRows)
        {
            file << row.g;

            for (std::size_t j = 0; j < maxPeakCount; ++j)
            {
                file << ' ';

                if (j < row.peaks.size())
                {
                    if (writeTemperatures)
                    {
                        file << row.peaks[j].temperature;
                    }
                    else
                    {
                        file << row.peaks[j].heatCapacity;
                    }
                }
                else
                {
                    file << "NaN";
                }
            }

            file << '\n';
        }

        file.close();

        if (!file)
        {
            throw std::runtime_error(
                "Error writing file: " + filename
            );
        }
    };

    const std::string suffix =
    config.seniorityZeroOnly ? "zero" : "all";

    writeWideFile(
        "peak_temperatures_vs_g_" + suffix + ".dat",
        true
    );

    writeWideFile(
        "peak_capacities_vs_g_" + suffix + ".dat",
        false
    );
    output.close();

    if (!output)
    {
        throw std::runtime_error(
            "Error while writing peak scan results."
        );
    }

    std::cout << "Peak scan saved to: "
              << filename << '\n';
}

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

        std::vector<SenioritySector> sectors;

        if (config.seniorityZeroOnly)
        {
            sectors = generateSeniorityZero(model);

            std::cout << "Seniority mode: zero\n";
        }
        else
        {
            SectorGenerator generator(model);
            sectors = generator.generate();

            std::cout << "Seniority mode: all\n";
        }

        std::cout << "Number of sectors: "
                  << sectors.size() << '\n';

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

        if (config.couplingScan.has_value())
        {
            scanHeatCapacityPeaks(config, sectors);
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

        std::cout << "Physical state count: ";
        if (config.seniorityZeroOnly)
        {
            const auto expected = countBoundedOccupations(
                sectors.front().effectivePairCapacity,
                sectors.front().pairCount
            );

            const auto actual = spectrum.physicalStateCount();

            if (actual != expected)
            {
                throw std::runtime_error(
                    "Incomplete seniority-zero spectrum."
                );
            }

            std::cout << "Seniority-zero dimension: "
                      << expected << '\n';
        }
        else
        {
            std::cout << spectrum.physicalStateCount() << '\n';
        }


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
        std::vector<ThermodynamicPoint> thermodynamicPoints;

        thermodynamicPoints.reserve(temperatures.size());
        for (double T : temperatures) {

            const ThermodynamicPoint result =
                thermo.evaluate(T);
            const auto point = thermo.evaluate(T);

            thermodynamicPoints.push_back(point);

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
        const auto peaks =
        findHeatCapacityPeaks(thermo, thermodynamicPoints);

        std::ofstream peakFile("heat_capacity_peaks.dat");

        if (!peakFile)
        {
            throw std::runtime_error(
                "Cannot open heat_capacity_peaks.dat"
            );
        }

        peakFile << std::scientific
                 << std::setprecision(
                        std::numeric_limits<long double>::max_digits10
                    );

        peakFile << "# T_peak C_V_peak\n";

        for (const auto& peak : peaks)
        {
            peakFile << peak.temperature << ' '
                     << peak.heatCapacity << '\n';
        }

        std::cout << "Detected heat capacity peaks: "
                  << peaks.size() << '\n';

        thermoFile.close();
        peakFile.close();

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
