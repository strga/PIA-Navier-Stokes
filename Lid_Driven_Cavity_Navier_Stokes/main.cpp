#include <iostream>
#include <chrono>
#include <limits>
#include "Simulation.h"

/*

mkdir build
cd build
cmake ..
make
./Lid_Driven_Cavity_Navier_Stokes

*/

int main() {

    int N;
    double domainLength, kinematicViscosity;


    // Input Validation for N, domainLength, and kinematicViscosity
    std::cout << "Enter the amount of points in cavity/domain (rec.: 100): ";
    while (!(std::cin >> N) || N <= 0) {
        std::cout << "Invalid input, enter a positive integer: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::cout << "Enter length of cavity/domain (rec.: 1): ";
    while (!(std::cin >> domainLength) || domainLength <= 0) {
        std::cout << "Invalid input, enter a positive number: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::cout << "Enter kinematic viscosity (rec.: 0.01): ";
    while (!(std::cin >> kinematicViscosity) || kinematicViscosity < 0) {
        std::cout << "Invalid input, enter a non-negative number: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::cout << "Calculating, please wait and make yourself a coffee...\n";

    auto start = std::chrono::high_resolution_clock::now();

    runSimulation(N, domainLength, kinematicViscosity);

    std::cout << "Thank you for waiting, here are the results, stored where this code is stored\n";

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by function: " << duration.count() << " microseconds\n";

    return 0;
}
