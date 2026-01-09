//  the field of computer science has provided us with many efficient ways to find codewords. this is NOT one of them.
//  however, simulated annealing is an interesting topic and a crude (and unsafe) version can be implemented as such
//  to find some codewords corredponding to an error correcting code.

#include <hyper-metropolis/metro.hpp>

#include <iostream>
#include <iomanip>



int main() {
    std::vector<std::vector<uint32_t>> ldpc32bit_chks = {
        {0, 2,  5,  8},     {1, 3,  6,  9},     {2, 4,  7,  10},    {0, 5,  11, 15},
        {1, 6,  12, 16},    {3, 7,  13, 17},    {4, 8,  14, 18},    {0, 9,  19, 23},
        {1, 10, 20, 24},    {2, 11, 21, 25},    {3, 12, 22, 26},    {4, 13, 23, 27},
        {5, 14, 24, 28},    {6, 15, 25, 29},    {7, 16, 26, 30},    {8, 17, 27, 31}
    };
    std::vector<uint32_t> initial_state = {
        0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 
        1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0
    };
    metro::hamiltonian::LDPC ldpc(ldpc32bit_chks, initial_state);

    float temp = 5;     //  temperature
    float cool = 0.97;  //  cooling multiplier per step
    uint32_t seed = 7;  //  mt seed

    std::mt19937 mt(seed);



    std::cout << "----- LDPC Simulated Annealing --------------------------------------------\n";
    std::cout << std::setw(10) << "Step" 
              << std::setw(10) << "Energy" 
              << std::setw(10) << "Temp" 
              << std::setw(10) << "State\n";
    std::cout << "---------------------------------------------------------------------------\n";
    


    uint32_t step = 0;
    while (ldpc.energy() != 0) {
        metro::algorithm::metropolis(ldpc, mt, temp);

        std::cout << std::setw(10) << step
                  << std::setw(10) << (int)ldpc.energy()
                  << std::setw(10) << std::fixed << std::setprecision(4) << temp
                  << "    ";

        auto state = ldpc.state();
        for (size_t i = 0; i < state.size(); ++i) {
            std::cout << state[i];
            if ((i+1) % 8 == 0 && i != state.size()-1) std::cout << " ";
        }
        std::cout << "\n";

        temp *= cool;
        step++;
    }



    std::cout << "----- Simulation Complete -------------------------------------------------\n";
    std::cout << "      Codeword found:             ";
    for (size_t i = 0; i < ldpc.state().size(); ++i) {
        std::cout << ldpc.state()[i];
        if ((i+1) % 8 == 0 && i != ldpc.state().size()-1) std::cout << " ";
    }
    std::cout << std::endl;



    return 0;
}