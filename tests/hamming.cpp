#include <hyper-metropolis/metro.hpp>

#include <vector>
#include <iostream>
#include <random>



int main() {

    std::vector<std::vector<uint32_t>> hamming74_chks = {
        {0, 2, 4, 6},
        {1, 2, 5, 6},
        {3, 4, 5, 6}
    };

    metro::hamiltonian::LDPC hamming(hamming74_chks);

    hamming.set_state({0, 0, 1, 1, 1, 1, 1});

    std::cout << "E0 = " << hamming.energy() << std::endl;

    std::mt19937 mt(0);

    while (hamming.energy() != 0) {
        metro::algorithm::metropolis(hamming, mt, 0.1);
        std::cout << "E = " << hamming.energy() << std::endl;

        auto state = hamming.state();
        for (int i = 0; i < state.size(); i++) {
            std::cout << (int)state[i] << ", ";
        }
        std::cout << std::endl;
    }

    return 0;
}