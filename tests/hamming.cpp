#include <hyper-metropolis/metro.hpp>

#include <vector>
#include <iostream>


int main() {

    std::vector<std::vector<uint32_t>> hamming74_chks = {
        {0, 2, 4, 6},
        {1, 2, 5, 6},
        {3, 4, 5, 6}
    };

    metro::hamiltonian::LDPC hamming(hamming74_chks);

    hamming.set_state({0, 0, 1, 1, 1, 1, 1});

    std::cout << hamming.get_energy() << std::endl;

    return 0;
}