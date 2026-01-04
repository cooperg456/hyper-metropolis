#include <hyper-metropolis/hamiltonian.hpp>
#include "helpers.hpp"



#pragma region LDPC

metro::hamiltonian::LDPC::LDPC(const std::vector<std::vector<uint32_t>> & parity_checks) 
: parchk(parity_checks), parchkT(helpers::transpose(parity_checks)) {}

float metro::hamiltonian::LDPC::compute_energy() {
	float E = 0;
	// Loop over each parity check
	for (uint32_t ci = 0; ci < parchk.size(); ci++) {
		uint8_t chk_parity = 0;
		for (uint32_t b : parchk[ci]) {
			// XOR bits affected by check
			chk_parity ^= state[b];
		}
		E += (float)chk_parity;
	}
	return E;
}

float metro::hamiltonian::LDPC::energy_delta(uint32_t i) const {
	float dE = 0;
	// Compute difference in energy if bit b were to flip
	// Loop over each parity check containing b
	for (uint32_t ci : parchkT[i]) {
		int8_t chk_parity = 0;
		for (int64_t bit : parchk[ci]) {
			// XOR bits affected by check
			chk_parity ^= state[bit];
		}
		dE += 1 - (float)(2 * chk_parity);
	}
	return dE;
}

void metro::hamiltonian::LDPC::flip_spin(uint32_t i) {
	state[i] ^= 1;  // flip bit
	state_energy += energy_delta(i);
}