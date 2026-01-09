#pragma once
#include <vector>
#include <cstdint>
#include <random>
#include <cmath>

#include <hyper-metropolis/hamiltonian.hpp>

namespace metro {
namespace algorithm {



	/// @brief step metropolis hastings algorithm
	/// @tparam hamiltonian derived of metro::hamiltonian::Base
	/// @param ham hamiltonian
	/// @param rng mersenne twister rng
	/// @param temp temperature value
	template<typename hamiltonian> 
	void metropolis(hamiltonian &ham, std::mt19937 &rng, float temp) {
		std::uniform_real_distribution<float> reals(0.0, 1.0);
		std::uniform_int_distribution<uint32_t> ints(0, ham.state().size() - 1);

		uint32_t idx = ints(rng);
		float dE = ham.energy_delta(idx);
		// Accept/reject flip
		if ((dE <= 0) || (reals(rng) < std::exp(-dE / temp))) {
			ham.flip_spin(idx);
		}
	}



}}	// namespace metro::algorithm