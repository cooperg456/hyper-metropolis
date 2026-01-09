#pragma once
#include <vector>
#include <cstdint>

namespace metro {
namespace hamiltonian {



	#pragma region Base

	/// @brief base class for constructing hamiltonians. contains all necessary objects to be interfaced by metro::algorithm
	template <typename Derived>
	class Base {
	public:

		//  ----------  ----------  virtuals  ----------  ----------  //

		///	@brief compute the energy of the current spin config
		///	@return float
		float compute_energy() const { return static_cast<const Derived*>(this)->compute_energy_impl(); }

		///	@brief compute the change in energy if the spin of site 'i' were to flip
		///	@param i index in spin configuration
		///	@return float 
		float energy_delta(uint32_t i) const { return static_cast<const Derived*>(this)->energy_delta_impl(i); }

		///	@brief flip spin of site 'i' and update config energy
		///	@param i index in spin configuration
		void flip_spin(uint32_t i) { return static_cast<Derived*>(this)->flip_spin_impl(i); }

		//  ----------  ----------  setters and getters  ----------  ----------  //

		///	@brief set the spin_config object and update config energy
		///	@param state 
		void set_state(const std::vector<uint32_t> &state) { 
			spin_config = state; 
			config_energy = compute_energy(); 
		}

		///	@brief return the current spin configuration
		///	@return const std::vector<uint32_t>& 
		const std::vector<uint32_t>& state() const { return spin_config; }

		///	@brief return the current config energy
		///	@return float 
		float energy() const { return config_energy; }

	protected:
		/// @brief construct a new base hamiltonian
		Base() = default;

		//  ----------  ----------  member vars  ----------  ----------  //

		std::vector<uint32_t> spin_config{};	///<	current spin configuration
		float config_energy{};					///<	energy of current spin configuration
	};

	#pragma endregion Base



	#pragma region LDPC

	/// @brief hamiltonian for classical LDPC codes
	class LDPC : public Base<LDPC> {
	public:
		
		///	@brief construct a new LDPC hamiltonian
		///	@param checks vector of parity checks
		LDPC(const std::vector<std::vector<uint32_t>> &checks) 
			: parchk(checks), nchks(static_cast<uint32_t>(parchk.size())), 
			parchkT(transpose(checks)), nbits(static_cast<uint32_t>(parchkT.size())) {}

		//  ----------  ----------  implementations  ----------  ----------  //

		/// @brief implementation for compute_energy()
		float compute_energy_impl() const {
			float E = 0.0;
			// loop over each parity check
			for (const std::vector<uint32_t> &chk : parchk) {
				uint32_t chk_parity = 0;
				for (uint32_t b : chk) {
					// xor bits affected by check
					chk_parity ^= spin_config[b];
				}
				E += static_cast<float>(chk_parity);
			}
			return E;
		}

		/// @brief implementation for energy_delta()
		float energy_delta_impl(uint32_t i) const {
			float dE = 0.0;
			// loop over each parity check containing 'i'
			for (uint32_t ci : parchkT[i]) {
				uint32_t chk_parity = 0;
				for (uint32_t b : parchk[ci]) {
					// xor bits affected by check
					chk_parity ^= spin_config[b];
				}
				dE += 1 - static_cast<float>(2 * chk_parity);
			}
			return dE;
		}

		/// @brief implementation for flip_spin()
		void flip_spin_impl(uint32_t i) {
			config_energy += energy_delta(i);
			spin_config[i] ^= 1;  // flip bit
		}

		//  ----------  ----------  member vars  ----------  ----------  //

		const std::vector<std::vector<uint32_t>> parchk;	///<	parity checks
		const std::vector<std::vector<uint32_t>> parchkT;	///<	parity checks transpose
		const uint32_t nbits;								///<	number of bits
		const uint32_t nchks;								///<	number of checks 

	private:	
	
		//  ----------  ----------  helper functions  ----------  ----------  //

		const std::vector<std::vector<uint32_t>> transpose(const std::vector<std::vector<uint32_t>> &mat) const {
			//	find num bits
			uint32_t max_bit = 0;
			for (const std::vector<uint32_t> &row : mat) {
				for (uint32_t b : row) {
					if (b > max_bit) {
						max_bit = b;
					}
				}
			}
			std::vector<std::vector<uint32_t>> trans(max_bit + 1);
			//	fill transpose
			for (uint32_t ci = 0; ci < mat.size(); ++ci) {
				for (uint32_t b : mat[ci]) {
					trans[b].push_back(ci);
				}
			}
			return trans;
		}
	
	};

	#pragma endregion LDPC


	
}}	// namespace metro::hamiltonian
