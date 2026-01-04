#pragma once

#include <vector>
#include <cstdint>



namespace metro {

	namespace hamiltonian {



		template<typename state_t> 
		class Base {
		public:

			virtual float compute_energy() = 0;

			virtual float energy_delta(uint32_t i) const = 0;

			virtual void flip_spin(uint32_t i) = 0;

			void set_state(const std::vector<state_t> & newState) { 
				state = newState; 
				state_energy = compute_energy(); 
			}

			const std::vector<state_t> & get_state() const { 
				return state; 
			}

			float get_energy() const { 
				return state_energy; 
			}

		protected:

			Base() {}

			std::vector<state_t> state{};

			float state_energy{};
		};



		class LDPC : public Base<uint8_t> {
		public:

			LDPC(const std::vector<std::vector<uint32_t>> & parity_checks);

			float compute_energy() override;

			float energy_delta(uint32_t i) const override;

			void flip_spin(uint32_t i) override;

		private:

			const std::vector<std::vector<uint32_t>> parchk;

			const std::vector<std::vector<uint32_t>> parchkT;		
		};

	}

}