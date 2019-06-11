#pragma once

#include <vector>

namespace DFT {

	class DFTAtom
	{
	public:
		static void Calculate(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
	protected:
		static void Normalize(std::vector<double>& Psi, double Rp, double deltaGrid);
	};

}

