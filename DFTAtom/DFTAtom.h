#pragma once

#include <vector>

namespace DFT {

	class DFTAtom
	{
	public:
		static void CalculateNonUniform(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
		static void CalculateUniform(int Z, int MultigridLevels, double alpha, double MaxR);

		static void NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid);
		static void NormalizeUniform(std::vector<double>& Psi, double h);
	};

}

