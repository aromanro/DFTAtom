#pragma once

#include <vector>
#include "Numerov.h"
#include "AufbauPrinciple.h"

namespace DFT {

	class DFTAtom
	{
	public:
		static const char orb[];

		static void CalculateNonUniform(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
		static void CalculateUniform(int Z, int MultigridLevels, double alpha, double MaxR);

	private:
		static void LoopOverLevels(Numerov<NumerovFunctionRegularGrid>& numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double MaxR, double h, bool& reallyConverged, double energyErr);
		static void LoopOverLevels(Numerov<NumerovFunctionNonUniformGrid>& numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double Rp, double deltaGrid, bool& reallyConverged, double energyErr);

		static void NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid);
		static void NormalizeUniform(std::vector<double>& Psi, double h);
	};

}

