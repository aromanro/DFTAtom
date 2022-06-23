#pragma once

#include <vector>
#include "Numerov.h"
#include "AufbauPrinciple.h"

namespace DFT {

	class DFTAtom
	{
	public:
		static const char orb[];

		static void CalculateNonUniformLDA(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
		static void CalculateUniformLDA(int Z, int MultigridLevels, double alpha, double MaxR);

		static void CalculateNonUniformLSDA(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
	private:
		static constexpr double fourM_PI = 4. * M_PI;

		static void LoopOverLevels(Numerov<NumerovFunctionRegularGrid>& numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double MaxR, double h, bool& reallyConverged, double energyErr, bool lda = true, bool isAlpha = true);
		static void LocateInterval(Numerov<NumerovFunctionRegularGrid>& numerov, double& TopEnergy, double& BottomEnergy, double MaxR, int L, int NumSteps, int NumNodes, double energyErr);

		static void LoopOverLevels(Numerov<NumerovFunctionNonUniformGrid>& numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double Rp, double deltaGrid, bool& reallyConverged, double energyErr, bool lda = true, bool isAlpha = true);
		static void LocateInterval(Numerov<NumerovFunctionNonUniformGrid>& numerov, double& TopEnergy, double& BottomEnergy, int L, int NumSteps, int NumNodes, double energyErr);

		static void NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid);
		static void NormalizeUniform(std::vector<double>& Psi, double h);
	};

}

