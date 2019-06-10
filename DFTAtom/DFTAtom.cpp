#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

#include "PoissonSolver.h"

#include "Integral.h"
#include "AufbauPrinciple.h"
#include "Numerov.h"
#include "VWNExcCor.h"

#include "DFTAtom.h"

namespace DFT {


	void DFTAtom::Calculate(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid)
	{
		static const char orb[] = { 's', 'p', 'd', 'f' };
		static const double energyErr = 1E-10;
		static const double derivErr = 1E-5;

		const double oneMinusAlpha = 1. - alpha;
		
		const int NumGridNodes = DFT::PoissonSolver::GetNumberOfNodes(MultigridLevels);

		const int NumSteps = NumGridNodes - 1;
		const double Rp = MaxR / (exp(NumSteps * deltaGrid) - 1.);


		std::vector<double> density(NumGridNodes);

		DFT::Potential potential;
		potential.m_potentialValues.resize(NumGridNodes);


		std::vector<DFT::Subshell> levels = DFT::AufbauPrinciple::GetSubshells(Z);
		std::sort(levels.begin(), levels.end());

		double Eold = 0;


		const double volume = 4. / 3. * M_PI * MaxR * MaxR * MaxR;
		const double constDens = Z / volume;

		density[0] = 0;
		for (int i = 1; i < density.size(); ++i)
			density[i] = constDens;


		DFT::PoissonSolver poissonSolver(MultigridLevels, deltaGrid);

		std::vector<double> UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
		std::vector<double> Vexc = DFT::VWNExchCor::Vexc(density);


		potential.m_potentialValues[0] = 0;
		for (int i = 1; i < NumGridNodes; ++i)
		{
			const double realPos = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = (-Z + UHartree[i]) / realPos + Vexc[i];
		}

		for (int sp = 0; sp < 500; ++sp)
		{
			std::cout << "Step: " << sp << std::endl;

			double Eelectronic = 0;

			std::vector<double> newDensity(density.size(), 0);

			DFT::Numerov<DFT::NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, MaxR, NumGridNodes);
			
			double BottomEnergy = - Z * Z - 1;

			bool reallyConverged = true;

			for (auto& level : levels)
			{
				const int NumNodes = level.m_N - level.m_L;

				double TopEnergy = 50;
				
				double toe = TopEnergy;
				double boe = BottomEnergy;
				double deltaEnergy = toe - boe;
				while (deltaEnergy > energyErr)
				{
					const double E = (toe + boe) / 2;

					int NumNodesCounted;
					numerov.SolveSchrodingerCountNodes(NumSteps, level.m_L, E, NumSteps, NumNodes, NumNodesCounted);

					if (NumNodesCounted > NumNodes)
						toe = E;
					else
						boe = E;

					deltaEnergy = toe - boe;
				}
				toe -= energyErr;

				TopEnergy = toe;

				boe = BottomEnergy;
				deltaEnergy = toe - boe;
				while (deltaEnergy > energyErr)
				{
					const double E = (toe + boe) / 2;

					int NumNodesCounted;
					numerov.SolveSchrodingerCountNodes(NumSteps, level.m_L, E, NumSteps, NumNodes, NumNodesCounted);

					if (NumNodesCounted < NumNodes)
						boe = E;
					else
						toe = E;

					deltaEnergy = toe - boe;
				}
				BottomEnergy = boe + energyErr;
				
				
				double topDelta = numerov.SolveSchrodingerMatch(NumSteps, level.m_L, TopEnergy, NumSteps);
				bool sgnTop = topDelta > 0;

				bool didNotConverge = true;
				for (int i = 0; i < 1000; ++i)
				{
					level.E = (TopEnergy + BottomEnergy) / 2;

					const double delta = numerov.SolveSchrodingerMatch(NumSteps, level.m_L, level.E, NumSteps);
					const double absdelta = abs(delta);

					if ((delta > 0) == sgnTop)
					{
						// the same sign as 'top'
						if (absdelta < abs(topDelta))
						{
							TopEnergy = level.E;
							topDelta = delta;
						}
						else
						{
							BottomEnergy = level.E;
						}
					}
					else
					{
						// different sign than 'top', the zero value is between them
						BottomEnergy = level.E;						
					}
					
					if (TopEnergy - BottomEnergy < energyErr && absdelta < derivErr && !isnan(absdelta))
					{
						didNotConverge = false;
						break;
					}						
				}

				if (didNotConverge) 
					reallyConverged = false;

				BottomEnergy = level.E - 3; // can happen sometimes to have it lower (see for example W, 4f is higher than 5s) 

				// now really solve it				
				std::vector<double> result = numerov.SolveSchrodingerMatchSolutionCompletely(NumSteps, level.m_L, level.E, NumSteps);

				// square the wavefunction
				// also integrate the square of wavefunction to get the normalization constant

				std::vector<double> result2(result.size());
				std::vector<double> result3(result.size());
				for (int i = 0; i < result.size(); ++i)
				{
					// if nonuniform, convert the function back!!!!!
					result[i] *= exp(i * deltaGrid * 0.5);

					result2[i] = result[i] * result[i];

					// if nonuniform, do the change for dr	
					const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
					result3[i] = result2[i] * cnst;
				}

				const double integralForSquare = DFT::Integral::SimpsonOneThird(1, result3); // for nonuniform case the step is 1

				std::cout << "Energy " << level.m_N + 1 << orb[level.m_L] << ": " << std::setprecision(12) << level.E << " Num nodes: " << NumNodes << std::endl;

				for (int i = 0; i < result.size() - 1; ++i)
					newDensity[i] += level.m_nrElectrons * result2[i] / integralForSquare;

				Eelectronic += level.m_nrElectrons * level.E;
			}


			for (int i = 1; i < density.size(); ++i)
			{
				const double position = Rp * (exp(i * deltaGrid) - 1.);

				// 4 * M_PI appears because we're in spherical coordinates
				// the actual integration for the 'true' wavefunction gives a 4 M_PI
				// the radial wavefunction is actually u / r, whence also the division by position * position to get the true density

				newDensity[i] /= 4. * M_PI * position * position;
				density[i] = alpha * density[i] + oneMinusAlpha * newDensity[i];
			}


			UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
			Vexc = DFT::VWNExchCor::Vexc(density);

			potential.m_potentialValues[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double position = Rp * (exp(i * deltaGrid) - 1.);

				potential.m_potentialValues[i] = (-Z + UHartree[i]) / position + Vexc[i];
			}

			// Nuclear energy:
			std::vector<double> nuclear(NumGridNodes);
			for (int i = 0; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				nuclear[i] = position * Z * density[i] * cnst;
			}
			const double Enuclear = -4 * M_PI * DFT::Integral::SimpsonOneThird(1, nuclear);

			// Exchange-correlation energy:
			std::vector<double> exccor(NumGridNodes);
			exccor[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				exccor[i] = position * position * density[i] * Vexc[i] * cnst;
			}
			double Exc = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, exccor);

			std::vector<double> eexcDeriv = DFT::VWNExchCor::eexcDif(density);
			exccor[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				exccor[i] = position * position * density[i] * eexcDeriv[i] * cnst;
			}
			const double eExcDif = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, exccor);
			Exc += eExcDif;

			// Hartree energy:
			std::vector<double> hartree(NumGridNodes);
			for (int i = 0; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				hartree[i] = position * density[i] * UHartree[i] * cnst;
			}
			const double Ehartree = -2 * M_PI * DFT::Integral::SimpsonOneThird(1, hartree);

			// potential energy:
			std::vector<double> potentiale(NumGridNodes);
			potentiale[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i) {
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);
				const double cnst = Rp * deltaGrid * expD;

				potentiale[i] = position * position * density[i] * potential.m_potentialValues[i] * cnst;
			}
			const double Epotential = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, potentiale);

			const double Ekinetic = Eelectronic - Epotential;
			const double Etotal = Eelectronic + Ehartree + eExcDif;

			std::cout << "Etotal = " << std::setprecision(12) << Etotal << " Ekin = " << std::setprecision(12) << Ekinetic << " Ecoul = " << std::setprecision(12) << -Ehartree << " Eenuc = " << std::setprecision(12) << Enuclear << " Exc = " << std::setprecision(12) << Exc << std::endl;

			if (abs((Eold - Etotal) / Etotal) < 1E-8 && reallyConverged)
			{
				std::cout << std::endl << "Finished!" << std::endl << std::endl;

				break;
			}
			Eold = Etotal;

			std::cout << "********************************************************************************" << std::endl;
		}

		// sort levels by energy, just in case the energy values for levels do not come up as in the expected order (there are exceptions to the aufbau principle, too)
		std::sort(levels.begin(), levels.end(), [](const auto& val1, const auto& val2) -> bool { return val1.E < val2.E; });

		for (const auto& level : levels)
			std::cout << level.m_N + 1 << orb[level.m_L] << level.m_nrElectrons << " ";
	}

}