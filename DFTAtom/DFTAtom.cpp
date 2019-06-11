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

	void DFTAtom::Normalize(std::vector<double>& Psi, double Rp, double deltaGrid)
	{
		std::vector<double> result2(Psi.size());
		for (int i = 0; i < result2.size(); ++i)
		{
			// if nonuniform, convert the function back!!!!!
			Psi[i] *= exp(i * deltaGrid * 0.5);

			result2[i] = Psi[i] * Psi[i];

			// if nonuniform, do the change for dr	
			const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
			result2[i] *= cnst;
		}

		const double integralForSquare = DFT::Integral::SimpsonOneThird(1, result2); // for nonuniform case the step is 1
		const double norm = sqrt(integralForSquare);

		for (int i = 0; i < Psi.size(); ++i)
			Psi[i] /= norm;
	}



	void DFTAtom::Calculate(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid)
	{
		static const char orb[] = { 's', 'p', 'd', 'f' };
		static const double energyErr = 1E-10;
		static const double derivErr = 1E-4; //set to bigger if matching the derivative in the match point

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
				//toe -= energyErr;
				toe -= deltaEnergy;

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
				//BottomEnergy = boe + energyErr;
				BottomEnergy = boe + deltaEnergy;
				
				// ***************************************************************************************************************************

				// what happens here is that I want to locate an interval where the difference between derivatives is one sign for 'top', the other one for 'bottom'
				// and locate the energy where the difference is closer to zero
				// but it doesn't always happen that way
				// in that case the choice would be just to pick the smallest value and that's it
				// I need to investigate this further, for now it's worse than the other method, although it kind of works...

				/*
				int matchPoint;
				std::vector<double> Psi = numerov.SolveSchrodingerMatchSolutionCompletely(NumSteps, level.m_L, TopEnergy, NumSteps, matchPoint);
				Normalize(Psi, Rp, deltaGrid);
				int nextPoint = matchPoint + 1;
				double deriv1 = (Psi[matchPoint] - Psi[matchPoint - 1]) / numerov.function.GetDerivativeStep(matchPoint, 1);
				double deriv2 = (Psi[nextPoint] - Psi[matchPoint]) / numerov.function.GetDerivativeStep(nextPoint, 1);
				double topDelta = deriv1 - deriv2;

				bool sgnTop = topDelta > 0;

				bool didNotConverge = true;
				for (int i = 0; i < 1000; ++i)
				{
					level.E = (TopEnergy + BottomEnergy) / 2;

					Psi = numerov.SolveSchrodingerMatchSolutionCompletely(NumSteps, level.m_L, level.E, NumSteps, matchPoint);
					Normalize(Psi, Rp, deltaGrid);
					nextPoint = matchPoint + 1;
					deriv1 = (Psi[matchPoint] - Psi[matchPoint - 1]) / numerov.function.GetDerivativeStep(matchPoint, 1);
					deriv2 = (Psi[nextPoint] - Psi[matchPoint]) / numerov.function.GetDerivativeStep(nextPoint, 1);
					const double delta = deriv1 - deriv2;

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
					
					const double energyDiff = TopEnergy - BottomEnergy;
					if (TopEnergy - BottomEnergy < energyErr && absdelta < derivErr && !isnan(absdelta))
					{
						didNotConverge = false;
						break;
					}
				}
				*/
				
				// ***************************************************************************************************************************

				double delta = numerov.SolveSchrodingerSolutionInZero(NumSteps, level.m_L, BottomEnergy, NumSteps);
				const bool sgnA = delta > 0;

				bool didNotConverge = true;
				for (int i = 0; i < 1000; ++i)
				{
					level.E = (TopEnergy + BottomEnergy) / 2;

					delta = numerov.SolveSchrodingerSolutionInZero(NumSteps, level.m_L, level.E, NumSteps);
					if ((delta > 0) == sgnA)
						BottomEnergy = level.E;
					else
						TopEnergy = level.E;

					const double absdelta = abs(delta);
					if (TopEnergy - BottomEnergy < energyErr && absdelta < derivErr && !isnan(absdelta))
					{
						didNotConverge = false;
						break;
					}						
				}

				// ***************************************************************************************************************************
											
				if (didNotConverge) 
					reallyConverged = false;

				BottomEnergy = level.E - 3; // can happen sometimes to have it lower (see for example W, 4f is higher than 5s) 

				// now really solve it	
				int matchPoint;
				std::vector<double> result = numerov.SolveSchrodingerMatchSolutionCompletely(NumSteps, level.m_L, level.E, NumSteps, matchPoint);
				//std::vector<double> result = numerov.SolveSchrodingerSolutionCompletely(NumSteps, level.m_L, level.E, NumSteps);
				Normalize(result, Rp, deltaGrid);

				std::cout << "Energy " << level.m_N + 1 << orb[level.m_L] << ": " << std::setprecision(12) << level.E << " Num nodes: " << NumNodes << std::endl;

				for (int i = 0; i < result.size() - 1; ++i)
					newDensity[i] += level.m_nrElectrons * result[i] * result[i];


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