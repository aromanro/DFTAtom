#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

#include "PoissonSolver.h"

#include "Integral.h"
#include "VWNExcCor.h"

#include "DFTAtom.h"

namespace DFT {

	const char DFTAtom::orb[] = { 's', 'p', 'd', 'f' };


	void DFTAtom::NormalizeUniform(std::vector<double>& Psi, double h)
	{
		std::vector<double> result2(Psi.size());
		for (int i = 0; i < result2.size(); ++i)
			result2[i] = Psi[i] * Psi[i];

		const double integralForSquare = DFT::Integral::Boole(h, result2);
		const double unorm = 1. / sqrt(integralForSquare);

		for (int i = 0; i < Psi.size(); ++i)
			Psi[i] *= unorm;
	}



	void DFTAtom::NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid)
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

		const double integralForSquare = DFT::Integral::Boole(1, result2); // for nonuniform case the step is 1
		const double unorm = 1. / sqrt(integralForSquare);

		for (int i = 0; i < Psi.size(); ++i)
			Psi[i] *= unorm;
	}



	void DFTAtom::CalculateUniform(int Z, int MultigridLevels, double alpha, double MaxR)
	{
		static const double energyErr = 1E-12;

		const double oneMinusAlpha = 1. - alpha;

		const int NumGridNodes = PoissonSolver::GetNumberOfNodes(MultigridLevels);

		const int NumSteps = NumGridNodes - 1;

		const double h = MaxR / NumSteps;

		std::vector<double> density(NumGridNodes);

		Potential potential;
		potential.m_potentialValues.resize(NumGridNodes);


		std::vector<Subshell> levels = AufbauPrinciple::GetSubshells(Z);
		std::sort(levels.begin(), levels.end());

		double Eold = 0;


		const double volume = 4. / 3. * M_PI * MaxR * MaxR * MaxR;
		const double constDens = Z / volume;

		density[0] = 0;
		for (int i = 1; i < density.size(); ++i)
			density[i] = constDens;


		PoissonSolver poissonSolver(MultigridLevels);

		std::vector<double> UHartree = poissonSolver.SolvePoissonUniform(Z, MaxR, density);
		std::vector<double> Vexc = VWNExchCor::Vexc(density);


		potential.m_potentialValues[0] = 0;
		for (int i = 1; i < NumGridNodes; ++i)
		{
			const double realPos = h * i;
			potential.m_potentialValues[i] = (-Z + UHartree[i]) / realPos + Vexc[i];
		}

		bool lastTimeConverged = false;

		for (int sp = 0; sp < 100; ++sp)
		{
			std::cout << "Step: " << sp << std::endl;

			double Eelectronic = 0;

			std::vector<double> newDensity(density.size(), 0);

			Numerov<NumerovFunctionRegularGrid> numerov(potential, 0, MaxR, NumGridNodes);

			bool reallyConverged = true;
			double BottomEnergy = -double(Z) * Z - 1.;

			LoopOverLevels(numerov, levels, newDensity, Eelectronic, BottomEnergy, NumSteps, MaxR, h, reallyConverged, energyErr);

			for (int i = 1; i < density.size(); ++i)
			{
				const double position = i * h;

				// 4 * M_PI appears because we're in spherical coordinates
				// the actual integration for the 'true' wavefunction gives a 4 M_PI
				// the radial wavefunction is actually u / r, whence also the division by position * position to get the true density

				newDensity[i] /= fourM_PI * position * position;
				density[i] = alpha * density[i] + oneMinusAlpha * newDensity[i];
			}


			UHartree = poissonSolver.SolvePoissonUniform(Z, MaxR, density);
			Vexc = VWNExchCor::Vexc(density);

			// Nuclear energy:
			std::vector<double> nuclear(NumGridNodes);

			// Exchange-correlation energy:
			std::vector<double> exccor(NumGridNodes);

			std::vector<double> eexcDeriv = VWNExchCor::eexcDif(density);

			// Hartree energy:
			std::vector<double> hartree(NumGridNodes);

			// potential energy:
			std::vector<double> potentiale(NumGridNodes);

			potential.m_potentialValues[0] = 0;
			nuclear[0] = 0;
			exccor[0] = 0;
			eexcDeriv[0] = 0;
			hartree[0] = 0;
			potentiale[0] = 0;

			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double position = i * h;

				potential.m_potentialValues[i] = (-Z + UHartree[i]) / position + Vexc[i];

				const double positiondensity = position * density[i];

				nuclear[i] = Z * positiondensity;

				const double position2density = position * position * density[i];

				exccor[i] = position2density * Vexc[i];
				eexcDeriv[i] = position2density * eexcDeriv[i];
				hartree[i] = positiondensity * UHartree[i];

				potentiale[i] = position2density * potential.m_potentialValues[i];
			}

			const double Enuclear = -fourM_PI * DFT::Integral::Boole(h, nuclear);

			double Exc = 4 * M_PI * DFT::Integral::Boole(h, exccor);

			const double eExcDif = fourM_PI * DFT::Integral::Boole(h, eexcDeriv);
			Exc += eExcDif;

			const double Ehartree = -2 * M_PI * DFT::Integral::Boole(h, hartree);

			const double Epotential = fourM_PI * DFT::Integral::Boole(h, potentiale);

			const double Ekinetic = Eelectronic - Epotential;
			const double Etotal = Eelectronic + Ehartree + eExcDif;

			std::cout << "Etotal = " << std::setprecision(12) << Etotal << " Ekin = " << std::setprecision(12) << Ekinetic << " Ecoul = " << std::setprecision(12) << -Ehartree << " Eenuc = " << std::setprecision(12) << Enuclear << " Exc = " << std::setprecision(12) << Exc << std::endl;

			if (abs((Eold - Etotal) / Etotal) < 1E-9 && reallyConverged && lastTimeConverged)
			{
				std::cout << std::endl << "Finished!" << std::endl << std::endl;

				break;
			}
			Eold = Etotal;
			lastTimeConverged = reallyConverged;

			std::cout << "********************************************************************************" << std::endl;
		}

		// sort levels by energy, just in case the energy values for levels do not come up as in the expected order (there are exceptions to the aufbau principle, too)
		std::sort(levels.begin(), levels.end(), [](const auto& val1, const auto& val2) -> bool { return val1.E < val2.E; });

		for (const auto& level : levels)
			std::cout << level.m_N + 1 << orb[level.m_L] << level.m_nrElectrons << " ";
	}


	void DFTAtom::LoopOverLevels(Numerov<NumerovFunctionRegularGrid> &numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double MaxR, double h, bool& reallyConverged, double energyErr)
	{
		for (auto& level : levels)
		{
			const int NumNodes = level.m_N - level.m_L;

			double TopEnergy = 50;

			// locate the interval to search into by using the number of nodes of the wavefunction

			LocateInterval(numerov, TopEnergy, BottomEnergy, MaxR, level.m_L, NumSteps, NumNodes, energyErr);

			// ***************************************************************************************************************************

			// locate the solution using the bisection method on the interval found above

			// it's the shooting method, it's supposed to shoot for 'zero' in the origin
			// sometimes it gets far away (whence the 1E15 comparison below)
			// the errors are too big, the line with the comment 'now really solve it' shoots from both 'infinity' and zero, matching the solutions in between				
			// there is another method that could be used, to shoot from both directions and do the match trying to have a fit for the derivative, for now I won't use it

			double delta = numerov.SolveSchrodingerSolutionInZero(MaxR, level.m_L, BottomEnergy, NumSteps);
			const bool sgnBottom = delta > 0;

			bool didNotConverge = true;
			for (int i = 0; i < 500; ++i)
			{
				level.E = (TopEnergy + BottomEnergy) / 2;

				delta = numerov.SolveSchrodingerSolutionInZero(MaxR, level.m_L, level.E, NumSteps);
				if ((delta > 0) == sgnBottom)
					BottomEnergy = level.E;
				else
					TopEnergy = level.E;

				const double absdelta = abs(delta);
				if (TopEnergy - BottomEnergy < energyErr && !isnan(absdelta) && absdelta < 1E15)
				{
					didNotConverge = false;
					break;
				}
			}
			level.E = BottomEnergy;

			// ***************************************************************************************************************************

			if (didNotConverge)
				reallyConverged = false;

			BottomEnergy = level.E - 3; // can happen sometimes to have it lower (see for example W, 4f is higher than 5s) 

			// now really solve it	
			long int matchPoint;
			std::vector<double> result = numerov.SolveSchrodingerMatchSolutionCompletely(MaxR, level.m_L, level.E, NumSteps, matchPoint);
			NormalizeUniform(result, h);

			std::cout << "Energy " << level.m_N + 1 << orb[level.m_L] << ": " << std::setprecision(12) << level.E << " Num nodes: " << NumNodes << std::endl;

			for (int i = 0; i < result.size() - 1; ++i)
				newDensity[i] += level.m_nrElectrons * result[i] * result[i];


			Eelectronic += level.m_nrElectrons * level.E;
		}
	}


	void DFTAtom::LocateInterval(Numerov<NumerovFunctionRegularGrid>& numerov, double& TopEnergy, double& BottomEnergy, double MaxR, int L, int NumSteps, int NumNodes, double energyErr)
	{
		double toe = TopEnergy;
		double boe = BottomEnergy;
		double deltaEnergy = toe - boe;
		while (deltaEnergy > energyErr)
		{
			const double E = (toe + boe) / 2;

			int NumNodesCounted;
			numerov.SolveSchrodingerCountNodes(MaxR, L, E, NumSteps, NumNodes, NumNodesCounted);

			if (NumNodesCounted > NumNodes)
				toe = E;
			else
				boe = E;

			deltaEnergy = toe - boe;
		}
		TopEnergy = toe;

		boe = BottomEnergy;
		deltaEnergy = toe - boe;
		while (deltaEnergy > energyErr)
		{
			const double E = (toe + boe) / 2;

			int NumNodesCounted;
			numerov.SolveSchrodingerCountNodes(MaxR, L, E, NumSteps, NumNodes, NumNodesCounted);

			if (NumNodesCounted < NumNodes)
				boe = E;
			else
				toe = E;

			deltaEnergy = toe - boe;
		}
		BottomEnergy = toe;
	}


	void DFTAtom::CalculateNonUniform(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid)
	{
		static const double energyErr = 1E-12; 

		const double oneMinusAlpha = 1. - alpha;
		
		const int NumGridNodes = PoissonSolver::GetNumberOfNodes(MultigridLevels);

		const int NumSteps = NumGridNodes - 1;
		const double Rp = MaxR / (exp(NumSteps * deltaGrid) - 1.);


		std::vector<double> density(NumGridNodes);

		Potential potential;
		potential.m_potentialValues.resize(NumGridNodes);


		std::vector<Subshell> levels = AufbauPrinciple::GetSubshells(Z);
		std::sort(levels.begin(), levels.end());

		double Eold = 0;

		const double volume = fourM_PI / 3. * MaxR * MaxR * MaxR;
		const double constDens = Z / volume;

		density[0] = 0;
		for (int i = 1; i < density.size(); ++i)
			density[i] = constDens;


		PoissonSolver poissonSolver(MultigridLevels, deltaGrid);

		std::vector<double> UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
		std::vector<double> Vexc = VWNExchCor::Vexc(density);


		potential.m_potentialValues[0] = 0;
		for (int i = 1; i < NumGridNodes; ++i)
		{
			const double realPos = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = (-Z + UHartree[i]) / realPos + Vexc[i];
		}

		bool lastTimeConverged = false;


		
		for (int sp = 0; sp < 100; ++sp)
		{
			std::cout << "Step: " << sp << std::endl;

			double Eelectronic = 0;

			std::vector<double> newDensity(density.size(), 0);

			Numerov<NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, MaxR, NumGridNodes);
			
			bool reallyConverged = true;
			double BottomEnergy = -double(Z) * Z - 1.;

			LoopOverLevels(numerov, levels, newDensity, Eelectronic, BottomEnergy, NumSteps, Rp, deltaGrid, reallyConverged, energyErr);

			for (int i = 1; i < density.size(); ++i)
			{
				const double position = Rp * (exp(i * deltaGrid) - 1.);

				// 4 * M_PI appears because we're in spherical coordinates
				// the actual integration for the 'true' wavefunction gives a 4 M_PI
				// the radial wavefunction is actually u / r, whence also the division by position * position to get the true density

				newDensity[i] /= fourM_PI * position * position;
				density[i] = alpha * density[i] + oneMinusAlpha * newDensity[i];
			}


			UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
			Vexc = DFT::VWNExchCor::Vexc(density);

			// Nuclear energy:
			std::vector<double> nuclear(NumGridNodes);
			
			// Exchange-correlation energy:
			std::vector<double> exccor(NumGridNodes);

			std::vector<double> eexcDeriv = DFT::VWNExchCor::eexcDif(density);

			// Hartree energy:
			std::vector<double> hartree(NumGridNodes);

			// potential energy:
			std::vector<double> potentiale(NumGridNodes);

			potential.m_potentialValues[0] = 0;
			nuclear[0] = 0;
			exccor[0] = 0;
			eexcDeriv[0] = 0;
			hartree[0] = 0;
			potentiale[0] = 0;

			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				potential.m_potentialValues[i] = (-Z + UHartree[i]) / position + Vexc[i];

				const double positiondensity = position * density[i] * cnst;

				nuclear[i] = Z * positiondensity;

				const double position2density = position * position * density[i] * cnst;

				exccor[i] = position2density * Vexc[i];
				eexcDeriv[i] = position2density * eexcDeriv[i];
				hartree[i] = positiondensity * UHartree[i];

				potentiale[i] = position2density * potential.m_potentialValues[i];
			}

			const double Enuclear = -fourM_PI * DFT::Integral::Boole(1, nuclear);
			double Exc = fourM_PI * DFT::Integral::Boole(1, exccor);

			const double eExcDif = fourM_PI * DFT::Integral::Boole(1, eexcDeriv);
			Exc += eExcDif;
			
			const double Ehartree = -2 * M_PI * DFT::Integral::Boole(1, hartree);

			const double Epotential = fourM_PI * DFT::Integral::Boole(1, potentiale);

			const double Ekinetic = Eelectronic - Epotential;
			const double Etotal = Eelectronic + Ehartree + eExcDif;

			std::cout << "Etotal = " << std::setprecision(12) << Etotal << " Ekin = " << std::setprecision(12) << Ekinetic << " Ecoul = " << std::setprecision(12) << -Ehartree << " Eenuc = " << std::setprecision(12) << Enuclear << " Exc = " << std::setprecision(12) << Exc << std::endl;

			if (abs((Eold - Etotal) / Etotal) < 1E-10 && reallyConverged && lastTimeConverged)
			{
				std::cout << std::endl << "Finished!" << std::endl << std::endl;

				break;
			}
			Eold = Etotal;
			lastTimeConverged = reallyConverged;

			std::cout << "********************************************************************************" << std::endl;
		}

		// sort levels by energy, just in case the energy values for levels do not come up as in the expected order (there are exceptions to the aufbau principle, too)
		std::sort(levels.begin(), levels.end(), [](const auto& val1, const auto& val2) -> bool { return val1.E < val2.E; });

		for (const auto& level : levels)
			std::cout << level.m_N + 1 << orb[level.m_L] << level.m_nrElectrons << " ";
	}

	void DFTAtom::LoopOverLevels(Numerov<NumerovFunctionNonUniformGrid>& numerov, std::vector<Subshell>& levels, std::vector<double>& newDensity, double& Eelectronic, double& BottomEnergy, int NumSteps, double Rp, double deltaGrid, bool& reallyConverged, double energyErr)
	{
		for (auto& level : levels)
		{
			const int NumNodes = level.m_N - level.m_L;

			double TopEnergy = 50;

			// locate the interval to search into by using the number of nodes of the wavefunction
			LocateInterval(numerov, TopEnergy, BottomEnergy, level.m_L, NumSteps, NumNodes, energyErr);

			// ***************************************************************************************************************************

			// locate the solution using the bisection method on the interval found above

			// it's the shooting method, it's supposed to shoot for 'zero' in the origin
			// sometimes it gets far away (whence the 1E15 comparison below)
			// the errors are too big, the line with the comment 'now really solve it' shoots from both 'infinity' and zero, matching the solutions in between				
			// there is another method that could be used, to shoot from both directions and do the match trying to have a fit for the derivative, for now I won't use it

			double delta = numerov.SolveSchrodingerSolutionInZero(NumSteps, level.m_L, BottomEnergy, NumSteps);
			const bool sgnBottom = delta > 0;

			bool didNotConverge = true;
			for (int i = 0; i < 500; ++i)
			{
				level.E = (TopEnergy + BottomEnergy) / 2;

				delta = numerov.SolveSchrodingerSolutionInZero(NumSteps, level.m_L, level.E, NumSteps);
				if ((delta > 0) == sgnBottom)
					BottomEnergy = level.E;
				else
					TopEnergy = level.E;

				const double absdelta = abs(delta);
				if (TopEnergy - BottomEnergy < energyErr && !isnan(absdelta) && absdelta < 1E15)
				{
					didNotConverge = false;
					break;
				}
			}
			level.E = BottomEnergy;

			// ***************************************************************************************************************************

			if (didNotConverge)
				reallyConverged = false;

			BottomEnergy = level.E - 3; // can happen sometimes to have it lower (see for example W, 4f is higher than 5s) 

			// now really solve it	
			long int matchPoint;
			std::vector<double> result = numerov.SolveSchrodingerMatchSolutionCompletely(NumSteps, level.m_L, level.E, NumSteps, matchPoint);
			NormalizeNonUniform(result, Rp, deltaGrid);

			std::cout << "Energy " << level.m_N + 1 << orb[level.m_L] << ": " << std::setprecision(12) << level.E << " Num nodes: " << NumNodes << std::endl;

			for (int i = 0; i < result.size() - 1; ++i)
				newDensity[i] += level.m_nrElectrons * result[i] * result[i];


			Eelectronic += level.m_nrElectrons * level.E;
		}
	}


	void DFTAtom::LocateInterval(Numerov<NumerovFunctionNonUniformGrid>& numerov, double& TopEnergy, double& BottomEnergy, int L, int NumSteps, int NumNodes, double energyErr)
	{
		double toe = TopEnergy;
		double boe = BottomEnergy;
		double deltaEnergy = toe - boe;
		while (deltaEnergy > energyErr)
		{
			const double E = (toe + boe) / 2;

			int NumNodesCounted;
			numerov.SolveSchrodingerCountNodes(NumSteps, L, E, NumSteps, NumNodes, NumNodesCounted);

			if (NumNodesCounted > NumNodes)
				toe = E;
			else
				boe = E;

			deltaEnergy = toe - boe;
		}
		TopEnergy = toe;

		boe = BottomEnergy;
		deltaEnergy = toe - boe;
		while (deltaEnergy > energyErr)
		{
			const double E = (toe + boe) / 2;

			int NumNodesCounted;
			numerov.SolveSchrodingerCountNodes(NumSteps, L, E, NumSteps, NumNodes, NumNodesCounted);

			if (NumNodesCounted < NumNodes)
				boe = E;
			else
				toe = E;

			deltaEnergy = toe - boe;
		}
		BottomEnergy = toe;
	}
}