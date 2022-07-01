#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

namespace DFT {

	class PoissonSolver
	{
	protected:
		static constexpr double fourM_PI = 4. * M_PI;

	public:
		PoissonSolver(int levels, double dGrid = 0, int Ncoarse = 3);

		// the radial Poisson eqn is L Uhartree = -4 * M_PI * r * Rho
		// minus is implied in the poisson solver, we deal here with electron charge

		std::vector<double> SolvePoissonUniform(int Z, double maxRadius, const std::vector<double>& density)
		{
			int coarsestLevel = static_cast<int>(SourceLevels.size() - 1);
			auto &Source = SourceLevels[0];

			PoissonSolver::FillR(Source, 0, maxRadius);

			const double delta = Source[1] - Source[0];
			const double delta2 = delta * delta;

			// the delta^2 multiplication is for convenience

			// Poisson equation also applies to the error and residual

			// so we go with the residual instead of the 'charge' and error instead of the real solution when computing at the other levels, whence the += in the code, the 'solution' gets corrected
			// the equation is Operator(error) = residual (for Poisson the operator is the Laplaceian, for non-uniform it also has a first derivative but still linear)

			// it happens that for our current 'guess' of the solution (zero everywhere except for known boundary values)		
			// we have the residual equal with the 'charge'

			const double delta2fourM_PI = delta2 * fourM_PI;
			for (int i = 0; i < Source.size(); ++i)
				Source[i] *= delta2fourM_PI * density[i];


			SetBoundaries(0, Z);

			FullCycle(1E-3, 1E-14);

			return PhiLevels[0];
		}

		std::vector<double> SolvePoissonNonUniform(int Z, double maxRadius, const std::vector<double>& density)
		{
			auto &Source = SourceLevels[0];

			PoissonSolver::FillRNonuniformR(Source, maxRadius, deltaGrid);


			// Poisson equation also applies to the error and residual

			// so we go with the residual instead of the 'charge' and error instead of the real solution when computing at the other levels, whence the += in the code, the 'solution' gets corrected
			// the equation is Operator(error) = residual (for Poisson the operator is the Laplaceian, for non-uniform it also has a first derivative but still linear)

			// it happens that for our current 'guess' of the solution (zero everywhere except for known boundary values)
			// we have the residual equal with the 'charge'
			const double Rp = maxRadius / (exp((density.size() - 1.) * deltaGrid) - 1.);
			const double delta2grid = deltaGrid * deltaGrid;
			const double Rp2delta2 = Rp * Rp * delta2grid;
			const double twodelta = 2. * deltaGrid;

			const double fourM_PIRp2delta2 = fourM_PI * Rp2delta2;
			for (int i = 1; i < Source.size() - 1; ++i)
				// step becomes 1 for the nonuniform grid, so multiplication with delta2 as above disappears from here
				Source[i] *= fourM_PIRp2delta2 * exp(i * twodelta) * density[i];

			SetBoundaries(0, Z);

			FullCycle(1E-3, 1E-14);

			return PhiLevels[0];
		}


		static void FillR(std::vector<double>& R, double firstR, double lastR);
		static void FillRNonuniformR(std::vector<double>& R, double lastR, double deltaGrid);

		void SetBoundaries(double lowBoundary, double highBoundary);

		double FullCycle(double errorMin = 0.001, double errorMinLast = 0.00001, double firstError = 0.3)
		{
			static const int numSweeps = 3;
			const int lastLevel = static_cast<int>(PhiLevels.size() - 1);

			Initialize(errorMin);

			// The following does this:

			// -  -  -  -  -  -  -  -  -  -  -  -  -  coarsest grid, highest index
			//  * *   *     *       
			//   * * * *   * *      
			//      *   * *   *      
			//           *     *      
			//                  *       
			// -------------------------------------- finest grid, index 0


			for (int i = static_cast<int>(PhiLevels.size() - 2); i > 0; --i)
			{
				Descend(lastLevel, i, errorMin, numSweeps); // this goes from the coarsest grid to level 'i'
				Ascend(i, lastLevel, errorMin, numSweeps); // this goes from level 'i' to the coarsest grid
			}

			Descend(lastLevel, 0, errorMinLast, numSweeps); // we are on the coarsest grid now, go to the finest grid

			// now full cycles until the error drops enough

			double err;
			for (int i = 0; i < 100; ++i)
			{
				err = VCycle(lastLevel, errorMinLast, numSweeps);
				if (err < errorMinLast) break;
			}

			return err;
		}

		// no need to improve it, it won't be called much
		static int GetNumberOfNodes(int levels, int Ncoarse = 3)
		{
			int size = Ncoarse;

			for (int i = 0; i < levels - 1; ++i)
				size = size * 2 - 1;

			return size;
		}


	protected:
		double GaussSeidel(int lvl);
		double IterateGaussSeidel(int lvl, double errorMin, int iterno);

		static void Prolong(const std::vector<double>& Phisrc, std::vector<double>& Phidst);
		void Restrict(int lvl);


		void Initialize(double errorMin);


		double Descend(int toLevel, double errorMin, int iterno);
		void Ascend(int toLevel, double errorMin, int iterno);

		double Descend(int fromLevel, int toLevel, double errorMin, int iterno);
		void Ascend(int fromLevel, int toLevel, double errorMin, int iterno);

		double VCycle(int level, double errorMin, int iterno)
		{
			Ascend(level, errorMin, iterno); // go to the coarser grid at 'level' from the finest grid
			return Descend(level, errorMin, iterno); // go back on the finest grid from the coarse grid at 'level'
		}

		// the finest grid is at index 0
		// the coarsest is at 'lastLevel'

		std::vector<std::vector<double>> PhiLevels;
		std::vector<std::vector<double>> SourceLevels;

		const double deltaGrid; //for nonuniform grid is non zero, for uniform, just let this be zero
		std::vector<double> deltaGridLevel;

		double m_lowBoundary;
		double m_highBoundary;
	};

}

