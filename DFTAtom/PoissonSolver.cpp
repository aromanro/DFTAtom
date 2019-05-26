#include "PoissonSolver.h"

#include <cassert>
#include <algorithm>

namespace DFT {

	PoissonSolver::PoissonSolver(int levels, double dGrid, int Ncoarse)
		: PhiLevels(levels), SourceLevels(levels), deltaGrid(dGrid), deltaGridLevel(levels), m_lowBoundary(0), m_highBoundary(0)
	{
		assert(Ncoarse >= 3);

		int size = Ncoarse;
		for (int i = levels - 1; i >= 0; --i)
		{
			PhiLevels[i].resize(size);
			SourceLevels[i].resize(size);
			size = size * 2 - 1;
		}

		double dGridLevel = deltaGrid;
		for (double& lvlDeltaGrid : deltaGridLevel)
		{
			lvlDeltaGrid = dGridLevel;
			dGridLevel *= 2;
		}
	}

	void PoissonSolver::SetBoundaries(double lowBoundary, double highBoundary)
	{
		m_lowBoundary = lowBoundary;
		m_highBoundary = highBoundary;
	}

	// Phi[i + 1] - Phi[i - 1] appears for non uniform grid from the first derivative
	// it's actually from ((Phi[i + 1] - Phi[i]) / h + (Phi[i] - Phi[i - 1]) / h) / 2 =
	// (h is actually 1) (Phi[i + 1] - Phi[i] + Phi[i] - Phi[i - 1]) * 0.5 = (Phi[i + 1] - Phi[i - 1]) * 0.5


	double PoissonSolver::GaussSeidel(int lvl)
	{
		const std::vector<double>& Source = SourceLevels[lvl];
		std::vector<double>& Phi = PhiLevels[lvl];
		double error2 = 0;

		const int limit = static_cast<int>(Phi.size() - 1);

		for (int i = 1; i < limit; ++i)
		{
			const double savePhi = Phi[i];

			// 0.5 comes from the 1 / diagonal term
			// Phi[i - 1] + Phi[i + 1] from off diagonal for laplaceian
			// the rest is for the first derivative out of the transform for non uniform grid

			Phi[i] = 0.5 * (Source[i] + Phi[i - 1] + Phi[i + 1]
				- deltaGridLevel[lvl] * (Phi[i + 1] - Phi[i - 1]) * 0.5); // this is for the non uniform grid

			const double dif = savePhi - Phi[i];
			error2 += dif * dif;
		}

		return sqrt(error2);
	}

	double PoissonSolver::IterateGaussSeidel(int lvl, double errorMin, int iterno)
	{
		double err = 1E10;
		for (int i = 0; i < iterno; ++i)
		{
			err = GaussSeidel(lvl);

			if (err < errorMin) break;
		}

		return err;
	}


	void PoissonSolver::Initialize(double errorMin)
	{
		std::fill(PhiLevels[0].begin(), PhiLevels[0].end(), 0);

		for (int i = 1; i < SourceLevels.size(); ++i)
		{
			const int im1 = i - 1;

			const int limit = static_cast<int>(SourceLevels[i].size() - 1);
			for (int p = 1; p < limit; ++p)
			{
				const int twop = 2 * p;

				SourceLevels[i][p] = 4 * SourceLevels[im1][twop];
				PhiLevels[i][p] = 0;
			}

			SourceLevels[i][0] = SourceLevels[i][limit] = 0;
			PhiLevels[i][0] = PhiLevels[i][limit] = 0;
		}

		const int coarsestIndex = static_cast<int>(PhiLevels.size() - 1);
		PhiLevels[coarsestIndex][0] = m_lowBoundary;
		PhiLevels[coarsestIndex][PhiLevels[coarsestIndex].size() - 1] = m_highBoundary;

		IterateGaussSeidel(coarsestIndex, errorMin, 15);
	}



	void PoissonSolver::Prolong(const std::vector<double>& Phisrc, std::vector<double>& Phidst)
	{
		assert(2 * (Phisrc.size() - 1) == Phidst.size() - 1);

		Phidst[0] += Phisrc[0];
		for (int i = 1; i < Phisrc.size(); ++i)
		{
			const int twoi = 2 * i;

			// correct the solution, use linear interpolation for 'in between' points
			Phidst[twoi] += Phisrc[i];
			Phidst[twoi - 1] += 0.5 * (Phisrc[i - 1] + Phisrc[i]);
		}
	}


	void PoissonSolver::Restrict(int lvl)
	{
		const int lvlm1 = lvl - 1;
		const std::vector<double>& Phisrc = PhiLevels[lvlm1];
		const std::vector<double>& Sourcesrc = SourceLevels[lvlm1];
		std::vector<double>& Phidst = PhiLevels[lvl];
		std::vector<double>& Sourcedst = SourceLevels[lvl];

		assert(2 * (Phidst.size() - 1) == Phisrc.size() - 1);
		assert(2 * (Sourcedst.size() - 1) == Sourcesrc.size() - 1);

		for (int i = 0; i < Phidst.size(); ++i)
			Phidst[i] = 0;

		for (int i = 1; i < Sourcedst.size() - 1; ++i)
		{
			const int twoi = 2 * i;

			// the 'source' is multiplied by delta^2 (delta = h or step, 1 for the non-uniform case)
			// but that means it needs adjusting the value with 4 when going to a coarser grid - for the source and second derivative
			// with 2 for the first derivative (2 * 0.5 = 1)
			const double adjresidual = 4. * (Sourcesrc[twoi] + Phisrc[twoi - 1] - 2. * Phisrc[twoi] + Phisrc[twoi + 1])
				- deltaGridLevel[lvl] * (Phisrc[twoi + 1] - Phisrc[twoi - 1]); // for non uniform grid

			Sourcedst[i] = adjresidual;
		}

		Sourcedst[0] = Sourcedst[Sourcedst.size() - 1] = 0; // no residual at boundaries, we know the solution exactly
	}




	void PoissonSolver::Ascend(int fromLevel, int toLevel, double errorMin, int iterno)
	{
		for (int i = fromLevel; i < toLevel; ++i)
		{
			IterateGaussSeidel(i, errorMin, iterno);
			Restrict(i + 1);
		}

		IterateGaussSeidel(toLevel, errorMin, iterno);
	}

	double PoissonSolver::Descend(int fromLevel, int toLevel, double errorMin, int iterno)
	{
		double err = 1E10;
		for (int i = fromLevel; i > toLevel; --i)
		{
			const int im1 = i - 1;
			Prolong(PhiLevels[i], PhiLevels[im1]);

			err = IterateGaussSeidel(im1, errorMin, iterno);
		}

		return err;
	}


	void PoissonSolver::Ascend(int toLevel, double errorMin, int iterno)
	{
		Ascend(0, toLevel, errorMin, iterno);
	}

	double PoissonSolver::Descend(int fromLevel, double errorMin, int iterno)
	{
		return Descend(fromLevel, 0, errorMin, iterno);
	}


	void PoissonSolver::FillR(std::vector<double>& R, double firstR, double lastR)
	{
		const int size = static_cast<int>(R.size());

		assert(size >= 3);

		const int N = size - 1;

		for (int i = 0; i < size; ++i)
			R[i] = (firstR * (N - i) + lastR * i) / N;
	}

	void PoissonSolver::FillRNonuniformR(std::vector<double>& R, double lastR, double deltaGrid)
	{
		const int size = static_cast<int>(R.size());

		assert(size >= 3);

		const int N = size - 1;
		const double Rp = lastR / (exp(N * deltaGrid) - 1.);

		for (int i = 0; i < size; ++i)
			R[i] = Rp * (exp(i * deltaGrid) - 1.);
	}

}


