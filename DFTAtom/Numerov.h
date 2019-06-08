#pragma once

#include <vector>

namespace DFT {

	class Potential
	{
	public:
		inline double operator()(int posIndex) const { return m_potentialValues[posIndex]; };

		std::vector<double> m_potentialValues;
	};


	class NumerovFunctionRegularGrid
	{
	public:
		NumerovFunctionRegularGrid(const Potential& pot, double /*delta*/, double /*Rmax*/, int /*numPoints*/) : m_pot(pot) {}

		inline double operator()(unsigned int l, double E, double position, int posIndex) const
		{
			const double effectivePotential = m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;

			return  2. * (effectivePotential - E);
		}

		inline double GetBoundaryValue(double position, double E) const
		{
			return exp(-position * sqrt(2. * abs(E))); //not tested!
		}

	protected:
		const Potential& m_pot;
	};


	class NumerovFunctionNonUniformGrid
	{
	public:
		NumerovFunctionNonUniformGrid(const Potential& pot, double delta, double Rmax, int numPoints)
			: m_pot(pot), m_delta(delta)
		{
			Rp = Rmax / (exp((numPoints - 1) * delta) - 1);
			const double Rp2 = Rp * Rp;

			twodelta = 2. * m_delta;
			const double delta2 = m_delta * m_delta;

			Rp2delta2 = Rp2 * delta2;
			delta2p4 = delta2 / 4.;
		}

		inline double operator()(unsigned int l, double E, double position, int posIndex) const
		{
			position = GetPosition(posIndex); // the passed value is ignored, use the real one

			const double effectivePotential = m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;

			return  2. * (effectivePotential - E) * Rp2delta2 * exp(posIndex * twodelta) + delta2p4;
		}

		inline double GetBoundaryValue(double position, double E) const
		{
			const double realPosition = GetPosition(static_cast<int>(position));

			return exp(-realPosition * sqrt(2. * abs(E)) - position * m_delta * 0.5);
		}

	protected:
		inline double GetPosition(int posIndex) const
		{
			return Rp * (exp(posIndex * m_delta) - 1.);
		}

		const Potential& m_pot;

		const double m_delta;

		double Rp;
		double twodelta;
		double delta2p4;
		double Rp2delta2;
	};


	template<class NumerovFunction> class Numerov
	{
	public:
		Numerov(const Potential& pot, double delta = 0, double Rmax = 0, int numPoints = 0) : function(pot, delta, Rmax, numPoints) {}


		inline void SolveSchrodingerCountNodes(double startPoint, unsigned int l, double E, unsigned int steps, int nodesLimit, int& nodesCount)
		{
			h = startPoint / steps;
			h2 = h * h;
			h2p12 = h2 / 12.;


			double position = startPoint;
			double solution = function.GetBoundaryValue(position, E);
		
			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			solution = function.GetBoundaryValue(position, E);
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;

			bool oldSgn = (solution >= 0);
			nodesCount = 0;

			for (int i = steps - 2; i > 0; --i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				prevSol = solution;

				solution = getU(w, funcVal);

				if (abs(solution) == std::numeric_limits<double>::infinity())
					return;

				const bool newSgn = (solution >= 0);
				if (newSgn != oldSgn)
				{
					++nodesCount;
					if (nodesCount > nodesLimit)
						return;

					oldSgn = newSgn;
				}
			}

			if (nodesCount <= nodesLimit)
			{
				solution = solution * (2 + h2 * funcVal) - prevSol;
				if ((solution >= 0) != oldSgn)
					++nodesCount;
			}
		}

		inline double SolveSchrodingerSolutionInZero(double startPoint, unsigned int l, double E, unsigned int steps)
		{
			h = startPoint / steps;
			h2 = h * h;
			h2p12 = h2 / 12.;

			double position = startPoint;
			double solution = function.GetBoundaryValue(position, E);
		
			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			solution = function.GetBoundaryValue(position, E);
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;
			
			for (int i = steps - 2; i > 0; --i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				prevSol = solution;
				solution = getU(w, funcVal);
			}

			solution = solution * (2 + h2 * funcVal) - prevSol;

			return solution;
		}


		inline std::vector<double> SolveSchrodingerSolutionCompletely(double startPoint, unsigned int l, double E, unsigned int steps)
		{
			h = startPoint / steps;
			h2 = h * h;
			h2p12 = h2 / 12.;

			std::vector<double> Psi(steps + 1);

			double position = startPoint;
			double solution = function.GetBoundaryValue(position, E);
			Psi[steps] = solution;
			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			solution = function.GetBoundaryValue(position, E);
			Psi[steps - 1] = solution;
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;
						
			for (int i = steps - 2; i > 0; --i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);				
				Psi[i] = solution = getU(w, funcVal);
			}

			Psi[0] = Psi[1] * (2 + h2 * funcVal) - Psi[2];

			return Psi;
		}


	protected:
		// 2.13
		inline double getU(double w, double funcVal) const
		{
			return w / (1. - h2p12 * funcVal);
		}

		NumerovFunction function;

		double h;
		double h2;
		double h2p12;
	};

}


