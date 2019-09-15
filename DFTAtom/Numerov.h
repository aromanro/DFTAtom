#pragma once

#include <vector>

namespace DFT {

	class Potential
	{
	public:
		inline double operator()(size_t posIndex) const { return m_potentialValues[posIndex]; };

		std::vector<double> m_potentialValues;
	};


	class NumerovFunctionRegularGrid
	{
	public:
		NumerovFunctionRegularGrid(const Potential& pot, double /*delta*/, double /*Rmax*/, size_t /*numPoints*/) : m_pot(pot) {}

		inline double GetEffectivePotential(unsigned int l, double position, size_t posIndex) const
		{
			return m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;
		}

		inline double operator()(unsigned int l, double E, double position, size_t posIndex) const
		{
			const double effectivePotential = GetEffectivePotential(l, position, posIndex);

			return  2. * (effectivePotential - E);
		}

		inline static double GetBoundaryValueFar(double position, double E)
		{
			return exp(-position * sqrt(2. * abs(E))); //not tested!
		}

		inline static double GetBoundaryValueZero(double position, unsigned int l)
		{
			return pow(position, static_cast<size_t>(l) + 1);
		}

		inline static double GetMaxRadiusIndex(double E, double stepSize)
		{
			return GetMaxRadius(E) / stepSize;
		}

		inline static double GetDerivativeStep(int posIndex, double h)
		{
			return h;
		}

	protected:
		inline static double GetMaxRadius(double E)
		{
			return 12. / sqrt(2. * abs(E));
		}

		const Potential& m_pot;
	};


	class NumerovFunctionNonUniformGrid
	{
	public:
		NumerovFunctionNonUniformGrid(const Potential& pot, double delta, double Rmax, size_t numPoints)
			: m_pot(pot), m_delta(delta)
		{
			Rp = Rmax / (exp((numPoints - 1) * delta) - 1);
			const double Rp2 = Rp * Rp;

			twodelta = 2. * m_delta;
			const double delta2 = m_delta * m_delta;

			Rp2delta2 = Rp2 * delta2;
			delta2p4 = delta2 / 4.;
		}

		inline double GetEffectivePotential(unsigned int l, double position, size_t posIndex) const
		{
			position = GetPosition(posIndex); // the passed value is ignored, use the real one

			return m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;
		}

		inline double operator()(unsigned int l, double E, double position, size_t posIndex) const
		{
			const double effectivePotential = GetEffectivePotential(l, position, posIndex);

			return  2. * (effectivePotential - E) * Rp2delta2 * exp(posIndex * twodelta) + delta2p4;
		}

		inline double GetBoundaryValueFar(double position, double E) const
		{
			const double realPosition = GetPosition(static_cast<int>(position));

			return exp(-realPosition * sqrt(2. * abs(E)));
		}

		inline double GetBoundaryValueZero(double position, unsigned int l) const
		{
			const int posInd = static_cast<int>(position);
			const double realPosition = GetPosition(posInd);

			return pow(realPosition, static_cast<size_t>(l) + 1) * exp(-position * m_delta * 0.5);
		}


		inline double GetMaxRadiusIndex(double E) const
		{
			const double maxRadius = GetMaxRadius(E);

			return log(maxRadius / Rp + 1.) / m_delta;
		}

		inline static double GetMaxRadius(double E)
		{
			return 12. / sqrt(2. * abs(E));
		}

		inline double GetDerivativeStep(int posIndex, double h) const
		{
			return Rp * exp(posIndex * m_delta) * (1. - exp(-m_delta));
		}

		inline double GetRp() const { return Rp; }
		inline double GetDelta() const { return m_delta; }
	protected:
		inline double GetPosition(size_t posIndex) const
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
		Numerov(const Potential& pot, double delta = 0, double Rmax = 0, size_t numPoints = 0) : function(pot, delta, Rmax, numPoints), h(1), h2(1), h2p12(1. / 12.) {}



		inline void SolveSchrodingerCountNodesFromNucleus(double endPoint, unsigned int l, double E, size_t steps, size_t nodesLimit, int& nodesCount)
		{
			if (endPoint == steps)
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadiusIndex(E));
				steps = static_cast<int>(endPoint);
			}
			else
			{
				h = endPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadius(E));
				steps = static_cast<int>(endPoint / h);
			}


			double position = 0;
			double prevSol = 0;
			double solution = 0;
			double wprev = 0;

			position += h;
			solution = function.GetBoundaryValueZero(position, l);
			double funcVal = function(l, E, position, 1);
			double w = (1 - h2p12 * funcVal) * solution;

			bool oldSgn = (solution > 0);
			nodesCount = 0;

			double effPotential = function.GetEffectivePotential(l, position, 1);
			bool firstClassicalReturnPoint = effPotential <= E;

			for (size_t i = 2; i <= steps; ++i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				solution = getU(w, funcVal);

				if (abs(solution) == std::numeric_limits<double>::infinity())
					return;

				const bool newSgn = (solution > 0);
				if (newSgn != oldSgn)
				{
					++nodesCount;
					if (nodesCount > nodesLimit)
						return;

					oldSgn = newSgn;
				}

				// bail out if hitting the classical turning point
				effPotential = function.GetEffectivePotential(l, position, i);
				if (effPotential <= E)
					firstClassicalReturnPoint = true;
				else if (firstClassicalReturnPoint && effPotential > E)
					return;
			}

		}

		inline void SolveSchrodingerCountNodes(double startPoint, unsigned int l, double E, size_t steps, size_t nodesLimit, int& nodesCount)
		{
			if (static_cast<size_t>(startPoint) == steps)
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadiusIndex(E));
				steps = static_cast<int>(startPoint);
			}
			else
			{
				h = startPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadius(E));
				steps = static_cast<int>(startPoint / h);
			}


			double position = startPoint;
			double solution = function.GetBoundaryValueFar(position, E);

			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			solution = function.GetBoundaryValueFar(position, E);
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;

			bool oldSgn = (solution > 0);
			nodesCount = 0;

			bool firstClassicalReturnPoint = false;
			for (size_t i = steps - 2; i > 0; --i)
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

				const bool newSgn = (solution > 0);
				if (newSgn != oldSgn)
				{
					++nodesCount;
					if (nodesCount > nodesLimit)
						return;

					oldSgn = newSgn;
				}

				// bail out if hitting the classical turning point
				const double effPotential = function.GetEffectivePotential(l, position, i);
				if (effPotential <= E)
					firstClassicalReturnPoint = true;
				else if (firstClassicalReturnPoint && effPotential > E)
					return;
			}

			if (nodesCount <= nodesLimit)
			{
				solution = solution * (2 + h2 * funcVal) - prevSol;
				if ((solution > 0) != oldSgn)
					++nodesCount;
			}
		}

		inline double SolveSchrodingerSolutionInZero(double startPoint, unsigned int l, double E, size_t steps)
		{
			if (static_cast<size_t>(startPoint) == steps)
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadiusIndex(E));
				steps = static_cast<int>(startPoint);
			}
			else
			{
				h = startPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadius(E));
				steps = static_cast<int>(startPoint / h);
			}

			double position = startPoint;
			double solution = function.GetBoundaryValueFar(position, E);

			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			solution = function.GetBoundaryValueFar(position, E);
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;

			for (size_t i = steps - 2; i > 0; --i)
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

		inline std::vector<double> SolveSchrodingerMatchSolutionCompletely(double startPoint, unsigned int l, double E, size_t steps, size_t& matchPoint)
		{
			const size_t highLimit = steps + 1;
			std::vector<double> Psi(highLimit);

			if (static_cast<size_t>(startPoint) == steps)
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadiusIndex(E));
				steps = static_cast<int>(startPoint);
			}
			else
			{
				h = startPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				startPoint = std::min(startPoint, function.GetMaxRadius(E));
				steps = static_cast<int>(startPoint / h);
			}
			for (size_t i = steps + 1; i < highLimit; ++i)
				Psi[i] = 0;

			h = startPoint / steps;
			h2 = h * h;
			h2p12 = h2 / 12.;

			const size_t size = steps + 1;


			double position = startPoint;
			double solution = function.GetBoundaryValueFar(position, E);
			Psi[steps] = solution;
			double prevSol = solution;
			double funcVal = function(l, E, position, steps);
			double wprev = (1 - h2p12 * funcVal) * solution;

			position -= h;
			Psi[steps - 1] = solution = function.GetBoundaryValueFar(position, E);
			funcVal = function(l, E, position, steps - 1);
			double w = (1 - h2p12 * funcVal) * solution;

			matchPoint = 2;
			for (size_t i = steps - 2; i > 0; --i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				Psi[i] = solution = getU(w, funcVal);

				//const double effPotential = function.GetEffectivePotential(l, position, i);				
				if (solution < Psi[i + 1] /*effPotential <= E*/ || abs(solution) > 1E15)
				{
					matchPoint = i;
					break;
				}

			}

			position = 0;
			prevSol = Psi[0] = solution = 0;
			wprev = 0;

			position += h;
			Psi[1] = solution = function.GetBoundaryValueZero(position, l);
			funcVal = function(l, E, position, 1);
			w = (1 - h2p12 * funcVal) * solution;

			for (size_t i = 2; i < matchPoint; ++i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				Psi[i] = solution = getU(w, funcVal);
			}

			w = 2. * w - wprev + h2 * solution * funcVal;
			position = h * matchPoint;
			funcVal = function(l, E, position, matchPoint);
			solution = getU(w, funcVal);

			const double factor = solution / Psi[matchPoint];

			Psi[matchPoint] = solution;
			for (size_t i = matchPoint + 1; i < size; ++i)
				Psi[i] *= factor;

			return std::move(Psi);
		}

		NumerovFunction function;

	protected:
		// 2.13
		inline double getU(double w, double funcVal) const
		{
			return w / (1. - h2p12 * funcVal);
		}

		double h;
		double h2;
		double h2p12;
	};

}


