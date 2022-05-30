#pragma once

#include <cassert>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

namespace DFT {

	// Vosko-Wilk-Nusair 
	// see Richard M. Martin, Electronic Structure, Basic Theory and Practical Methods
	// also here: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-exchange-term
	// and here: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-6-3

	class VWNExchCor
	{
	protected:
		static constexpr double fourM_PI = 4. * M_PI;

		// values for 'paramagnetic' variant (used for LDA)
		static constexpr double AP = 0.0310907; // actually 0.5 * A
		static constexpr double y0P = -0.10498;
		static constexpr double bP = 3.72744;
		static constexpr double cP = 12.93532;
		static constexpr double Y0P = y0P * y0P + bP * y0P + cP;

		// values for 'feromagnetic' variant (useful for LSDA)
		static constexpr double AF = 0.01554535; // actually 0.5 * A
		static constexpr double y0F = -0.325;
		static constexpr double bF = 7.06042;
		static constexpr double cF = 18.0578;
		static constexpr double Y0F = y0F * y0F + bF * y0F + cF;

		// values for 'spin stiffness' (useful for LSDA)
		static constexpr double Aalpha = -1. / (6. * M_PI * M_PI); // actually 0.5 * A
		static constexpr double y0Falpha = -0.0047584;
		static constexpr double bFalpha = 1.13107;
		static constexpr double cFalpha = 13.0045;
		static constexpr double Y0Falpha = y0Falpha * y0Falpha + bFalpha * y0Falpha + cFalpha;

		inline static double F(double y /*sqrt(rs)*/, double dify /*y - y0*/, double A, double y0, double b, double c, double Y0, double Y)
		{
			const double Q = sqrt(4 * c - b * b);
			const double atanQ = atan(Q / (2. * y + b));

			return A * (log(y * y / Y) + 2. * b / Q * atanQ - b * y0 / Y0 * (log(dify * dify / Y) + 2. * (b + 2. * y0) / Q * atanQ)); // B.5
		}

		inline static double ecDif(double y, double dify, double A, double y0, double b, double c, double Y)
		{
			return A * (c * dify - b * y0 * y) / (dify * Y);
		}

		inline static double f(double zeta)
		{
			return (pow(1. + zeta, 4./3.) + pow(1. - zeta, 4. / 3.) - 2.) / (2. * (pow(2., 1. / 3.) - 1.)); // eq 5 from NIST
		}

	public:
		static std::vector<double> Vexc(const std::vector<double>& n)
		{
			static const double
				X1 = pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient
							
			std::vector<double> res(n.size());

			for (int i = 0; i < n.size(); ++i)
			{
				const double ro = n[i];
				
				const double rs = pow(3. / (fourM_PI*ro), 1. / 3.);
				const double y = sqrt(rs);
				const double Y = y * y + bP * y + cP;
				const double dify = y - y0P;

				res[i] = -X1 / rs // exchange term
					//the following make the Vc as in B.1
					+ F(y, dify, AP, y0P, bP, cP, Y0P, Y) // B.5
					- 1. / 3. * ecDif(y, dify, AP, y0P, bP, cP, Y); // B.6
			}

			return res;
		}



		static std::vector<double> eexcDif(const std::vector<double>& n)
		{
			static const double X1 = 0.25 * pow(3. / (2.*M_PI), 2. / 3.);  // Exchange energy coefficient
				
			std::vector<double> res(n.size());

			for (int i = 0; i < n.size(); ++i)
			{
				const double ro = n[i];
				
				const double rs = pow(fourM_PI / 3.*ro, -1. / 3.);
				const double y = sqrt(rs);
				const double Y = y * y + bP * y + cP;
				const double dify = y - y0P;

				res[i] = X1 / rs // exchange term
					+ 1./ 3. * ecDif(y, dify, AP, y0P, bP, cP, Y); // B.6
			}

			return res;
		}
	};

}