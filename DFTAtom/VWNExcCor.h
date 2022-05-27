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
		static constexpr double A = 0.0310907; // actually 0.5 * A
		static constexpr double y0 = -0.10498;
		static constexpr double b = 3.72744;
		static constexpr double c = 12.93532;
		static constexpr double Y0 = y0 * y0 + b * y0 + c;

		// values for 'feromagnetic' variant (useful for LSDA)
		static constexpr double AF = 0.01554535; // actually 0.5 * A
		static constexpr double y0F = -0.325;
		static constexpr double bF = 7.06042;
		static constexpr double cF = 18.0578;
		static constexpr double Y0F = y0F * y0F + bF * y0F + cF;

		// values for 'spin stiffness' (useful for LSDA)
		static constexpr double Aalpha = -1. / (6. * M_PI * M_PI); // actually 0.5 * A
		static constexpr double y0F = -0.0047584;
		static constexpr double bF = 1.13107;
		static constexpr double cF = 13.0045;
		static constexpr double Y0F = y0F * y0F + bF * y0F + cF;

	public:
		static std::vector<double> Vexc(const std::vector<double>& n)
		{
			static const double
				X1 = pow(3. / (2.*M_PI), 2. / 3.),  // Exchange energy coefficient
				Q = sqrt(4 * c - b * b);
				

			std::vector<double> res(n.size());

			for (int i = 0; i < n.size(); ++i)
			{
				const double ro = n[i];
				
				const double rs = pow(3. / (fourM_PI*ro), 1. / 3.);

				const double y = sqrt(rs);
				const double Y = y * y + b * y + c;

				const double atanQ = atan(Q / (2.*y + b));
				const double dify = y - y0;

				res[i] = -X1 / rs // exchange term
					//the following make the Vc as in B.1
					    + A * (log(y*y / Y) + 2.*b / Q * atanQ - b * y0 / Y0 * (log(dify*dify / Y) + 2.*(b + 2.*y0) / Q * atanQ)) // B.5
					    - A / 3. * (c * dify - b * y0 * y) / (dify * Y); // B.6
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
				const double Y = y * y + b * y + c;
				const double dify = y - y0;

				res[i] = X1 / rs // exchange term
					+ A / 3. * (c * dify - b * y0 * y) / (dify * Y); // B.6
			}

			return res;
		}
	};

}