#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>

namespace DFT {


	class  ChachiyoExchCor
	{
	protected:
		static constexpr double a = (M_LN2 - 1.) / (2. * M_PI * M_PI);
		static constexpr double b = 20.4562557;

		static constexpr double fourM_PI = 4. * M_PI;
		static constexpr double threeDivM_PI = 3. / M_PI;

	public:
		static std::vector<double> exc(const std::vector<double>& n)
		{
			size_t sz = n.size();
			std::vector<double> res(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double ro = n[i];

				// exchange
				res[i] = -3. / 4. * pow(threeDivM_PI * ro, 1. / 3.); // Dirac exchange

				// correlation
				const double rs = pow(3. / (fourM_PI * ro), 1. / 3.);
				const double bprs = b / rs;

				res[i] += a * log(1. + bprs + bprs / rs);
			}

			return res;
		}


		static std::vector<double> excDeriv(const std::vector<double>& n)
		{
			size_t sz = n.size();
			std::vector<double> res(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double ro = n[i];

				// exchange
				res[i] = -1. / 4. * pow(threeDivM_PI, 1. / 3.) * pow(ro, -2. / 3.);

				// correlation
				const double rs = pow(3. / (fourM_PI * ro), 1. / 3.);
				const double bprs = b / rs;
				const double bprs2 = bprs / rs;
				const double bprs3 = bprs2 / rs;

				res[i] += a * (bprs2 + 2. * bprs3) / (1. + bprs + bprs2) * rs / (3. * ro);
			}

			return res;
		}
	};

}
