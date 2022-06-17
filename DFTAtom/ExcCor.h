#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>

namespace DFT {

	// see https://aip.scitation.org/doi/10.1063/1.4958669
	class ChachiyoExchCorParam
	{
	public:
		static constexpr double b = 20.4562557;
		static constexpr double b1 = 27.4203609;
	};

	// see https://aip.scitation.org/doi/10.1063/1.4964758
	class ChachiyoExchCorImprovedParam
	{
	public:
		static constexpr double b = 21.7392245;
		static constexpr double b1 = 28.3559732;
	};

	template<class Params> class ChachiyoExchCor
	{
	protected:
		static constexpr double a = (M_LN2 - 1.) / (2. * M_PI * M_PI);
		static constexpr double b = Params::b;

		static constexpr double a1 = (M_LN2 - 1.) / (4. * M_PI * M_PI);
		static constexpr double b1 = Params::b1;

		static constexpr double fourM_PI = 4. * M_PI;
		static constexpr double threeDivM_PI = 3. / M_PI;

	public:

		static std::vector<double> Vexc(const std::vector<double>& n)
		{
			static const double	X1 = pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient

			size_t sz = n.size();
			std::vector<double> res(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double ro = n[i];
				const double rs = pow(3. / (fourM_PI * ro), 1. / 3.);

				// exchange
				res[i] = -X1 / rs;

				// correlation
				const double bprs = b / rs;

				res[i] += a * log(1. + bprs + bprs / rs);
			}

			return res;
		}


		static std::vector<double> eexcDif(const std::vector<double>& n)
		{
			static const double X1 = 0.25 * pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient

			size_t sz = n.size();
			std::vector<double> res(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double ro = n[i];
				const double rs = pow(3. / (fourM_PI * ro), 1. / 3.);

				// exchange
				res[i] = X1 / rs;

				// correlation
				const double bprs = b / rs;
				const double bprs2 = bprs / rs;

				res[i] += a / (1. + bprs + bprs2) * (bprs + 2. * bprs2) * rs / 3.;
			}

			return res;
		}

	protected:
		inline static double f(double zeta)
		{
			static const double div = 2. * (pow(2., 1. / 3.) - 1.);

			return (pow(1. + zeta, 4. / 3.) + pow(1. - zeta, 4. / 3.) - 2.) / div; // eq 5 from NIST
		}
	};

}
