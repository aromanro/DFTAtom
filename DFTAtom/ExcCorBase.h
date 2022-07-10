#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

namespace DFT
{

	class ExcCorBase
	{
	protected:
		static constexpr double aThird = 1. / 3.;

		inline static double f(double zeta)
		{
			static const double mul = 1. / (2. * (pow(2., aThird) - 1.));

			return mul * (pow(1. + zeta, 4. * aThird) + pow(1. - zeta, 4. * aThird) - 2.); // eq 5 from NIST
		}

		inline static double df(double zeta)
		{
			static const double mul = 2. / (3. * (pow(2., aThird) - 1.));

			return mul * (pow(1. + zeta, aThird) - pow(1. - zeta, aThird)); // eq 5 from NIST
		}
	};

}

