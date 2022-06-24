#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

namespace DFT
{

	class ExcCorBase
	{
	protected:
		inline static double f(double zeta)
		{
			static const double mul = 1. / (2. * (pow(2., 1. / 3.) - 1.));

			return mul * (pow(1. + zeta, 4. / 3.) + pow(1. - zeta, 4. / 3.) - 2.); // eq 5 from NIST
		}

		inline static double df(double zeta)
		{
			static const double mul = 4. / (6. * (pow(2., 1. / 3.) - 1.));

			return mul * (pow(1. + zeta, 1. / 3.) - pow(1. - zeta, 1. / 3.)); // eq 5 from NIST
		}
	};

}

