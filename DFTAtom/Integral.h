#pragma once

#include <vector>
#include <cassert>

namespace DFT {

	class Integral
	{
	public:

		static double SimpsonOneThird(const double delta, const std::vector<double>& values)
		{
			assert(values.size() >= 5);
			assert(values.size() % 2);

			double sum = values.front() + values.back();

			double sum4 = 0;
			double sum2 = 0;

			const int szm1 = static_cast<int>(values.size()) - 1;

			for (int i = 1; i < szm1; ++i)
			{
				sum4 += values[i++];
				if (i < szm1) sum2 += values[i];
			}

			sum += 4. * sum4 + 2. * sum2;


			return sum * delta / 3.;
		}

	};

}




