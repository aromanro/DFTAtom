#pragma once

#include <vector>
#include <cassert>

namespace DFT {

	class Integral
	{
	public:

		static double Trapezoid(const double delta, const std::vector<double>& values)
		{
			assert(values.size() >= 2);

			double sum = 0.5 * (values.front() + values.back());
		
			const int szm1 = static_cast<int>(values.size()) - 1;

			for (int i = 1; i < szm1; ++i)
				sum += values[i];

			return sum * delta;
		}

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

		static double Boole(const double delta, const std::vector<double>& values)
		{
			assert(values.size() > 4);
			assert(values.size() % 4 == 1);

			double sum = 7. * (values.front() + values.back());

			double sum32 = 0;
			double sum12 = 0;
			double sum14 = 0;

			const int szm = static_cast<int>(values.size() - 1);
			
			for (int i = 1; i < szm; ++i)
			{
				sum32 += values[i++];

				if (i < szm)
				{
					if (i % 4 == 0) sum14 += values[i];
					else sum12 += values[i];
				}
			}

			sum += 32. * sum32 + 12. * sum12 + 14 * sum14;

			return sum * delta * 2. / 45.;
		}

	};

}




