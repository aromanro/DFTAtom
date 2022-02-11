#pragma once

#include <vector>
#include <cassert>

namespace DFT {

	class Integral
	{
	public:
		template<typename T> static T Trapezoid(const double delta, const std::vector<T>& values)
		{
			assert(values.size() >= 2);

			T sum = 0.5 * (values.front() + values.back());

			const int szm1 = static_cast<int>(values.size()) - 1;

			for (int i = 1; i < szm1; ++i)
				sum += values[i];

			return sum * delta;
		}

		template<typename T> static T SimpsonOneThird(const double delta, const std::vector<T>& values)
		{
			assert(values.size() >= 5);
			assert(values.size() % 2);

			T sum = values.front() + values.back();

			T sum4 = 0;
			T sum2 = 0;

			const int szm1 = static_cast<int>(values.size()) - 1;

			for (int i = 1; i < szm1; ++i)
			{
				sum4 += values[i++];
				if (i < szm1) sum2 += values[i];
			}

			sum += 4. * sum4 + 2. * sum2;


			return sum * delta / 3.;
		}

		template<typename T> static T Boole(const double delta, const std::vector<T>& values)
		{
			assert(values.size() > 4);
			assert(values.size() % 4 == 1);

			T sum = 7. * (values.front() + values.back());

			T sum32 = 0;
			T sum12 = 0;
			T sum14 = 0;

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

			sum += 32. * sum32 + 12. * sum12 + 14. * sum14;

			return sum * delta * 2. / 45.;
		}
	};

}




