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

			constexpr double coef = 1. / 3.;

			return sum * delta * coef;
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

			constexpr double coef = 2. / 45.;

			return sum * delta * coef;
		}

		template<typename T> static T Romberg(const double delta, const std::vector<T>& values, const double err = 1E-18, const int minSteps = 3)
		{
			const int sz = static_cast<int>(values.size());

			assert(sz % 2);

			const int numPoints = sz - 1;

			int n = numPoints;
			int cnt = 0;
			while (n)
			{
				++cnt;
				n >>= 1;
			}

			std::vector<T> Rprev(cnt), Rcur(cnt);
			double h = delta * numPoints;
			Rprev[0] = 0.5 * h * (values[0] + values[numPoints]); // R(0,0) - trapezoid

			n = numPoints;
			for (int i = 1;i < cnt;++i)
			{
				const int oldStep = n;
				n >>= 1;
				//if (n == 0) break;

				T sum = 0;
				for (int j = n; j < numPoints; j += oldStep) 
					sum += values[j];
			
				h *= 0.5;
				Rcur[0] = 0.5 * Rprev[0] + h * sum; // R(n, 0) - still trapezoid

				T nk = 1;
				for (int m = 1; m <= i; ++m) 
				{
					nk *= 4;
					Rcur[m] = Rcur[m - 1] + (Rcur[m - 1] - Rprev[m - 1]) / (nk - 1); // R(n, m)
				}

				if (i >= minSteps && fabs(Rcur[i] - Rcur[i - 1]) < err)
					return Rcur[i];

				Rcur.swap(Rprev);
			}

			return Rprev.back();
		}
	};

}




