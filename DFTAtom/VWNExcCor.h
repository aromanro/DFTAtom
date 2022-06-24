#pragma once

#include <cassert>
#include <vector>

#include "ExcCorBase.h"

namespace DFT {

	// Vosko-Wilk-Nusair 
	// see Richard M. Martin, Electronic Structure, Basic Theory and Practical Methods
	// also here: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-exchange-term
	// and here: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-6-3

	class VWNExchCor : public ExcCorBase
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
		static constexpr double y0alpha = -0.0047584;
		static constexpr double balpha = 1.13107;
		static constexpr double calpha = 13.0045;
		static constexpr double Y0alpha = y0alpha * y0alpha + balpha * y0alpha + calpha;

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

	public:
		static std::vector<double> Vexc(const std::vector<double>& n)
		{
			static const double	X1 = pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient
							
			std::vector<double> res(n.size());

			for (int i = 0; i < n.size(); ++i)
			{
				const double ro = n[i];
				if (ro < 1E-18)
				{
					res[i] = 0;
					continue;
				}
				
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
				if (ro < 1E-18)
				{
					res[i] = 0;
					continue;
				}

				const double rs = pow(3. / (fourM_PI * ro), 1. / 3.);
				const double y = sqrt(rs);
				const double Y = y * y + bP * y + cP;
				const double dify = y - y0P;

				res[i] = X1 / rs // exchange term
					+ 1./ 3. * ecDif(y, dify, AP, y0P, bP, cP, Y); // B.6
			}

			return res;
		}

		static std::vector<double> Vexc(const std::vector<double>& na, const std::vector<double>& nb, std::vector<double>& va, std::vector<double>& vb)
		{
			int sz = static_cast<int>(na.size());
			if (sz != nb.size()) return {};

			static const double	X1 = pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient
			static const double X2 = pow(2., 1. / 3.);
			static const double fdd = 1. / (9. * (pow(2., 1. / 3.) - 1.));

			std::vector<double> res(sz);

			va.resize(sz);
			vb.resize(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double roa = na[i];
				const double rob = nb[i];
				const double n = roa + rob;
				if (n < 1E-18)
				{
					res[i] = 0;
					continue;
				}

				const double rs = pow(3. / (fourM_PI * n), 1. / 3.);

				const double ep = -X1 / rs;
				const double ef = X2 * ep;

				const double zeta = (roa - rob) / n;

				const double fval = f(zeta);
				const double dfval = df(zeta);

				const double y = sqrt(rs);

				const double YP = y * y + bP * y + cP;
				const double difyP = y - y0P;


				const double ecp = F(y, difyP, AP, y0P, bP, cP, Y0P, YP); // B.5

				const double YF = y * y + bF * y + cF;
				const double difyF = y - y0F;

				const double ecf = F(y, difyF, AF, y0F, bF, cF, Y0F, YF); // B.5

				const double YA = y * y + balpha * y + calpha;
				const double difyA = y - y0alpha;

				const double eca = F(y, difyA, Aalpha, y0alpha, balpha, calpha, Y0alpha, YA); // B.5

				const double ecpd = ecDif(y, difyP, AP, y0P, bP, cP, YP);
				const double ecfd = ecDif(y, difyF, AF, y0F, bF, cF, YF);
				const double ecad = ecDif(y, difyA, Aalpha, y0alpha, balpha, calpha, YA);

				const double deltaec = ecf - ecp;
				const double deltaecd = ecfd - ecpd;
				const double beta = fdd * deltaec / eca - 1.;
				const double betad = fdd * fdd * (deltaecd * eca - ecad * deltaec) / (eca * eca);
				const double interp = fval / fdd * (1 + beta * pow(zeta, 4.));
				const double interpd = fval / fdd * pow(zeta, 4.) * betad;

				res[i] = ep + (ef - ep) * fval // exchange term
					// paramagnetic part:
					+ ecp
					// the rest of it:
					+ eca * interp;

				// TODO: something is wrong here, fix it!

				// now the derivative
				const double deriv = 1. / 3. * (ecpd // paramagnetic part
					+ ecad * interp);// +eca * interpd);

				va[i] = res[i] - (1. + zeta) * deriv;
				vb[i] = res[i] - (1. - zeta) * deriv;

				res[i] -= deriv;
			}

			return res;
		}


		static std::vector<double> eexcDif(const std::vector<double>& na, const std::vector<double>& nb)
		{
			int sz = static_cast<int>(na.size());
			if (sz != nb.size()) return {};

			static const double	X1 = 0.25 * pow(3. / (2. * M_PI), 2. / 3.);  // Exchange energy coefficient
			static const double X2 = pow(2., 1. / 3.);
			static const double fdd = 1. / (9. * (pow(2., 1. / 3.) - 1.));

			std::vector<double> res(sz);

			for (int i = 0; i < sz; ++i)
			{
				const double roa = na[i];
				const double rob = nb[i];
				const double n = roa + rob;
				if (n < 1E-18)
				{
					res[i] = 0;
					continue;
				}


				const double rs = pow(3. / (fourM_PI * n), 1. / 3.);

				const double ep = X1 / rs;
				const double ef = X2 * ep;

				const double zeta = (roa - rob) / n;

				const double fval = f(zeta);
				const double dfval = df(zeta);

				const double y = sqrt(rs);

				const double YP = y * y + bP * y + cP;
				const double difyP = y - y0P;


				const double ecp = F(y, difyP, AP, y0P, bP, cP, Y0P, YP); // B.5

				const double YF = y * y + bF * y + cF;
				const double difyF = y - y0F;

				const double ecf = F(y, difyF, AF, y0F, bF, cF, Y0F, YF); // B.5

				const double YA = y * y + balpha * y + calpha;
				const double difyA = y - y0alpha;

				const double eca = F(y, difyA, Aalpha, y0alpha, balpha, calpha, Y0alpha, YA); // B.5

				const double ecpd = ecDif(y, difyP, AP, y0P, bP, cP, YP);
				const double ecfd = ecDif(y, difyF, AF, y0F, bF, cF, YF);
				const double ecad = ecDif(y, difyA, Aalpha, y0alpha, balpha, calpha, YA);

				const double deltaec = ecf - ecp;
				const double deltaecd = ecfd - ecpd;
				const double beta = fdd * deltaec / eca - 1.;
				const double betad = fdd * fdd * (deltaecd * eca - ecad * deltaec) / (eca * eca);
				const double interp = fval / fdd * (1 + beta * pow(zeta, 4.));
				const double interpd = fval / fdd * pow(zeta, 4.) * betad;

				// TODO: something is wrong here, fix it!

				const double deriv = 1. / 3. * (ecpd // paramagnetic part
					+ ecad * interp);// +eca * interpd);

				res[i] = ep + (ef - ep) * fval // exchange term
					// correlation term:
					+ deriv;
			}

			return res;
		}
	};
}