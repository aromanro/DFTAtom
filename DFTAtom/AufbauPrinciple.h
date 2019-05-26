#pragma once

namespace DFT {

	struct Subshell
	{
		Subshell(int N = 0, int L = 0, int nrElectrons = 0) : m_N(N), m_L(L), m_nrElectrons(nrElectrons), E(0) {}


		bool operator<(const Subshell& o) const
		{
			return m_N < o.m_N || (m_N == o.m_N && m_L < o.m_L);
		}

		int m_N;
		int m_L;
		int m_nrElectrons;

		double E;
	};


	class AufbauPrinciple
	{
	public:

		static std::vector<Subshell> GetSubshells(int Z)
		{
			std::vector<Subshell> levels;
			bool exitLoops = false;

			int electronCount = 0;

			for (int NplusL = 0; !exitLoops && NplusL < 10; ++NplusL)
				for (int N = 0; N <= NplusL; ++N)
				{
					const int L = NplusL - N;

					if (L <= N)
					{
						int nrElectrons = 2 * (2 * L + 1);

						if (Z - electronCount < nrElectrons)
							nrElectrons = Z - electronCount;

						electronCount += nrElectrons;

						levels.emplace_back(Subshell(N, L, nrElectrons));

						if (electronCount == Z)
						{
							exitLoops = true;
							break;
						}
					}
				}

			std::cout << std::endl;

			return levels;
		}

	};

}




