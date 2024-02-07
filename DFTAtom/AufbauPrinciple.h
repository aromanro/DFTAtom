#pragma once

namespace DFT {

	struct Subshell
	{
		Subshell(int N = 0, int L = 0, int nrElectrons = 0) : m_N(N), m_L(L), m_nrElectrons(nrElectrons) {}


		bool operator<(const Subshell& o) const
		{
			return m_N < o.m_N || (m_N == o.m_N && m_L < o.m_L);
		}

		int m_N;
		int m_L;
		int m_nrElectrons;

		double E = 0;
	};


	class AufbauPrinciple
	{
	public:
		inline static int getMaxNrAlphaElectrons(int L)
		{
			return 2 * L + 1;
		}

		inline static int getMaxNrElectrons(int L)
		{
			return 2 * getMaxNrAlphaElectrons(L);
		}

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
						int nrElectrons = getMaxNrElectrons(L);

						AdjustForLanthanidesAndActinides(nrElectrons, Z, N, L);

						if (Z - electronCount < nrElectrons)
							nrElectrons = Z - electronCount;

						// exceptions lanthanides and actinides
						AdjustForLanthanidesAndActinides(nrElectrons, Z, N, L);

						if (nrElectrons > 0)
						{
							electronCount += nrElectrons;
							levels.emplace_back(Subshell(N, L, nrElectrons));
						}

						if (electronCount == Z)
						{
							exitLoops = true;
							break;
						}
					}
				}

			return levels;
		}

	private:
		static void AdjustForTransitionMetals(int& nrElectrons, int Z, int N, int L)
		{
			// exceptions transition metals
			if (TransitionMetalsException(Z, L)) // the electron is obtained from 's'
			{
				if (Z <= 29) // 3d gets 1 electron
				{
					if (3 == N) // from 4s
						--nrElectrons;
				}
				else if (Z <= 47) // 4d gets 1 electron
				{
					if (4 == N) // from 5s
						--nrElectrons;
				}
				else // 5d gets 1 electron
					if (5 == N) // from 6s
						--nrElectrons;
			}
			else if (PdException(Z, N, L)) // Pd is an exception to the above, 4d gets 2 electrons
				nrElectrons -= 2;
		}

		static void AdjustForLanthanidesAndActinides(int& nrElectrons, int Z, int N, int L)
		{
			if (3 == L) // f
			{
				if (LaCeGaException(Z, N)) // La, Ce and Ga, 4f loses one electron, it will go on 5d
					--nrElectrons;
				else if (4 == N)
				{
					if (89 == Z || 90 == Z) // Ac, 5f loses one electron, it will go in 6d, Th, 5f loses two electrons, they will go in 6d
						nrElectrons = 0;
					else if (PaUNpCdException(Z)) // Pa, 5f loses one electron (2 still remain, so not set to 0), it goes on 6d, U, similarly, loses 1, 3 remain, Np and Cm as for Uranium
						--nrElectrons;
				}
			}
			else if (LrException(Z, N, L)) // Lr, 6d loses the electron, goes into 7p
				nrElectrons = 0;
		}

		static bool TransitionMetalsException(int Z, int L)
		{
			return (24 == Z || 29 == Z || 41 == Z || 42 == Z || 44 == Z || 45 == Z || 47 == Z || 78 == Z || 79 == Z) && 0 == L;
		}

		static bool PdException(int Z, int N, int L)
		{
			return 46 == Z && 4 == N && 0 == L;
		}

		static bool LaCeGaException(int Z, int N)
		{
			return (57 == Z || 58 == Z || 64 == Z) && 3 == N;
		}

		static bool PaUNpCdException(int Z)
		{
			return 91 == Z || 92 == Z || 93 == Z || 96 == Z;
		}

		static bool LrException(int Z, int N, int L)
		{
			return 103 == Z && 5 == N && 2 == L;
		}
	};
}
