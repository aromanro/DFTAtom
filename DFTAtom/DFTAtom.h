#pragma once

namespace DFT {

	class DFTAtom
	{
	public:
		static void Calculate(int Z, int MultigridLevels, double alpha, double MaxR, double deltaGrid);
	};

}

