#include "DFTAtomFrame.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>

#include "PoissonSolver.h"

#include "Integral.h"
#include "AufbauPrinciple.h"
#include "Numerov.h"
#include "VWNExcCor.h"


#include "wx/aboutdlg.h"
#include "wx/statline.h"
#include "wx/generic/aboutdlgg.h"

#include "DFTAtomApp.h"
#include "OptionsFrame.h"


#define ID_TIMER   101
#define ID_EXECUTE 102

wxDECLARE_APP(DFTAtomApp);

wxBEGIN_EVENT_TABLE(DFTAtomFrame, wxFrame)
EVT_MENU(wxID_EXIT, DFTAtomFrame::OnExit)
EVT_MENU(wxID_PREFERENCES, DFTAtomFrame::OnOptions)
EVT_MENU(wxID_ABOUT, DFTAtomFrame::OnAbout)
EVT_MENU(ID_EXECUTE, DFTAtomFrame::OnExecute)
EVT_UPDATE_UI(ID_EXECUTE, DFTAtomFrame::OnUpdateExecute)
EVT_TIMER(ID_TIMER, DFTAtomFrame::OnTimer)
wxEND_EVENT_TABLE()


class MyStream : public std::ostream, std::streambuf
{
public:
	MyStream(std::string& bufStr, std::mutex& strMutex) : std::ostream(this), m_bufStr(bufStr), m_strMutex(strMutex) {}

	int overflow(int c)
	{
		log(c);

		return 0;
	}

private:
	void log(char c)
	{
		std::lock_guard<std::mutex> lock(m_strMutex);
		m_bufStr += c;
	}

	std::string& m_bufStr;
	std::mutex& m_strMutex;
};


class RedirectStream
{
public:
	RedirectStream(std::ostream& old, std::ostream& dst)
		: s(old), backupbuf(old.rdbuf())
	{
		s.rdbuf(dst.rdbuf());
	}

	~RedirectStream()
	{
		s.rdbuf(backupbuf);
	}

protected:
	std::ostream& s;
	std::streambuf* const backupbuf;
};





DFTAtomFrame::DFTAtomFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame(NULL, wxID_ANY, title, pos, size), inExecution(false), timer(this, ID_TIMER)
{
	wxMenu *menuFile = new wxMenu;

	menuFile->Append(ID_EXECUTE, "Execute");
	menuFile->Append(wxID_EXIT);

	wxMenu *menuView = new wxMenu;
	menuView->Append(wxID_PREFERENCES);

	wxMenu *menuHelp = new wxMenu;
	menuHelp->Append(wxID_ABOUT);

	wxMenuBar *menuBar = new wxMenuBar;
	menuBar->Append(menuFile, "&File");
	menuBar->Append(menuView, "&View");
	menuBar->Append(menuHelp, "&Help");

	SetMenuBar(menuBar);

	CreateStatusBar();
	SetStatusText("Welcome to DFTAtom!");

	richTextCtrl = new wxRichTextCtrl(this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxVSCROLL | wxHSCROLL | wxBORDER_NONE /*| wxWANTS_CHARS*/ | wxTE_MULTILINE | wxTE_READONLY);

	wxFont font(14, wxROMAN, wxNORMAL, wxNORMAL);

	richTextCtrl->SetFont(font);
	richTextCtrl->GetCaret()->Hide();

	Layout();	
}


DFTAtomFrame::~DFTAtomFrame()
{
}


void DFTAtomFrame::OnOptions(wxCommandEvent& WXUNUSED(event))
{
	DFTAtomApp& app = wxGetApp();
	OptionsFrame optionsFrame(app.options, "Options", this);
	
	if (wxID_OK == optionsFrame.ShowModal())
	{
		app.options = optionsFrame.options;
		app.options.Save();
	}
}


void DFTAtomFrame::OnExit(wxCommandEvent& WXUNUSED(event))
{
	Close(true);
}

void DFTAtomFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
	wxAboutDialogInfo info;

	info.SetName("DFTAtom");

	static const int majorVer = 1;
	static const int minorVer = 0;
	wxString verStr = wxString::Format("%d.%d", majorVer, minorVer);
	info.SetVersion(verStr,	wxString::Format("Version %s", verStr));

	info.SetDescription("   Density Functional Theory Atom Application   ");
	info.SetLicense("GNU GPL v3.0, see LICENSE file for details");

	info.AddDeveloper("Adrian Roman");

	info.SetWebSite("https://github.com/aromanro/DFTAtom", "GitHub repository");

	wxAboutBox(info, this);	
}

void DFTAtomFrame::OnTimer(wxTimerEvent& WXUNUSED(event))
{
	if (!inExecution)
	{
		timer.Stop();
		if (wxIsBusy()) wxEndBusyCursor();
	}

	{
		std::lock_guard<std::mutex> lock(bufferStrMutex);
		if (!bufferStr.empty())
		{
			richTextCtrl->WriteText(wxString(bufferStr));
			bufferStr.clear();
		}
	}
	richTextCtrl->ScrollIntoView(richTextCtrl->GetLastPosition(), WXK_END);
}

void DFTAtomFrame::OnExecute(wxCommandEvent& WXUNUSED(event))
{
	if (inExecution.exchange(true)) return;
	
	wxBeginBusyCursor();
	richTextCtrl->Clear();
	timer.Start(100);
	DFTAtomApp& app = wxGetApp();
	Options options = app.options;

	std::thread([this, options]()
	{
		MyStream myStream(bufferStr, bufferStrMutex);
		RedirectStream redirect(std::cout, myStream);

		const int Z = options.Z;
		const int MultigridLevels = options.MultigridLevels;

		const double alpha = options.alpha;
		const double oneMinusAlpha = 1. - alpha;

		const double MaxR = options.MaxR;
		const double deltaGrid = options.deltaGrid;

		const double energyErr = 1E-10;


		const int NumGridNodes = DFT::PoissonSolver::GetNumberOfNodes(MultigridLevels);

		const int NumSteps = NumGridNodes - 1;
		const double Rp = MaxR / (exp(NumSteps * deltaGrid) - 1.);


		std::vector<double> density(NumGridNodes);

		DFT::Potential potential;
		potential.m_potentialValues.resize(NumGridNodes);


		std::vector<DFT::Subshell> levels = DFT::AufbauPrinciple::GetSubshells(Z);
		std::sort(levels.begin(), levels.end());

		double Eold = 0;


		const double volume = 4. / 3. * M_PI * MaxR * MaxR * MaxR;
		const double constDens = Z / volume;

		density[0] = 0;
		for (int i = 1; i < density.size(); ++i)
			density[i] = constDens;


		DFT::PoissonSolver poissonSolver(MultigridLevels, deltaGrid);

		std::vector<double> UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
		std::vector<double> Vexc = DFT::VWNExchCor::Vexc(density);


		potential.m_potentialValues[0] = 0;
		for (int i = 1; i < NumGridNodes; ++i)
		{
			const double realPos = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = (-Z + UHartree[i]) / realPos + Vexc[i];
		}

		for (int sp = 0; sp < 500; ++sp)
		{
			std::cout << "Step: " << sp << std::endl;

			double Eelectronic = 0;

			std::vector<double> newDensity(density.size(), 0);

			DFT::Numerov<DFT::NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, MaxR, NumGridNodes);
			
			double BottomEnergy = -0.5 * Z * Z - 1;

			for (auto& level : levels)
			{
				const int NumNodes = level.m_N - level.m_L;

				double TopEnergy = 50;
				

				double toe = TopEnergy;
				double boe = BottomEnergy;
				double deltaEnergy = toe - boe;
				while (deltaEnergy > energyErr)
				{
					const double E = (toe + boe) / 2;

					int NumNodesCounted;
					numerov.SolveSchrodingerCountNodes(Z, NumSteps, level.m_L, E, NumSteps, NumNodes, NumNodesCounted);

					if (NumNodesCounted > NumNodes)
						toe = E;
					else
						boe = E;

					deltaEnergy = toe - boe;
				}
				toe -= energyErr;
				TopEnergy = toe;

				boe = BottomEnergy;
				deltaEnergy = toe - boe;
				while (deltaEnergy > energyErr)
				{
					const double E = (toe + boe) / 2;

					int NumNodesCounted;
					numerov.SolveSchrodingerCountNodes(Z, NumSteps, level.m_L, E, NumSteps, NumNodes, NumNodesCounted);

					if (NumNodesCounted < NumNodes)
						boe = E;
					else
						toe = E;

					deltaEnergy = toe - boe;
				}
				BottomEnergy = boe + energyErr;

				double delta = numerov.SolveSchrodingerSolutionInZero(Z, NumSteps, level.m_L, BottomEnergy, NumSteps);
				const bool sgnA = delta > 0;

				for (int i = 0; i < 100; ++i)
				{
					level.E = (TopEnergy + BottomEnergy) / 2;

					delta = numerov.SolveSchrodingerSolutionInZero(Z, NumSteps, level.m_L, level.E, NumSteps);

					if ((delta > 0) == sgnA)
						BottomEnergy = level.E;
					else
						TopEnergy = level.E;

					if (TopEnergy - BottomEnergy < energyErr) break;
				}

				// now really solve it
				level.E = (TopEnergy + BottomEnergy) / 2;
				BottomEnergy = level.E;
				std::vector<double> result = numerov.SolveSchrodingerSolutionCompletely(Z, NumSteps, level.m_L, level.E, NumSteps);

				// square the wavefunction
				// also integrate the square of wavefunction to get the normalization constant

				std::vector<double> result2(result.size());
				std::vector<double> result3(result.size());
				for (int i = 0; i < result.size(); ++i)
				{
					// if nonuniform, convert the function back!!!!!
					result[i] *= exp(i * deltaGrid * 0.5);

					result2[i] = result[i] * result[i];

					// if nonuniform, do the change for dr	
					const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
					result3[i] = result2[i] * cnst;
				}

				const double integralForSquare = DFT::Integral::SimpsonOneThird(1, result3); // for nonuniform case the step is 1

				std::cout << "Energy: " << std::setprecision(12) << level.E << " Num nodes: " << NumNodes << std::endl;

				for (int i = 0; i < result.size() - 1; ++i)
					newDensity[i] += level.m_nrElectrons * result2[i] / integralForSquare;

				Eelectronic += level.m_nrElectrons * level.E;
			}


			for (int i = 1; i < density.size(); ++i)
			{
				const double position = Rp * (exp(i * deltaGrid) - 1.);

				// 4 * M_PI appears because we're in spherical coordinates
				// the actual integration for the 'true' wavefunction gives a 4 M_PI
				// the radial wavefunction is actually u / r, whence also the division by position * position to get the true density

				newDensity[i] /= 4. * M_PI * position * position;
				density[i] = alpha * density[i] + oneMinusAlpha * newDensity[i];
			}


			UHartree = poissonSolver.SolvePoissonNonUniform(Z, MaxR, density);
			Vexc = DFT::VWNExchCor::Vexc(density);

			potential.m_potentialValues[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double position = Rp * (exp(i * deltaGrid) - 1.);

				potential.m_potentialValues[i] = (-Z + UHartree[i]) / position + Vexc[i];
			}

			// Nuclear energy:
			std::vector<double> nuclear(NumGridNodes);
			for (int i = 0; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				nuclear[i] = position * Z * density[i] * cnst;
			}
			const double Enuclear = -4 * M_PI * DFT::Integral::SimpsonOneThird(1, nuclear);

			// Exchange-correlation energy:
			std::vector<double> exccor(NumGridNodes);
			exccor[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				exccor[i] = position * position * density[i] * Vexc[i] * cnst;
			}
			double Exc = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, exccor);

			std::vector<double> eexcDeriv = DFT::VWNExchCor::eexcDif(density);
			exccor[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				exccor[i] = position * position * density[i] * eexcDeriv[i] * cnst;
			}
			const double eExcDif = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, exccor);
			Exc += eExcDif;

			// Hartree energy:
			std::vector<double> hartree(NumGridNodes);
			for (int i = 0; i < NumGridNodes; ++i)
			{
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);

				const double cnst = Rp * deltaGrid * expD;

				hartree[i] = position * density[i] * UHartree[i] * cnst;
			}
			const double Ehartree = -2 * M_PI * DFT::Integral::SimpsonOneThird(1, hartree);

			// potential energy:
			std::vector<double> potentiale(NumGridNodes);
			potentiale[0] = 0;
			for (int i = 1; i < NumGridNodes; ++i) {
				const double expD = exp(deltaGrid * i);
				const double position = Rp * (expD - 1.);
				const double cnst = Rp * deltaGrid * expD;

				potentiale[i] = position * position * density[i] * potential.m_potentialValues[i] * cnst;
			}
			const double Epotential = 4 * M_PI * DFT::Integral::SimpsonOneThird(1, potentiale);

			const double Ekinetic = Eelectronic - Epotential;
			const double Etotal = Eelectronic + Ehartree + eExcDif;

			std::cout << "Etotal = " << std::setprecision(12) << Etotal << " Ekin = " << std::setprecision(12) << Ekinetic << " Ecoul = " << std::setprecision(12) << -Ehartree << " Eenuc = " << std::setprecision(12) << Enuclear << " Exc = " << std::setprecision(12) << Exc << std::endl;

			if (abs(Eold - Etotal) < 1E-9)
			{
				std::cout << std::endl << "Finished!" << std::endl << std::endl;

				break;
			}
			Eold = Etotal;

			std::cout << "********************************************************************************" << std::endl;
		}

		const char orb[] = { 's', 'p', 'd', 'f' };
		for (const auto& level : levels)
			std::cout << level.m_N + 1 << orb[level.m_L] << level.m_nrElectrons << " ";

		inExecution = false;
	}).detach();
}

void DFTAtomFrame::OnUpdateExecute(wxUpdateUIEvent& event)
{
	event.Enable(!inExecution);
}

