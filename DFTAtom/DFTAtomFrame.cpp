#include "DFTAtomFrame.h"

#include <string>
#include <thread>

#include "wx/aboutdlg.h"
#include "wx/statline.h"
#include "wx/generic/aboutdlgg.h"

#include "DFTAtomApp.h"
#include "OptionsFrame.h"

#include "DFTAtom.h"


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

	virtual int overflow(int c) override
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

private:
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
		app.options.Close();
		app.options = optionsFrame.options;
		app.options.Open();
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
		
		//DFT::DFTAtom::CalculateUniform(options.Z, options.MultigridLevels, options.alpha, options.MaxR);
		DFT::DFTAtom::CalculateNonUniform(options.Z, options.MultigridLevels, options.alpha, options.MaxR, options.deltaGrid);

		inExecution = false;
	}).detach();
}

void DFTAtomFrame::OnUpdateExecute(wxUpdateUIEvent& event)
{
	event.Enable(!inExecution);
}

