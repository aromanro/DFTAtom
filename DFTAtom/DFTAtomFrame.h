#pragma once

#define wxNEEDS_DECL_BEFORE_TEMPLATE
#define _MATH_DEFINES_DEFINED

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
	#pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers
#ifndef WX_PRECOMP
	#include "wx/wx.h"
#endif


#include <wx/richtext/richtextctrl.h>

#include <atomic>
#include <mutex>

class DFTAtomFrame : public wxFrame
{
public:
	DFTAtomFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
	~DFTAtomFrame();

private:
	std::atomic_bool inExecution;

	std::string bufferStr;
	std::mutex bufferStrMutex;
	
	wxTimer timer;
	
	void OnExit(wxCommandEvent& event);
	void OnOptions(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnExecute(wxCommandEvent& event);
	void OnUpdateExecute(wxUpdateUIEvent& event);
	void OnTimer(wxTimerEvent& event);

	wxRichTextCtrl* richTextCtrl;

	wxDECLARE_EVENT_TABLE();
};

