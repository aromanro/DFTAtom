#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/statline.h>
#include <wx/valgen.h>
#include <wx/valnum.h>
#include <wx/panel.h>
#include <wx/bookctrl.h>

#include "DFTAtomApp.h"
#include "OptionsFrame.h"
#include "DFTAtomFrame.h"

#define ID_Z 101
#define ID_MULTIGRID 102
#define ID_R 103
#define ID_DELTA 104
#define ID_ALPHA 105

wxDECLARE_APP(DFTAtomApp);


wxIMPLEMENT_CLASS(OptionsFrame, wxDialog);

wxBEGIN_EVENT_TABLE(OptionsFrame, wxDialog)
EVT_CLOSE(OptionsFrame::OnClose)
wxEND_EVENT_TABLE()

OptionsFrame::OptionsFrame(const Options& opt, const wxString & title, wxWindow* parent)
	: wxDialog(parent, wxID_ANY, title, wxDefaultPosition, wxSize(400, 250)) 
{
	options = opt;
	CreateControls();
	Layout();
	Centre();
}

bool OptionsFrame::TransferDataFromWindow()
{
	if (!wxDialog::TransferDataFromWindow()) return false;
	
	wxTextCtrl* MultilevelCtrl = (wxTextCtrl*)FindWindow(ID_MULTIGRID);
	wxString str = MultilevelCtrl->GetValue();
	long int val = 0;
	if (!str.ToLong(&val)) return false;
	if (val < 10 || val > 20)
	{
		wxMessageBox("Please enter between 10 and 20 levels", "Validation", wxOK | wxICON_INFORMATION, this);

		return false; 
	}

	return true;
}

void OptionsFrame::OnClose(wxCloseEvent& event)
{
	event.Skip();
}


void OptionsFrame::CreateControls()
{
	// box to contain them all
	wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);
	SetSizer(vbox);	

	// box with margin to contain option controls
	wxBoxSizer* boxSizer = new wxBoxSizer(wxVERTICAL);
	vbox->Add(boxSizer, 0, wxALIGN_CENTER_HORIZONTAL| wxGROW | wxALL, 5);

	// *****************************************************************
	// Controls

	wxBoxSizer* box = new wxBoxSizer(wxHORIZONTAL);
	boxSizer->Add(box, 0, wxGROW|wxALL, 5);

	wxStaticText* label = new wxStaticText(this, wxID_STATIC, "&Z:", wxDefaultPosition, wxSize(60, -1), wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL, 5);
	
	wxString str = wxString::Format(wxT("%i"), options.Z);
	wxTextCtrl* ZCtrl = new wxTextCtrl(this, ID_Z, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(ZCtrl, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	box->Add(5, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5); // pushes to the right

	label = new wxStaticText(this, wxID_STATIC, "&Multigrid levels:", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

	str = wxString::Format(wxT("%i"), options.MultigridLevels);
	wxTextCtrl* MultigridCtrl = new wxTextCtrl(this, ID_MULTIGRID, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(MultigridCtrl, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);


	box = new wxBoxSizer(wxHORIZONTAL);
	boxSizer->Add(box, 0, wxGROW|wxALL, 5);

	label = new wxStaticText(this, wxID_STATIC, "Max &R:", wxDefaultPosition, wxSize(60, -1), wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

	str = wxString::Format(wxT("%f"), options.MaxR);
	wxTextCtrl* RCtrl = new wxTextCtrl(this, ID_R, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(RCtrl, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	box->Add(5, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5); // pushes to the right

	label = new wxStaticText(this, wxID_STATIC, "&Delta grid:", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

	str = wxString::Format(wxT("%f"), options.deltaGrid);
	wxTextCtrl* DeltaCtrl = new wxTextCtrl(this, ID_DELTA, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(DeltaCtrl, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);


	box = new wxBoxSizer(wxHORIZONTAL);
	boxSizer->Add(box, 0, wxGROW|wxALL, 5);

	label = new wxStaticText(this, wxID_STATIC, "Mi&xing:", wxDefaultPosition, wxSize(60, -1), wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

	str = wxString::Format(wxT("%f"), options.alpha);
	wxTextCtrl* AlphaCtrl = new wxTextCtrl(this, ID_ALPHA, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(AlphaCtrl, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	// *****************************************************************
	// Validators

	wxIntegerValidator<int> val1(&options.Z, wxNUM_VAL_DEFAULT);
	val1.SetRange(1, 118);
	ZCtrl->SetValidator(val1);

	wxIntegerValidator<int> val2(&options.MultigridLevels, wxNUM_VAL_DEFAULT);
	val2.SetRange(1, 20);
	MultigridCtrl->SetValidator(val2);

	wxFloatingPointValidator<double> dblVal1(&options.MaxR);
	dblVal1.SetRange(1, 30);
	dblVal1.SetPrecision(2);
	RCtrl->SetValidator(dblVal1);

	wxFloatingPointValidator<double> dblVal2(&options.deltaGrid);
	dblVal2.SetRange(0, 1);
	dblVal2.SetPrecision(8);
	DeltaCtrl->SetValidator(dblVal2);

	wxFloatingPointValidator<double> dblVal3(&options.alpha);
	dblVal3.SetRange(0, 1);
	dblVal3.SetPrecision(2);
	AlphaCtrl->SetValidator(dblVal3);


	// *****************************************************************
	
	wxStaticLine* line = new wxStaticLine(this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL);
	boxSizer->Add(line, 0, wxGROW|wxALL, 5);

	// bottom box with ok & cancel buttons
	wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
	
	wxButton *okButton = new wxButton(this, wxID_OK, "Ok", wxDefaultPosition, wxSize(70, 30));
	wxButton *closeButton = new wxButton(this, wxID_CANCEL, "Cancel", wxDefaultPosition, wxSize(70, 30));

	hbox->Add(okButton, 1);
	hbox->Add(closeButton, 1, wxLEFT, 5);

	vbox->Add(hbox, 0, wxALIGN_CENTER | wxTOP | wxBOTTOM, 10);
}
