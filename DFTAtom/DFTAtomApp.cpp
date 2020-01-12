#include "DFTAtomApp.h"
#include "DFTAtomFrame.h"


wxIMPLEMENT_APP(DFTAtomApp);

bool DFTAtomApp::OnInit()
{
	if (!wxApp::OnInit())
		return false;
	
	options.Load();

	frame = new DFTAtomFrame("DFTAtom", wxPoint(50, 50), wxSize(1024, 800));
	frame->Show(true);

	return true;
}
