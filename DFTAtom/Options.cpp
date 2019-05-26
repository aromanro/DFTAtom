#include "Options.h"

#include <wx/stdpaths.h> 

Options::Options()
	: Z(36), MultigridLevels(12), MaxR(10.), deltaGrid(0.001), alpha(0.5),
	m_fileconfig(nullptr)
{
}


void Options::Open()
{
	if (m_fileconfig) return;

	wxString dir = wxStandardPaths::Get().GetConfigDir() + wxFileName::GetPathSeparator();

	if(!wxFileName::DirExists(dir))
		wxFileName::Mkdir(dir, 0777, wxPATH_MKDIR_FULL);

	wxString iniFilePath = dir + "DFTAtom.ini";

	m_fileconfig = new wxFileConfig("DFTAtom", wxEmptyString, iniFilePath);

	wxConfigBase::Set(m_fileconfig);
}


void Options::Close()
{
	delete m_fileconfig;
	m_fileconfig = NULL;
	wxConfigBase::Set(NULL);
}

void Options::Load()
{
	wxConfigBase *conf = wxConfigBase::Get(false);
	if (conf)
	{		
		Z = conf->ReadLong("/Z", 36);
		MultigridLevels = conf->ReadLong("/MultigridLevels", 12);

		MaxR = conf->ReadDouble("/MaxR", 10.);
		deltaGrid = conf->ReadDouble("/deltaGrid", 0.001);
		alpha = conf->ReadDouble("/alpha", 0.5);
	}
}

void Options::Save()
{
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{		
		conf->Write("/Z", static_cast<long int>(Z));
		conf->Write("/MultigridLevels", static_cast<long int>(MultigridLevels));
		
		conf->Write("/MaxR", MaxR);
		conf->Write("/deltaGrid", deltaGrid);
		conf->Write("/alpha", alpha);
	}

	if (m_fileconfig)
		m_fileconfig->Flush();
}
