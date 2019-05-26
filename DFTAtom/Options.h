#pragma once


#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/fileconf.h>

class Options
{
public:
	Options();

	void Load();
	void Save();

	void Open();
	void Close();

	int Z;
	int MultigridLevels;
	double MaxR;
	double deltaGrid;
	double alpha;

protected:
	wxFileConfig *m_fileconfig;
};

