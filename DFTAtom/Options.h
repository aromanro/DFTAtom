#pragma once


#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/fileconf.h>

class Options
{
public:
	Options();

	// avoid double deletion of m_fileconfig at destruction if copied
	Options(const Options& other)
		:
		Z(other.Z),
		MultigridLevels(other.MultigridLevels),
		MaxR(other.MaxR),
		deltaGrid(other.deltaGrid),
		alpha(other.alpha),
		m_fileconfig(nullptr)
	{
	}

	Options& operator=(const Options& other)
	{
		Z = other.Z;
		MultigridLevels = other.MultigridLevels;
		MaxR = other.MaxR;
		deltaGrid = other.deltaGrid;
		alpha = other.alpha;
		m_fileconfig = nullptr;

		return *this;
	}


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

