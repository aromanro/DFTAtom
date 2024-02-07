#pragma once


#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/fileconf.h>

class Options
{
public:
	Options();

	~Options()
	{
		delete m_fileconfig;
	}

	// avoid double deletion of m_fileconfig at destruction if copied
	Options(const Options& other)
		:
		Z(other.Z),
		MultigridLevels(other.MultigridLevels),
		MaxR(other.MaxR),
		deltaGrid(other.deltaGrid),
		alpha(other.alpha),
		method(other.method)
	{
		m_fileconfig = nullptr;
	}

	Options& operator=(const Options& other)
	{
		Z = other.Z;
		MultigridLevels = other.MultigridLevels;
		MaxR = other.MaxR;
		deltaGrid = other.deltaGrid;
		alpha = other.alpha;
		method = other.method;
		m_fileconfig = nullptr;

		return *this;
	}


	void Load();
	void Save();

	int Z;
	int MultigridLevels;
	double MaxR;
	double deltaGrid;
	double alpha;

	int method; // 0 LDA, 1 LSDA

protected:
	void Open();
	void Close();

	wxFileConfig *m_fileconfig = nullptr;
};

