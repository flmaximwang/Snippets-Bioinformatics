/*
----------------------------------------------------------------------------
This file is part of MASTER.

MASTER is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

MASTER is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MASTER.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014 Jianfu Zhou, Gevorg Grigoryan
----------------------------------------------------------------------------
*/

#include "ParseOptions.h"

void ParseOptions::setPdbFile(const string & fn)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("output PDB files already exist");
	}
	_pdbfnames.push_back(fn);
}

void ParseOptions::setPdbFiles(const string & list)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("output PDB files already exist");
	}
	file2array(list, _pdbfnames);
}

void ParseOptions::setPdsFile(const string & fn)
{
	if (_pdsfnames.size() != 0)
	{
		usage();
		error("input PDS files already exist");
	}
	_pdsfnames.push_back(fn);
}

void ParseOptions::setPdsFiles(const string & list)
{
	if (_pdsfnames.size() != 0)
	{
		usage();
		error("input PDS files already exist");
	}
	file2array(list, _pdsfnames);
}

void ParseOptions::usage()
{
	int w = 100, p1 = 3, p2 = p1 + 12; // based on the length of the longest option name
    string description = 
        "Parses and extracts various protein data from PDS file(s). If you use this tool for your research, please cite:"
        "\n\nZhou, G. Grigoryan, Protein Science, 24(4): 508-524, 2015";
    cout << endl << optionUsage("", suiteName(), w, 0, 0) << endl << endl;
    cout << optionUsage("", description, w, 0, 0) << "\n\n";
    cout << optionUsage("Input options. Either --pds or --pdsList must be given:", "", w, 0, 0) << "\n";
	cout << optionUsage("--pds", "input PDS file.", w, p1, p2) << endl;
	cout << optionUsage("--pdsList", "list of input PDS files, one per line.", w, p1, p2) << endl;

    cout << optionUsage("", "\nOutput options. At least one of --pdb, --pdbList, --seq, --brk, and --nr must be given.", w, 0, 0) << "\n";
	cout << optionUsage("--pdb", "write PDB information into this file name; use this if --pds was given.", w, p1, p2) << endl;
	cout << optionUsage("--pdbList", "a file with a list of output PDB file names (one per input PDS file); use this if --pdsList was given.", w, p1, p2) << endl;
	cout << optionUsage("--seq", "file name for storing the sequences of input PDS file(s); one per line. If 'stdout' is given, will print to screen instead of a file.", w, p1, p2) << endl;
	cout << optionUsage("--brk", "file name for storing the locations of chain breaks of input PDS file(s); one per line. Chain break locations are given as "
	                             "zero-initiated indices corresponding to the residue right before the break. If is 'stdout', will print to screen instead of a file.", w, p1, p2) << endl;
	cout << optionUsage("--nr", "file name for writing the number of residues of input PDS files (one per line). If is 'stdout', will print to screen instead of file.", w, p1, p2) << endl;
}
