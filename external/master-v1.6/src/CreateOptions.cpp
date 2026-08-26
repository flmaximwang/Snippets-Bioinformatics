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

#include "CreateOptions.h"

CreateOptions::CreateOptions()
{
	_bin = true;
	_contRmsdThresh = 2.0; //--gRMSD
	_dcut = 25.0;
	_dstep = 5.0;
	_ext = ".pds";
	_legalAA.clear();
	_nr = false; // --nr
	_nrPdb = false; // --nrPDB
	_optRmsdThresh = 1.0; //--lRMSD
	_phistep = 10.0;
	_psistep = 10.0;
	_sep = "";
	_seqIdenThresh = 0.9; // --seqID
	_ter = "";
}

void CreateOptions::setContRmsdThresh(const char* s)
{
	if (!((sscanf(s, "%lf", &(_contRmsdThresh)) == 1) && (_contRmsdThresh >= 0.0)))
	{
		usage();
		error("bad contact RMSD threshold value");
	}
}

void CreateOptions::setDistCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_dcut)) == 1) && (_dcut > 0.0)))
	{
		usage();
		error("bad distance cutoff value");
	}
}

void CreateOptions::setDistStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_dstep)) == 1) && (_dstep > 0.0)))
	{
		usage();
		error("bad distance step value");
	}
}

void CreateOptions::setLegalAA(const bool & onlyNat)
{
	// legal residue names that are considered as "protein" here
	_legalAA.clear();
	_legalAA.push_back("ALA");
	_legalAA.push_back("ARG");
	_legalAA.push_back("ASN");
	_legalAA.push_back("ASP");
	_legalAA.push_back("CYS");
	_legalAA.push_back("GLN");
	_legalAA.push_back("GLU");
	_legalAA.push_back("GLY");
	_legalAA.push_back("HIS"); _legalAA.push_back("HSC"); _legalAA.push_back("HSD"); _legalAA.push_back("HSE"); _legalAA.push_back("HSP");
	_legalAA.push_back("ILE");
	_legalAA.push_back("LEU");
	_legalAA.push_back("LYS");
	_legalAA.push_back("MET"); _legalAA.push_back("MSE");
	_legalAA.push_back("PHE");
  	_legalAA.push_back("PRO");
	_legalAA.push_back("SER");
	_legalAA.push_back("THR");
	_legalAA.push_back("TRP");
	_legalAA.push_back("TYR");
	_legalAA.push_back("VAL");

	if (!onlyNat) {
		_legalAA.push_back("CSO");
		_legalAA.push_back("HIP");
		_legalAA.push_back("PTR");
		_legalAA.push_back("SEC");
		_legalAA.push_back("SEP");
		_legalAA.push_back("TPO");

#if defined(DEBUG_ONLYNAT)
		cout << "Unnatural amino acids included.\n";
		exit(-1);
#endif
	}

#if defined(DEBUG_ONLYNAT)
	cout << "Only natural amino acids.\n";
	exit(-1);
#endif
}

void CreateOptions::setOptRmsdThresh(const char * s)
{
	if (!((sscanf(s, "%lf", &(_optRmsdThresh)) == 1) && (_optRmsdThresh >= 0.0)))
	{
		usage();
		error("bad optimal RMSD threshold value");
	}
}

void CreateOptions::setPdbFile(const string & fn)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("input PDB files already exist");
	}
	_pdbfnames.push_back(fn);
}

void CreateOptions::setPdbFiles(const string & list)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("input PDB files already exist");
	}
	file2array(list, _pdbfnames);
}

void CreateOptions::setPdsFile(const string & fn, const bool & chk)
{
	if (chk)
	{
		if (_pdsfnames.size() != 0)
		{
			usage();
			error("output PDS files already exist");
		}
	}
	_pdsfnames.push_back(fn);
}

void CreateOptions::setPdsFiles(const string & list)
{
	if (_pdsfnames.size() != 0)
	{
		usage();
		error("output PDS files already exist");
	}
	file2array(list, _pdsfnames);
}

void CreateOptions::setPdsType(const string & t)
{
	_type = t;
	if (!((_type == "query") || (_type == "target")))
	{
		usage();
		error("bad output PDS type value");
	}
}

void CreateOptions::setPhiStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_phistep)) == 1) && (_phistep > 0.0)))
	{
		usage();
		error("bad phi step value");
	}
}

void CreateOptions::setPostPdbFile(const string & fn)
{
	if (_opdbfnames.size() != 0)
	{
		usage();
		error("output post-processed PDB files already exist");
	}
	_opdbfnames.push_back(fn);
}

void CreateOptions::setPostPdbFiles(const string & list)
{
	if (_opdbfnames.size() != 0)
	{
		usage();
		error("output post-processed PDB files already exist");
	}
	file2array(list, _opdbfnames);
}

void CreateOptions::setPsiStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_psistep)) == 1) && (_psistep > 0.0)))
	{
		usage();
		error("bad psi step value");
	}
}

void CreateOptions::setSeqIdenThresh(const char * s)
{
	if (!((sscanf(s, "%lf", &(_seqIdenThresh)) == 1) && (_seqIdenThresh >= 0.0) && (_seqIdenThresh <= 1.0)))
	{
		usage();
		error("bad sequence identity threshold value");
	}
}

void CreateOptions::usage()
{
	int w = 100, p1 = 3, p2 = p1 + 16; // the length of the longest option name + 3
	string description = 
	    "This program converts PDB file(s) to PDS file(s). Note, PDS file formats differ based on whether the structures are intended as queries "
	    "(i.e., --type query) or database entries to be searched against (--type target). If you use this tool for your research, please cite:"
        "\n\nZhou, G. Grigoryan, Protein Science, 24(4): 508-524, 2015";
    cout << endl << optionUsage("", suiteName(), w, 0, 0) << endl << endl;
    cout << optionUsage("", description, w, 0, 0) << "\n\n";
    cout << optionUsage("Input/output options:", "", w, 0, 0) << "\n";
	cout << optionUsage("--type", "whether input PDB structures are intended as future queries (query) or database structures (target).", w, p1, p2) << endl;
	cout << optionUsage("--pdb", "input PDB file.", w, p1, p2) << endl;
	cout << optionUsage("--pdbList", "file with a list of paths to input PDB files, one per line. Either --pdb or --pdbList must be given.", w, p1, p2) << endl;
	cout << optionUsage("--pds", "optional: output PDS file name. By default, takes the base name of the input PDB file and appends \".pds\" to it.", w, p1, p2) << endl;
	cout << optionUsage("--pdsList", "optional: a file with a list of output PDS file names (one per input PDB file).", w, p1, p2) << endl;
	cout << optionUsage("--opdb", "optional: name of output post-processed PDB file (useful for keeping track of how the various PDB weirdnesses got parsed).", w, p1, p2) << endl;
	cout << optionUsage("--opdbList", "optional: a file with a list of file names for post-processed PDBs (one per input PDB file).", w, p1, p2) << endl;
    cout << optionUsage("--cleanPDB", "optional: admit only amino acids HSC, HSD, HSE, HSP, MSE, CSO, HIP, PTR, SEC, SEP, and TPO in addition to the 20 natural ones, "
									  "and skip any residue missing backbone atom(s). By default, admit all residues and stop upon encountering any residue missing backbone atom(s).", w, p1, p2) << endl;

    cout << optionUsage("\nStructural redundancy removal options:", "", w, 0, 0) << "\n";
	cout << optionUsage("--nr", "optional: remove redundancy in the input structure prior to writing the PDS. The non-redundant part of the structure is defined "
	                            "in terms of non-redundant chains and inter-chain interfaces; see paper. Briefly, two chains are considered redundant if they "
	                            "are nearly identical in sequence (see --seqID below) and in structure (see --lRMSD below), and their chain neighborhoods (i.e., "
	                            "other chains making contacts with the two being compared) are also nearly identical in sequence and structure, both individually "
	                            "and as a group (see --gRMSD below).", w, p1, p2) << endl;
    cout << optionUsage("--seqID", "chain-to-chain sequence identity threshold for considering two chains to be redundant with --nr; default is 0.9 (i.e., 90%).", w, p1, p2) << endl;
    cout << optionUsage("--lRMSD", "chain-to-chain RMSD threshold for considering two chains to be redundant with --nr; default is 1.0.", w, p1, p2) << endl;
    cout << optionUsage("--gRMSD", "RMSD threshold for neighboring chains upon superimposing the central chain, used with --nr; default is 2.0.", w, p1, p2) << endl;

    cout << optionUsage("\nExpert data-structure options (changes not recommended):", "", w, 0, 0) << "\n";
    cout << optionUsage("--dCut", "upper limit for tabulating CA-CA distances (must be identical between queries and targets); default is 25.0.", w, p1, p2) << endl;
	cout << optionUsage("--dStep", "bin size for tabulating CA-CA distances (only matters for target structures); default is 5.0.", w, p1, p2) << endl;
	cout << optionUsage("--phiStep", "bin size for tabulating phi angles (only matters for target structures); default is 10.0.", w, p1, p2) << endl;
	cout << optionUsage("--psiStep", "bin size for tabulating psi angles (only matters for target structures); default is 10.0.", w, p1, p2) << endl;
//	cout << optionUsage("--nrPDB", "if specified, dump PDB file of each unique neighborhood for large structures. Only works when --nr has been given.", w, p1, p2) << endl;
	
//	cout << optionUsage("--binary", "output binary file? 1 for binary (default), 0 for text.", w, p1, p2) << endl;	  
	cout << endl;
}

