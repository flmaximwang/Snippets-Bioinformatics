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

#include "Common.h"
#include "CreateOptions.h"
#include "QueryStruct.h"
#include "TargetStruct.h"

void parseCmdLine(int, char**, CreateOptions &);

int main(int argc, char *argv[])
{
	int i;
	fstream ofs;
	FILE* ifp;
	CreateOptions copts;
	parseCmdLine(argc, argv, copts);

	if ("query" == copts.getPdsType())
	{
		for (i = 0; i < copts.getPdbFiles().size(); i++)
		{
			QueryStruct qs(copts.getPdbFiles()[i], copts.getLegalAA(), copts.getDistCut());
			openFileCPP(ofs, copts.getPdsFiles()[i], copts.getBinary() ? (ios::out | ios::trunc | ios::binary) : (ios::out | ios::trunc));
			qs.writeQueryStructFile(ofs, copts);
			ofs.close();
			openFileC(ifp, copts.getPdsFiles()[i], copts.getBinary() ? "rb" : "r");
			qs.checkQueryStructFile(ifp, copts.getPdsFiles()[i], copts);
			fclose(ifp);
			if (copts.getPostPdbFiles().size() > 0)
			{
				Structure& sys = qs.getProteinSys();
				sys.writePDB(copts.getPostPdbFiles()[i]);
			}
		}
	}
	else
	{
		if ("target" == copts.getPdsType())
		{
			for (i = 0; i < copts.getPdbFiles().size(); i++)
			{
				TargetStruct ts(copts.getPdbFiles()[i], copts.getLegalAA(), copts.getNonRed(), copts.getSeqIdenThresh(), copts.getOptRmsdThresh(), copts.getContRmsdThresh(), copts.getNonRedPdb(), copts.getPhiStep(), copts.getPsiStep(), copts.getDistCut(), copts.getDistStep());
				openFileCPP(ofs, copts.getPdsFiles()[i], copts.getBinary() ? (ios::out | ios::trunc | ios::binary) : (ios::out | ios::trunc));
				ts.writeTargetStructFile(ofs, copts);
				ofs.close();
				openFileC(ifp, copts.getPdsFiles()[i], copts.getBinary() ? "rb" : "r");
				ts.checkTargetStructFile(ifp, copts.getPdsFiles()[i], copts);
				fclose(ifp);
				if (copts.getPostPdbFiles().size() > 0)
				{
					Structure& sys = ts.getProteinSys();
					sys.writePDB(copts.getPostPdbFiles()[i]);
				}
			}
		}
		else
		{
			copts.usage();
			error("bad output PDS type value");
		}
	}

	return 0;
}

void parseCmdLine(int argc, char** argv, CreateOptions & copts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
    	static struct option opts[] =
    	{
			{"pdb", 1, 0, 1},
			{"pdbList", 1, 0, 2},
			{"pds", 1, 0, 3},
			{"pdsList", 1, 0, 4},
			{"phiStep", 1, 0, 5},
	  		{"psiStep", 1, 0, 6},
	  		{"dCut", 1, 0, 7},
	  		{"dStep", 1, 0, 8},
			{"type", 1, 0, 9},
			{"opdb", 1, 0, 10},
			{"opdbList", 1, 0, 11},
			{"nr", 0, 0, 12},
			{"gRMSD", 1, 0, 13},
			{"lRMSD", 1, 0, 14},
			{"seqID", 1, 0, 15},
			{"nrPDB", 0, 0, 16},
			{"cleanPDB", 0, 0, 17},
			{"onlyNat", 0, 0, 18},
			{0, 0, 0, 0}
    	};

    	int c = getopt_long (argc, argv, "", opts, &oind);
    	if (c == -1)
    	{
			break;
    	}

    	switch (c) {
			case 1:				
				copts.setPdbFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 2:
				copts.setPdbFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 3:
				copts.setPdsFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 4:
				copts.setPdsFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;
			
      		case 5:
			  	copts.setPhiStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 6:
		        copts.setPsiStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 7:
		        copts.setDistCut(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 8:
			  	copts.setDistStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 9:
				copts.setPdsType(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;
			
			case 10:
				copts.setPostPdbFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 11:
				copts.setPostPdbFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 12:
				copts.setNonRed(true);
				spec[string(opts[oind].name)] = true;
				break;

			case 13:
		        copts.setContRmsdThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 14:
		        copts.setOptRmsdThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 15:
		        copts.setSeqIdenThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 16:
				copts.setNonRedPdb(true);
				spec[string(opts[oind].name)] = true;
				break;

			case 17:
				spec[string(opts[oind].name)] = true;
				break;

			case 18:
				spec[string(opts[oind].name)] = true;
				break;

      		case '?':
				copts.usage();
				exit(-1);

      		default:
        		printf ("?? getopt returned character code %d ??\n", c);
				copts.usage();
				exit(-1);
    	}
	}

  	// make sure all required options have been specified
  	if (!(((spec.find(string("pdb")) != spec.end()) || (spec.find(string("pdbList")) != spec.end()))
		&& (spec.find(string("type")) != spec.end())))
  	{
		copts.usage();
		error("not all required options specified");
  	}

	if (spec.find(string("cleanPDB")) != spec.end()) {
		copts.setLegalAA(spec.find(string("onlyNat")) != spec.end());
	} else if (spec.find(string("onlyNat")) != spec.end()) {
		warning("'--onlyNat' does not work without '--cleanPDB'");
	}

	if (optind < argc)
	{
		copts.usage();
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
    	{
			printf ("%s ", argv[optind++]);
    	}
		printf ("\n");
		exit(-1);
	}

  	int i;
	if (0 == copts.getPdsFiles().size())
	{
		for (i = 0; i < copts.getPdbFiles().size(); i++)
		{
			copts.setPdsFile(fileName(copts.getPdbFiles()[i], false, true) + copts.getFileExt(), false);
		}
	}
	ASSERT(copts.getPdsFiles().size() == copts.getPdbFiles().size(), "numbers of input PDB files and output PDS files not consistent");
	if (copts.getPostPdbFiles().size() > 0)
  	{
  		ASSERT(copts.getPostPdbFiles().size() == copts.getPdbFiles().size(), "numbers of input PDB files and output post-processed PDB files not consistent");
		for (i = 0; i < copts.getPostPdbFiles().size(); i++)
		{
			copts.getPostPdbFiles()[i] = fileName(copts.getPostPdbFiles()[i], false, true) + ".pdb";
		}
	}

#if defined(DEBUG_CMDLINE)
	cout << "output binary file? " << copts.getBinary() << endl;
	cout << "contact RMSD threshold: " << copts.getContRmsdThresh() << endl;
	cout << "distance cutoff: " << copts.getDistCut() << endl;
	cout << "distance bin size: " << copts.getDistStep() << endl;
	cout << "output PDS file extension: " << copts.getFileExt() << endl;
	cout << "optimal RMSD threshold: " << copts.getOptRmsdThresh() << endl;
	cout << "phi bin size: " << copts.getPhiStep() << endl;
	cout << "psi bin size: " << copts.getPsiStep() << endl;
	cout << "separator: " << copts.getWordSep() << endl;
	cout << "sequence identity threshold: " << copts.getSeqIdenThresh() << endl;	
	cout << "non-redundant? " << copts.getNonRed() << endl;
	cout << "output non-redundant neighbor PDBs? " << copts.getNonRedPdb() << endl;
	cout << "terminator: " << copts.getWordTer() << endl;
	cout << "output PDS type: " << copts.getPdsType() << endl;
	exit(-1);
#endif
}
