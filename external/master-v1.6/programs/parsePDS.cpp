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
along with MaDCaT.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014 Jianfu Zhou, Gevorg Grigoryan             
----------------------------------------------------------------------------
*/

#include "Common.h"
#include "Parse.h"
#include "ParseOptions.h"

void parseCmdLine(int, char**, ParseOptions &);

int main(int argc, char *argv[])
{	
	ParseOptions popts;
	int i, numres;
	FILE *ifp = NULL;
	fstream ofsBrk, ofsNumRes, ofsSeq;
	string pdbfname;
	
	parseCmdLine(argc, argv, popts);

	if ((!popts.getBrkFile().empty()) && (popts.getBrkFile() != "stdout"))
	{
		openFileCPP(ofsBrk, popts.getBrkFile(), (ios::out | ios::trunc));
	}

	if ((!popts.getNumResFile().empty()) && (popts.getNumResFile() != "stdout"))
	{
		openFileCPP(ofsNumRes, popts.getNumResFile(), (ios::out | ios::trunc));
	}

	if ((!popts.getSeqFile().empty()) && (popts.getSeqFile() != "stdout"))
	{
		openFileCPP(ofsSeq, popts.getSeqFile(), (ios::out | ios::trunc));
	}

	for (i = 0; i < popts.getPdsFiles().size(); i++)
	{
		openFileC(ifp, popts.getPdsFiles()[i], "rb");
		Parse::extractNumRes(ifp, numres, popts.getPdsFiles()[i]);

		if ((!popts.getPdbFiles().empty()) || (!popts.getBrkFile().empty()))
		{
			if (!popts.getPdbFiles().empty())
			{
				pdbfname = popts.getPdbFiles()[i];
			}
			else
			{
				pdbfname = fileName(popts.getPdsFiles()[i], false, true) + ".tmp.pdb";
			}
			Parse::outPdb(ifp, numres, pdbfname, popts.getPdsFiles()[i]);

			if (!popts.getBrkFile().empty())
			{
				Parse::outBrk(ofsBrk, pdbfname);

				if (popts.getPdbFiles().empty())
				{
					remove(pdbfname.c_str());
				}
			}
		}

		if (!popts.getNumResFile().empty())
		{
			Parse::outNumRes(ofsNumRes, numres);
		}

		if (!popts.getSeqFile().empty())
		{
			Parse::outSeq(ifp, numres, ofsSeq, popts.getPdsFiles()[i]);
		}

		fclose(ifp);
	}

	if ((!popts.getBrkFile().empty()) && (popts.getBrkFile() != "stdout"))
	{
		ofsBrk.close();
	}

	if ((!popts.getNumResFile().empty()) && (popts.getNumResFile() != "stdout"))
	{
		ofsNumRes.close();
	}

	if ((!popts.getSeqFile().empty()) && (popts.getSeqFile() != "stdout"))
	{
		ofsSeq.close();
	}

	return 0;
}

void parseCmdLine(int argc, char** argv, ParseOptions & popts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
    	static struct option opts[] =
    	{
			{"pds", 1, 0, 1},
	  		{"pdsList", 1, 0, 2},
			{"pdb", 1, 0, 3},
			{"pdbList", 1, 0, 4},
			{"seq", 1, 0, 5},
			{"brk", 1, 0, 6},
			{"nr", 1, 0, 7},
      		{0, 0, 0, 0}
    	};

    	int c = getopt_long (argc, argv, "", opts, &oind);
    	if (c == -1)
    	{
			break;
    	}

    	switch (c) {
			case 1:
				popts.setPdsFile(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 2:
				popts.setPdsFiles(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 3:
				popts.setPdbFile(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 4:
				popts.setPdbFiles(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 5:
				popts.setSeqFile(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 6:
				popts.setBrkFile(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 7:
				popts.setNumResFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

      		case '?':
				popts.usage();
				exit(-1);

      		default:
        		printf ("?? getopt returned character code %d ??\n", c);
				popts.usage();
				exit(-1);
    	}
	}

  	// make sure all required options have been specified
  	if (!(((spec.find("pds") != spec.end()) || (spec.find("pdsList") != spec.end())) && 
		((spec.find("pdb") != spec.end()) || (spec.find("pdbList") != spec.end()) || 
		(spec.find("seq") != spec.end()) || (spec.find("brk") != spec.end()) || 
		(spec.find("nr") != spec.end()))))
  	{
		popts.usage();
	    error("not all required options specified");
	}

	if (optind < argc)
	{
		popts.usage();
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
    	{
			printf ("%s ", argv[optind++]);
    	}
		printf ("\n");
		exit(-1);
	}

	if (!popts.getPdbFiles().empty())
	{
		ASSERT(popts.getPdbFiles().size() == popts.getPdsFiles().size(), "numbers of input PDS files and output PDB files not consistent");
	}
}
