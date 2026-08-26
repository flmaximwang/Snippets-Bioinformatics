/* example_search.cpp */

#include "master.h"

int main(int argc, char* argv[])
{
	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2)
	{
		cout << "USAGE:\nexample_search <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	bool bbrmsd;
	unsigned int i;
	SearchResults mlist;
	string qpdbfname = (string(argv[1]) + "/" + "frag1WP700173.pdb");
	QueryStruct qs;
	double rthresh;
	Structure Sq, *St = new Structure[3];
	string *tpdbfname = new string[3];
	TargetStruct *ts = new TargetStruct[3];

	tpdbfname[0] = (string(argv[1]) + "/" + "1g2c.pdb");
	tpdbfname[1] = (string(argv[1]) + "/" + "1wp7.pdb");
	tpdbfname[2] = (string(argv[1]) + "/" + "2d4u.pdb");

	// parse PDB files using MSL System and fill query and target structures
	Sq.readPDB(qpdbfname);
	// set query structure
	// where the 2nd argument is an identifier for query structure, and users can set whatever they want
	// the 3rd argument specifies whether or not to clean the PDB (true: Yes, false: No)
	qs.fillQueryStruct(Sq, qpdbfname, true);
	for (i = 0; i < 3; i++)
	{
		St[i].readPDB(tpdbfname[i]);
		// set target structure
		// where the 2nd argument is an identifier for target structure, and users can set whatever they want
		// the 3rd argument specifies whether or not to clean the PDB (true: Yes, false: No)
		// the 4th argument specifies whether or not to remove redundant chains (true: Yes, false: No)
		ts[i].fillTargetStruct(St[i], tpdbfname[i], true, true);
	}

	// search for all the matches within 1.0A full-BB RMSD cutoff
	bbrmsd = true; // choose full-BB RMSD based search
	rthresh = 1.0; // use 1.0A as RMSD cutoff
	qs.setSearchConstraints(rthresh); // set RMSD cutoff
	mlist.setQs(qs); // initialize match list

	for (i = 0; i < 3; i++) // search procedure
	{
		masterSearch(qs, ts[i], mlist, bbrmsd); // search
	}

	cout << mlist; // print out match list

	mlist.clear(); // clear match list for next round search

	// search for only top 100 matches within 2.0A full-BB RMSD cutoff
	rthresh = 2.0; // use 2.0A as RMSD cutoff
	qs.setSearchConstraints(rthresh); // set RMSD cutoff
	mlist.setQs(qs); // initialize match list
	mlist.setMatchMost(100); // keep at most 100 matches with the lowest RMSDs

	for (i = 0; i < 3; i++) // search procedure
	{
		masterSearch(qs, ts[i], mlist, bbrmsd); // search
	}

	cout << mlist; // print out match list

	delete [] St;
	delete [] tpdbfname;
	delete [] ts;

	return 0;
}

