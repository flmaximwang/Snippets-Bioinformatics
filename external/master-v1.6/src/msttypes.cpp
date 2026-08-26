#include "msttypes.h"

using namespace MST;

/* --------- Structure --------- */
Structure::Structure() {
  numResidues = numAtoms = 0;
}

Structure::Structure(string pdbFile, string options) {
  name = pdbFile;
  numResidues = numAtoms = 0;
  readPDB(pdbFile, options);
}

Structure::Structure(istream& is, string options) {
  name = "";
  numResidues = numAtoms = 0;
  readPDB(is, options);
}

Structure::Structure(const Structure& S) {
  copy(S);
}

void Structure::copy(const Structure& S) {
  name = S.name;
  numResidues = S.numResidues;
  numAtoms = S.numAtoms;
  for (int i = 0; i < S.chainSize(); i++) {
    chains.push_back(new Chain(S[i]));
    chains.back()->setParent(this);
    chainsByID[S[i].getID()] = chains.back();
    chainsBySegID[S[i].getSegID()] = chains.back();
  }
}

Structure::Structure(Chain& C) {
  numResidues = numAtoms = 0;
  appendChain(new Chain(C));
}

Structure::Structure(Residue& R) {
  numResidues = numAtoms = 0;
  Chain* newChain = appendChain("A", true);
  newChain->appendResidue(new Residue(R));
}

Structure::Structure(const vector<Atom*>& atoms) {
  numResidues = numAtoms = 0;
  addAtoms(atoms);
}

Structure::Structure(const vector<Residue*>& residues) {
  numResidues = numAtoms = 0;
  for (int i = 0; i < residues.size(); i++) addResidue(residues[i]);
}

/* The assumption is that if a Structure is deleted, all
 * of its children objects are no longer needed and should
 * go away. If a user needs to hold on to these, they
 * should generate copies as needed via copy constructors. */
Structure::~Structure() {
  deletePointers();
}

void Structure::deletePointers() {
  for (int i = 0; i < chains.size(); i++) delete(chains[i]);
}

void Structure::reset() {
  deletePointers();
  chains.resize(0);
  chainsByID.clear();
  chainsBySegID.clear();
  name = "";
  numResidues = numAtoms = 0;
}

Structure& Structure::operator=(const Structure& A) {
  reset();
  copy(A);
  return *this;
}

void Structure::readPDB(const string& pdbFile, string options) {
  name = pdbFile;
  fstream ifh; MstUtils::openFile(ifh, pdbFile, fstream::in, "Structure::readPDB");
  readPDB(ifh, options);
  ifh.close();
}

void Structure::readPDB(istream& is, string options) {
  int lastresnum = -999999;
  string lastresname = "XXXXXX";
  string lasticode = "";
  string lastchainID = "";
  string lastalt = " ";
  Chain* chain = NULL;
  Residue* residue = NULL;

  // various parsing options (the wonders of dealing with the good-old PDB format)
  bool ter = true;                   // flag to indicate that chain terminus was reached. Initialize to true so as to create a new chain upon reading the first atom.
  bool usesegid = false;             // use segment IDs to name chains instead of chain IDs? (useful when the latter are absent OR when too many chains, so need multi-letter names)
  bool skipHetero = false;           // skip hetero-atoms?
  bool charmmFormat = false;         // the PDB file was written by CHARMM (slightly different format)
  bool charmm19Format = false;       // upon reading, convert from all-hydrogen topology (param22 and higher) to the CHARMM19 united-atom topology (matters for HIS protonation states)
  bool fixIleCD1 = true;             // rename CD1 in ILE to CD (as is standard in MM packages)
  bool iCodesAsSepResidues = true;   // consequtive residues that differ only in their insertion code will be treated as separate residues
  bool uniqChainIDs = true;          // make sure chain IDs are unique, even if they are not unique in the read file
  bool ignoreTER = false;            // if true, will not pay attention to TER lines in deciding when chains end/begin
  bool verbose = true;               // report various warnings when weird things are found and fixed?

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("USESEGID") != string::npos) usesegid = true;
  if (options.find("SKIPHETERO") != string::npos) skipHetero = true;
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("ALLOW DUPLICATE CIDS") != string::npos) uniqChainIDs = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("IGNORE-ICODES") != string::npos) iCodesAsSepResidues = true;
  if (options.find("IGNORE-TER") != string::npos) ignoreTER = true;
  if (options.find("QUIET") != string::npos) verbose = false;

  // read line by line
  string line;
  while (getline(is, line)) {
    if (line.find("END") == 0) break;
    if ((line.find("TER") == 0) && !ignoreTER) { ter = true; continue; }
    if ((skipHetero && (line.find("ATOM") != 0)) || (!skipHetero && (line.find("ATOM") != 0) && (line.find("HETATM") != 0))) continue;

    /* Now read atom record. Sometimes PDB lines are too short (if they do not contain some
     * of the last optional columns). We don't want to read past the end of the string! */
    line += string(100, ' ');
    int atominx = MstUtils::toInt(line.substr(6, 5));
    string atomname = MstUtils::trim(line.substr(12, 4));
    string alt = line.substr(16, 1);
    string resname = MstUtils::trim(line.substr(17, 4));
    string chainID = MstUtils::trim(line.substr(21, 1));
    int resnum = charmmFormat ? MstUtils::toInt(line.substr(23, 4)) : MstUtils::toInt(line.substr(22, 4));
    string icode = charmmFormat ? " " : line.substr(26, 1);
    mstreal x = MstUtils::toReal(line.substr(30, 8));
    mstreal y = MstUtils::toReal(line.substr(38, 8));
    mstreal z = MstUtils::toReal(line.substr(46, 8));
    string segID = MstUtils::trim(line.substr(72, 4));
    mstreal B = MstUtils::toReal(line.substr(60, 6), false);
    mstreal occ = MstUtils::toReal(line.substr(54, 6), false);
    bool het = (line.find("HETATM") == 0);

    // use segment ID's instead of chain ID's?
    if (usesegid) {
      chainID = segID;
    } else if ((chainID.compare("") == 0) && (segID.size() > 0) && (isalnum(segID[0]))) {
      // use first character of segment name if no chain name is specified, a segment ID is specified, and the latter starts with an alphanumeric character
      chainID = segID.substr(0, 1);
    }

    // create a new chain object, if necessary
    if ((chainID.compare(lastchainID) != 0) || ter) {
      chain = new Chain(chainID, segID);
      this->appendChain(chain, uniqChainIDs);
      // non-unique chains will be automatically renamed (unless the user specified not to rename chains), BUT we need to
      // remember the name that was actually read, since this name is what will be used to determine when the next chain comes
      if (verbose && chainID.compare(chain->getID())) {
        MstUtils::warn("chain name '" + chainID + "' was repeated in '" + name + "', renaming the chain to '" + chain->getID() + "'", "Structure::readPDB");
      }

      // start to count residue numbers in this chain
      lastresnum = -999999;
      lastresname = "";
      ter = false;
    }

    if (charmm19Format) {
      if (resname.compare("HSE") == 0) resname = "HSD";   // neutral HIS, proton on ND1
      if (resname.compare("HSD") == 0) resname = "HIS";   // neutral HIS, proton on NE2
      if (resname.compare("HSC") == 0) resname = "HSP";   // doubley-protonated +1 HIS
    }
    // many PDB files in the Protein Data Bank call the delta carbon of isoleucine CD1, but
    // the convention in basically all MM packages is to call it CD, since there is only one
    if (fixIleCD1 && atomname.compare("CD1") == 0) atomname = "CD";

    // if necessary, make a new residue
    bool reallyNewAtom = true; // is this a truely new atom, as opposed to an alternative position?
    if ((resnum != lastresnum) || resname.compare(lastresname) || (iCodesAsSepResidues && (icode.compare(lasticode))))  {
      // this corresponds to a case, where the alternative location flag is being used to
      // designate two (or more) different possible amino acids at a particular position
      // (e.g., where the density is not clear to assign one). In this case, we shall keep
      // only the first option, because we don't know any better.
      if ((resnum == lastresnum) && resname.compare(lastresname) && (alt != lastalt)) {
        continue;
      }
      residue = new Residue(resname, resnum, icode[0]);
      chain->appendResidue(residue);
    } else if (alt.compare(" ")) {
      // if this is not a new residue AND the alternative location flag is specified,
      // figure out if another location for this atom has already been given. If not,
      // then treat this as the "primary" location, and whatever other locations
      // are specified will be treated as alternatives.
      Atom* a = residue->findAtom(atomname, false);
      if (a) {
        reallyNewAtom = false;
        a->addAlternative(x, y, z, B, occ, alt[0]);
      }
    }
    // if necessary, make a new atom
    if (reallyNewAtom) {
      residue->appendAtom(new Atom(atominx, atomname, x, y, z, B, occ, het, alt[0]));
    }

    // remember previous values for determining whether something interesting happens next
    lastresnum = resnum;
    lasticode = icode;
    lastresname = resname;
    lastchainID = chainID;
    lastalt = alt;
  }
}

void Structure::writePDB(const string& pdbFile, string options) const {
  fstream ofs; MstUtils::openFile(ofs, pdbFile, fstream::out, "Structure::writePDB(string, string)");
  writePDB(ofs, options);
  ofs.close();
}

void Structure::writePDB(ostream& ofs, string options) const {
  options = MstUtils::uc(options);

///  my $chainstr = shift; // probably want to implement this eventually. Or maybe some more generic selection mechanism based on regular expressions applied onto full atom strings.

  // various formating options (the wonders of dealing with the good-old PDB format)
  bool charmmFormat = false;         // the PDB file is intended for use in CHARMM or some other MM package
  bool charmm19Format = false;       // upon writing, convert from all-hydrogen topology (param 22 and higher) to CHARMM19 united-atom topology (matters for HIS protonation states)
  bool charmm22Format = false;       // upon writing, convert from CHARMM19 united-atom topology to all-hydrogen param 22 topology (matters for HIS protonation states). Also works for converting generic PDB files downloaded from the PDB.
  bool genericFormat = false;        // upon writing, convert to a generic PDB naming convention (no protonation state specified for HIS)
  bool fixIleCD1 = true;             // rename CD1 in ILE to CD in the output file (as is standard in MM packages)
  bool renumber = false;             // upon writing, renumber residue and atom names to start from 1 and go in order
  bool noend = false;                // do not write END at the end of the PDB file (e.g., useful for concatenating chains from several structures)
  bool noter = false;                // do not demark the end of each chain with TER (this is not _really_ necessary, assuming chain names are unique, and it is sometimes nice not to have extra lines other than atoms)

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("CHARMM22") != string::npos) charmm22Format = true;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("RENUMBER") != string::npos) renumber = true;
  if (options.find("NOEND") != string::npos) noend = true;
  if (options.find("NOTER") != string::npos) noter = true;
  if (charmm19Format && charmm22Format) MstUtils::error("CHARMM 19 and 22 formatting options cannot be specified together", "Structure::writePDB");

  int atomIndex = 0;
  for (int ci = 0; ci < this->chainSize(); ci++) {
    Chain& chain = (*this)[ci];
    for (int ri = 0; ri < chain.residueSize(); ri++) {
      Residue residue = chain[ri]; // NOTE: using a copy constructor here, in case residue details will be changed for formatting reasons upon writing
      residue.setParent(&chain);   // by default, objects are copied as being disembodied (so as not to create inconsistent states)
      for (int ai = 0; ai < residue.atomSize(); ai++) {
        Atom& atom = residue[ai]; // NOTE: no need to copy atom, since residue copying does a deep copy
        atomIndex++;
        // dirty details of formating for MM purposes converting
        if (charmmFormat) {
          if (residue.isNamed("ILE") && atom.isNamed("CD1")) atom.setName("CD");
          if (atom.isNamed("O") && (ri == chain.residueSize() - 1)) atom.setName("OT1");
          if (atom.isNamed("OXT") && (ri == chain.residueSize() - 1)) atom.setName("OT2");
          if (residue.isNamed("HOH")) residue.setName("TIP3");
        }
        if (charmm19Format) {
          if (residue.isNamed("HSD")) residue.setName("HIS"); // neutral HIS, proton on ND1
          if (residue.isNamed("HSE")) residue.setName("HSD"); // neutral HIS, proton on NE2
          if (residue.isNamed("HSC")) residue.setName("HSP"); // doubley-protonated +1 HIS
        } else if (charmm22Format) {
          /* this will convert from CHARMM19 to CHARMM22 as well as from a generic downlodaded
           * PDB file to one ready for use in CHARMM22. The latter is because in the all-hydrogen
           * topology, HIS protonation state must be explicitly specified, so there is no HIS per se.
           * Whereas in typical downloaded PDB files HIS is used for all histidines (usually, one
           * does not even really know the protonation state). Whether sometimes people do specify it
           * nevertheless, and what naming format they use to do so, I am not sure (welcome to the
           * PDB file format). But certainly almost always it is just HIS. Below HIS is renamed to
           * HSD, the neutral form with proton on ND1. This is an assumption; not a perfect one, but
           * something needs to be assumed. Doing this renaming will make the PDB file work in MM
           * packages with the all-hydrogen model. */
          if (residue.isNamed("HSD")) residue.setName("HSE"); // neutral HIS, proton on NE2
          if (residue.isNamed("HIS")) residue.setName("HSD"); // neutral HIS, proton on ND1
          if (residue.isNamed("HSP")) residue.setName("HSC"); // doubley-protonated +1 HIS
        } else if (genericFormat) {
          if (residue.isNamed("HSD")) residue.setName("HIS");
          if (residue.isNamed("HSP")) residue.setName("HIS");
          if (residue.isNamed("HSE")) residue.setName("HIS");
          if (residue.isNamed("HSC")) residue.setName("HIS");
          if (residue.isNamed("ILE") && atom.isNamed("CD")) atom.setName("CD1");
        }

        // write the atom line
        if (renumber) {
          ofs << atom.pdbLine(ri+1, atomIndex) << endl;
        } else {
          ofs << atom.pdbLine() << endl;
        }

      }
      if (!noter && (ri == chain.residueSize() - 1)) {
        ofs << "TER" << endl;
      }
    }
    if (!noend && (ci == this->chainSize() - 1)) {
      ofs << "END" << endl;
    }
  }
}

bool Structure::appendChain(Chain* C, bool allowRename) {
  chains.push_back(C);
  bool cidUnique = (chainsByID.find(C->getID()) == chainsByID.end());

  // if allowed to rename and there is a name clash, try to pick a unique chain name
  if (allowRename && !cidUnique) {
    string goodNames = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    int f = -1;
    for (int i = 0; i < goodNames.size(); i++) {
      if (chainsByID.find(goodNames.substr(i, 1)) == chainsByID.end()) { f = i; break; }
    }
    if (f < 0) {
      MstUtils::warn("ran out of reasonable single-letter chain names, will use more than one character (PDB sctructure may be repeating chain IDs upon writing, but should still have unique segment IDs)!", "Structure::appendChain");
      string longName; bool found = false;
      for (int i = 0; i < goodNames.size(); i++) {
        for (int k = 0; k < 999; k++) {
          longName = goodNames.substr(i, 1) + MstUtils::toString(k);
          if (chainsByID.find(longName) == chainsByID.end()) { found = true; break; }
        }
        if (found) break;
      }
      MstUtils::assert(found, "ran out of even multi-character chain names -- your PDB structure really has more than 36,000 chains???", "Structure::appendChain");
      C->setID((string) longName);
      C->setSegID(C->getID());
    } else {
      C->setID(goodNames.substr(f, 1));
      cidUnique = true;
    }
  }

  chainsByID[C->getID()] = C;
  chainsBySegID[C->getSegID()] = C;
  C->setParent(this);
  numAtoms += C->atomSize();
  numResidues += C->residueSize();
  return cidUnique;
}

Chain* Structure::appendChain(string cid, bool allowRename) {
  Chain* newChain = new Chain(cid, cid);
  this->appendChain(newChain, allowRename);
  return newChain;
}

void Structure::deleteChain(Chain* chain) {
  // make sure the chain is from this Structure
  int i = 0;
  for (i = 0; i < chains.size(); i++) {
    if (chains[i] == chain) break;
  }
  if (i == chains.size()) MstUtils::error("chain not from this structure", "Structure::deleteChain");

  chains.erase(chains.begin() + i);
  numResidues -= chain->residueSize();
  numAtoms -= chain->atomSize();
  if (chainsByID.find(chain->getID()) != chainsByID.end()) chainsByID.erase(chain->getID());
  if (chainsBySegID.find(chain->getSegID()) != chainsBySegID.end()) chainsBySegID.erase(chain->getSegID());
  delete chain;
}

Residue& Structure::getResidue(int i) const {
  if ((i < 0) || (i >= residueSize()))
    MstUtils::error("residue index " + MstUtils::toString(i) + " out of range for Structure", "Structure::getResidue(int)");
  int io = i;
  for (int ci = 0; ci < chainSize(); ci++) {
    Chain& chain = getChain(ci);
    if (i >= chain.residueSize()) {
      i -= chain.residueSize();
    } else {
      return chain[i];
    }
  }
  MstUtils::error("something strange happened when fetching residue " + MstUtils::toString(io) + " from Structure object that reports " + MstUtils::toString(this->residueSize()) + " residues; most likely, various counters are inconsistent in Structure object", "Structure::getResidue(int)");
  return *(new Residue()); // just to make the compiler happy and not throw a warning; this is never reached
}

vector<Atom*> Structure::getAtoms() const {
  vector<Atom*> vec(this->atomSize());
  int ii = 0;
  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    for (int j = 0; j < c.residueSize(); j++) {
      Residue& r = c[j];
      for (int k = 0; k < r.atomSize(); k++) {
        vec[ii] = &(r[k]);
        ii++;
      }
    }
  }

  return vec;
}

vector<Residue*> Structure::getResidues() const {
  vector<Residue*> vec;

  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    vector<Residue*> chainResidues = c.getResidues();
    vec.insert(vec.end(), chainResidues.begin(), chainResidues.end());
  }

  return vec;
}

void Structure::renumber(int startResNum, int startAtomIndex) {
  int index = startAtomIndex;
  for (int i = 0; i < chains.size(); i++) {
    Chain& chain = (*this)[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      res.setNum(startResNum+j);
      for (int k = 0; k < res.atomSize(); k++) {
        Atom& a = res[k];
        a.setIndex(index);
        index++;
      }
    }
  }
}

void Structure::addAtom(Atom* A) {
  if ((A->getParent() == NULL) || (A->getParent()->getParent() == NULL)) MstUtils::error("cannot add a disembodied Atom", "Structure::addAtom");
  Residue* oldResidue = A->getParent();
  Chain* oldChain = oldResidue->getParent();
  Chain* newChain; Residue* newResidue; Atom* newAtom;

  // is there a chain matching the Atom's parent chain? If not, create one.
  newChain = getChainByID(oldChain->getID());
  if (newChain == NULL) {
    newChain = new Chain(oldChain->getID(), oldChain->getSegID());
    this->appendChain(newChain);
  }

  // is there a residue matching the Atoms parent residue? If not, create one.
  newResidue = newChain->findResidue(oldResidue->getName(), oldResidue->getNum(), oldResidue->getIcode());
  if (newResidue == NULL) {
    newResidue = new Residue(oldResidue->getName(), oldResidue->getNum(), oldResidue->getIcode());
    newChain->appendResidue(newResidue);
  }

  // finally, insert the atom into the correct residue
  newAtom = new Atom(*A);
  newResidue->appendAtom(newAtom);
}

void Structure::addAtoms(vector<Atom*>* atoms) {
  for (int i = 0; i < atoms->size(); i++) addAtom((*atoms)[i]);
}

Residue* Structure::addResidue(Residue* res) {
  if (res->getParent() == NULL) MstUtils::error("cannot add a disembodied Residue", "Structure::addResidue");
  Chain* oldChain = res->getParent();

  // is there a chain matching the Residue's parent chain? If not, create one.
  Chain* newChain = getChainByID(oldChain->getID());
  if (newChain == NULL) {
    newChain = new Chain(oldChain->getID(), oldChain->getSegID());
    this->appendChain(newChain);
  }

  // append a copy of the given residue into the correct chain
  Residue* newResidue = new Residue(*res);
  newChain->appendResidue(newResidue);
  return newResidue;
}

int Structure::getResidueIndex(Residue* res) {
  Chain* parentChain = res->getParent();
  if (parentChain == NULL)
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(res) + "'", "Structure::getResidueIndex(Residue*)");

  int n = 0;
  bool found = false;
  for (int i = 0; i < chainSize(); i++) {
    Chain& chain = (*this)[i];
    if (&chain  == parentChain) {
      n += chain.getResidueIndex(res);
      found = true;
      break;
    }
    n += chain.residueSize();
  }
  if (!found) MstUtils::error("residue not from Structure '" + MstUtils::toString(res) + "'", "Structure::getResidueIndex(Residue*)");

  return n;
}

/* --------- Chain --------- */
Chain::Chain() {
  numAtoms = 0;
  parent = NULL;
}

Chain::Chain(const Chain& C) {
  numAtoms = C.numAtoms;
  parent = NULL;
  for (int i = 0; i < C.residueSize(); i++) {
    residues.push_back(new Residue(C[i]));
    residues.back()->setParent(this);
    residueIndexInChain[residues.back()] = i;
  }
  cid = C.cid;
  sid = C.sid;
}

Chain::Chain(const string& chainID, const string& segID) {
  numAtoms = 0;
  parent = NULL;
  cid = chainID;
  sid = segID;
}

Chain::Chain(const string& chainID, const string& segID, const vector<Residue*>& residues) {
	numAtoms = 0;
	parent = NULL;
	cid = chainID;
	sid = segID;
	appendResidueCopies(residues);
}

Chain::~Chain() {
  for (int i = 0; i < residues.size(); i++) delete(residues[i]);
}

vector<Atom*> Chain::getAtoms() {
  vector<Atom*> vec;

  for (int i = 0; i < residues.size(); i++) {
    Residue& r = getResidue(i);
    vector<Atom*> resAtoms = r.getAtoms();
    vec.insert(vec.end(), resAtoms.begin(), resAtoms.end());
  }

  return vec;
}

int Chain::getResidueIndex(const Residue* res) {
  if (residueIndexInChain.find((Residue*) res) == residueIndexInChain.end())
    MstUtils::error("passed residue does not appear in chain's index map", "Chain::residueIndexInChain");
  return residueIndexInChain[(Residue*) res];
}

int Chain::getIndex() const {
  const Structure* P = getParent();
  if (P == NULL) MstUtils::error("cannot get index of disembodied chain", "Chain::getIndex()");
  int idx = -1;
  for (int i = 0; i < P->chainSize(); i++) {
    if (&(P->getChain(i)) == this) return i;
  }
  MstUtils::error("strange error: Chain does not appear to belong to its parent Structure", "Chain::getIndex()");
  return -1;
}

void Chain::incrementNumAtoms(int delta) {
  numAtoms += delta;
  if (parent != NULL) {
    parent->incrementNumAtoms(delta);
  }
}

void Chain::appendResidue(Residue* R) {
  incrementNumAtoms(R->atomSize());
  if (parent != NULL) {
    parent->incrementNumResidues();
  }
  residues.push_back(R);
  R->setParent(this);
  residueIndexInChain[R] = residues.size() - 1;
}

void Chain::appendResidueCopies(const vector<Residue*>& residues) {
	for (int i = 0; i < residues.size(); i++) {
		insertResidueCopy(residues[i], -1);
	}
}

void Chain::insertResidue(Residue* R, int index) {
  if ((index < 0) || (index > residues.size()))
    MstUtils::error("residue index '" + MstUtils::toString(index) + "' out of range", "Chain::insertResidue(Residue*, int)");
  incrementNumAtoms(R->atomSize());
  if (parent != NULL) {
    parent->incrementNumResidues();
  }
  residues.insert(residues.begin() + index, R);
  R->setParent(this);
  residueIndexInChain[R] = index;
  for (int i = index + 1; i < residues.size(); i++) {
    residueIndexInChain[residues[i]]++;
  }
}
Residue* Chain::insertResidueCopy(Residue* R, int index) {
  Residue* newRes = new Residue(*R);
  if (index == -1) this->appendResidue(newRes);
  else this->insertResidue(newRes, index);
  return newRes;
}

Residue* Chain::insertResidueCopy(Residue& R, int index) {
  Residue* newRes = new Residue(R);
  if (index == -1) this->appendResidue(newRes);
  else this->insertResidue(newRes, index);
  return newRes;
}

Residue* Chain::findResidue(string resname, int resnum) {
  for (int i = 0; i < residueSize(); i++) {
    Residue& res = (*this)[i];
    if ((res.getNum() == resnum) && (res.isNamed(resname))) return &res;
  }
  return NULL;
}

Residue* Chain::findResidue(string resname, int resnum, char icode) {
  for (int i = 0; i < residueSize(); i++) {
    Residue& res = (*this)[i];
    if ((res.getNum() == resnum) && (res.isNamed(resname)) && (res.getIcode() == icode)) return &res;
  }
  return NULL;
}

/* --------- Residue --------- */
Residue::Residue() {
  resname = "UNK";
  resnum = 1;
  parent = NULL;
  icode = ' ';
}

Residue::Residue(const Residue& R, bool copyAlt) {
  parent = NULL;
  for (int i = 0; i < R.atomSize(); i++) {
    atoms.push_back(new Atom(R[i], copyAlt));
    atoms.back()->setParent(this);
  }
  resnum = R.resnum;
  resname = R.resname;
  icode = R.icode;
}

Residue::Residue(string _resname, int _resnum, char _icode) {
  resname = _resname;
  resnum = _resnum;
  parent = NULL;
  icode = _icode;
}

Residue::~Residue() {
  for (int i = 0; i < atoms.size(); i++) delete(atoms[i]);
}

void Residue::appendAtom(Atom* A) {
  atoms.push_back(A);
  A->setParent(this);
  if (parent != NULL) {
    parent->incrementNumAtoms();
  }
}

void Residue::appendAtoms(const vector<Atom*>& A) {
  for (int i = 0; i < A.size(); i++) appendAtom(A[i]);
}

void Residue::deleteAtom(int i) {
  if ((i < 0) || (i >= atoms.size())) {
    MstUtils::error("index out of range of atom vector in residue", "Residue::deleteAtom");
  }
  delete atoms[i];
  atoms.erase(atoms.begin() + i);
  if (parent != NULL) {
    parent->incrementNumAtoms(-1);
  }
}

void Residue::copyAtoms(Residue& R, bool copyAlt) {
  copyAtoms(R.getAtoms(), copyAlt);
}

void Residue::copyAtoms(const vector<Atom*>& _atoms, bool copyAlt) {
  deleteAtoms();
  for (int i = 0; i < _atoms.size(); i++) appendAtom(new Atom(_atoms[i], copyAlt));
}

void Residue::makeAlternativeMain(int altInd) {
  for (int i = 0; i < atomSize(); i++) {
    (*this)[i].makeAlternativeMain(altInd);
  }
}

void Residue::deleteAtoms() {
  if (parent != NULL) {
    parent->incrementNumAtoms(-atoms.size());
  }
  for (int i = 0; i < atoms.size(); i++) {
    delete atoms[i];
  }
  atoms.resize(0);
}

void Residue::replaceAtoms(const vector<Atom*>& newAtoms, vector<int>* toRemove) {
  bool delAll = (toRemove == NULL);
  int N = delAll ? atoms.size() : toRemove->size();
  if (parent != NULL) {
    parent->incrementNumAtoms(newAtoms.size() - N);
  }
  Structure* S = getStructure();

  // delete those atoms needing deletion
  for (int i = 0; i < N; i++) {
    int ai = delAll ? i : (*toRemove)[i];
    if ((ai < 0) || (ai >= atoms.size())) {
      MstUtils::error("index out of range of atom vector in residue", "Residue::replaceAtoms");
    }
    delete atoms[ai];
    atoms[ai] = NULL;
  }

  // create a new atom vector without them
  vector<Atom*> oldAtoms = atoms;
  atoms.resize(atoms.size() - N + newAtoms.size());
  int k = 0;
  for (int i = 0; i < oldAtoms.size(); i++) {
    if (oldAtoms[i] == NULL) continue;
    atoms[k] = oldAtoms[i];
    k++;
  }

  // and now append the new atoms
  for (int i = 0; i < newAtoms.size(); i++) {
    atoms[k] = newAtoms[i];
    atoms[k]->setParent(this);
    k++;
  }
}

void Residue::replaceAtoms(const vector<Atom*>& newAtoms, const vector<Atom*>& oldAtoms) {
  vector<int> oldAtomIndices;
  for (int i = 0; i < oldAtoms.size(); i++) {
    for (int k = 0; k < atoms.size(); k++) {
      if (atoms[k] == oldAtoms[i]) {
        oldAtomIndices.push_back(k); break;
      }
    }
  }
  replaceAtoms(newAtoms, &oldAtomIndices);
}

Atom* Residue::findAtom(string _name, bool strict) const {
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i]->isNamed(_name)) return atoms[i];
  }
  if (strict) MstUtils::error("could not find atom named '" + _name + "' in residue " + MstUtils::toString(*this), "Residue::findAtom");
  return NULL;
}

string Residue::getChainID(bool strict) {
  if (parent == NULL) {
    if (strict) MstUtils::error("residue has no parent", "Residue::getChainID");
    return "";
  }
  return parent->getID();
}

int Residue::getResidueIndex() const {
  Chain* parentChain = NULL; Structure* parentStructure = NULL;
  parentChain = getParent();
  if (parentChain != NULL) parentStructure = parentChain->getParent();
  if ((parentChain == NULL) || (parentStructure == NULL))
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(this) + "'", "Residue::getResidueIndex()");

  int n = 0;
  bool found = false;
  for (int i = 0; i < parentStructure->chainSize(); i++) {
    Chain& chain = (*parentStructure)[i];
    if (&chain == parentChain) {
      n += chain.getResidueIndex(this);
      found = true;
      break;
    }
    n += chain.residueSize();
  }
  if (!found) MstUtils::error("residue not in its parent Structure '" + MstUtils::toString(this) + "'", "Residue::getResidueIndex()");

  return n;
}

int Residue::getResidueIndexInChain() const {
  if (getParent() == NULL)
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(this) + "'", "Residue::getResidueIndexInChain()");
  return getParent()->getResidueIndex(this);
}

/* --------- Atom --------- */
Atom::atomInfo::atomInfo() {
  parent = NULL;
  het = false;
  name = MstUtils::copyStringC("UNK");
  alternatives = NULL;
  alt = ' ';
  index = 0;
  occ = B = 0;
}

Atom::atomInfo::atomInfo(const atomInfo& other, bool copyAlt) {
  index = other.index;
  name = NULL;
  setName(other.name);
  B = other.B;
  occ = other.occ;
  het = other.het;
  alt = other.alt;
  parent = other.parent;
  if (copyAlt && (other.alternatives != NULL)) {
    alternatives = new vector<altInfo>(other.alternatives->size());
    for (int i = 0; i < other.alternatives->size(); i++) {
      (*alternatives)[i] = (*(other.alternatives))[i];
    }
  } else {
    alternatives = NULL;
  }
}

Atom::atomInfo::atomInfo(int _index, const string& _name, mstreal _B, mstreal _occ, bool _het, char _alt, Residue* _parent) {
  index = _index;
  name = NULL;
  setName(_name);
  B = _B;
  occ = _occ;
  het = _het;
  alt = _alt;
  parent = _parent;
  alternatives = NULL;
}

Atom::atomInfo::~atomInfo() {
  if (name != NULL) delete[] name;
  if (alternatives != NULL) delete alternatives;
}

Atom::Atom() {
  x = y = z = 0;
  info = new atomInfo();
}

Atom::Atom(const Atom& A, bool copyAlt) {
  x = A.x;
  y = A.y;
  z = A.z;
  info = NULL;
  if (A.hasInfo()) info = new atomInfo(*(A.info));
}

Atom::Atom(int _index, const string& _name, mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, bool _het, char _alt, Residue* _parent) {
  x = _x;
  y = _y;
  z = _z;
  info = new atomInfo(_index, _name, _B, _occ, _het, _alt, _parent);
}

Atom::~Atom() {
  if (info != NULL) delete info;
}

mstreal& Atom::operator[](int i) {
  switch(i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      MstUtils::error("invalid coordinate index " + MstUtils::toString(i), "Atom::operator[](int)");
  }
  return x; // just to silence the warning from some compilres; in reality, this is never reached
}

const mstreal& Atom::operator[](int i) const {
  switch(i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      MstUtils::error("invalid coordinate index " + MstUtils::toString(i), "Atom::operator[](int) const");
  }
  return x; // just to silence the warning from some compilres; in reality, this is never reached
}

CartesianPoint Atom::getCoor() const {
  CartesianPoint coor(x, y, z); return coor;
}

void Atom::atomInfo::setName(const char* _name) {
  if (name != NULL) delete[] name;
  name = new char[strlen(_name)+1];
  strcpy(name, _name);
}

CartesianPoint Atom::atomInfo::getAltCoor(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltCoor");
  altInfo& targ = (*alternatives)[altInd];
  CartesianPoint coor(targ.x, targ.y, targ.z);
  return coor;
}
CartesianPoint Atom::getAltCoor(int altInd) const { return info->getAltCoor(altInd); }

mstreal Atom::atomInfo::getAltB(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltB");
  return (*alternatives)[altInd].B;
}

mstreal Atom::atomInfo::getAltOcc(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltOcc");
  return (*alternatives)[altInd].occ;
}

char Atom::atomInfo::getAltLocID(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltLocID");
  return (*alternatives)[altInd].alt;
}

void Atom::setCoor(const CartesianPoint& xyz) {
  x = xyz[0]; y = xyz[1]; z = xyz[2];
}

void Atom::setCoor(const Atom& a) {
  x = a.getX(); y = a.getY(); z = a.getZ();
}

void Atom::atomInfo::setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z) {
  altInfo& A = (*alternatives)[ai];
  A.x = _x; A.y = _y; A.z = _z;
}

void Atom::swapWithAlternative(int altInd) {
  if (!hasInfo() || (altInd >= numAlternatives())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(numAlternatives()) + " alternatives available)", "Atom::swapWithAlternative");
  atomInfo::altInfo& targ = (*(info->alternatives))[altInd];
  atomInfo::altInfo temp = targ;
  targ.x = x; targ.y = y; targ.z = z; targ.occ = getOcc(); targ.B = getB(); targ.alt = getAlt();
  x = temp.x; y = temp.y; z = temp.z; setOcc(temp.occ); setB(temp.B); setAlt(temp.alt);
}

void Atom::makeAlternativeMain(int altInd) {
  if (!hasInfo() || (altInd >= numAlternatives())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(numAlternatives()) + " alternatives available)", "Atom::makeAlternativeMain");
  atomInfo::altInfo& targ = (*(info->alternatives))[altInd];
  x = targ.x; y = targ.y; z = targ.z; setOcc(targ.occ); setB(targ.B); setAlt(targ.alt);
}

void Atom::atomInfo::addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt) {
  if (alternatives == NULL) {
    alternatives = new vector<altInfo>(0);
  }
  alternatives->push_back(altInfo(_x, _y, _z, _B, _occ, _alt));
}

void Atom::atomInfo::removeLastAlternative() {
  alternatives->pop_back();
}

void Atom::atomInfo::removeAlternative(int i) {
  alternatives->erase(alternatives->begin() + i);
}

void Atom::atomInfo::clearAlternatives() {
  if (alternatives != NULL) { delete alternatives; alternatives = NULL; }
}

string Atom::pdbLine(int resIndex, int atomIndex) {
  char line[100]; // a PDB line is at most 80 characters, so this is plenty
  string resname = "UNK"; string chainID = "?"; string segID = "?";
  int resnum = 1; char icode = ' ';

  // chain and residue info
  Residue* parent = getParent();
  if (parent != NULL) {
    resname = parent->getName();
    if (resname.length() > 4) resname = resname.substr(0, 4);
    resnum = parent->getNum();
    icode = parent->getIcode();
    Chain* chain = parent->getParent();
    if (chain != NULL) {
      chainID = chain->getID();
      if (chainID.length() > 1) chainID = chainID.substr(0, 1);
      segID = chain->getSegID();
      if (segID.length() > 4) segID = segID.substr(0, 4);
    }
  }

  // atom name placement is different when it is 4 characters long
  char atomname[5];
  if (strlen(getNameC()) < 4) { sprintf(atomname, " %-.3s", getNameC()); }
  else { sprintf(atomname, "%.4s", getNameC()); }

  // moduli are used to make sure numbers do not go over prescribe field widths (this is not enforced by sprintf like with strings)
  sprintf(line, "%6s%5d %-4s%c%-4s%.1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %.4s",
          isHetero() ? "HETATM" : "ATOM  ", atomIndex % 100000, atomname, getAlt(), resname.c_str(), chainID.c_str(),
          resnum % 10000, icode, x, y, z, getOcc(), getB(), segID.c_str());

  return (string) line;
}

mstreal Atom::distance(const Atom& another) const {
  return sqrt((x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z));
}

mstreal Atom::distance2(const Atom& another) const {
  return (x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z);
}

ostream& MST::operator<<(ostream &_os, const Atom& _atom) {
  _os << _atom.getName() << _atom.getAlt() << " " << _atom.getIndex() << " " << (_atom.isHetero() ? "HETERO" : "");
  _os << _atom.getX() << " " << _atom.getY() << " " << _atom.getZ() << " : " << _atom.getOcc() << " " << _atom.getB();
  return _os;
}

void Atom::stripInfo() {
  if (info != NULL) {
    delete info;
    info = NULL;
  }
}

/* --------- AtomPointerVector --------- */

void AtomPointerVector::copyCoordinates(const AtomPointerVector& other) {
  if (size() != other.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::copyCoordinates");
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor(other[i]);
  }
}

void AtomPointerVector::push_back(const Residue& R) {
  int sz = this->size();
  this->resize(sz + R.atomSize());
  for (int i = 0; i < R.atomSize(); i++) {
    (*this)[i + sz] = &(R[i]);
  }
}

CartesianPoint AtomPointerVector::getGeometricCenter() {
  CartesianPoint C(3, 0);
  for (int i = 0; i < this->size(); i++) {
    C += CartesianPoint(*((*this)[i]));
  }
  C /= this->size();
  return C;
}

void AtomPointerVector::getGeometricCenter(mstreal& xc, mstreal& yc, mstreal& zc) {
  xc = 0; yc = 0; zc = 0;
  Atom* a;
  for (int i = 0; i < this->size(); i++) {
    a = (*this)[i];
    xc += (*a)[0]; yc += (*a)[1]; zc += (*a)[2];
  }
  xc /= size(); yc /= size(); zc /= size();
}

void AtomPointerVector::center() {
  CartesianPoint C = getGeometricCenter();
  for (int i = 0; i < this->size(); i++) {
    Atom& a = *((*this)[i]);
    for (int k = 0; k < 3; k++) a[k] -= C[k];
  }
}

mstreal AtomPointerVector::radiusOfGyration() {
  CartesianPoint center = getGeometricCenter();
  mstreal s = 0;
  for (int i = 0; i < size(); i++) {
    s += CartesianPoint((*this)[i]).distance2(center);
  }
  return sqrt(s / size());
}

mstreal AtomPointerVector::boundingSphereRadiusCent() {
  CartesianPoint center = getGeometricCenter();
  mstreal r = 0;
  for (int i = 0; i < size(); i++) {
    mstreal dist = CartesianPoint((*this)[i]).distance(center);
    if (dist > r) r = dist;
  }
  return r;
}

void AtomPointerVector::deletePointers() {
  for (int i = 0; i < size(); i++) delete (*this)[i];
  resize(0);
}

AtomPointerVector AtomPointerVector::clone() const {
  AtomPointerVector into;
  clone(into);
  return into;
}

AtomPointerVector AtomPointerVector::subvector(int beg, int end) {
  return AtomPointerVector(vector<Atom*>(begin() + beg, begin() + end));
}

void AtomPointerVector::clone(AtomPointerVector& into) const {
  int L = into.size();
  into.resize(L + size());
  for (int i = 0; i < size(); i++) {
    Atom* newAtom = new Atom(*((*this)[i]));
    newAtom->setParent(NULL);
    into[L + i] = newAtom;
  }
}

AtomPointerVector& AtomPointerVector::operator+=(const AtomPointerVector& rhs) {
  if (size() != rhs.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::operator+=");
  for (int i = 0; i < size(); i++) {
    Atom& A = *((*this)[i]);
    const Atom& B = *(rhs[i]);
    A.setCoor(A.getX() + B.getX(), A.getY() + B.getY(), A.getZ() + B.getZ());
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator-=(const AtomPointerVector& rhs) {
  if (size() != rhs.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::operator+=");
  for (int i = 0; i < size(); i++) {
    Atom& A = *((*this)[i]);
    Atom& B = *(rhs[i]);
    A.setCoor(A.getX() - B.getX(), A.getY() - B.getY(), A.getZ() - B.getZ());
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator/=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor((*this)[i]->getX()/s, (*this)[i]->getY()/s, (*this)[i]->getZ()/s);
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator*=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor((*this)[i]->getX()*s, (*this)[i]->getY()*s, (*this)[i]->getZ()*s);
  }
  return *this;
}


ostream& MST::operator<<(ostream &_os, const AtomPointerVector& _atoms) {
  for (int i = 0; i < _atoms.size(); i++) {
    _os << *(_atoms[i]) << endl;
  }
  return _os;
}

/* --------- CartesianPoint --------- */

CartesianPoint::CartesianPoint(const Atom& A) {
  this->resize(3, 0);
  (*this)[0] = A.getX();
  (*this)[1] = A.getY();
  (*this)[2] = A.getZ();
}

CartesianPoint& CartesianPoint::operator+=(const CartesianPoint &rhs) {
  if (size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator+=");
  for (int i = 0; i < size(); i++) {
    (*this)[i] += rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator-=(const CartesianPoint &rhs) {
  if (size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator-=");
  for (int i = 0; i < size(); i++) {
    (*this)[i] -= rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator*=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i] *= s;
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator/=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i] /= s;
  }
  return *this;
}

const CartesianPoint CartesianPoint::operator+(const CartesianPoint &other) const {
  CartesianPoint result = *this;
  result += other;
  return result;
}

const CartesianPoint CartesianPoint::operator-(const CartesianPoint &other) const {
  CartesianPoint result = *this;
  result -= other;
  return result;
}

const CartesianPoint CartesianPoint::operator-() const {
  return CartesianPoint(size(), 0) - *this;
}

const CartesianPoint CartesianPoint::operator*(const mstreal& s) const {
  CartesianPoint result = *this;
  result *= s;
  return result;
}

const CartesianPoint CartesianPoint::operator/(const mstreal& s) const {
  CartesianPoint result = *this;
  result /= s;
  return result;
}

mstreal CartesianPoint::norm() const {
  mstreal n = 0;
  for (int i = 0; i < size(); i++) n += (*this)[i]*(*this)[i];
  return sqrt(n);
}

mstreal CartesianPoint::norm2() const {
  mstreal n = 0;
  for (int i = 0; i < size(); i++) n += (*this)[i]*(*this)[i];
  return n;
}

mstreal CartesianPoint::mean() const {
  return sum()/size();
}

mstreal CartesianPoint::stdev() const {
  return sqrt(var());
}

mstreal CartesianPoint::var() const {
  mstreal m = mean();
  return norm2()/size() - m*m;
}

mstreal CartesianPoint::sum() const {
  mstreal s = 0;
  for (int i = 0; i < size(); i++) s += (*this)[i];
  return s;
}

mstreal CartesianPoint::median() const {
  if (size() == 0) return 0; // could also throw an error in this case
  vector<mstreal> vec = *this;
  sort(vec.begin(), vec.end());
  if (vec.size() % 2 == 0) return (vec[vec.size() / 2] + vec[vec.size() / 2 - 1])/2;
  return vec[vec.size() / 2];
}

CartesianPoint CartesianPoint::cross(const CartesianPoint& other) const {
  if (size() != 3) MstUtils::error("don't know how to compute cross produces for dimensions other than 3", "CartesianPoint::cross");
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::cross");

  CartesianPoint C(3, 0);
  C[0] = getY()*other.getZ() - getZ()*other.getY();
  C[1] = getZ()*other.getX() - getX()*other.getZ();
  C[2] = getX()*other.getY() - getY()*other.getX();

  return C;
}

CartesianPoint CartesianPoint::elemProd(const CartesianPoint& other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::elemProd");
  CartesianPoint P(size());
  for (int i = 0; i < size(); i++) P[i] = (*this)[i] * other[i];
  return P;
}

mstreal CartesianPoint::dot(const CartesianPoint& other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::dot");

  mstreal d = 0;
  for (int i = 0; i < size(); i++) {
    d += (*this)[i] * other[i];
  }

  return d;
}

mstreal CartesianPoint::distance(const CartesianPoint& another) const {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance(CartesianPoint&)");
  mstreal d = 0;
  for (int i = 0; i < this->size(); i++) {
    d += ((*this)[i] - another[i])*((*this)[i] - another[i]);
  }
  return sqrt(d);
}

mstreal CartesianPoint::distancenc(const CartesianPoint& another) const {
  mstreal d = 0, dd;
  const CartesianPoint& p = *this;
  for (int i = 0; i < p.size(); i++) {
    dd = p[i] - another[i];
    d += dd*dd;
  }
  return sqrt(d);
}

mstreal CartesianPoint::distance2(const CartesianPoint& another) const {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance2(CartesianPoint&)");
  mstreal d = 0, dd;
  for (int i = 0; i < this->size(); i++) {
    dd = (*this)[i] - another[i];
    d += dd*dd;
  }
  return d;
}

mstreal CartesianPoint::distance2nc(const CartesianPoint& another) const {
  mstreal d = 0, dd;
  const CartesianPoint& p = *this;
  for (int i = 0; i < p.size(); i++) {
    dd = p[i] - another[i];
    d += dd*dd;
  }
  return d;
}

/* --------- CartesianGeometry --------- */
mstreal CartesianGeometry::dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, bool radians) {
  CartesianPoint AB = _p1 - _p2;
  CartesianPoint CB = _p3 - _p2;
  CartesianPoint DC = _p4 - _p3;

  if (AB.norm() == 0.0 || CB.norm() == 0.0 || DC.norm() == 0.0) MstUtils::error("some points coincide in dihedral calculation", "CartesianGeometry::dihedralRadians");

  CartesianPoint ABxCB = AB.cross(CB).getUnit();
  CartesianPoint DCxCB = DC.cross(CB).getUnit();

  // the following is necessary for values very close to 1 but just above
  double dotp = ABxCB * DCxCB;
  if (dotp > 1.0) {
    dotp = 1.0;
  } else if (dotp < -1.0) {
    dotp = -1.0;
  }

  double angle = acos(dotp);
  if (ABxCB * DC > 0) {
    angle *= -1;
  }
  if (!radians) angle *= 180/M_PI;
  return angle;
}

mstreal CartesianGeometry::dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4, bool radians) {
  return dihedral(*_p1, *_p2, *_p3, *_p4, radians);
}

/* --------- MstUtils --------- */
void MstUtils::openFile (fstream& fs, string filename, ios_base::openmode mode, string from) {
  fs.open(filename.c_str(), mode);
  if (!fs.is_open()) {
    if (!from.empty()) from += " -> ";
    MstUtils::error("could not open file '" + filename + "'", from + "MstUtils::openFile");
  }
}

string MstUtils::uc(const string& str){
  string ret = str;
  for (int i = 0; i < ret.length(); i++) {
    ret[i] = toupper(ret[i]);
  }
  return ret;
}

string MstUtils::lc(const string& str){
  string ret = str;
  for (int i = 0; i < ret.length(); i++) {
    ret[i] = tolower(ret[i]);
  }
  return ret;
}

string MstUtils::trim(const string& str, string delimiters) {
  int i = str.find_first_not_of(delimiters);
  if (i == string::npos) return "";
  int j = str.find_last_not_of(delimiters);
  return str.substr(i, j - i + 1);
}
vector<string> MstUtils::trim(const vector<string>& strings, string delimiters) {
  vector<string> ret = strings;
  for (int i = 0; i < strings.size(); i++) ret[i] = MstUtils::trim(strings[i], delimiters);
  return ret;
}


void MstUtils::warn(const string& message, string from) {
  string head = from.empty() ? "Warning: " : "Warning in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;
}

void MstUtils::error(const string& message, string from, int code) {
  string head = from.empty() ? "Error: " : "Error in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;

  // print backtrace
  void *array[100];
  size_t size = backtrace(array, 100);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  exit(code);
}

void MstUtils::assert(bool condition, string message, string from, int exitCode) {
  if(!condition) {
    MstUtils::error(message, from, exitCode);
  }
}

string MstUtils::wrapText(const string& message, int width, int leftSkip, int startingOffset) {
  string text;
  int b = 0, e, off = startingOffset, word = 0;
  while (b < message.size()) {
    // find the end of the next word
    e = message.find_first_of(" ", b);
    if (e == string::npos) e = message.size();
    int n = e - b;
    // if including the next word on this line will go over the width, start a new line
    if ((off + n >= width) && (word > 0)) {
      text += "\n" + string(leftSkip, ' ');
      off = leftSkip;
      word = 0;
    }
    text += message.substr(b, n) + " ";
    off += n + 1;
    b = e+1;
    word++;
  }
  return text;
}

char* MstUtils::copyStringC(const char* str) {
  char* copy = new char[strlen(str) + 1];
  strcpy(copy, str);
  return copy;
}

int MstUtils::toInt(const string& num, bool strict) {
  int ret = 0;
  try { ret = stoi(num); }
  catch (...) {
    if (strict) MstUtils::error("failed to convert '" + num + "' to integer", "MstUtils::toInt");
  }
  return ret;
}

bool MstUtils::isReal(const string& num) {
  double ret;
  return (sscanf(num.c_str(), "%lf", &ret) == 1);
}

MST::mstreal MstUtils::toReal(const string& num, bool strict) {
  double ret = 0.0;
  if ((sscanf(num.c_str(), "%lf", &ret) != 1) && strict) MstUtils::error("failed to convert '" + num + "' to mstreal", "MstUtils::toReal");
  return (mstreal) ret;
}
