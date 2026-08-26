#ifndef _MSTTYPES_H
#define _MSTTYPES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <locale>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <limits>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <chrono>
#include <execinfo.h>
#include <signal.h>
#undef assert

using namespace std;

namespace MST {

// forward declarations
class Chain;
class Residue;
class Atom;
class Structure;
class CartesianPoint;

typedef double mstreal;

class Structure {
  friend class Chain;

  public:
    Structure();
    Structure(string pdbFile, string options = "");
    Structure(istream& is, string options = "");
    Structure(const Structure& S);
    Structure(Chain& C);
    Structure(Residue& R);
    Structure(const vector<Atom*>& atoms);
    Structure(const vector<Residue*>& residues);
    ~Structure();

    void readPDB(const string& pdbFile, string options = "");
    void readPDB(istream& is, string options = "");
    void writePDB(const string& pdbFile, string options = "") const;
    void writePDB(ostream& ofs, string options = "") const;
    void reset();
    Structure& operator=(const Structure& A);

    int chainSize() const { return chains.size(); }
    int residueSize() const { return numResidues; }
    int positionSize() const { return residueSize(); }  // for interchangability with MSL
    int atomSize() const { return numAtoms; }
    Chain* getChainByID(string id) { return (chainsByID.find(id) != chainsByID.end()) ? chainsByID[id] : NULL; }
    Chain* getChainBySegID(string id) { return (chainsBySegID.find(id) != chainsBySegID.end()) ? chainsBySegID[id] : NULL; }
    Chain& getChain(int i) const { return (*this)[i]; }
    Residue& getResidue(int i) const;
    Chain& operator[](int i) const { return *(chains[i]); }
    vector<Atom*> getAtoms() const;
    vector<Residue*> getResidues() const;
    void setName(const string& _name) { name = _name; }
    string getName() const { return name; }
    void renumber(int startResNum=1, int startAtomIndex=1); // make residue numbering consequitive in each chain and atom index consequitive throughout

    /* ----- functions that grow/shrink structure ----- */
    /* returns false if the chain name collides with existing chain names and no suitable single-letter
     * chain name was found as replacement OR if renaiming was not allowed. This could still mean that
     * a multi-character name is picked that is unique, but that's not technically correct for output,
     * so false will still be returned. Note that in the latter case, the segment ID will be renamed as
     * well to be the same multi-character name. So, although in the output chain names will repeat,
     * segment names will still be unique. If it fails to find an even multi-character name, errors. */
    bool appendChain(Chain* C, bool allowRename = true);
    Chain* appendChain(string cid, bool allowRename = true);
    void deleteChain(Chain* chain);

    /* makes a copy of the atom, then decides where it is supposed to go and inserts it
     * into the Structure, creating a new Chain and/or Residue as needed. */
    void addAtom(Atom* A);
    void addAtom(Atom& A) { addAtom(&A); }
    void addAtoms(vector<Atom*> atoms) { addAtoms(&atoms); }
    void addAtoms(vector<Atom*>* atoms);

    /* makes a copy of the residue, then decides where it is supposed to go and
     * inserts it into the Structure, creating a new Chain as needed. */
    Residue* addResidue(Residue* res);

    /* ----- functions that grow/shrink structure ----- */

    int getResidueIndex(Residue* res);

    /* == and != operators are needed to convert vector<Structure> into python
     * lists via boost.python. This is because python lists are quite a bit more
     * powerful than C++ vectors, enabling, for example, contains queries. */
    bool operator==(const Structure& other) { return (this == &other); }
    bool operator!=(const Structure& other) { return (this != &other); }

  protected:
    void incrementNumAtoms(int delta = 1) { numAtoms += delta; }
    void incrementNumResidues(int delta = 1) { numResidues += delta; }
    void deletePointers();
    void copy(const Structure& S);

  private:
    vector<Chain*> chains;
    string name;
    int numResidues, numAtoms;
    // NOTE: thse two maps are maintained for convenience and will not guarantee the lack of collisions. That is,
    // if more than one chain use the same ID or segment ID, these maps will only store the last one added.
    map<string, Chain*> chainsByID;
    map<string, Chain*> chainsBySegID;
};

class Chain {
  friend class Residue;
  friend class Structure;

  public:
    Chain();
    Chain(const Chain& C);
    Chain(const string& chainID, const string& segID);
    Chain(const string& chainID, const string& segID, const vector<Residue*>& residues);
    ~Chain();

    int residueSize() const { return residues.size(); }
    int positionSize() const { return residueSize(); }  // for interchangability with MSL
    int atomSize() const { return numAtoms; }
    Residue& operator[](int i) const { return *(residues[i]); }
    Residue& getResidue(int i) const { return (*this)[i]; }
    vector<Residue*> getResidues() { return residues; }
    vector<Atom*> getAtoms();
    string getID() const { return cid; }
    string getSegID() const { return sid; }
    Structure* getParent() const { return parent; }
    Structure* getStructure() const { return getParent(); }
    int getResidueIndex(const Residue* res); // this in the in-chain residue index!
    int getIndex() const; // index of this chain within the larger structure

    /* convenience functoins, not efficient (linear search). If you need to do this a lot,
     * call getResidues() and construct your own data structure (e.g., a map<>) for fast lookups. */
    Residue* findResidue(string resname, int resnum);
    Residue* findResidue(string resname, int resnum, char icode);

    void setID(string _cid) { cid = _cid; }
    void setSegID(string _sid) { sid = _sid; }

    /* ----- functions that grow/shrink structure ----- */
    void appendResidue(Residue* R);
    void insertResidue(Residue* R, int index); // insert the Residue in such a way that it ends up being at index i
    void appendResidueCopies(const vector<Residue*>& residues);
    Residue* insertResidueCopy(Residue* R, int index = -1); // same, but copies the residue first
    Residue* insertResidueCopy(Residue& R, int index = -1); // same, but copies the residue first
    /* ----- functions that grow/shrink structure ----- */

  protected:
    void setParent(Structure* p) { parent = p; } // will not itself update residue/atom counts in parent
    void incrementNumAtoms(int delta = 1);

  private:
    vector<Residue*> residues;
    map<Residue*, int> residueIndexInChain; // to enable quick look-ups of up/down-stream residues
    Structure* parent;
    int numAtoms;
    string cid, sid;
};

class Residue {
  friend class Structure;
  friend class Chain;
  friend class Atom;

  public:
    Residue();
    Residue(const Residue& R, bool copyAlt = true);
    Residue(string _resname, int _resnum, char _icode = ' ');
    ~Residue();

    int atomSize() const { return atoms.size(); }
    vector<Atom*> getAtoms() { return atoms; }
    Atom& operator[](int i) const { return *(atoms[i]); }
    Atom& getAtom(int i) const { return *(atoms[i]); }
    Chain* getChain() const { return parent; }
    string getChainID(bool strict = true);
    string getName() const { return resname; }
    int getNum() const { return resnum; }
    char getIcode() const { return icode; }
    bool isNamed(const string& _name) const { return (resname.compare(_name) == 0); }
    bool isNamed(const char* _name) const { return (strcmp(resname.c_str(), _name) == 0); }
    Atom* findAtom(string _name, bool strict = true) const; // returns NULL if not found and if strict is false
    bool atomExists(string _name) { return (findAtom(_name, false) != NULL); } // mostly for interchangeability with MSL, better to use findAtom and check for NULL
    Chain* getParent() const { return parent; }
    Structure* getStructure() const { return (parent == NULL) ? NULL : parent->getParent(); }

    void setName(const char* _name) { resname = (string) _name; }
    void setName(const string& _name) { resname = _name; }
    void setIcode(char _icode) { icode = _icode; }
    void setNum(int num) { resnum = num; }
    void copyAtoms(Residue& R, bool copyAlt = true);
    void copyAtoms(const vector<Atom*>& _atoms, bool copyAlt = true);
    void compactify() { atoms.shrink_to_fit(); } // saves memory by adjusting capacity to match length

    /* for all atoms in this residue, overwrite the main coordinate set with the
     * coordinate set from the alternative with the specified index. */
    void makeAlternativeMain(int altInd);

    /* ----- functions that grow/shrink structure ----- */
    void appendAtom(Atom* A);
    void appendAtoms(const vector<Atom*>& A);
    void deleteAtoms();
    void deleteAtom(int ind);

    /* replaces the residue's atom vector with the given vector of atoms. By default,
     * all old atoms are deleted (i.e., removed from the residue's atom vector and
     * destructed). However, if the second argument is passed, will only delete the
     * atoms that were at the specified indices in the old atom array. Note that
     * this function is flexible enough to do things like erase a set of one or more
     * atoms, insert a set of one or more atoms, both, replace the entire
     * set of atoms with a new set, destroying the old ones, or even simply change
     * the permutation of the existing atoms. The order of atoms in the new vector
     * will be: any old ones that survived, in their initial order, followed by any
     * newly added atoms, in the specified order. */
    void replaceAtoms(const vector<Atom*>& newAtoms, vector<int>* oldAtoms = NULL);
    // same as above, but old Atoms are identified via Atom pointers
    void replaceAtoms(const vector<Atom*>& newAtoms, const vector<Atom*>& oldAtoms);
    /* ----- end functions that grow/shrink structure -- */

    int getResidueIndex() const;
    int getResidueIndexInChain() const;

    friend ostream & operator<<(ostream &_os, const Residue& _res) {
      if (_res.getParent() != NULL) {
        _os << _res.getParent()->getID() << ",";
      }
      _os << _res.getNum() << " " << _res.getName();
      return _os;
    }
    // so that sets of residues, for example, are sorted in the right order
    friend bool operator<(const Residue& r1, const Residue& r2) {
      return r1.getResidueIndex() < r2.getResidueIndex();
    }

  protected:
    void setParent(Chain* _parent) { parent = _parent; } // will not itself update residue/atom counts in parents

  private:
    vector<Atom*> atoms;
    Chain* parent;
    int resnum;
    string resname;
    char icode;
};

class Atom {
  friend class Structure;
  friend class Chain;
  friend class Residue;
  friend class AtomPointerVector;

  public:
    Atom();
    Atom(const Atom& A, bool copyAlt = true);
    Atom(const Atom* A, bool copyAlt = true) : Atom(*A, copyAlt) {}
    Atom(int _index, const string& _name, mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
    ~Atom();

    mstreal getX() const { return x; }
    mstreal getY() const{ return y; }
    mstreal getZ() const{ return z; }
    bool hasInfo() const { return info != NULL; }
    mstreal& operator[](int i);
    const mstreal& operator[](int i) const;
    CartesianPoint getCoor() const;
    CartesianPoint getAltCoor(int altInd) const;
    mstreal getAltB(int altInd) const { return info->getAltB(altInd); }
    mstreal getAltOcc(int altInd) const { return info->getAltOcc(altInd); }
    char getAltLocID(int altInd) const  { return info->getAltLocID(altInd); }
    mstreal getB() const { return info->B; }
    mstreal getOcc() const { return info->occ; }
    string getName() const { return string(info->name); }
    char* getNameC() { return info->name; }
    bool isHetero() const { return info->het; }
    int getIndex() const { return info->index; }
    char getAlt() const { return info->alt; }
    void setAlt(char a) const { info->alt = a; }
    bool isNamed(const char* _name) const { return (strcmp(info->name, _name) == 0); }
    bool isNamed(const string& _name) const { return isNamed(_name.c_str()); }
    bool hasAlternatives() const { return (info->alternatives != NULL); }
    int numAlternatives() const { return (info->alternatives == NULL) ? 0 : info->alternatives->size(); }
    Residue* getParent() const { return info->parent; }
    Residue* getResidue() const { return info->parent; }
    Chain* getChain() const { return (info->parent == NULL) ? NULL : info->parent->getParent(); }
    Structure* getStructure() { Chain* chain = getChain(); return (chain == NULL) ? NULL : chain->getParent(); }

    void setName(const char* _name) { info->setName(_name); }
    void setName(const string& _name) { info->setName(_name); }
    void setX(mstreal _x) { x = _x; }
    void setY(mstreal _y) { y = _y; }
    void setZ(mstreal _z) { z = _z; }
    void setCoor(mstreal _x, mstreal _y, mstreal _z) { x = _x; y = _y; z = _z; }
    void setCoor(const CartesianPoint& xyz);
    void setCoor(const Atom& a);
    void setCoor(const Atom* a) { setCoor(*a); }
    void setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z) { info->setAltCoor(ai, _x, _y, _z); }
    void setOcc(mstreal _occ) { info->occ = _occ; }
    void setB(mstreal _B) { info->B = _B; }
    void seetHet(bool _het) { info->het = _het; }
    void setIndex(int _index) { info->index = _index; }

    /* make the alternative with the specified index the main one, making the current
     * main position the alternative with the specified index. Calling this twice with
     * the same index will return things back to the way they were originally. */
    void swapWithAlternative(int altInd);

    /* overwrite the main coordinate set with the coordinate set from the alternative
     * with the specified index. */
    void makeAlternativeMain(int altInd);

    void addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt = ' ') { info->addAlternative(_x, _y, _z, _B, _occ, _alt); }
    void addAlternative(const Atom& a) { addAlternative(a.getX(), a.getY(), a.getZ(), a.getB(), a.getOcc()); }
    void removeLastAlternative() { info->removeLastAlternative(); }
    void removeAlternative(int i) { info->removeAlternative(i); }
    void clearAlternatives() { info->clearAlternatives(); }

    string pdbLine() { return pdbLine((info->parent == NULL) ? 1 : info->parent->getNum(), getIndex()); }
    string pdbLine(int resIndex, int atomIndex);

    mstreal distance(const Atom& another) const;
    mstreal distance(const Atom* another) const { return distance(*another); }
    mstreal distance2(const Atom& another) const;
    mstreal distance2(const Atom* another) const { return distance2(*another); }
    mstreal angle(const Atom& A, const Atom& B, bool radians = false) const;
    mstreal angle(const Atom* A, const Atom* B, bool radians = false) const;
    mstreal dihedral(const Atom& A, const Atom& B, const Atom& C, bool radians = false) const;
    mstreal dihedral(const Atom* A, const Atom* B, const Atom* C, bool radians = false) const;

    /* Sets the coordinates of the atom based on internal coordinates relative
     * to three other atoms: thA, anA, diA, A (this atom). The internal
     * coordinates are distance diA-A (di), angle anA-diA-A (an), and dihedral
     * angle thA - anA - diA - A (th). */
    bool build(const Atom& diA, const Atom& anA, const Atom& thA, mstreal di, mstreal an, mstreal th, bool radians = false);
    bool build(const Atom* diA, const Atom* anA, const Atom* thA, mstreal di, mstreal an, mstreal th, bool radians = false) {
      return build(*diA, *anA, *thA, di, an, th, radians);
    }

    friend ostream & operator<<(ostream &_os, const Atom& _atom);

    /* Strips the Atom of all of its additional information, besides the 3D coor-
     * dinate. This can be useful when memory needs to be preserved and all the
     * other info is not necessary. NOTE: after calling this function, attempts
     * to access any properties other than coordinates will lead to undefined
     * behavior, including segmentation faults. For efficiency, checking for
     * whether the information exists or has been stripepd is not performed. */
    void stripInfo();

  protected:
    void setParent(Residue* _parent) { info->parent = _parent; } // will not itself update residue/atom counts in parents

  private:
    mstreal x, y, z;
    class atomInfo {
      public:
        // data structure for storing information about alternative atom locations
        class altInfo {
          public:
            altInfo() { x = y = z = occ = B = 0; alt = ' '; }
            altInfo(const altInfo& A) { x = A.x; y = A.y; z = A.z; B = A.B; occ = A.occ; alt = A.alt; }
            altInfo(mstreal _x, mstreal _y, mstreal _z, mstreal _occ, mstreal _B, char _alt) { x = _x; y = _y; z = _z; B = _B; occ = _occ; alt = _alt; }
            mstreal x, y, z, occ, B;
            char alt;
        };

        atomInfo();
        atomInfo(const atomInfo& other, bool copyAlt = true);
        atomInfo(int _index, const string& _name, mstreal _B, mstreal _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
        ~atomInfo();

        void setName(const char* _name);
        void setName(const string& _name) { setName(_name.c_str()); }
        CartesianPoint getAltCoor(int altInd) const;
        mstreal getAltB(int altInd) const;
        mstreal getAltOcc(int altInd) const;
        char getAltLocID(int altInd) const;
        void setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z);
        void addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt);
        void removeLastAlternative();
        void removeAlternative(int i);
        void clearAlternatives();

        mstreal occ, B;
        char *name, alt;
        Residue* parent;
        bool het;
        int index;
        vector<altInfo>* alternatives; /* since this is a pointer, and will be NULL for most atoms, it's fine
                                        * to use vector here in terms of memory, but very convenient for use */
    };
    atomInfo* info;
};
ostream & operator<<(ostream &_os, const Atom& _atom); // this just to silence a silly compiler warning

/* The following several classes look and feel like MSL classes, BUT (importantly) their
 * use is absolutely optional, and none of the basic MST datastructures use them. On
 * the other hand, they can be created from those basic types for a similar use as
 * in MSL when needed (and only when needed). For example, CartesianPoint knows how
 * to construct itself from atom, both via a constructor and assignment operator.
 * Similarly, AtomPointerVector knows how to construct itself from vector<Atom*>. */
class CartesianPoint : public vector<mstreal> {
  /* this class it no limited to 3D vectors, though some of the functions will only
   * work with 3D vectors. The intention is to make it general, such that if a 3D
   * vector is required (or another dimension), appropriate assertions are made. */
  public:
    // inherit a bunch of useful constructors from vector
    CartesianPoint() : vector<mstreal>() { }
    CartesianPoint(size_t sz) : vector<mstreal>(sz) { }
    CartesianPoint(size_t sz, mstreal val) : vector<mstreal>(sz, val) { }
    CartesianPoint(const CartesianPoint& other) : vector<mstreal>(other) { }
    CartesianPoint(const vector<mstreal>& other) : vector<mstreal>(other) { }
    CartesianPoint(mstreal x, mstreal y, mstreal z) : vector<mstreal>(3, 0) { (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; }

    CartesianPoint(const Atom& A);
    CartesianPoint(const Atom* A) : CartesianPoint(*A) {}
    CartesianPoint& operator+=(const CartesianPoint& rhs);
    CartesianPoint& operator-=(const CartesianPoint& rhs);
    CartesianPoint& operator/=(const mstreal& s);
    CartesianPoint& operator*=(const mstreal& s);
    const CartesianPoint operator+(const CartesianPoint &other) const;
    const CartesianPoint operator-(const CartesianPoint &other) const;
    const CartesianPoint operator*(const mstreal& s) const;
    const CartesianPoint operator/(const mstreal& s) const;
    const CartesianPoint operator-() const;
    // CartesianPoint& operator=(const Atom& A);
    const double operator*(const CartesianPoint& other) const { return this->dot(other); }

    mstreal norm() const;
    mstreal norm2() const;
    mstreal mean() const;
    mstreal stdev() const;
    mstreal var() const;
    mstreal median() const;
    mstreal sum() const;
    CartesianPoint cross(const CartesianPoint& other) const;
    mstreal dot(const CartesianPoint& other) const;
    CartesianPoint getUnit() const { double L = norm(); return (*this/L); };
    CartesianPoint elemProd(const CartesianPoint& other) const; // element-wise product

    // a few special access operations
    mstreal getX() const { return (*this)[0]; }
    mstreal getY() const { return (*this)[1]; }
    mstreal getZ() const { return (*this)[2]; }

    mstreal distance(const CartesianPoint& another) const;
    mstreal distance(const CartesianPoint* another) const { return distance(*another); }
    mstreal distancenc(const CartesianPoint& another) const; // no size check (for speed)
    mstreal distancenc(const CartesianPoint* another) const { return distancenc(*another); } // no size check (for speed)
    mstreal distance2(const CartesianPoint& another) const;
    mstreal distance2(const CartesianPoint* another) const { return distance2(*another); }
    mstreal distance2nc(const CartesianPoint& another) const; // no size check (for speed)
    mstreal distance2nc(const CartesianPoint* another) const { return distance2nc(*another); } // no size check (for speed)

    friend ostream & operator<<(ostream &_os, const CartesianPoint& _p) {
      for (int i = 0; i < _p.size(); i++) {
        _os << _p[i];
        if (i != _p.size() - 1) _os << " ";
      }
      return _os;
    }
};

class CartesianGeometry {
  public:
    static mstreal dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, bool radians = false);
    static mstreal dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4, bool radians = false);
};

class AtomPointerVector : public vector<Atom*> {
  public:
    // inherit a bunch of useful constructors from vector
    AtomPointerVector() : vector<Atom*>() { }
    AtomPointerVector(size_t sz, Atom* val = NULL) : vector<Atom*>(sz, val) { }
    AtomPointerVector(const AtomPointerVector& other) : vector<Atom*>(other) { }
    AtomPointerVector(const vector<Atom*>& other) : vector<Atom*>(other) { }
    void copyCoordinates(const AtomPointerVector& other);

    using vector<Atom*>::push_back;    // base push_back of vector class
    void push_back(const Residue& R);  // overloaded push_back for Residues
    void push_back(const Residue* R) { push_back(*R); }

    CartesianPoint getGeometricCenter();
    void getGeometricCenter(mstreal& xc, mstreal& yc, mstreal& zc);
    void center();
    mstreal radiusOfGyration();
    // computes the smallest radius of a sphere centered at the centroid of the
    // atom vector, which encloses the set of all atoms.
    mstreal boundingSphereRadiusCent();
    void deletePointers();

    AtomPointerVector clone() const;
    void clone(AtomPointerVector& into) const;
    AtomPointerVector subvector(int beg, int end); // returns the range [beg, end)

    AtomPointerVector& operator+=(const AtomPointerVector& rhs);
    AtomPointerVector& operator-=(const AtomPointerVector& rhs);
    AtomPointerVector& operator/=(const mstreal& s);
    AtomPointerVector& operator*=(const mstreal& s);

    friend ostream & operator<<(ostream &_os, const AtomPointerVector& _atoms);
};
ostream & operator<<(ostream &_os, const AtomPointerVector& _atoms); // this just to silence a silly compiler warning
}

/* Utilities class, with a bunch of useful static functions, is defined outside of the MST namespace because:
 * 1) it really represents a different beast, not an MST type
 * 2) some of its functions (like assert) are likely to clash with function names in other project
 */
class MstUtils {
  public:
    static char* copyStringC(const char* str);
    static void openFile(fstream& fs, string filename, ios_base::openmode mode = ios_base::in, string from = "");
    static vector<string> trim(const vector<string>& strings, string delimiters = " \t\n\v\f\r");
    static string trim(const string& str, string delimiters = " \t\n\v\f\r");
    static void warn(const string& message, string from = "");
    static void error(const string& message, string from = "", int code = -1);
    static void assert(bool condition, string message = "error: assertion failed", string from = "", int exitCode = -1);
    static string uc(const string& str);                        // returns an upper-case copy of the input string
    static string lc(const string& str);                        // returns an lower-case copy of the input string
    static string wrapText(const string& message, int width, int leftSkip = 0, int startingOffset = 0);
    static int toInt(const string& num, bool strict = true);
    static MST::mstreal toReal(const string& num, bool strict = true);
    static bool isReal(const string& num);

    template <class T>
    static string toString(const T& obj) { return toString(&obj); }
    template <class T>
    static string toString(const T* obj);
};

template <class T>
string MstUtils::toString(const T* obj) {
  stringstream ss;
  ss << *obj;
  return ss.str();
}

#endif
