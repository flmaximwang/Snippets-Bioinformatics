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

void error(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	fprintf(stdout, "Error: %s.\n", wrapVsprintf(format, args).c_str());
	va_end(args);
	exit(-1);
}

bool exist(const char *path)
{
	struct stat buffer;
	if (stat(path, &buffer) == 0)
	{
		return true;
	}
	return false;
}

void file2array(string _filename, vector<string> & lines)
{
	FILE* ifp;
	int maxline = 1000;
	char *line, *tline;
	line = (char*) malloc(sizeof(char)*maxline);

	ifp = fopen(_filename.c_str(), "r");
	if (ifp == NULL)
	{
		error("unable to open file %s", _filename.c_str());
	}

	while (fgets(line, maxline, ifp) != NULL)
	{
		if (line[strlen(line)-1] != '\n')
		{
			error("lines in file %s are over %d long. Please increase max line limit and recompile", _filename.c_str(), maxline);
		}
		tline = trim(line);
		if (strlen(tline) > 0)
		{
			lines.push_back(line);
		}
		if (feof(ifp))
		{
			break;
		}
	}
	if (ferror(ifp))
	{
		error("an error occurred while reading file %s", _filename.c_str());
	}
	fclose(ifp);
	free(line);
}

string fileName(const string & fname, const bool & ext, const bool & path)
{
	if (ext && path)
	{
		return fname;
	}
	else
	{
		if (path)
		{
			if (fname.find_last_of(".") == string::npos)
			{
				return fname;
			}
			else
			{
				return fname.substr(0, fname.find_last_of("."));
			}
		}
		else
		{
			if (fname.find_last_of("/") == string::npos)
			{
				if (fname.find_last_of(".") == string::npos)
				{
					return fname;
				}
				else
				{
					if (ext)
					{
						return fname;
					}
					else
					{
						return fname.substr(0, fname.find_last_of("."));
					}
				}
			}
			else
			{
				if (fname.find_last_of(".") == string::npos)
				{
					return fname.substr(fname.find_last_of("/") + 1);
				}
				else
				{
					if (ext)
					{
						return fname.substr(fname.find_last_of("/") + 1);
					}
					else
					{
						return fname.substr(fname.find_last_of("/") + 1, 
											fname.find_last_of(".") - 1 - fname.find_last_of("/"));
					}
				}
			}
		}
	}
}

void findBrk(Structure & sys, vector<int> & brk)
{
	brk.clear();
	
	vector<Residue*> residues = sys.getResidues();
	for (int i = 0; i < (residues.size() - 1); i++)
	{
		if (hasBreak(*(residues[i]), *(residues[i + 1])))
		{
			brk.push_back(i);
		}
	}
}

bool hasBreak(Residue & cr, Residue & nr)
{
	Atom *A, *B;
	if ((A = cr.findAtom("C", false)) != NULL)
	{
		if ((B = nr.findAtom("N", false)) != NULL)
		{
			if (A->distance(B) > 2.0)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if ((B = nr.findAtom("NT", false)) != NULL)
		{
			if (A->distance(B) > 2.0)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	if (!((cr.getChainID() == nr.getChainID()) && ((cr.getNum() == nr.getNum() - 1) || (cr.getNum() == nr.getNum()))))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool hasFullBackbone(Residue & r)
{
	if (((r.atomExists("N")) || (r.atomExists("NT"))) && (r.atomExists("CA")) && (r.atomExists("C")) && 
		((r.atomExists("O")) || (r.atomExists("OT1")) || (r.atomExists("OT2")) || (r.atomExists("OXT"))))
	{
		return true;
	}
	else
	{
		return false;
	}
}

int isBackboneAtom(const char* an)
{
	if (strcmp(an, "N") == 0)
	{
		return 0;
	}
	if (strcmp(an, "NT") == 0)
	{
		return 10;
	}
	if (strcmp(an, "CA") == 0)
	{
		return 1;
	}
	if (strcmp(an, "C") == 0)
	{
		return 2;
	}
	if (strcmp(an, "O") == 0)
	{
		return 3;
	}
	if (strcmp(an, "OT1") == 0)
	{
		return 13;
	}
	if (strcmp(an, "OT2") == 0)
	{
		return 23;
	}
	if (strcmp(an, "OXT") == 0)
	{
		return 33;
	}

	return -1;
}

bool isDigit(const string & str)
{
	if (str.length() <= 0)
	{
		return false;
	}	
	int i;
	for (i = 0; i < str.length(); i++)
	{
		if (isdigit(str[i]))
		{
			continue;
		}
		return false;
	}
	return true;
}

bool isProtein(const vector<string> & legalAA, const string & resname)
{
	for(int i = 0; i < legalAA.size(); i++)
	{
		if (0 == legalAA[i].compare(resname))
		{
			return true;
		}
	}
	return false;
}

void openFileC (FILE* & fp, const string & fname, const char* mode)
{
	fp = fopen(fname.c_str(), mode);
	ASSERT(fp != NULL, "could not open file %s", fname.c_str());
}


void openFileCPP(fstream & fs, const string & fname, const ios_base::openmode & mode)
{
	fs.open(fname.c_str(), mode);
	ASSERT(fs.is_open(), "could not open file %s", fname.c_str());
}

string optionUsage(string opt, string mes, int w, int p1, int p2) {
	// first print the name of the option
	string text(p1, ' ');
	text += opt;
	if (p2 > text.size()) text += string(p2 - text.size(), ' ');

	// next print the description text
	int i = 0, k, L = text.size(), n, kk;
	bool newLine = false;
	while (i < mes.size()) {
		k = mes.find_first_of(" ", i);
		if (k == string::npos) k = mes.size();
		kk = mes.find_first_of("\n", i);
		if ((kk != string::npos) && (kk < k)) {
			mes[kk] = ' ';
			k = kk;
			newLine = true;
		}
		n = k - i;
		if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
		text += mes.substr(i, n) + " ";
		L += n + 1;
		i = k+1;
		if (newLine) {
			newLine = false;
			text += "\n" + string(p2, ' '); L = p2;
		}
	}
	return text;
}

string pad(string str, int len)
{
	if (str.size() == len)
  	{
		return str;
  	}
  	else
  	{
		if (str.size() < len)
		{
    		for (int i = 0; i < len - str.size(); i++)
    		{
				str += " ";
    		}
    		return str;
  		}
		else
		{
    		return str.substr(0, len);
		}
  	}
}

string suiteName() {
	return (string) "MASTER (Method of Accelerated Search for Tertiary Ensemble Representatives), version " + string(MASTER_VERSION);
}

vector<string> tokenize(const std::string & _input, const std::string & _delimiter, bool _allowEmtpy){
	vector<string> results;

	if (_input == "") {
		if (_allowEmtpy) {
			results.push_back(_input);
		}
		return results;
	}
	
	if (_allowEmtpy) {
		size_t prePos = 0;
		size_t pos  = _input.find(_delimiter);
		unsigned int delimiterSize = _delimiter.size();
		string left = _input, right;

		while (pos != std::string::npos) {
			results.push_back(left.substr(prePos, pos));
			if( pos + delimiterSize <= left.size() ) {
				left = left.substr(pos + delimiterSize, left.size() );
			} else {
				left = "";
			}
			pos  = left.find(_delimiter);
		}

		results.push_back(left);
	} else {
		int start  = _input.find_first_not_of(_delimiter);
		int end    = 0;
		string cur = _input;


		while (start != std::string::npos){
			end    = _input.find_first_of(_delimiter, start);
			results.push_back(_input.substr(start, end-start));
			start  = _input.find_first_not_of(_delimiter, end);
		}
	}
	return results;

}

char* trim(char* str)
{
	char* nptr;
	int i;
	for (i = 0; i < strlen(str); i++)
	{
		if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' '))
		{
			break;
		}
	}
	nptr = str + i;
	if (i == strlen(str))
	{
		return nptr;
	}

	for (i = strlen(str) - 1; i >= 0; i--)
	{
		if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' '))
		{
			str[i + 1] = '\0';
			break;
		}
	}
	return nptr;
}

// places a string end character after the last non-space character and returns a new
// point that points to the first non-space character
/*
char* trimSpace(char* str) {
  char* nptr; int i;
  for (i = 0; i < strlen(str); i++) {
    if (str[i] != ' ') { break; }
  }
  nptr = str + i;
  if (i == strlen(str)) { return nptr; }

  for (i = strlen(str)-1; i >= 0; i--) {
    if (str[i] != ' ') {
      str[i+1] = '\0';
      break;
    }
  }
  return nptr;
}
*/
char* trimSpace(char* str, const size_t & sz)
{
	char* nptr;
	int i;
	for (i = 0; i < sz; i++)
	{
		if (str[i] != ' ')
		{
			break;
		}
	}
	nptr = str + i;

	for (i = sz - 1; i >= 0; i--)
	{
		if (str[i] != ' ')
		{
			str[i + 1] = '\0';
			break;
		}
	}

	return nptr;
}

void warning(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	fprintf(stdout, "Warning: %s.\n", wrapVsprintf(format, args).c_str());
	va_end(args);
}

string wrapVsprintf(const char* format, va_list args)
{
	const unsigned int lenStr = 256;
	char* buffer = new char[lenStr];

	vsprintf(buffer, format, args);
	if ((strlen(buffer) + 1) > lenStr)
	{
		delete [] buffer;
		const unsigned int lenStrNew = (strlen(buffer) + 1);
		buffer = new char[lenStrNew];

		vsprintf(buffer, format, args);
	}

	string s = string(buffer);
	delete [] buffer;

	return s;
}

void writeString (fstream & ofs, const string & str, const bool & bin)
{
	if (bin)
	{
		ofs.write(str.c_str(), sizeof(char) * str.size());
	}
	else
	{
		ofs << str;
	}
}
