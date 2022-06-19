#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "classes/protein.h"

using namespace std;

enum VarType
{
	SV_INT,
	SV_FLOAT,
	SV_POINT,
	SV_STRING,
};

struct ScriptVar
{
	string name;
	Star value;
	VarType vt;
};

std::vector<string> script_lines;
ScriptVar script_var[256];
int vars = 0;
int program_counter = 0;


VarType type_from_name(const char* varname)
{
	switch (varname[0])
	{
		case '%':	return SV_INT;
		case '&':	return SV_FLOAT;
		case '@':	return SV_POINT;
		case '$':	return SV_STRING;
		default:	throw 0xbadfa12;
	}	
}

int find_var_index(const char* varname)
{
	int i;
	for (i=0; i<vars; i++)
	{
		if (!strcmp(script_var[i].name.c_str(), varname)) return i;
	}
	
	return -1;
}

Point interpret_Cartesian_literal(const char* param)
{
	char const* next;
	Point pt;
	
	pt.x = atof(&param[1]);
	next = strchr(&param[1], ',');
	if (next)
	{
		next = &next[1];
		pt.y = atof(next);
		
		next = strchr(next, ',');
		if (next)
		{
			pt.z = atof(next);
		}
	}
	
	return pt;
}

float interpret_single_float(const char* param)
{
	int n;
	Point pt;
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return 0;
		return script_var[n].value.n;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return 0;
		return script_var[n].value.f;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return 0;
		return script_var[n].value.ppt->magnitude();
		
		case '$':
		n = find_var_index(param);
		if (n<0) return 0;
		return atof( script_var[n].value.psz );
		
		case '[':
		pt = interpret_Cartesian_literal(param);
		return pt.magnitude();
		
		case '"':
		return atof(&param[1]);
		
		default:
		return atof(param);
	}
	
}

int interpret_single_int(const char* param)
{
	return round(interpret_single_float(param));
}

Point interpret_single_point(const char* param)
{
	int n;
	Point pt(0,0,0);
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return pt;
		pt.x = script_var[n].value.n;
		return pt;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return pt;
		pt.x = script_var[n].value.f;
		return pt;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return pt;
		return *(script_var[n].value.ppt);
		
		case '$':
		n = find_var_index(param);
		if (n<0) return pt;
		return interpret_Cartesian_literal( script_var[n].value.psz );
		
		case '[':
		return interpret_Cartesian_literal(param);
		
		case '"':
		return interpret_Cartesian_literal(&param[1]);
		
		default:
		return pt;
	}
}

char* interpret_single_string(const char* param)
{
	int n;
	char* buffer = new char[256];
	for (n=0; n<256; n++) buffer[n] = 0;
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return buffer;
		sprintf(buffer, "%d", script_var[n].value.n);
		return buffer;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return buffer;
		sprintf(buffer, "%f", script_var[n].value.f);
		return buffer;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return buffer;
		strcpy(buffer, script_var[n].value.ppt->printable().c_str());
		return buffer;
		
		case '$':
		n = find_var_index(param);
		if (n<0) return buffer;
		strcpy(buffer, script_var[n].value.psz);
		return buffer;
		
		case '[':
		strcpy(buffer, interpret_Cartesian_literal(param).printable().c_str());
		return buffer;
		
		case '"':
		strcpy(buffer, &param[1]);
		if (buffer[strlen(buffer)-1] == '"') buffer[strlen(buffer)-1] = '\0';
		return buffer;
		
		default:
		return buffer;
	}
}

int main(int argc, char** argv)
{
	int i, j, k, l, m, n;
	float f;
	char* psz;
	string script_fname = "";
	string PDB_fname = "";
	FILE* pf;
	
	for (i=1; i<argc; i++)
	{
		if (argv[i][0] == '-')
		{
			// TODO:
		}
		else
		{
			if (pf = fopen(argv[i], "rb"))
			{
				fclose(pf);
			
				char argvil[256];
				strcpy(argvil, argv[i]);
				for (j=0; argvil[j]; j++) if (argvil[j] >= 'A' && argvil[j] <= 'Z') argvil[j] |= 0x20;
				
				if (strstr(argvil, ".pdb"))
					PDB_fname = argv[i];
				else
					script_fname = argv[i];
			}
			else
			{
				cout << "Unrecognized argument or file not found: " << argv[i] << endl;
				return 0xbadf12e;
			}
		}
	}
	
	if (!PDB_fname.length() || !script_fname.length())
	{
		cout << "Usage:" << endl << "interpreter protein.pdb script_filename" << endl;
		cout << "interpreter script_filename protein.pdb" << endl;
		cout << endl;
		return 0;
	}
		
	for (i=0; i<256; i++)
	{
		script_var[i].name = "";
		script_var[i].value.n = 0;
	}
	vars = 0;
	
	char pdbname[1024];
	strcpy(pdbname, PDB_fname.c_str());
	char* dot = strchr(pdbname, '.');
	char* slash = strrchr(pdbname, '/');
	if (!slash) slash = strrchr(pdbname, '\\');
	
	script_var[vars].name = "$PDB";
	script_var[vars].value.psz = new char[strlen(pdbname)+4];
	strcpy(script_var[vars].value.psz, pdbname);
	vars++;
	
	if (dot) dot[0] = '\0';
	Protein p(&slash[1]);
	script_var[vars].name = "$PROTEIN";
	script_var[vars].value.psz = new char[strlen(slash)+4];
	strcpy(script_var[vars].value.psz, &slash[1]);
	vars++;
	if (dot) dot[0] = '.';
	
	pf = fopen(pdbname, "rb");
	p.load_pdb(pf);
	fclose(pf);
	
	if (pf = fopen(script_fname.c_str(), "rb"))
	{
		while (!feof(pf))
		{
			char buffer[1024];
			buffer[0] = '\0';
			fgets(buffer, 1023, pf);
			if (buffer[0] && buffer[0] != '#')
			{
				while (strlen(buffer) && buffer[strlen(buffer)-1] <= ' ') buffer[strlen(buffer)-1] = '\0';
        		script_lines.push_back(buffer);
			}
		}
		
		fclose(pf);
	}
	
	// TODO: Set some automatic variables.
	
	while (program_counter < script_lines.size())
	{
		char buffer[1024];
		char buffer1[1024];
		for (m=0; m<1024; m++) buffer[m] = buffer1[m] = '\0';
		strcpy(buffer, script_lines[program_counter].c_str());
		char** fields = chop_spaced_fields(buffer);
		if (fields && fields[0] && fields[0][0] && fields[0][1])
		{
			// Debug cout.
			/*cout << endl << "***Line: " << buffer << endl << "***Fields: ";
			for (n=0; fields[n]; n++) cout << (n?"|":"") << fields[n];
			cout << endl;*/
			
			if (fields[0][strlen(fields[0])-1] == ':') continue;
			
			// Interpret the script.
			if (!strcmp(fields[0], "HELIX"))
			{
				float phi, psi;
				l = 2;
				
				if (!strcmp(fields[1], "ALPHA")) { phi = ALPHA_PHI; psi = ALPHA_PSI; }
				else if (!strcmp(fields[1], "PI")) { phi = PI_PHI; psi = PI_PSI; }
				else if (!strcmp(fields[1], "3.10")) { phi = _310_PHI; psi = _310_PSI; }
				else if (!strcmp(fields[1], "PPRO1")) { phi = POLYPRO1_PHI; psi = POLYPRO1_PSI; }
				else if (!strcmp(fields[1], "PPRO2")) { phi = POLYPRO2_PHI; psi = POLYPRO2_PSI; }
				
				// Not technically helices, but included for consistency.
				else if (!strcmp(fields[1], "BETA")) { phi = BETA_PHI; psi = BETA_PSI; }
				else if (!strcmp(fields[1], "STRAIGHT")) { phi = M_PI; psi = M_PI; }
				
				else
				{
					phi = interpret_single_float(fields[1]);
					psi = interpret_single_float(fields[2]);
					l++;
				}
				
				int sr, er;
				sr = interpret_single_int(fields[l]);
				er = interpret_single_int(fields[l+1]);
				
				p.make_helix(sr, er, phi, psi);
				
			}	// HELIX
			
			else if (!strcmp(fields[0], "LET"))
			{
				n = find_var_index(fields[1]);
				if (n<0)n = vars++;
				script_var[n].name = fields[1];
				script_var[n].vt = type_from_name(fields[1]);
				
				switch (script_var[n].vt)
				{
					case SV_INT:
					if (!strcmp(fields[2], "=")) script_var[n].value.n = interpret_single_int(fields[3]);
					else if (!strcmp(fields[2], "+=")) script_var[n].value.n += interpret_single_int(fields[3]);
					else if (!strcmp(fields[2], "-=")) script_var[n].value.n -= interpret_single_int(fields[3]);
					else if (!strcmp(fields[2], "*=")) script_var[n].value.n *= interpret_single_int(fields[3]);
					else if (!strcmp(fields[2], "/=")) script_var[n].value.n /= interpret_single_int(fields[3]);
					else if (!strcmp(fields[2], "++")) script_var[n].value.n++;
					else if (!strcmp(fields[2], "--")) script_var[n].value.n--;
					else
					{
						cout << "Unimplemented operator " << fields[2] << " for int assignment." << endl;
						return 0x51974c5;
					}
					break;
					
					case SV_FLOAT:
					if (!strcmp(fields[2], "=")) script_var[n].value.f = interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "+=")) script_var[n].value.f += interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "-=")) script_var[n].value.f -= interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "*=")) script_var[n].value.f *= interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "/=")) script_var[n].value.f /= interpret_single_float(fields[3]);
					else
					{
						cout << "Unimplemented operator " << fields[2] << " for float assignment." << endl;
						return 0x51974c5;
					}
					break;
					
					case SV_POINT:
					script_var[n].value.ppt = new Point();
					if (!strcmp(fields[2], "=")) *(script_var[n].value.ppt) = interpret_single_point(fields[3]);
					else if (!strcmp(fields[2], "+=")) *(script_var[n].value.ppt) = script_var[n].value.ppt->add(interpret_single_point(fields[3]));
					else if (!strcmp(fields[2], "-=")) *(script_var[n].value.ppt) = script_var[n].value.ppt->subtract(interpret_single_point(fields[3]));
					else if (!strcmp(fields[2], "*="))
					{
						f = interpret_single_float(fields[3]);
						script_var[n].value.ppt->x *= f;
						script_var[n].value.ppt->y *= f;
						script_var[n].value.ppt->z *= f;
					}
					else if (!strcmp(fields[2], "/="))
					{
						f = interpret_single_float(fields[3]);
						script_var[n].value.ppt->x /= f;
						script_var[n].value.ppt->y /= f;
						script_var[n].value.ppt->z /= f;
					}
					else
					{
						cout << "Unimplemented operator " << fields[2] << " for Cartesian assignment." << endl;
						return 0x51974c5;
					}
					
					break;
					
					case SV_STRING:
					script_var[n].value.psz = new char(255);
					psz = interpret_single_string(fields[3]);
					if (!strcmp(fields[2], "=")) strcpy(script_var[n].value.psz, psz);
					else
					{
						delete[] psz;
						cout << "Unimplemented operator " << fields[2] << " for string assignment." << endl;
						return 0x51974c5;		// If you use your imagination, that says "syntax".
					}
					delete[] psz;
					break;
					
					default:
					;
				}
				
				
			}	// LET
			
			else if (!strcmp(fields[0], "ECHO"))
			{
				for (l=1; fields[l]; l++)
				{
					if (!strcmp(fields[l], "~")) goto _no_newline_on_echo;
					else
					{
						psz = interpret_single_string(fields[l]);
						cout << psz;
						delete[] psz;
					}
				}
				
				cout << endl;
				_no_newline_on_echo:
				;
			}	// ECHO
			
			else if (!strcmp(fields[0], "GOTO"))
			{
				sprintf(buffer1, "%s:", fields[1]);
				for (n=0; n<script_lines.size(); n++)
				{
					// debug cout
					psz = new char[256];
					strcpy(psz, script_lines[n].c_str());
					if (!strcmp(psz, buffer1))
					{
						delete[] psz;
						program_counter = n+1;
						goto _found_goto_target;
					}
					delete[] psz;
				}
				cout << "Label not found: \"" << buffer1 << '"' << endl;
				return 0x51974c5;
				
				_found_goto_target:
				continue;
			}	// GOTO
			
			else if (!strcmp(fields[0], "END") || !strcmp(fields[0], "EXIT") || !strcmp(fields[0], "QUIT"))
			{
				return 0;
			}	// END
			
			else
			{
				cout << "Unimplemented command: \"" << fields[0] << '"' << endl;
				return 0x51974c5;
			}
		}
		
		program_counter++;
	}
	
	return 0;
}










