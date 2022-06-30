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

#define _HAS_DOT 0x8000
#define _VARNUM_MASK 0xff
#define _VARFLAGS_MASK 0xff00

enum VarType
{
	SV_NONE,
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

string script_fname = "";
std::vector<string> script_lines;
ScriptVar script_var[256];
int vars = 0;
int program_counter = 0;

Protein p("TheProtein");

void raise_error(std::string message)
{
	cout << script_fname << ":" << program_counter+1 << " Error: ";
	cout << message << endl;
	throw 0xbadc0de;
}

VarType type_from_name(const char* varname)
{
	switch (varname[0])
	{
		case '%':	return SV_INT;
		case '&':	return SV_FLOAT;
		case '@':	return SV_POINT;
		case '$':	return SV_STRING;
		default:	raise_error("Variable names may only start with %, &, @, or $.");
	}
	return SV_NONE;	
}

int find_var_index(const char* varname)
{
	int i;
	
	for (i=0; i<vars; i++)
	{
		if (!strcmp(script_var[i].name.c_str(), varname)) return i;
	}
	
	char buffer[256];
	char* c;
	int flags = 0;
	strcpy(buffer, varname);
	
	if (c=strchr(buffer,',')) *c=0;
	else if (c=strchr(buffer,'.')) { flags |= _HAS_DOT; *c=0; }
	else if (c=strchr(buffer,':')) *c=0;
	
	for (i=0; i<vars; i++)
	{
		if (!strcmp(script_var[i].name.c_str(), buffer)) return i|flags;
	}
	
	return -1;
}

float interpret_single_float(const char* param);

Point interpret_Cartesian_literal(const char* param)
{
	char const* next;
	Point pt;
	
	// pt.x = atof(&param[1]);
	pt.x = interpret_single_float(&param[1]);
	next = strchr(&param[1], ',');
	if (next)
	{
		next = &next[1];
		// pt.y = atof(next);
		pt.y = interpret_single_float(next);
		
		next = strchr(next, ',');
		if (next)
		{
			// pt.z = atof(&next[1]);
			pt.z = interpret_single_float(&next[1]);
		}
	}
	
	return pt;
}

float interpret_single_float(const char* param)
{
	int n;
	Point pt;
	
	if ((param[0] >= '0' && param[0] <= '9')
		||
		param[0] == '-'
	   )
	{
		return atof(param);
	}
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return 0;
		n &= _VARNUM_MASK;
		return script_var[n].value.n;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return 0;
		n &= _VARNUM_MASK;
		return script_var[n].value.f;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return 0;
		if (n & _HAS_DOT)
		{
			param = strchr(param, '.');
			if (!param) return 0;
			param++;
			if (!strcmp(param, "x")) return script_var[n&_VARNUM_MASK].value.ppt->x;
			else if (!strcmp(param, "X")) return script_var[n&_VARNUM_MASK].value.ppt->x;
			else if (!strcmp(param, "y")) return script_var[n&_VARNUM_MASK].value.ppt->y;
			else if (!strcmp(param, "Y")) return script_var[n&_VARNUM_MASK].value.ppt->y;
			else if (!strcmp(param, "z")) return script_var[n&_VARNUM_MASK].value.ppt->z;
			else if (!strcmp(param, "Z")) return script_var[n&_VARNUM_MASK].value.ppt->z;
			else
			{
				raise_error((std::string)"Cartesian has no member named " + (std::string)param);
				return 0;
			}
		}
		n &= _VARNUM_MASK;
		return script_var[n].value.ppt->magnitude();
		
		case '$':
		n = find_var_index(param);
		if (n<0) return 0;
		n &= _VARNUM_MASK;
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
	AminoAcid* aa;
	
	if (param[0] >= '0' && param[0] <= '9')
	{
		aa = p.get_residue(atoi(param));
		if (aa)
		{
			pt = aa->get_CA_location();
			return pt;
		}
	}
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return pt;
		n &= _VARNUM_MASK;
		pt.x = script_var[n].value.n;
		aa = p.get_residue(script_var[n].value.n);
		if (aa) pt = aa->get_CA_location();
		return pt;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return pt;
		n &= _VARNUM_MASK;
		pt.x = script_var[n].value.f;
		return pt;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return pt;
		if (n &  _HAS_DOT)
		{
			param = strchr(param, '.');
			if (!param) return 0;
			param++;
			if (!strcmp(param, "x")) 		pt = Point(script_var[n&_VARNUM_MASK].value.ppt->x, 0, 0);
			else if (!strcmp(param, "X")) 	pt = Point(script_var[n&_VARNUM_MASK].value.ppt->x, 0, 0);
			else if (!strcmp(param, "y")) 	pt = Point(0, script_var[n&_VARNUM_MASK].value.ppt->y, 0);
			else if (!strcmp(param, "Y")) 	pt = Point(0, script_var[n&_VARNUM_MASK].value.ppt->y, 0);
			else if (!strcmp(param, "z")) 	pt = Point(0, 0, script_var[n&_VARNUM_MASK].value.ppt->z);
			else if (!strcmp(param, "Z")) 	pt = Point(0, 0, script_var[n&_VARNUM_MASK].value.ppt->z);
			else
			{
				raise_error((std::string)"Cartesian has no member named " + (std::string)param);
				return Point(0,0,0);
			}
			return pt;

		}
		n &= _VARNUM_MASK;
		return *(script_var[n].value.ppt);
		
		case '$':
		n = find_var_index(param);
		if (n<0) return pt;
		n &= _VARNUM_MASK;
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
	char* buffer = new char[65536];
	for (n=0; n<65536; n++) buffer[n] = 0;
	
	switch (param[0])
	{
		case '%':
		n = find_var_index(param);
		if (n<0) return buffer;
		n &= _VARNUM_MASK;
		sprintf(buffer, "%d", script_var[n].value.n);
		return buffer;
		
		case '&':
		n = find_var_index(param);
		if (n<0) return buffer;
		n &= _VARNUM_MASK;
		sprintf(buffer, "%f", script_var[n].value.f);
		return buffer;
		
		case '@':
		n = find_var_index(param);
		if (n<0) return buffer;
		if (n &  _HAS_DOT)
		{
			sprintf(buffer, "%f", interpret_single_float(param));
			return buffer;
		}
		n &= _VARNUM_MASK;
		strcpy(buffer, script_var[n].value.ppt->printable().c_str());
		return buffer;
		
		case '$':
		n = find_var_index(param);
		if (n<0) return buffer;
		n &= _VARNUM_MASK;
		if (!script_var[n].value.psz) return buffer;
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
	Point pt;
	std::string builder;
	string PDB_fname = "";
	FILE* pf;
		
	for (i=0; i<256; i++)
	{
		script_var[i].name = "";
		script_var[i].value.n = 0;
	}
	vars = 0;
	
	for (i=1; i<argc; i++)
	{
		script_var[vars].name = (string)"$arg" + std::to_string(i);
		script_var[vars].vt = SV_STRING;
		script_var[vars].value.psz = new char[strlen(argv[i])+4];
		strcpy(script_var[vars].value.psz, argv[i]);
		vars++;
		
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
				
				/*if (strstr(argvil, ".pdb"))
					PDB_fname = argv[i];
				else*/
				if (!script_fname.length())
					script_fname = argv[i];
			}
			else
			{
				cout << "Unrecognized argument or file not found: " << argv[i] << endl;
				return 0xbadf12e;
			}
		}
	}
	
	if (/*!PDB_fname.length() ||*/ !script_fname.length())
	{
		/*cout << "Usage:" << endl << "interpreter protein.pdb script_filename" << endl;
		cout << "interpreter script_filename protein.pdb" << endl;*/
		cout << "Error: no script filename supplied." << endl;
		cout << endl;
		return 0;
	}
	
	
	if (pf = fopen(script_fname.c_str(), "rb"))
	{
		while (!feof(pf))
		{
			char buffer[1024];
			buffer[0] = '\0';
			fgets(buffer, 1023, pf);
			if (buffer[0]) // && buffer[0] != '#')
			{
				while (strlen(buffer) && buffer[strlen(buffer)-1] <= ' ') buffer[strlen(buffer)-1] = '\0';
        		script_lines.push_back(buffer);
			}
		}
		
		fclose(pf);
	}
	
	
	while (program_counter < script_lines.size())
	{
		// cout << endl << "Process line " << program_counter << "... " << flush;
		char buffer[1024];
		char buffer1[1024];
		for (m=0; m<1024; m++) buffer[m] = buffer1[m] = '\0';
		strcpy(buffer, script_lines[program_counter].c_str());
		char** fields = chop_spaced_fields(buffer);
		if (fields && fields[0] && fields[0][0] && fields[0][1])
		{
			for (k=0; fields[k]; k++)
			{
				if (fields[k][0] == '#')
				{
					// cout << "Ignoring comment " << fields[k] << endl << flush;
					fields[k] = nullptr;
					break;
				}
			}
			// cout << "command " << fields[0] << endl << flush;
			
			// Debug cout.
			/*cout << endl << "***Line: " << buffer << endl << "***Fields: ";
			for (n=0; fields[n]; n++) cout << (n?"|":"") << fields[n];
			cout << endl;*/
			
			if (fields[0][strlen(fields[0])-1] == ':') goto _pc_continue;
			
			// Interpret the script.
			if (!strcmp(fields[0], "HELIX"))
			{
				float phi, psi;
				l = 2;
				
				if (!fields[1]) raise_error("No parameters given for HELIX.");
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
					if (!fields[2]) raise_error("Insufficient parameters given for HELIX.");
					psi = interpret_single_float(fields[2]);
					l++;
				}
				
				int sr, er;
				if (!fields[l]) raise_error("Insufficient parameters given for HELIX.");
				sr = interpret_single_int(fields[l]);
				if (!fields[l+1]) raise_error("Insufficient parameters given for HELIX.");
				er = interpret_single_int(fields[l+1]);
				if (fields[l+2]) raise_error("Too many parameters given for HELIX.");
				
				p.make_helix(sr, er, phi, psi);
				
			}	// HELIX
			
			else if (!strcmp(fields[0], "SEARCH"))
			{
				l = 1;
				int sr, er, esr;
				if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
				sr = interpret_single_int(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
				er = interpret_single_int(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
				psz = interpret_single_string(fields[l]);
				esr = er - strlen(psz);
				
				n = 0; k = 0;
				for (i=sr; i<esr; i++)
				{
					m = 0;
					for (j=0; psz[j]; j++)
					{
						m += p.get_residue(i+j)->similarity_to(psz[j]);
					}
					
					if (m > n)
					{
						k = i;
						n = m;
					}
				}
				
				delete[] psz;
				l++;
				
				n = find_var_index(fields[l]);
				if (n<0) n = vars++;
				n &= _VARNUM_MASK;
				if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
				if (fields[l+1]) raise_error("Too many parameters given for SEARCH.");
				script_var[n].name = fields[l];
				script_var[n].vt = type_from_name(fields[l]);
				
				switch (script_var[n].vt)
				{
					case SV_INT: script_var[n].value.n = k; break;
					case SV_FLOAT: script_var[n].value.f = k; break;
					case SV_STRING:
					builder = k;
					strcpy(script_var[n].value.psz, builder.c_str());
					break;
					
					default:
					raise_error("Bad destination variable type for SEARCH.");
					return 0xbadfa12;
				}
				
			}	// SEARCH
			
			else if (!strcmp(fields[0], "ALIGN"))
			{
				int sr, er, eachend;
				Point sp, ep;
				l = 1;
				if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
				sr = interpret_single_int(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
				er = interpret_single_int(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
				eachend = interpret_single_int(fields[l++]);
				if (!eachend) eachend = 1;
				if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
				sp = interpret_single_point(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
				ep = interpret_single_point(fields[l++]);
				if (fields[l]) raise_error("Too many parameters given for ALIGN.");
				
				// From start and end residues inwards for a total of eachend, average the CA locations.
				Point sl(0,0,0);
				for (i=0; i<eachend; i++)	// It's actually easier to average manually than to screw around with object arrays.
				{
					AminoAcid* aa = p.get_residue(sr+i);
					if (aa) sl = sl.add(aa->get_CA_location());
				}
				
				if (eachend > 1) sl.scale(sl.magnitude()/eachend);
				
				// Translate the range so that the starting average moves to the target start point.
				SCoord motion = sp.subtract(sl);
				for (i=sr; i<=er; i++)
				{
					AminoAcid* aa = p.get_residue(i);
					if (aa)
					{
						MovabilityType fmov = aa->movability;
						aa->movability = MOV_ALL;
						aa->aamove(motion);
						aa->movability = fmov;
					}
				}
				
				
				Point el(0,0,0);
				
				for (i=0; i<eachend; i++)
				{
					AminoAcid* aa = p.get_residue(er-i);
					if (aa) el = el.add(aa->get_CA_location());
				}
				
				if (eachend > 1) el.scale(el.magnitude()/eachend);
				
				// Rotate the range about the start point so the ending average moves to the target end point.
				Rotation rot = align_points_3d(&el, &ep, &sp);
				LocatedVector lv = rot.v;
				lv.origin = sp;
				
				// TODO: Offer the option to align middle residues to point toward some outside landmark.
				for (i=sr; i<=er; i++)
				{
					AminoAcid* aa = p.get_residue(i);
					if (aa)
					{
						MovabilityType fmov = aa->movability;
						aa->movability = MOV_ALL;
						aa->rotate(lv, rot.a);
						aa->movability = fmov;
					}
				}
				
			}	// ALIGN
			
			else if (!strcmp(fields[0], "CONNECT"))
			{
				l=1;
				int sr, er, ct, iters=250;
				if (!fields[l]) raise_error("Insufficient parameters given for CONNECT.");
				sr = interpret_single_int(fields[l++]);
				if (!fields[l]) raise_error("Insufficient parameters given for CONNECT.");
				ct = interpret_single_int(fields[l++]);
				if (fields[l]) iters = interpret_single_int(fields[l++]);
				if (fields[l]) raise_error("Too many parameters given for CONNECT.");
				er = ct - sgn(ct-sr);
				
				Atom *a1, *a2, *a3;
				AminoAcid *cta = p.get_residue(ct), *era = p.get_residue(er);
				
				if (!era || !cta) goto _no_connect;
				
				Point* pt3;
				if (ct > sr)
				{
					// Line up the C and O of er to the expected prevaa C and O of ct.
					a1 = era->get_atom("C");
					a2 = era->get_atom("O");
					pt3 = cta->predict_previous_COCA();
				}
				else
				{
					// Line up the N and HN (or substitute) of er to the expected nextaa N and HN of ct.
					a1 = era->get_atom("N");
					a2 = era->HN_or_substitute();
					pt3 = cta->predict_next_NHCA();
				}
				a3 = era->get_atom("CA");
				
				p.conform_backbone(sr, er, a1, pt3[0], a2, pt3[1], iters);
				delete[] pt3;
				p.backconnect(sr, er);
				
				_no_connect:
				goto _pc_continue;
			}	// CONNECT
			
			else if (!strcmp(fields[0], "GEN"))
			{
				if (!fields[1]) raise_error("No sequence given for GEN.");
				psz = interpret_single_string(fields[1]);
				
				p.add_sequence(psz);
				p.conform_backbone(1, p.get_seq_length(), 50);
				goto _prot_deets;
			} // GEN
			
			else if (!strcmp(fields[0], "DELETE"))
			{
				if (!fields[1]) raise_error("Insufficient parameters given for DELETE.");
				int sr = interpret_single_int(fields[1]), er=0;
				if (fields[2]) er = interpret_single_int(fields[2]);
				if (fields[3]) raise_error("Too many parameters given for DELETE.");
				
				if (er) p.delete_residues(sr, er);
				else p.delete_residue(sr);
			} // DELETE
			
			else if (!strcmp(fields[0], "LOAD"))
			{
				if (!fields[1]) raise_error("Insufficient parameters given for LOAD.");
				psz = interpret_single_string(fields[1]);
				if (fields[2]) raise_error("Too many parameters given for LOAD.");
				
				pf = fopen(psz, "rb");
				if (!pf)
				{
					raise_error( (std::string)"Failed to open " + (std::string)psz + (std::string)" for reading.");
					return 0xbadf12e;
				}
				p.load_pdb(pf);
				
				fclose(pf);
				
				_prot_deets:
				
				script_var[vars].name = "$PDB";
				script_var[vars].value.psz = new char[strlen(psz)+4];
				script_var[vars].vt = SV_STRING;
				strcpy(script_var[vars].value.psz, psz);
				vars++;
				
				const char* pname = p.get_name().c_str();
				script_var[vars].name = "$PROTEIN";
				script_var[vars].value.psz = new char[strlen(pname)+4];
				script_var[vars].vt = SV_STRING;
				strcpy(script_var[vars].value.psz, pname);
				vars++;
				
				delete[] psz;
	
				int seqlen = p.get_seq_length();
				script_var[vars].name = "%SEQLEN";
				script_var[vars].vt = SV_INT;
				script_var[vars].value.n = seqlen;
				vars++;
				
				script_var[vars].name = "$SEQUENCE";
				script_var[vars].value.psz = new char[seqlen+4];
				script_var[vars].vt = SV_STRING;
				strcpy(script_var[vars].value.psz, p.get_sequence().c_str());
				vars++;
				
				std::vector<std::string> rem_hx = p.get_remarks("650 HELIX");
				for (l=0; l<rem_hx.size(); l++)
				{
					char buffer[1024];
					char buffer1[1024];
					strcpy(buffer, rem_hx[l].c_str());
					char** fields = chop_spaced_fields(buffer);
					
					sprintf(buffer1, "%c%s.s", '%', fields[3]);
					script_var[vars].name = buffer1;
					script_var[vars].vt = SV_INT;
					script_var[vars].value.n = atoi(fields[4]);
					vars++;
					
					sprintf(buffer1, "%c%s.e", '%', fields[3]);
					script_var[vars].name = buffer1;
					script_var[vars].vt = SV_INT;
					script_var[vars].value.n = atoi(fields[5]);
					vars++;
					
					delete[] fields;
				}
			}
			
			else if (!strcmp(fields[0], "SAVE"))
			{
				if (!fields[1]) raise_error("Insufficient parameters given for SAVE.");
				psz = interpret_single_string(fields[1]);
				if (fields[3]) raise_error("Too many parameters given for SAVE.");
				
				pf = fopen(psz, "wb");
				if (!pf)
				{
					raise_error( (std::string)"Failed to open " + (std::string)psz + (std::string)" for writing.");
					return 0xbadf12e;
				}
				p.save_pdb(pf);
				p.end_pdb(pf);
				
				cout << "Wrote " << psz << "." << endl;
				
				fclose(pf);
				delete[] psz;
				
				if (fields[2])
				{
					if ( !strcmp("QUIT", fields[2]) || !strcmp("EXIT", fields[2]) || !strcmp("END", fields[2]) )
						return 0;
					else raise_error("Too many parameters given for SAVE.");
				}
			}
			
			else if (!strcmp(fields[0], "LET"))
			{
				if (!fields[1]) raise_error("No parameters given for LET.");
				n = find_var_index(fields[1]);
				if (n<0)
				{
					n = vars++;
					script_var[n].name = fields[1];
					script_var[n].vt = type_from_name(fields[1]);
				}
				int flags = n & _VARFLAGS_MASK;
				n &= _VARNUM_MASK;
				
				if (!fields[2]) raise_error("No operator given for LET.");
				if (!fields[3] && strcmp(fields[2], "++") && strcmp(fields[2], "--") ) raise_error("No rvalue given for LET.");
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
						raise_error( (std::string)"Unimplemented operator " + (std::string)fields[2] + (std::string)" for int assignment.");
						return 0x51974c5;
					}
					l=0;
					break;
					
					case SV_FLOAT:
					if (!strcmp(fields[2], "=")) script_var[n].value.f = interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "+=")) script_var[n].value.f += interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "-=")) script_var[n].value.f -= interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "*=")) script_var[n].value.f *= interpret_single_float(fields[3]);
					else if (!strcmp(fields[2], "/=")) script_var[n].value.f /= interpret_single_float(fields[3]);
					else
					{
						raise_error( (std::string)"Unimplemented operator " + (std::string)fields[2] + (std::string)" for float assignment.");
						return 0x51974c5;
					}
					l=0;
					break;
					
					case SV_POINT:
					if (flags & _HAS_DOT)
					{
						float* ff = nullptr;
						
						char* param = strchr(fields[1], '.');
						if (!param) raise_error( (std::string)"Missing member after dot.");
						param++;
						if (!strcmp(param, "x")) ff = &(script_var[n].value.ppt->x);
						else if (!strcmp(param, "X")) ff = &(script_var[n].value.ppt->x);
						else if (!strcmp(param, "y")) ff = &(script_var[n].value.ppt->y);
						else if (!strcmp(param, "Y")) ff = &(script_var[n].value.ppt->y);
						else if (!strcmp(param, "z")) ff = &(script_var[n].value.ppt->z);
						else if (!strcmp(param, "Z")) ff = &(script_var[n].value.ppt->z);
						else
						{
							raise_error((std::string)"Cartesian has no member named " + (std::string)param);
							return 0xbadd07;
						}
						
						if (ff)
						{
							if (!strcmp(fields[2], "=")) *ff = interpret_single_float(fields[3]);
							else if (!strcmp(fields[2], "+=")) *ff += interpret_single_float(fields[3]);
							else if (!strcmp(fields[2], "-=")) *ff -= interpret_single_float(fields[3]);
							else if (!strcmp(fields[2], "*=")) *ff *= interpret_single_float(fields[3]);
							else if (!strcmp(fields[2], "/=")) *ff /= interpret_single_float(fields[3]);
							else
							{
								raise_error( (std::string) "Unimplemented operator " + (std::string)fields[2] + (std::string)" for float assignment.");
								return 0x51974c5;
							}
							l = 0;
							n = -1;
						}
						break;
					}
					else
					{
						if (!strcmp(fields[2], "="))
						{
							script_var[n].value.ppt = new Point();
							*(script_var[n].value.ppt) = interpret_single_point(fields[3]);
						}
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
							raise_error( (std::string)"Unimplemented operator " + (std::string)fields[2] + (std::string)" for Cartesian assignment.");
							return 0x51974c5;
						}
					}
					
					l=0;
					
					break;
					
					case SV_STRING:
					
					psz = interpret_single_string(fields[3]);
					
					l = m = 0;
					if (fields[4] && !strcmp(fields[4], "FROM"))
					{
						l+=2;
						if (!fields[5]) raise_error("Insufficient parameters given for LET.");
						m = interpret_single_int(fields[5]);
						if (m < 0) m = 0;
						if (m)
						{
							if (m > strlen(psz)) psz[0] = 0;
							else
							{
								m--;
								psz += m;
								if (fields[6] && !strcmp(fields[6], "FOR"))
								{
									l+=2;
									if (!fields[7]) raise_error("Insufficient parameters given for LET.");
									k = interpret_single_int(fields[7]);
									if (k >= 0 && k < strlen(psz)) psz[k] = 0;
								}
							}
						}
					}
					
					if (!strcmp(fields[2], "="))
					{
						script_var[n].value.psz = new char[65536];
						strcpy(script_var[n].value.psz, psz);
					}
					else if (!strcmp(fields[2], "+="))
					{
						builder = script_var[n].value.psz;
						builder.append(psz);
						script_var[n].value.psz = new char[65536];
						strcpy(script_var[n].value.psz, builder.c_str());
					}
					else
					{
						psz -= m;
						delete[] psz;
						raise_error( (std::string)"Unimplemented operator " + (std::string)fields[2] + (std::string)" for string assignment.");
						return 0x51974c5;		// If you use your imagination, that says "syntax".
					}
					psz -= m;
					delete[] psz;
					break;
					
					default:
					;
				}	// switch (script_var[n].vt)
				
				while (n >= 0 && fields[3+l] && fields[4+l] && fields[5+l])
				{
					if (!fields[5+l]) raise_error("Insufficient parameters given for LET.");
					switch (script_var[n].vt)
					{
						case SV_INT:
						m = interpret_single_int(fields[5+l]);
						if (!strcmp(fields[4+l], "+")) script_var[n].value.n += m;
						else if (!strcmp(fields[4+l], "-")) script_var[n].value.n -= m;
						else if (!strcmp(fields[4+l], "*")) script_var[n].value.n *= m;
						else if (!strcmp(fields[4+l], "/")) script_var[n].value.n /= m;
						else if (!strcmp(fields[4+l], "^")) script_var[n].value.n = pow(script_var[n].value.n, m);
						else
						{
							raise_error( (std::string)"Bad operator " + (std::string)fields[4+l] + (std::string)" for int.");
							return 0x51974c5;
						}
						break;
						
						case SV_FLOAT:
						f = interpret_single_float(fields[5+l]);
						if (!strcmp(fields[4+l], "+")) script_var[n].value.f += f;
						else if (!strcmp(fields[4+l], "-")) script_var[n].value.f -= f;
						else if (!strcmp(fields[4+l], "*")) script_var[n].value.f *= f;
						else if (!strcmp(fields[4+l], "/")) script_var[n].value.f /= f;
						else if (!strcmp(fields[4+l], "^")) script_var[n].value.f = pow(script_var[n].value.f, f);
						else
						{
							// cout << "Bad operator " << fields[4+l] << " for float." << endl;
							raise_error( (std::string)"Bad operator " + (std::string)fields[4+l] + (std::string)" for float.");
							return 0x51974c5;
						}
						break;
						
						case SV_POINT:
						pt = interpret_single_point(fields[5+l]);
						if (!strcmp(fields[4+l], "+")) *script_var[n].value.ppt = script_var[n].value.ppt->add(pt);
						else if (!strcmp(fields[4+l], "-")) *script_var[n].value.ppt = script_var[n].value.ppt->subtract(pt);
						else if (!strcmp(fields[4+l], "*")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() * pt.magnitude());
						else if (!strcmp(fields[4+l], "/")) script_var[n].value.ppt->scale(pt.magnitude());
						else
						{
							// cout << "Bad operator " << fields[4+l] << " for point." << endl;
							raise_error( (std::string)"Bad operator " + (std::string)fields[4+l] + (std::string)" for point.");
							return 0x51974c5;
						}
						break;
						
						case SV_STRING:
						builder = script_var[n].value.psz;
						psz = interpret_single_string(fields[5+l]);
						if (!strcmp(fields[4+l], "+"))
						{
							builder.append(psz);
							script_var[n].value.psz = new char[65536];
							strcpy(script_var[n].value.psz, builder.c_str());
						}
						else
						{
							// cout << "Bad operator " << fields[4+l] << " for string." << endl;
							raise_error( (std::string)"Bad operator " + (std::string)fields[4+l] + (std::string)" for string.");
							return 0x51974c5;
						}
						break;
						
						default:
						;
					}
					
					l += 2;
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
				
				cout << endl << flush;
				_no_newline_on_echo:
				;
			}	// ECHO
			
			else if (!strcmp(fields[0], "DUMP"))
			{
				for (j=0; j<vars; j++)
				{
					cout << j << ": " << script_var[j].name << " ";
					switch (script_var[j].vt)
					{
						case SV_INT:	cout << "int " << script_var[j].value.n << endl; break;
						case SV_FLOAT:	cout << "float " << script_var[j].value.f << endl; break;
						case SV_POINT:	cout << "point " << *script_var[j].value.ppt << endl; break;
						case SV_STRING:	cout << "string " << script_var[j].value.psz << endl; break;
						default:		cout << "??? " << hex << script_var[j].value.n << dec << endl; break;
					}
				}
			}	// DUMP
			
			else if (!strcmp(fields[0], "GOTO"))
			{
				if (!fields[1]) raise_error("Insufficient parameters given for GOTO.");
				sprintf(buffer1, "%s:", fields[1]);
				for (n=0; n<script_lines.size(); n++)
				{
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
				raise_error( (std::string)"Label not found: \"" + (std::string)buffer1 + (std::string)"\"");
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
				raise_error( (std::string)"Unimplemented command: \"" + (std::string)fields[0] + (std::string)"\"");
				return 0x51974c5;
			}
		}
		
		_pc_continue:
		program_counter++;
	}
	
	return 0;
}










