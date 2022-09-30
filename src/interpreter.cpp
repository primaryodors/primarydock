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
    case '%':
        return SV_INT;
    case '&':
        return SV_FLOAT;
    case '@':
        return SV_POINT;
    case '$':
        return SV_STRING;
    default:
        raise_error("Variable names may only start with %, &, @, or $.");
    }
    return SV_NONE;
}

int interpret_single_int(const char* param);

float interpret_single_float(const char* param);

char* interpret_single_string(const char* param);

int find_var_index(const char* varname, char** out_varname = nullptr)
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

    char* c1 = strchr(buffer+1,'%');
    char* c2 = strchr(buffer+1,'$');

    if (!c1 && !c2) c = c1;
    if (!c1 &&  c2) c = c2;
    if ( c1 && !c2) c = c1;
    if ( c1 &&  c2) c = min(c1, c2);

    if (c > buffer)
    {
        char buffer1[256];
        strcpy(buffer1, c);

        while (c[1])
        {
            Star s;
            if (*c == '%') s.n = interpret_single_int(buffer1);
            if (*c == '$') s.psz = interpret_single_string(buffer1);

            if (s.n)
            {
                char buffer2[512];
                if (*c == '%')
                {
                    *c = 0;
                    sprintf(buffer2, "%s%lld%s", buffer, s.n, &c[strlen(buffer1)]);
                }
                if (*c == '$')
                {
                    *c = 0;
                    sprintf(buffer2, "%s%s%s", buffer, s.psz, &c[strlen(buffer1)]);
                }
                strcpy(buffer, buffer2);

                if (out_varname)
                {
                    // delete[] *out_varname;
                    *out_varname = new char[strlen(buffer)+4];
                    strcpy(*out_varname, buffer);
                }

                break;
            }

            buffer1[strlen(buffer1)-1] = 0;
        }
    }

    if (c=strchr(buffer,',')) *c=0;
    else if (buffer[0] != '%') if (c=strchr(buffer,'.'))
        {
            flags |= _HAS_DOT;
            *c=0;
        }
        else if (c=strchr(buffer,':')) *c=0;

    for (i=0; i<vars; i++)
    {
        if (!strcmp(script_var[i].name.c_str(), buffer)) return i|flags;
    }

    return -1;
}

int set_variable(char* vname, Star vvalue)
{
    int n = find_var_index(vname);
    if (n<0)
    {
        n = vars;
        vars++;
    }

    script_var[n].name = vname;
    script_var[n].vt = type_from_name(vname);
    script_var[n].value = vvalue;

    return n;
}

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
        sprintf(buffer, "%lld", script_var[n].value.n);
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
        strcpy(buffer, param);
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
        /*cout << "Usage:" << endl << "peptiditor protein.pdb script_filename" << endl;
        cout << "peptiditor script_filename protein.pdb" << endl;*/
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
        char** ofields = fields;
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

            if (!fields[0]) goto _pc_continue;

            if (fields[0][strlen(fields[0])-1] == ':') goto _pc_continue;

        _interpret_command:

            // Interpret the script.
            if (!strcmp(fields[0], "HELIX"))
            {
                float phi, psi;
                l = 2;

                if (!fields[1]) raise_error("No parameters given for HELIX.");
                if (!strcmp(fields[1], "ALPHA"))
                {
                    phi = ALPHA_PHI;
                    psi = ALPHA_PSI;
                }
                else if (!strcmp(fields[1], "PI"))
                {
                    phi = PI_PHI;
                    psi = PI_PSI;
                }
                else if (!strcmp(fields[1], "3.10"))
                {
                    phi = _310_PHI;
                    psi = _310_PSI;
                }
                else if (!strcmp(fields[1], "PPRO1"))
                {
                    phi = POLYPRO1_PHI;
                    psi = POLYPRO1_PSI;
                }
                else if (!strcmp(fields[1], "PPRO2"))
                {
                    phi = POLYPRO2_PHI;
                    psi = POLYPRO2_PSI;
                }

                // Not technically helices, but included for consistency.
                else if (!strcmp(fields[1], "BETA"))
                {
                    phi = BETA_PHI;
                    psi = BETA_PSI;
                }
                else if (!strcmp(fields[1], "STRAIGHT"))
                {
                    phi = M_PI;
                    psi = M_PI;
                }

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

            else if (!strcmp(fields[0], "REGION"))
            {
                int sr, er;
                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for REGION.");
                psz = interpret_single_string(fields[l++]);
                if (!psz[0]) psz = fields[l-1];
                if (!fields[l]) raise_error("Insufficient parameters given for REGION.");
                sr = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for REGION.");
                er = interpret_single_int(fields[l++]);
                if (fields[l]) raise_error("Too many parameters given for REGION.");

                p.set_region(psz, sr, er);

                char lbuffer[256];
                sprintf(lbuffer, "%s%s.s", "%", psz);
                Star s;
                s.n = sr;
                set_variable(lbuffer, s);
                sprintf(lbuffer, "%s%s.e", "%", psz);
                s.n = er;
                set_variable(lbuffer, s);
            }	// REGION

            else if (!strcmp(fields[0], "CENTER"))
            {
                l = 1;
                Point newcen(0,0,0);
                if (fields[l]) newcen = interpret_single_point(fields[l++]);
                if (fields[l]) raise_error("Too many parameters given for CENTER.");
                p.move_piece(1, 9999, newcen);
            }	// CENTER

            else if (!strcmp(fields[0], "UPRIGHT"))
            {
                l = 1;
                if (fields[l]) raise_error("Too many parameters given for UPRIGHT.");

                int seql = p.get_seq_length();
                p.move_piece(1, 9999, Point(0,0,0));

                Point extracellular[256], cytoplasmic[256];
                int exr_n=0, cyt_n=0;

                for (i=1; i<=7; i++)
                {
                    int sr = p.get_region_start((std::string)"TMR" + std::to_string(i));
                    if (!sr) continue;
                    int er = p.get_region_end((std::string)"TMR" + std::to_string(i));

                    for (j=0; j<4; j++)
                    {
                        if (i % 1)			// TMR1, TMR3, TMR5, TMR7 begin on the extracellular side and descend.
                        {
                            extracellular[exr_n++] = p.get_atom_location(sr+j, "CA");
                            cytoplasmic[cyt_n++] = p.get_atom_location(er-j, "CA");
                        }
                        else				// TMR2, TMR4, TMR6 ascend from the cytoplasmic side.
                        {
                            cytoplasmic[cyt_n++] = p.get_atom_location(sr+j, "CA");
                            extracellular[exr_n++] = p.get_atom_location(er-j, "CA");
                        }
                    }
                }

                if (!exr_n || !cyt_n) raise_error("Cannot UPRIGHT protein without transmembrane regions named TMR{n}.");

                Point exrdir = average_of_points(extracellular, exr_n);
                Point cytdir = average_of_points(cytoplasmic, cyt_n);

                Rotation rot = align_points_3d(&exrdir, new Point(0,1e9,0), &cytdir);

                p.rotate_piece(1, 9999, rot, 0);

                // TODO: Rotate to place TMR4 +Z to TMR1.
                int sr = p.get_region_start("TMR4");
                if (sr)
                {
                    int er = p.get_region_end("TMR4");

                    Point tmr1[64], tmr4[64];
                    int tmr1_n=0, tmr4_n=0;

                    for (i=sr; i<=er; i++)
                    {
                        tmr4[tmr4_n++] = p.get_atom_location(i, "CA");
                    }

                    sr = p.get_region_start("TMR1");
                    er = p.get_region_end("TMR1");

                    for (i=sr; i<=er; i++)
                    {
                        tmr1[tmr1_n++] = p.get_atom_location(i, "CA");
                    }

                    Point tmr1dir = average_of_points(tmr1, tmr1_n);
                    Point tmr4dir = average_of_points(tmr4, tmr4_n);


                    tmr1dir.y = tmr4dir.y = 0;

                    rot = align_points_3d(&tmr4dir, new Point(0,0,1e9), &tmr1dir);

                    p.rotate_piece(1, 9999, rot, 0);

                }
            }	// UPRIGHT

            else if (!strcmp(fields[0], "RENUMBER"))
            {
				l = 1;
                int sr, er, nsr;
                if (!fields[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)fields[0] + (std::string)".");
                sr = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)fields[0] + (std::string)".");
                er = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)fields[0] + (std::string)".");
                nsr = interpret_single_int(fields[l++]);
                if (fields[l]) raise_error((std::string)"Too many parameters given for " + (std::string)fields[0] + (std::string)".");

				p.renumber_residues(sr, er, nsr);
            }	// RENUMBER

            else if (!strcmp(fields[0], "SEARCH"))
            {
                l = 1;
                int sr, er, esr, sim;
                if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
                sr = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
                er = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
                psz = interpret_single_string(fields[l++]);
                esr = er - strlen(psz);

                int threshold = -1;
                int num_eq;

                if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
                // if (fields[l+1]) threshold = interpret_single_int(fields[l++]);
                if (!strcmp(fields[l], "TH"))
                {
                    l++;
                    threshold = interpret_single_int(fields[l++]);
                }

                n = 0;
                k = 0;
                for (i=sr; i<esr; i++)
                {
                    m = num_eq = 0;
                    for (j=0; psz[j]; j++)
                    {
                        char c = psz[j], aac = p.get_residue(i+j)->get_letter();
                        if (c == 'X') c = aac;

                        if (c == aac) num_eq++;

                        sim = p.get_residue(i+j)->similarity_to(c);
                        // cout << c << "/" << aac << " " << sim << "  ";

                        m += sim;
                    }
                    // cout << "___ m: " << m << ", n: " << n << endl;

                    if (m > n && num_eq >= threshold)
                    {
                        k = i;
                        n = m;
                    }
                }
                sim = n;

                delete[] psz;

                n = find_var_index(fields[l]);
                if (n<0) n = vars++;
                n &= _VARNUM_MASK;
                if (!fields[l]) raise_error("Insufficient parameters given for SEARCH.");
                if (fields[l+1] && fields[l+2]) raise_error("Too many parameters given for SEARCH.");
                script_var[n].name = fields[l];
                script_var[n].vt = type_from_name(fields[l]);

                switch (script_var[n].vt)
                {
                case SV_INT:
                    script_var[n].value.n = k;
                    break;
                case SV_FLOAT:
                    script_var[n].value.f = k;
                    break;
                case SV_STRING:
                    builder = std::to_string(k);
                    script_var[n].value.psz = new char[builder.length()+2];
                    strcpy(script_var[n].value.psz, builder.c_str());
                    break;

                default:
                    raise_error("Bad destination variable type for SEARCH.");
                    return 0xbadfa12;
                }

                l++;
                if (fields[l])
                {
                    n = find_var_index(fields[l]);
                    if (n<0) n = vars++;
                    n &= _VARNUM_MASK;
                    script_var[n].name = fields[l];
                    script_var[n].vt = type_from_name(fields[l]);

                    switch (script_var[n].vt)
                    {
                    case SV_INT:
                        script_var[n].value.n = sim;
                        break;
                    case SV_FLOAT:
                        script_var[n].value.f = sim;
                        break;
                    case SV_STRING:
                        builder = std::to_string(sim);
                        script_var[n].value.psz = new char[builder.length()+2];
                        strcpy(script_var[n].value.psz, builder.c_str());
                        break;

                    default:
                        raise_error("Bad destination variable type for SEARCH.");
                        return 0xbadfa12;
                    }
                }

            }	// SEARCH

            else if (!strcmp(fields[0], "BEND"))
            {
                int sr, er;
                bb_rot_dir bbrotd;
                float theta;

                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for BEND.");
                sr = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for BEND.");
                er = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for BEND.");
                char* tmp = interpret_single_string(fields[l++]);
                if (!strcmp(tmp, "N-CA")) bbrotd = N_asc;
                else if (!strcmp(tmp, "CA-C")) bbrotd = CA_asc;
                else if (!strcmp(tmp, "CA-N")) bbrotd = CA_desc;
                else if (!strcmp(tmp, "C-CA")) bbrotd = C_desc;
                else raise_error("Unknown direction parameter given for BEND.");

                if ((bbrotd == N_asc || bbrotd == CA_asc) && er < sr) raise_error("Cannot rotate ascending bond in the descending direction.");
                if ((bbrotd == CA_desc || bbrotd == C_desc) && er > sr) raise_error("Cannot rotate descending bond in the ascending direction.");

                if (!fields[l]) raise_error("Insufficient parameters given for BEND.");
                theta = interpret_single_float(fields[l++]) * fiftyseventh;
                if (fields[l]) raise_error("Too many parameters given for BEND.");

                p.rotate_backbone_partial(sr, er, bbrotd, theta);

            }	// BEND

            else if (!strcmp(fields[0], "BRIDGE"))
            {
                int r1, r2, iters=50;

                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for BRIDGE.");
                r1 = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for BRIDGE.");
                r2 = interpret_single_int(fields[l++]);
                if (fields[l]) iters = interpret_single_int(fields[l++]);
                if (fields[l]) raise_error("Too many parameters given for BRIDGE.");

                Star a1, a2;
                a1.paa = p.get_residue(r1);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r1) + (std::string)" not found in protein.");
                a2.paa = p.get_residue(r2);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r2) + (std::string)" not found in protein.");

                a1.pmol->movability = MOV_FLEXONLY;
                a2.pmol->movability = MOV_FLEXONLY;

                Molecule* mm[5];
                for (i=0; i<5; i++) mm[i] = nullptr;
                mm[0] = a1.pmol;
                mm[1] = a2.pmol;

                Molecule::multimol_conform(mm, iters);

            }	// BRIDGE

            else if (!strcmp(fields[0], "BENERG"))
            {
                // Read non-covalent binding strength between two side chains into a float var.
                int r1, r2;

                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for BENERG.");
                r1 = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for BENERG.");
                r2 = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for BENERG.");
                n = find_var_index(fields[l]);
                if (n<0)
                {
                    n = vars++;
                    script_var[n].name = fields[l];
                    script_var[n].vt = type_from_name(fields[l]);
                }
                int flags = n & _VARFLAGS_MASK;
                n &= _VARNUM_MASK;
                l++;
                if (fields[l]) raise_error("Too many parameters given for BENERG.");

                if (script_var[n].vt != SV_FLOAT) raise_error("Wrong variable type given for BENERG; float required.");

                Star a1, a2;
                a1.paa = p.get_residue(r1);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r1) + (std::string)" not found in protein.");
                a2.paa = p.get_residue(r2);
                if (!a2.n) raise_error((std::string)"Residue " + to_string(r2) + (std::string)" not found in protein.");

                script_var[n].value.f = -(a1.pmol->get_intermol_binding(a2.pmol));
            }

            else if (!strcmp(fields[0], "PTALIGN"))
            {
                Point point, align, center;
                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for PTALIGN.");
                point = interpret_single_point(fields[l++]);

                if (!fields[l]) raise_error("Insufficient parameters given for PTALIGN.");
                align = interpret_single_point(fields[l++]);

                if (!fields[l]) raise_error("Insufficient parameters given for PTALIGN.");
                center = interpret_single_point(fields[l++]);
                
                Rotation rot = align_points_3d(&point, &align, &center);

                if (!fields[l]) raise_error("Insufficient parameters given for PTALIGN.");
                Star s;
                rot.v.r = 1;
                Point p = rot.v;
                s.ppt = &p;
                set_variable(fields[l++], s);

                if (!fields[l]) raise_error("Insufficient parameters given for PTALIGN.");
                s.f = rot.a * fiftyseven;
                set_variable(fields[l++], s);

                if (fields[l]) raise_error("Too many parameters given for PTALIGN.");
            }

            else if (!strcmp(fields[0], "PTROTATE"))
            {
                Point point, origin, axis;
                float theta;
                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for PTROTATE.");
                point = interpret_single_point(fields[l++]);

                if (!fields[l]) raise_error("Insufficient parameters given for PTROTATE.");
                origin = interpret_single_point(fields[l++]);

                if (!fields[l]) raise_error("Insufficient parameters given for PTROTATE.");
                axis = interpret_single_point(fields[l++]);

                if (!fields[l]) raise_error("Insufficient parameters given for PTROTATE.");
                theta = interpret_single_float(fields[l++]) * fiftyseventh;

                Point result = rotate3D(point, origin, axis, theta);

                if (!fields[l]) raise_error("Insufficient parameters given for PTROTATE.");
                Star s;
                s.ppt = &result;
                set_variable(fields[l++], s);

                if (fields[l]) raise_error("Too many parameters given for PTROTATE.");
            }

            else if (!strcmp(fields[0], "ALIGN"))
            {
                int sr, er, asr, aer, eachend;
                Point sp, ep;
                l = 1;
                if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                sr = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                er = interpret_single_int(fields[l++]);
                if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                eachend = interpret_single_int(fields[l++]);
                if (!eachend) eachend = 1;
                if (fields[l+2] && fields[l+3])
                {                   
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    asr = interpret_single_int(fields[l++]);                   
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    sp = interpret_single_point(fields[l++]);
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    aer = interpret_single_int(fields[l++]);   
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    ep = interpret_single_point(fields[l++]);
                }
                else
                {
                    asr = sr;
                    aer = er;                    
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    sp = interpret_single_point(fields[l++]);
                    if (!fields[l]) raise_error("Insufficient parameters given for ALIGN.");
                    ep = interpret_single_point(fields[l++]);
                }
                if (fields[l]) raise_error("Too many parameters given for ALIGN.");

                // From start and end residues inwards for a total of eachend, average the CA locations.
                Point sl(0,0,0);
                for (i=0; i<eachend; i++)	// It's actually easier to average manually than to screw around with object arrays.
                {
                    AminoAcid* aa = p.get_residue(asr+i);
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
                    AminoAcid* aa = p.get_residue(aer-i);
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

				Star sv;

                int seqlen = p.get_seq_length();
                strcpy(buffer, "%SEQLEN");
                sv.n = seqlen;
				set_variable(buffer, sv);

                strcpy(buffer, "$SEQUENCE");
                sv.psz = new char[seqlen+4];
                strcpy(sv.psz, p.get_sequence().c_str());
				set_variable(buffer, sv);
            } // DELETE

            else if (!strcmp(fields[0], "MCOORD"))
            {
                l = 1;
                Atom* ma;
                string elem_sym;
                int elem_charge=0;
                int ncr=0;						// number of coordinating residues.
                int resnos[13];					// more than the task will ever conceivably require.
                std::vector<string> cratoms;	// coordinating residue atoms.																						oh god I hooked up with this guy once who used kratom and after we both finished goddamn all he did was yipyapyipyapyipyap all night until finally I said babe, I n**d my sleep, and I musta finally dozed off at something like 5am, yeah kinda like the song. Then it's morning and he goes home and a hot drummer/songwriter is coming over to audition for my rock band that day and I'm running on 3 hours of sleep, what a way to make a first impression, thank flying spaghetti monster our guitarist was my roommate because at least someone in the house was awake. Ahhh, my 30s, when I still had (false) hope. Good times, those.
                bool force_tyrosine_O = false;
                bool thiolate = false;
                Atom** ba = nullptr;

            _yes_I_used_goto_for_this:
                if (!strcmp(fields[l], "YO"))
                {
                    force_tyrosine_O = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(fields[l], "YAr"))
                {
                    force_tyrosine_O = false;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(fields[l], "Th8"))
                {
                    thiolate = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }

                if (!fields[l]) raise_error("Insufficient parameters given for MCOORD.");
                elem_sym = fields[l++];
                if (!fields[l]) raise_error("Insufficient parameters given for MCOORD.");
                elem_charge = interpret_single_int(fields[l++]);
                Point pt;
                ma = new Atom(elem_sym.c_str(), &pt, elem_charge);
                string mname = elem_sym.append("1");
                ma->name = new char[8];
                strcpy(ma->name, mname.c_str());
                strcpy(ma->aa3let, "MTL");
                ma->residue = 0;

                if (!fields[l]) raise_error("Insufficient parameters given for MCOORD.");
                for (; fields[l]; l++)
                {
                    bool local_O = force_tyrosine_O;

                _another_goto:
                    if (!strcmp(fields[l], "YO"))
                    {
                        local_O = true;
                        l++;
                        goto _another_goto;
                    }
                    else if (!strcmp(fields[l], "YAr"))
                    {
                        local_O = false;
                        l++;
                        goto _another_goto;
                    }

                    if (!fields[l]) raise_error("Insufficient parameters given for MCOORD.");
                    k = interpret_single_int(fields[l]);
                    AminoAcid* aa = p.get_residue(k);
                    resnos[ncr] = k;

                    // cout << aa->get_3letter() << k << " is " << (aa->is_tyrosine_like() ? "" : "not ") << "tyrosine-like." << endl;
                    if (aa->is_tyrosine_like() && !local_O)
                    {
                        Ring* rr = aa->get_most_distal_arom_ring();
                        if (rr)
                        {
                            n = rr->get_atom_count();
                            // Get members 1 and 1+floor(n/2).
                            cratoms.push_back(rr->get_atom(1)->name);
                            cout << "Found metal coord atom " << *aa << ":" << cratoms[ncr] << endl;
                            ncr++;

                            resnos[ncr] = k;
                            j = 1 + (n/2);
                            cratoms.push_back(rr->get_atom(j)->name);
                            cout << "Found metal coord atom " << *aa << ":" << cratoms[ncr] << endl;
                            ncr++;

                            goto _found_coord_atom;
                        }
                    }

                    ba = aa->get_most_bindable(1, ma);

                    if (ba && ba[0])
                    {
                        cratoms.push_back(ba[0]->name);
                        cout << "Found metal coord atom " << *aa << ":" << cratoms[ncr] << endl;
                        ncr++;
                    }
                    else raise_error((std::string)"No metal coordination atom found for " + (std::string)aa->get_3letter() + to_string(k));

                _found_coord_atom:
                    ;

                }

                if (ncr < 3) raise_error("MCOORD requires at least 3 coordinating atoms.");
                p.coordinate_metal(ma, ncr, resnos, cratoms);

            } // MCOORD

            else if (!strcmp(fields[0], "LOAD"))
            {
                if (!fields[1]) raise_error("Insufficient parameters given for LOAD.");
                psz = interpret_single_string(fields[1]);
				n = 0;
				// if (fields[2]) n = atoi(fields[2]);
                if (fields[2] /*&& fields[3]*/) raise_error("Too many parameters given for LOAD.");

                pf = fopen(psz, "rb");
                if (!pf)
                {
                    raise_error( (std::string)"Failed to open " + (std::string)psz + (std::string)" for reading.");
                    return 0xbadf12e;
                }
                p.load_pdb(pf, n);

                fclose(pf);

            _prot_deets:

				char buffer[1024];
				char buffer1[1024];
				Star sv;

				strcpy(buffer, "$PDB");
				sv.psz = new char[strlen(psz)+4];
                strcpy(sv.psz, psz);
				set_variable(buffer, sv);

                const char* pname = p.get_name().c_str();
                strcpy(buffer, "$PROTEIN");
                sv.psz = new char[strlen(pname)+4];
                strcpy(sv.psz, pname);
				set_variable(buffer, sv);

                delete[] psz;

                int seqlen = p.get_seq_length();
                strcpy(buffer, "%SEQLEN");
                sv.n = seqlen;
				set_variable(buffer, sv);

                strcpy(buffer, "$SEQUENCE");
                sv.psz = new char[seqlen+4];
                strcpy(sv.psz, p.get_sequence().c_str());
				set_variable(buffer, sv);

                std::vector<std::string> rem_hx = p.get_remarks("650 HELIX");
                for (l=0; l<rem_hx.size(); l++)
                {
                    strcpy(buffer, rem_hx[l].c_str());
                    char** fields = chop_spaced_fields(buffer);

					if (!fields[3] || !fields[4] || !fields[5]) continue;

                    sprintf(buffer1, "%c%s.s", '%', fields[3]);
                    sv.n = atoi(fields[4]);
                    set_variable(buffer1, sv);

                    sprintf(buffer1, "%c%s.e", '%', fields[3]);
                    sv.n = atoi(fields[5]);
                    set_variable(buffer1, sv);

					p.set_region(fields[3], atoi(fields[4]), atoi(fields[5]));

                    delete[] fields;
                }
            }

            else if (!strcmp(fields[0], "SAVE"))
            {
                if (!fields[1]) raise_error("Insufficient parameters given for SAVE.");
                psz = interpret_single_string(fields[1]);
                if (fields[2] && fields[3]) raise_error("Too many parameters given for SAVE.");

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
                    else raise_error((std::string)"Too many parameters given for SAVE: " + (std::string)fields[2]);
                }
            }

            else if (!strcmp(fields[0], "LET"))
            {
                if (!fields[1]) raise_error("No parameters given for LET.");
                n = find_var_index(fields[1], &fields[1]);
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
                    else if (!strcmp(fields[2], "&=")) script_var[n].value.n &= interpret_single_int(fields[3]);
                    else if (!strcmp(fields[2], "|=")) script_var[n].value.n |= interpret_single_int(fields[3]);
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
                        else if (!strcmp(fields[4+l], "&")) script_var[n].value.n &= m;
                        else if (!strcmp(fields[4+l], "|")) script_var[n].value.n |= m;
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
                        else if (!strcmp(fields[4+l], "/")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() / pt.magnitude());
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
                        else if (!strcmp(fields[4+l], "FOR"))
                        {
                            l += 2;
                            continue;
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
                    case SV_INT:
                        cout << "int " << script_var[j].value.n << endl;
                        break;
                    case SV_FLOAT:
                        cout << "float " << script_var[j].value.f << endl;
                        break;
                    case SV_POINT:
                        cout << "point " << *script_var[j].value.ppt << endl;
                        break;
                    case SV_STRING:
                        cout << "string " << script_var[j].value.psz << endl;
                        break;
                    default:
                        cout << "??? " << hex << script_var[j].value.n << dec << endl;
                        break;
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

            else if (!strcmp(fields[0], "IF"))
            {
                if (!fields[1]) raise_error("Insufficient parameters given for IF.");
                if (!fields[2]) raise_error("Insufficient parameters given for IF.");

                // TODO: NOT operator and IF %var GOTO blah.
                l = 2;

                // If the operator is =, and both l-value and r-value are strings, do a direct comparison.
                if (!fields[l]) raise_error("Insufficient parameters given for IF.");
                if (!strcmp(fields[l], "THEN"))
                {
                    l--;
                    if (interpret_single_float(fields[l])) goto _evaluated_true;
                    else goto _evaluated_false;
                }
                if (!strcmp(fields[l], "="))
                {
                    if (!fields[l+1]) raise_error("Insufficient parameters given for IF.");
                    if (fields[l-1][0] == fields[l+1][0]
                            ||
                            fields[l-1][0] == '$' && fields[l+1][0] == '"'
                            ||
                            fields[l-1][0] == '"' && fields[l+1][0] == '$'
                       )
                    {
                        char *lvalue = interpret_single_string(fields[l-1]),
                              *rvalue = interpret_single_string(fields[l+1]);
                        if (strcmp(lvalue, rvalue)) goto _evaluated_false;
                        else goto _evaluated_true;
                    }
                    else goto _just_interpret_floats;
                }
                else
                {
                    // Otherwise, interpret both values as floats.
                _just_interpret_floats:
                    float lvalue = interpret_single_float(fields[l-1]),
                          rvalue = interpret_single_float(fields[l+1]);

                    if (!strcmp(fields[l], "="))
                    {
                        if (lvalue == rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(fields[l], "!="))
                    {
                        if (lvalue != rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(fields[l], ">"))
                    {
                        if (lvalue > rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(fields[l], "<"))
                    {
                        if (lvalue < rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(fields[l], ">="))
                    {
                        if (lvalue >=rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(fields[l], "<="))
                    {
                        if (lvalue <= rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else raise_error( (std::string)"Unknown operator " + (std::string)fields[l] + (std::string)" for comparison.");
                }

            _evaluated_true:
                l += 2;
                if (!strcmp(fields[l], "THEN")) l++;
                fields = &fields[l];

                if (!fields[0]) goto _pc_continue;

                while (!strcmp(fields[0], "OR"))
                {
                    fields = &fields[4];
                    if (!fields[0]) goto _pc_continue;
                }

                if (!strcmp(fields[0], "AND")) strcpy(fields[0], "IF");
                goto _interpret_command;

            _evaluated_false:
                if (fields[l+2] && !strcmp(fields[l+2], "OR"))
                {
                    l += 2;
                    fields = &fields[l];
                    strcpy(fields[0], "IF");
                    goto _interpret_command;
                }

                program_counter++;
                strcpy(buffer, script_lines[program_counter].c_str());
                fields = chop_spaced_fields(buffer);
                if (!fields || !fields[0]) goto _pc_continue;
                if (strcmp(fields[0], "ELSE")) continue;

                fields = &fields[1];
                goto _interpret_command;
            }	// IF

            else if (!strcmp(fields[0], "ELSE")) goto _pc_continue;

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
        delete[] ofields;
        program_counter++;
    }

    return 0;
}










