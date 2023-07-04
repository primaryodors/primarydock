#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
#if _WIN32
#include <Windows.h>
#endif
#include "classes/protein.h"
#include "classes/group.h"

using namespace std;

#define _HAS_DOT 0x800000
#define _VARNUM_MASK 0xffff
#define _VARFLAGS_MASK 0xff0000

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
ScriptVar script_var[_VARNUM_MASK+1];
int vars = 0;
int program_counter = 0;

// Protein p("TheProtein");
Protein* strands[26];
Protein* working = nullptr;
char g_chain = 'A';

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

bool file_exists(std::string fname)
{
    struct stat s;
    if (stat(fname.c_str(), &s) == 0) return true;
    else return false;
}

bool download_file(std::string url, std::string destination)
{
    #if _WIN32
    URLDownloadToFile(NULL, url.c_str(), destination.c_str(), 0, NULL);

    #elif defined(__linux__) || defined(__sun) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__APPLE__)
    pid_t child_pid = fork();
    if (child_pid == -1)
    {
        perror("fork");
    }
    else if (child_pid == 0)
    {
        execlp("wget", "wget", "-O", destination.c_str(), url.c_str(), (char *)NULL);
    }

    int status = 0;
    waitpid(child_pid, &status, 0);
    if (status) return false;

    #else
    #error It appears your operating system is not supported yet - would you be willing to add it and submit a pull request?
    #endif

    return file_exists(destination);
}

float contact_energy(Protein* a, Protein* b)
{
    std::vector<AminoAcid*> vca = a->get_contact_residues(b);
    float f=0;

    int i, n = vca.size();

    Molecule* mols[n+4];
    for (i=0; i<n; i++) mols[i] = (Molecule*)vca[i];
    mols[i] = nullptr;

    Molecule::conform_molecules(mols, 20);

    for (i=0; i<n; i++)
    {
        f -= mols[i]->get_intermol_binding(mols);
    }

    return f;
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
    char* c = nullptr;
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

int set_variable(const char* vname, const Star vvalue)
{
    int n = find_var_index(vname);
    if (n<0)
    {
        n = vars;
        vars++;
    }

    script_var[n].name = vname;
    script_var[n].vt = type_from_name(vname);
    if (script_var[n].vt == SV_STRING)
    {
        if (script_var[n].value.psz) delete script_var[n].value.psz;
        script_var[n].value.psz = new char[strlen(vvalue.psz)+4];
        strcpy(script_var[n].value.psz, vvalue.psz);
    }
    else script_var[n].value = vvalue;

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

Point interpret_single_point(const char* param, Point old_value = Point(0,0,0))
{
    int n;
    Point pt = old_value;
    AminoAcid* aa;

    if (param[0] >= '0' && param[0] <= '9')
    {
        aa = working->get_residue(atoi(param));
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
        aa = working->get_residue(script_var[n].value.n);
        if (aa) pt = aa->get_CA_location();
        return pt;

    case '&':
        n = find_var_index(param);
        if (n<0) return pt;
        n &= _VARNUM_MASK;
        if (pt.magnitude()) pt.scale(script_var[n].value.f);
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

    for (i=0; i<26; i++) strands[i] = nullptr;

    strands[0] = new Protein("TheProt");
    working = strands[0];
    g_chain = 'A';

    for (i=0; i<=_VARNUM_MASK; i++)
    {
        script_var[i].name = "";
        script_var[i].value.n = 0;
    }
    vars = 0;

    bool script_loaded = false;
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
        else if (!script_loaded)
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
                
                script_loaded = true;
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
        /*cout << "Usage:" << endl << "pepteditor protein.pdb script_filename" << endl;
        cout << "pepteditor script_filename protein.pdb" << endl;*/
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
        char** words = chop_spaced_words(buffer);
        char** owords = words;
        char chain = 'A';
        char new_chain;
        if (words && words[0] && words[0][0] && words[0][1])
        {
            for (k=0; words[k]; k++)
            {
                if (words[k][0] == '#')
                {
                    // cout << "Ignoring comment " << words[k] << endl << flush;
                    words[k] = nullptr;
                    break;
                }
            }
            // cout << "command " << words[0] << endl << flush;

            // Debug cout.
            /*cout << endl << "***Line: " << buffer << endl << "***words: ";
            for (n=0; words[n]; n++) cout << (n?"|":"") << words[n];
            cout << endl;*/

            if (!words[0]) goto _pc_continue;

            if (words[0][strlen(words[0])-1] == ':') goto _pc_continue;

        _interpret_command:

            // Interpret the script.
            if (0)
            {
                ;
            }
            else if (!strcmp(words[0], "ALIGN"))
            {
                int sr, er, asr, aer, eachend;
                Point sp, ep;
                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                er = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                eachend = interpret_single_int(words[l++]);
                if (!eachend) eachend = 1;
                if (words[l+2] && words[l+3])
                {                   
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    asr = interpret_single_int(words[l++]);                   
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    sp = interpret_single_point(words[l++]);
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    aer = interpret_single_int(words[l++]);   
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    ep = interpret_single_point(words[l++]);
                }
                else
                {
                    asr = sr;
                    aer = er;                    
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    sp = interpret_single_point(words[l++]);
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    ep = interpret_single_point(words[l++]);
                }
                if (words[l]) raise_error("Too many parameters given for ALIGN.");

                // From start and end residues inwards for a total of eachend, average the CA locations.
                Point sl(0,0,0);
                for (i=0; i<eachend; i++)	// It's actually easier to average manually than to screw around with object arrays.
                {
                    AminoAcid* aa = working->get_residue(asr+i);
                    if (aa) sl = sl.add(aa->get_CA_location());
                }

                if (eachend > 1) sl.scale(sl.magnitude()/eachend);

                // Translate the range so that the starting average moves to the target start point.
                SCoord motion = sp.subtract(sl);
                for (i=sr; i<=er; i++)
                {
                    AminoAcid* aa = working->get_residue(i);
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
                    AminoAcid* aa = working->get_residue(aer-i);
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
                    AminoAcid* aa = working->get_residue(i);
                    if (aa)
                    {
                        MovabilityType fmov = aa->movability;
                        aa->movability = MOV_ALL;
                        aa->rotate(lv, rot.a);
                        aa->movability = fmov;
                    }
                }

            }	// ALIGN

            else if (!strcmp(words[0], "BEND"))
            {
                int sr, er;
                bb_rot_dir bbrotd;
                float theta;

                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for BEND.");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for BEND.");
                er = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for BEND.");
                char* tmp = interpret_single_string(words[l++]);
                if (!strcmp(tmp, "N-CA")) bbrotd = N_asc;
                else if (!strcmp(tmp, "CA-C")) bbrotd = CA_asc;
                else if (!strcmp(tmp, "CA-N")) bbrotd = CA_desc;
                else if (!strcmp(tmp, "C-CA")) bbrotd = C_desc;
                else raise_error("Unknown direction parameter given for BEND.");

                if ((bbrotd == N_asc || bbrotd == CA_asc) && er < sr) raise_error("Cannot rotate ascending bond in the descending direction.");
                if ((bbrotd == CA_desc || bbrotd == C_desc) && er > sr) raise_error("Cannot rotate descending bond in the ascending direction.");

                if (!words[l]) raise_error("Insufficient parameters given for BEND.");
                theta = interpret_single_float(words[l++]) * fiftyseventh;
                if (words[l]) raise_error("Too many parameters given for BEND.");

                working->rotate_backbone_partial(sr, er, bbrotd, theta);

            }	// BEND

            else if (!strcmp(words[0], "BENERG"))
            {
                // Read non-covalent binding strength between two side chains into a float var.
                int r1, r2;

                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for BENERG.");
                r1 = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for BENERG.");
                r2 = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for BENERG.");
                n = find_var_index(words[l]);
                if (n<0)
                {
                    n = vars++;
                    script_var[n].name = words[l];
                    script_var[n].vt = type_from_name(words[l]);
                }
                int flags = n & _VARFLAGS_MASK;
                n &= _VARNUM_MASK;
                l++;
                if (words[l]) raise_error("Too many parameters given for BENERG.");

                if (script_var[n].vt != SV_FLOAT) raise_error("Wrong variable type given for BENERG; float required.");

                Star a1, a2;
                a1.paa = working->get_residue(r1);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r1) + (std::string)" not found in protein.");
                a2.paa = working->get_residue(r2);
                if (!a2.n) raise_error((std::string)"Residue " + to_string(r2) + (std::string)" not found in protein.");

                script_var[n].value.f = -(a1.pmol->get_intermol_binding(a2.pmol));
            }

            else if (!strcmp(words[0], "BRIDGE"))
            {
                int r1, r2, iters=50;

                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for BRIDGE.");
                r1 = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for BRIDGE.");
                r2 = interpret_single_int(words[l++]);
                if (words[l]) iters = interpret_single_int(words[l++]);
                if (words[l]) raise_error("Too many parameters given for BRIDGE.");

                Star a1, a2;
                a1.paa = working->get_residue(r1);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r1) + (std::string)" not found in protein.");
                a2.paa = working->get_residue(r2);
                if (!a1.n) raise_error((std::string)"Residue " + to_string(r2) + (std::string)" not found in protein.");

                a1.pmol->movability = MOV_FLEXONLY;
                a2.pmol->movability = MOV_FLEXONLY;

                Molecule* mm[5];
                for (i=0; i<5; i++) mm[i] = nullptr;
                mm[0] = a1.pmol;
                mm[1] = a2.pmol;

                Molecule::conform_molecules(mm, iters);

            }	// BRIDGE

            else if (!strcmp(words[0], "CENTER"))
            {
                l = 1;
                Point newcen(0,0,0);
                if (words[l]) newcen = interpret_single_point(words[l++]);
                if (words[l]) raise_error("Too many parameters given for CENTER.");
                working->move_piece(1, 9999, newcen);
            }	// CENTER

            else if (!strcmp(words[0], "CONNECT"))
            {
                l=1;
                int sr, er, ct, iters=250;
                if (!words[l]) raise_error("Insufficient parameters given for CONNECT.");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for CONNECT.");
                ct = interpret_single_int(words[l++]);
                if (words[l]) iters = interpret_single_int(words[l++]);
                if (words[l]) raise_error("Too many parameters given for CONNECT.");
                er = ct - sgn(ct-sr);

                Atom *a1, *a2, *a3;
                AminoAcid *cta = working->get_residue(ct), *era = working->get_residue(er);

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

                working->conform_backbone(sr, er, a1, pt3[0], a2, pt3[1], iters);
                delete[] pt3;
                working->backconnect(sr, er);

            _no_connect:
                goto _pc_continue;
            }	// CONNECT

            else if (!strcmp(words[0], "CTNRG"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for CTNRG.");
                char* sa = interpret_single_string(words[1]);
                if (!words[2]) raise_error("Insufficient parameters given for CTNRG.");
                char* sb = interpret_single_string(words[2]);
                if (!words[3]) raise_error("Insufficient parameters given for CTNRG.");
                SCoord dir(0,0,0);
                if (words[4])
                {
                    dir = interpret_single_point(words[4]);
                    if (words[5]) raise_error("Too many parameters given for CTNRG.");
                }

                sa[0] -= 65;
                sb[0] -= 65;

                if (dir.r)
                {
                    dir.r = 1;
                    float former, latter;

                    former = contact_energy(strands[sa[0]], strands[sb[0]]);
                    for (i=0; i<200; i++)
                    {
                        strands[sb[0]]->move_piece(1, strands[sb[0]]->get_end_resno(), dir);
                        latter = contact_energy(strands[sa[0]], strands[sb[0]]);

                        if (latter > former)
                        {
                            strands[sb[0]]->undo();
                            dir.r *= -0.75;
                        }
                        else
                        {
                            former = latter;
                            dir.r *= 1.1;
                        }

                        if (fabs(dir.r) < 0.05)
                        {
                            cout << "Got it in " << i << endl;
                            break;
                        }
                    }
                }

                Star s;
                s.f = contact_energy(strands[sa[0]], strands[sb[0]]);

                set_variable(words[3], s);
            }   // CTNRG

            else if (!strcmp(words[0], "DELETE"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for DELETE.");
                int sr = interpret_single_int(words[1]), er=0;
                if (words[2]) er = interpret_single_int(words[2]);
                if (words[3]) raise_error("Too many parameters given for DELETE.");

                if (er) working->delete_residues(sr, er);
                else working->delete_residue(sr);

				Star sv;

                int seqlen = working->get_seq_length();
                strcpy(buffer, "%SEQLEN");
                sv.n = seqlen;
				set_variable(buffer, sv);

                strcpy(buffer, "$SEQUENCE");
                sv.psz = new char[seqlen+4];
                strcpy(sv.psz, working->get_sequence().c_str());
				set_variable(buffer, sv);
            } // DELETE

            else if (!strcmp(words[0], "DISULF"))
            {
                if (words[1])
                {
                    if (!words[2]) raise_error("DISULF takes either zero or two parameters.");
                    if ( words[3]) raise_error("DISULF takes either zero or two parameters.");
                    j = interpret_single_int(words[1]);
                    k = interpret_single_int(words[2]);
                    bool result = working->disulfide_bond(j, k);
                    if (!result) cout << "WARNING: Failed to form disulfide bond between " << j << " and " << k << "." << endl;
                }
                else
                {
                    n = working->get_end_resno();
                    for (j=1; j<n; j++)
                    {
                        AminoAcid* resj = working->get_residue(j);
                        if (!resj->is_thiol()) continue;
                        for (k=j+1; k<=n; k++)
                        {
                            AminoAcid* resk = working->get_residue(k);
                            if (resk->is_thiol())
                            {
                                if (working->disulfide_bond(j, k)) cout << "Disulfide bonded " << *resj << "-" << *resk << "." << endl;
                            }
                        }
                    }
                }
            } // DISULF

            else if (!strcmp(words[0], "DOWNLOAD"))
            {
                if (!words[1]) raise_error("Insufficient parameters for DOWNLOAD.");
                if (!words[2]) raise_error("Insufficient parameters for DOWNLOAD.");

                std::string url = "", retfmt = "", destfn = "";
                bool skip_if_exists = false;

                if (words[3])
                {
                    if (!strcmp(words[3], "ONCE")) skip_if_exists = true;
                    else destfn = interpret_single_string(words[3]);

                    if (!strcmp(words[4], "ONCE")) skip_if_exists = true;

                    if (words[5]) raise_error("Too many parameters for DOWNLOAD.");
                }

                if (!skip_if_exists || !file_exists(destfn))
                {
                    FILE* fp = fopen("data/dlsrc.dat", "rb");
                    if (!fp) raise_error("Please ensure data/dlsrc.dat file exists.");
                    while (!feof(fp))
                    {
                        fgets(buffer1, 1022, fp);
                        if (buffer1[0] == '#') continue;
                        char** dls = chop_spaced_words(buffer1);

                        if (!strcmp(dls[0], words[1]))
                        {
                            url = dls[1];
                            retfmt = dls[2];
                        }

                        delete[] dls;
                    }

                    if (!url.size()) raise_error("Download source not found in data file.");
                    j = url.find('%');
                    url.replace(j, 1, interpret_single_string(words[2]));

                    if (retfmt == "PDBDATA")
                    {
                        if (!download_file(url, destfn)) raise_error("Download failed.");
                    }
                    else if (retfmt.substr(0, 5) == "JSON:")
                    {
                        //
                    }
                    else raise_error("Unimplemented return format.");
                }
            }   // DOWNLOAD

            else if (!strcmp(words[0], "DUMP"))
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
                        if (script_var[j].value.ppt) cout << "point " << *script_var[j].value.ppt << endl;
                        else cout << "(empty)" << endl;
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

            else if (!strcmp(words[0], "ECHO"))
            {
                for (l=1; words[l]; l++)
                {
                    if (!strcmp(words[l], "~")) goto _no_newline_on_echo;
                    else
                    {
                        psz = interpret_single_string(words[l]);
                        cout << psz << flush;
                        delete[] psz;
                    }
                }

                cout << endl << flush;
            _no_newline_on_echo:
                ;
            }	// ECHO

            else if (!strcmp(words[0], "ELSE")) goto _pc_continue;

            else if (!strcmp(words[0], "END") || !strcmp(words[0], "EXIT") || !strcmp(words[0], "QUIT") || !strcmp(words[0], "DIE"))
            {
                if (words[1])
                {
                    l = interpret_single_int(words[1]);
                    if (l) return l;
                    else psz = interpret_single_string(words[1]);
                    if (strlen(psz)) cout << psz << endl;
                }
                return 0;
            }	// END

            else if (!strcmp(words[0], "GEN"))
            {
                if (!words[1]) raise_error("No sequence given for GEN.");
                psz = interpret_single_string(words[1]);

                working->add_sequence(psz);
                // working->conform_backbone(1, working->get_seq_length(), 50); // Takes too long.
                working->make_helix(1, working->get_seq_length(), M_PI, M_PI);
                goto _prot_deets;
            } // GEN

            else if (!strcmp(words[0], "GOTO"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for GOTO.");
                sprintf(buffer1, "%s:", words[1]);
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

            else if (!strcmp(words[0], "HELICES"))
            {
                int sr, er, hc, hf = 0;

                er = working->get_end_resno();
                sr = -1;
                for (hc = working->get_start_resno(); hc <= er; hc++)
                {
                    AminoAcid *aa = working->get_residue(hc), *aan = working->get_residue(hc+1);
                    if (!aa) continue;

                    bool is_helix = aa->is_alpha_helix();

                    // Allow single outliers.
                    if (aan) is_helix |= aan->is_alpha_helix();

                    if (is_helix)
                    {
                        if (sr < 0) sr = hc;
                    }
                    else
                    {
                        if (sr > 0 && sr < hc-4)
                        {
                            hf++;
                            std::string hxname = (std::string)"\%helix" + std::to_string(hf) + (std::string)".start";
                            Star s;
                            s.n = sr;
                            set_variable(hxname.c_str(), s);

                            hxname = (std::string)"\%helix" + std::to_string(hf) + (std::string)".end";
                            s.n = hc - 1;
                            set_variable(hxname.c_str(), s);

                        }
                        
                        sr = -1;
                    }
                }
            }   // HELICES

            else if (!strcmp(words[0], "HELIX"))
            {
                float phi, psi;
                l = 2;

                if (!words[1]) raise_error("No parameters given for HELIX.");
                if (!strcmp(words[1], "ALPHA"))
                {
                    phi = ALPHA_PHI;
                    psi = ALPHA_PSI;
                }
                else if (!strcmp(words[1], "PI"))
                {
                    phi = PI_PHI;
                    psi = PI_PSI;
                }
                else if (!strcmp(words[1], "3.10"))
                {
                    phi = _310_PHI;
                    psi = _310_PSI;
                }
                else if (!strcmp(words[1], "PPRO1"))
                {
                    phi = POLYPRO1_PHI;
                    psi = POLYPRO1_PSI;
                }
                else if (!strcmp(words[1], "PPRO2"))
                {
                    phi = POLYPRO2_PHI;
                    psi = POLYPRO2_PSI;
                }

                // Not technically helices, but included for consistency.
                else if (!strcmp(words[1], "BETA"))
                {
                    phi = BETA_PHI;
                    psi = BETA_PSI;
                }
                else if (!strcmp(words[1], "STRAIGHT"))
                {
                    phi = M_PI;
                    psi = M_PI;
                }

                else
                {
                    phi = interpret_single_float(words[1])*fiftyseventh;
                    if (!words[2]) raise_error("Insufficient parameters given for HELIX.");
                    psi = interpret_single_float(words[2])*fiftyseventh;
                    l++;
                }

                int sr, er, sa;
                if (!words[l]) raise_error("Insufficient parameters given for HELIX.");
                sr = interpret_single_int(words[l]);
                if (!words[l+1]) raise_error("Insufficient parameters given for HELIX.");
                sa = er = interpret_single_int(words[l+1]);
                if (words[l+2]) sa = interpret_single_int(words[l+2]);
                if (words[l+2] && words[l+3]) raise_error("Too many parameters given for HELIX.");

                working->make_helix(sr, er, sa, phi, psi);

            }	// HELIX

            else if (!strcmp(words[0], "HOMOLOGY"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for HOMOLOGY.");
                Protein tpl("template");
                
                words[1] = interpret_single_string(words[1]);
                FILE* fp = fopen(words[1], "rb");
                if (!fp) raise_error((std::string)"File not found: " + (std::string)words[1]);
                else
                {
                    tpl.load_pdb(fp);
                    fclose(fp);
                    if (words[2])
                    {
                        words[2] = interpret_single_string(words[2]);
                        fp = fopen(words[2], "rb");
                        if (!fp) raise_error((std::string)"File not found: " + (std::string)words[2]);
                        else
                        {
                            Protein ref("reference");
                            ref.load_pdb(fp);
                            fclose(fp);
                            working->homology_conform(&tpl, &ref);
                        }
                    }
                    else working->homology_conform(&tpl, working);
                }
            }   // HOMOLOGY

            else if (!strcmp(words[0], "HYDRO"))
            {
                int resno, endres = working->get_end_resno();
                cout << "Hydrogenating...";
                for (resno=1; resno<=endres; resno++)
                {
                    AminoAcid* res = working->get_residue(resno);
                    if (res)
                    {
                        res->hydrogenate();

                        AminoAcid** rcc = working->get_residues_can_clash(resno);
                        if (!rcc) continue;
                        for (n=0; rcc[n]; n++);
                        Molecule* mols[n+4];
                        for (i=0; i<=n; i++) mols[i] = (Molecule*)rcc[i];
                        float f = res->get_intermol_clashes(mols);
                        if (f >= 5)
                        {
                            working->minimize_residue_clashes(resno);
                        }
                    }
                    cout << "." << flush;
                }
                cout << endl;
            }   // HYDRO

            else if (!strcmp(words[0], "IF"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for IF.");
                if (!words[2]) raise_error("Insufficient parameters given for IF.");

                // TODO: NOT operator and IF %var GOTO blah.
                l = 2;

                // If the operator is =, and both l-value and r-value are strings, do a direct comparison.
                if (!words[l]) raise_error("Insufficient parameters given for IF.");
                if (!strcmp(words[l], "THEN"))
                {
                    l--;
                    if (interpret_single_float(words[l])) goto _evaluated_true;
                    else if (strlen(interpret_single_string(words[l]))) goto _evaluated_true;
                    else goto _evaluated_false;
                }

                if (!strcmp(words[l], "==")) words[l][1] = 0;

                if (!strcmp(words[l], "=") || !strcmp(words[l], "!=") || !strcmp(words[l], "=*"))
                {
                    if (!words[l+1]) raise_error("Insufficient parameters given for IF.");
                    if (words[l-1][0] == words[l+1][0]
                            ||
                            words[l-1][0] == '$' && words[l+1][0] == '"'
                            ||
                            words[l-1][0] == '"' && words[l+1][0] == '$'
                       )
                    {
                        char *lvalue = interpret_single_string(words[l-1]),
                              *rvalue = interpret_single_string(words[l+1]);
                        if (!strcmp(words[l], "="))
                        {
                            if (strcmp(lvalue, rvalue)) goto _evaluated_false;
                            else goto _evaluated_true;
                        }
                        if (!strcmp(words[l], "!="))
                        {
                            if (!strcmp(lvalue, rvalue)) goto _evaluated_false;
                            else goto _evaluated_true;
                        }
                        if (!strcmp(words[l], "=*"))
                        {
                            if (!strstr(lvalue, rvalue)) goto _evaluated_false;
                            else goto _evaluated_true;
                        }
                    }
                    else goto _just_interpret_floats;
                }
                else
                {
                    // Otherwise, interpret both values as floats.
                _just_interpret_floats:
                    float lvalue = interpret_single_float(words[l-1]),
                          rvalue = interpret_single_float(words[l+1]);

                    if (!strcmp(words[l], "="))
                    {
                        if (lvalue == rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(words[l], "!="))
                    {
                        if (lvalue != rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(words[l], ">"))
                    {
                        if (lvalue > rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(words[l], "<"))
                    {
                        if (lvalue < rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(words[l], ">="))
                    {
                        if (lvalue >=rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else if (!strcmp(words[l], "<="))
                    {
                        if (lvalue <= rvalue) goto _evaluated_true;
                        else goto _evaluated_false;
                    }
                    else raise_error( (std::string)"Unknown operator " + (std::string)words[l] + (std::string)" for comparison.");
                }

            _evaluated_true:
                l += 2;
                if (!strcmp(words[l], "THEN")) l++;
                words = &words[l];

                if (!words[0]) goto _pc_continue;

                while (!strcmp(words[0], "OR"))
                {
                    words = &words[4];
                    if (!words[0]) goto _pc_continue;
                }

                if (!strcmp(words[0], "AND")) strcpy(words[0], "IF");
                goto _interpret_command;

            _evaluated_false:
                if (words[l+2] && !strcmp(words[l+2], "OR"))
                {
                    l += 2;
                    words = &words[l];
                    strcpy(words[0], "IF");
                    goto _interpret_command;
                }

                program_counter++;
                if (!script_lines[program_counter].c_str()) break;
                strcpy(buffer, script_lines[program_counter].c_str());
                words = chop_spaced_words(buffer);
                if (!words || !words[0]) goto _pc_continue;
                if (strcmp(words[0], "ELSE")) continue;

                words = &words[1];
                goto _interpret_command;
            }	// IF

            else if (!strcmp(words[0], "LET"))
            {
                if (!words[1]) raise_error("No parameters given for LET.");
                n = find_var_index(words[1], &words[1]);
                if (n<0)
                {
                    n = vars++;
                    script_var[n].name = words[1];
                    script_var[n].vt = type_from_name(words[1]);
                }
                int flags = n & _VARFLAGS_MASK;
                n &= _VARNUM_MASK;

                if (!words[2]) raise_error("No operator given for LET.");
                if (!words[3] && strcmp(words[2], "++") && strcmp(words[2], "--") ) raise_error("No rvalue given for LET.");
                switch (script_var[n].vt)
                {
                case SV_INT:
                    if (!strcmp(words[2], "=")) script_var[n].value.n = interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "+=")) script_var[n].value.n += interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "-=")) script_var[n].value.n -= interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "*=")) script_var[n].value.n *= interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "/=")) script_var[n].value.n /= interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "&=")) script_var[n].value.n &= interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "|=")) script_var[n].value.n |= interpret_single_int(words[3]);
                    else if (!strcmp(words[2], "++")) script_var[n].value.n++;
                    else if (!strcmp(words[2], "--")) script_var[n].value.n--;
                    else
                    {
                        raise_error( (std::string)"Unimplemented operator " + (std::string)words[2] + (std::string)" for int assignment.");
                        return 0x51974c5;
                    }
                    l=0;
                    break;

                case SV_FLOAT:
                    if (!strcmp(words[2], "=")) script_var[n].value.f = interpret_single_float(words[3]);
                    else if (!strcmp(words[2], "+=")) script_var[n].value.f += interpret_single_float(words[3]);
                    else if (!strcmp(words[2], "-=")) script_var[n].value.f -= interpret_single_float(words[3]);
                    else if (!strcmp(words[2], "*=")) script_var[n].value.f *= interpret_single_float(words[3]);
                    else if (!strcmp(words[2], "/=")) script_var[n].value.f /= interpret_single_float(words[3]);
                    else
                    {
                        raise_error( (std::string)"Unimplemented operator " + (std::string)words[2] + (std::string)" for float assignment.");
                        return 0x51974c5;
                    }
                    l=0;
                    break;

                case SV_POINT:
                    if (flags & _HAS_DOT)
                    {
                        double* ff = nullptr;

                        char* param = strchr(words[1], '.');
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
                            if (!strcmp(words[2], "=")) *ff = interpret_single_float(words[3]);
                            else if (!strcmp(words[2], "+=")) *ff += interpret_single_float(words[3]);
                            else if (!strcmp(words[2], "-=")) *ff -= interpret_single_float(words[3]);
                            else if (!strcmp(words[2], "*=")) *ff *= interpret_single_float(words[3]);
                            else if (!strcmp(words[2], "/=")) *ff /= interpret_single_float(words[3]);
                            else
                            {
                                raise_error( (std::string) "Unimplemented operator " + (std::string)words[2] + (std::string)" for float assignment.");
                                return 0x51974c5;
                            }
                            l = 0;
                            n = -1;
                        }
                        break;
                    }
                    else
                    {
                        if (!strcmp(words[2], "="))
                        {
                            if (!script_var[n].value.ppt) script_var[n].value.ppt = new Point();
                            *(script_var[n].value.ppt) = interpret_single_point(words[3], *script_var[n].value.ppt);
                        }
                        else if (!strcmp(words[2], "+=")) *(script_var[n].value.ppt) = script_var[n].value.ppt->add(interpret_single_point(words[3]));
                        else if (!strcmp(words[2], "-=")) *(script_var[n].value.ppt) = script_var[n].value.ppt->subtract(interpret_single_point(words[3]));
                        else if (!strcmp(words[2], "*="))
                        {
                            f = interpret_single_float(words[3]);
                            script_var[n].value.ppt->x *= f;
                            script_var[n].value.ppt->y *= f;
                            script_var[n].value.ppt->z *= f;
                        }
                        else if (!strcmp(words[2], "/="))
                        {
                            f = interpret_single_float(words[3]);
                            script_var[n].value.ppt->x /= f;
                            script_var[n].value.ppt->y /= f;
                            script_var[n].value.ppt->z /= f;
                        }
                        else
                        {
                            raise_error( (std::string)"Unimplemented operator " + (std::string)words[2] + (std::string)" for Cartesian assignment.");
                            return 0x51974c5;
                        }
                    }

                    l=0;

                    break;

                case SV_STRING:

                    psz = interpret_single_string(words[3]);

                    l = m = 0;
                    if (words[4] && !strcmp(words[4], "FROM"))
                    {
                        l+=2;
                        if (!words[5]) raise_error("Insufficient parameters given for LET.");
                        m = interpret_single_int(words[5]);
                        if (m < 0) m = 0;
                        if (m)
                        {
                            if (m > strlen(psz)) psz[0] = 0;
                            else
                            {
                                m--;
                                psz += m;
                                if (words[6] && !strcmp(words[6], "FOR"))
                                {
                                    l+=2;
                                    if (!words[7]) raise_error("Insufficient parameters given for LET.");
                                    k = interpret_single_int(words[7]);
                                    if (k >= 0 && k < strlen(psz)) psz[k] = 0;
                                }
                            }
                        }
                    }

                    if (!strcmp(words[2], "="))
                    {
                        script_var[n].value.psz = new char[65536];
                        strcpy(script_var[n].value.psz, psz);
                    }
                    else if (!strcmp(words[2], "+="))
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
                        raise_error( (std::string)"Unimplemented operator " + (std::string)words[2] + (std::string)" for string assignment.");
                        return 0x51974c5;		// If you use your imagination, that says "syntax".
                    }
                    psz -= m;
                    delete[] psz;
                    break;

                default:
                    ;
                }	// switch (script_var[n].vt)

                while (n >= 0 && words[3+l] && words[4+l] && words[5+l])
                {
                    if (!words[5+l]) raise_error("Insufficient parameters given for LET.");
                    switch (script_var[n].vt)
                    {
                    case SV_INT:
                        m = interpret_single_int(words[5+l]);
                        if (!strcmp(words[4+l], "+")) script_var[n].value.n += m;
                        else if (!strcmp(words[4+l], "-")) script_var[n].value.n -= m;
                        else if (!strcmp(words[4+l], "*")) script_var[n].value.n *= m;
                        else if (!strcmp(words[4+l], "/")) script_var[n].value.n /= m;
                        else if (!strcmp(words[4+l], "^")) script_var[n].value.n = pow(script_var[n].value.n, m);
                        else if (!strcmp(words[4+l], "&")) script_var[n].value.n &= m;
                        else if (!strcmp(words[4+l], "|")) script_var[n].value.n |= m;
                        else
                        {
                            raise_error( (std::string)"Bad operator " + (std::string)words[4+l] + (std::string)" for int.");
                            return 0x51974c5;
                        }
                        break;

                    case SV_FLOAT:
                        f = interpret_single_float(words[5+l]);
                        if (!strcmp(words[4+l], "+")) script_var[n].value.f += f;
                        else if (!strcmp(words[4+l], "-")) script_var[n].value.f -= f;
                        else if (!strcmp(words[4+l], "*")) script_var[n].value.f *= f;
                        else if (!strcmp(words[4+l], "/")) script_var[n].value.f /= f;
                        else if (!strcmp(words[4+l], "^")) script_var[n].value.f = pow(script_var[n].value.f, f);
                        else
                        {
                            // cout << "Bad operator " << words[4+l] << " for float." << endl;
                            raise_error( (std::string)"Bad operator " + (std::string)words[4+l] + (std::string)" for float.");
                            return 0x51974c5;
                        }
                        break;

                    case SV_POINT:
                        pt = interpret_single_point(words[5+l]);
                        if (!strcmp(words[4+l], "+")) *script_var[n].value.ppt = script_var[n].value.ppt->add(pt);
                        else if (!strcmp(words[4+l], "-")) *script_var[n].value.ppt = script_var[n].value.ppt->subtract(pt);
                        else if (!strcmp(words[4+l], "*")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() * pt.magnitude());
                        else if (!strcmp(words[4+l], "/")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() / pt.magnitude());
                        else
                        {
                            // cout << "Bad operator " << words[4+l] << " for point." << endl;
                            raise_error( (std::string)"Bad operator " + (std::string)words[4+l] + (std::string)" for point.");
                            return 0x51974c5;
                        }
                        break;

                    case SV_STRING:
                        builder = script_var[n].value.psz;
                        psz = interpret_single_string(words[5+l]);
                        if (!strcmp(words[4+l], "+"))
                        {
                            builder.append(psz);
                            script_var[n].value.psz = new char[65536];
                            strcpy(script_var[n].value.psz, builder.c_str());
                        }
                        else if (!strcmp(words[4+l], "FOR"))
                        {
                            l += 2;
                            continue;
                        }
                        else
                        {
                            // cout << "Bad operator " << words[4+l] << " for string." << endl;
                            raise_error( (std::string)"Bad operator " + (std::string)words[4+l] + (std::string)" for string.");
                            return 0x51974c5;
                        }
                        break;

                    default:
                        ;
                    }

                    l += 2;
                }
            }	// LET

            else if (!strcmp(words[0], "LOAD"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for LOAD.");
                psz = interpret_single_string(words[1]);
				n = 0;
                chain = 'A';
                
				if (words[2]) chain = words[2][0];
                if (words[2] && words[3])
                {
                    new_chain = words[3][0];
                    if (words[4]) raise_error("Too many parameters given for LOAD.");

                    if (!strands[new_chain - 65]) strands[new_chain - 65] = new Protein(psz);
                    working = strands[new_chain - 65];
                    g_chain = new_chain;
                }

                pf = fopen(psz, "rb");
                if (!pf)
                {
                    raise_error( (std::string)"Failed to open " + (std::string)psz + (std::string)" for reading.");
                    return 0xbadf12e;
                }
                working->load_pdb(pf, n, chain);

                fclose(pf);

            _prot_deets:

				char buffer[1024];
				char buffer1[1024];
				Star sv;

				// strcpy(buffer, "$PDB");
                sprintf(buffer, "$PDB%c", g_chain);
				sv.psz = new char[strlen(psz)+4];
                strcpy(sv.psz, psz);
				set_variable(buffer, sv);

                const char* pname = working->get_name().c_str();
                // strcpy(buffer, "$PROTEIN");
                sprintf(buffer, "$PROTEIN%c", g_chain);
                sv.psz = new char[strlen(pname)+4];
                strcpy(sv.psz, pname);
				set_variable(buffer, sv);

                delete[] psz;

                int seqlen = working->get_seq_length();
                // strcpy(buffer, "%SEQLEN");
                sprintf(buffer, "%cSEQLEN%c", '%', g_chain);
                sv.n = seqlen;
				set_variable(buffer, sv);

                // strcpy(buffer, "$SEQUENCE");
                sprintf(buffer, "$SEQUENCE%c", g_chain);
                sv.psz = new char[seqlen+4];
                strcpy(sv.psz, working->get_sequence().c_str());
				set_variable(buffer, sv);

                std::vector<std::string> rem_hx = working->get_remarks("650 HELIX");
                for (l=0; l<rem_hx.size(); l++)
                {
                    strcpy(buffer, rem_hx[l].c_str());
                    char** words = chop_spaced_words(buffer);

					if (!words[3] || !words[4] || !words[5]) continue;

                    sprintf(buffer1, "%c%c.%s.s", '%', g_chain, words[3]);
                    sv.n = atoi(words[4]);
                    set_variable(buffer1, sv);

                    sprintf(buffer1, "%c%c.%s.e", '%', g_chain, words[3]);
                    sv.n = atoi(words[5]);
                    set_variable(buffer1, sv);

					working->set_region(words[3], atoi(words[4]), atoi(words[5]));

                    delete[] words;
                }

                for (l=1; l<=7; l++)
                {
                    int bw50 = working->get_bw50(l);
                    if (bw50 > 0)
                    {
                        sprintf(buffer1, "%c%c.%d.50", '%', g_chain, l);
                        sv.n = bw50;
                        set_variable(buffer1, sv);
                    }
                }

                for (l=12; l<=67; l+=11)
                {
                    int bw50 = working->get_bw50(l);
                    if (bw50 > 0)
                    {
                        sprintf(buffer1, "%c%c.%d.50", '%', g_chain, l);
                        sv.n = bw50;
                        set_variable(buffer1, sv);
                    }
                }
            }

            else if (!strcmp(words[0], "MCOORD"))
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
                if (!strcmp(words[l], "YO"))
                {
                    force_tyrosine_O = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(words[l], "YAr"))
                {
                    force_tyrosine_O = false;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(words[l], "Th8"))
                {
                    thiolate = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }

                if (!words[l]) raise_error("Insufficient parameters given for MCOORD.");
                elem_sym = words[l++];
                if (!words[l]) raise_error("Insufficient parameters given for MCOORD.");
                elem_charge = interpret_single_int(words[l++]);
                Point pt;
                ma = new Atom(elem_sym.c_str(), &pt, elem_charge);
                n = working->get_metals_count() + 1;
                string mname = elem_sym + std::to_string(n);
                ma->name = new char[8];
                strcpy(ma->name, mname.c_str());
                strcpy(ma->aa3let, "MTL");
                ma->residue = 0;

                if (!words[l]) raise_error("Insufficient parameters given for MCOORD.");
                for (; words[l]; l++)
                {
                    bool local_O = force_tyrosine_O;

                    _another_goto:
                    if (!strcmp(words[l], "YO"))
                    {
                        local_O = true;
                        l++;
                        goto _another_goto;
                    }
                    else if (!strcmp(words[l], "YAr"))
                    {
                        local_O = false;
                        l++;
                        goto _another_goto;
                    }

                    if (!words[l]) raise_error("Insufficient parameters given for MCOORD.");
                    k = interpret_single_int(words[l]);
                    AminoAcid* aa = working->get_residue(k);
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
                working->coordinate_metal(ma, ncr, resnos, cratoms);

            } // MCOORD

            else if (!strcmp(words[0], "MOVE"))
            {
                l = 1;
                Point newcen(0,0,0);
                int sr, er;
                if (words[l]) sr = interpret_single_int(words[l++]);
                else raise_error("Not enough parameters given for MOVE.");
                if (words[l]) er = interpret_single_int(words[l++]);
                else raise_error("Not enough parameters given for MOVE.");
                if (words[l]) newcen = interpret_single_point(words[l++]);
                if (words[l]) raise_error("Too many parameters given for MOVE.");
                working->move_piece(sr, er, newcen);
            }	// MOVE

            else if (!strcmp(words[0], "PTALIGN"))
            {
                Point point, align, center;
                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for PTALIGN.");
                point = interpret_single_point(words[l++]);

                if (!words[l]) raise_error("Insufficient parameters given for PTALIGN.");
                align = interpret_single_point(words[l++]);

                if (!words[l]) raise_error("Insufficient parameters given for PTALIGN.");
                center = interpret_single_point(words[l++]);
                
                Rotation rot = align_points_3d(&point, &align, &center);

                if (!words[l]) raise_error("Insufficient parameters given for PTALIGN.");
                Star s;
                rot.v.r = 1;
                s.ppt = new Point(rot.v);
                set_variable(words[l++], s);

                if (!words[l]) raise_error("Insufficient parameters given for PTALIGN.");
                s.f = rot.a * fiftyseven;
                set_variable(words[l++], s);

                if (words[l]) raise_error("Too many parameters given for PTALIGN.");
            }

            else if (!strcmp(words[0], "PTROTATE"))
            {
                Point point, origin, axis;
                float theta;
                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for PTROTATE.");
                point = interpret_single_point(words[l++]);

                if (!words[l]) raise_error("Insufficient parameters given for PTROTATE.");
                origin = interpret_single_point(words[l++]);

                if (!words[l]) raise_error("Insufficient parameters given for PTROTATE.");
                axis = interpret_single_point(words[l++]);

                if (!words[l]) raise_error("Insufficient parameters given for PTROTATE.");
                theta = interpret_single_float(words[l++]) * fiftyseventh;

                Point result = rotate3D(point, origin, axis, theta);

                if (!words[l]) raise_error("Insufficient parameters given for PTROTATE.");
                Star s;
                s.ppt = &result;
                set_variable(words[l++], s);

                if (words[l]) raise_error("Too many parameters given for PTROTATE.");
            }

            else if (!strcmp(words[0], "REGION"))
            {
                int sr, er;
                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for REGION.");
                psz = interpret_single_string(words[l++]);
                if (!psz[0]) psz = words[l-1];
                if (!words[l]) raise_error("Insufficient parameters given for REGION.");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for REGION.");
                er = interpret_single_int(words[l++]);
                if (words[l]) raise_error("Too many parameters given for REGION.");

                working->set_region(psz, sr, er);

                char lbuffer[256];
                sprintf(lbuffer, "%s%s.s", "%", psz);
                Star s;
                s.n = sr;
                set_variable(lbuffer, s);
                sprintf(lbuffer, "%s%s.e", "%", psz);
                s.n = er;
                set_variable(lbuffer, s);
            }	// REGION

            else if (!strcmp(words[0], "REMARK"))
            {
                sprintf(buffer1, "%s\n", script_lines[program_counter].c_str());
                for (l=0; buffer1[l] != 'R'; l++);
                working->add_remark(buffer1+l);
            }   // REMARK

            else if (!strcmp(words[0], "RENUMBER"))
            {
				l = 1;
                int sr, er, nsr;
                if (!words[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)words[0] + (std::string)".");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)words[0] + (std::string)".");
                er = interpret_single_int(words[l++]);
                if (!words[l]) raise_error((std::string)"Insufficient parameters given for " + (std::string)words[0] + (std::string)".");
                nsr = interpret_single_int(words[l++]);
                if (words[l]) raise_error((std::string)"Too many parameters given for " + (std::string)words[0] + (std::string)".");

				working->renumber_residues(sr, er, nsr);
            }	// RENUMBER

            else if (!strcmp(words[0], "ROTATE"))
            {
				if (!words[1]) raise_error("Insufficient parameters given for ROTATE.");
                SCoord axis = interpret_single_point(words[1]);
				if (!words[2]) raise_error("Insufficient parameters given for ROTATE.");
                float theta = interpret_single_float(words[2]) * fiftyseventh;
                if (words[3]) raise_error("Too many parameters given for ROTATE.");

                LocatedVector lv = axis;
                cout << "Axis is " << (Point)axis << endl;
                int sr = working->get_start_resno(), er = working->get_end_resno();
                lv.origin = working->get_region_center(sr, er);

                for (i=sr; i<=er; i++)
                {
                    AminoAcid* aa = working->get_residue(i);
                    if (aa)
                    {
                        MovabilityType fmov = aa->movability;
                        aa->movability = MOV_ALL;
                        aa->rotate(lv, theta);
                        aa->movability = fmov;
                    }
                }
            }

            else if (!strcmp(words[0], "ROTBOND"))
            {
				l = 1;
                raise_error((std::string)"Unimplemented command " + (std::string)words[0] + (std::string)".");
            }

            else if (!strcmp(words[0], "SAVE"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for SAVE.");
                psz = interpret_single_string(words[1]);
                if (words[2] && words[3]) raise_error("Too many parameters given for SAVE.");

                pf = fopen(psz, "wb");
                if (!pf)
                {
                    raise_error( (std::string)"Failed to open " + (std::string)psz + (std::string)" for writing.");
                    return 0xbadf12e;
                }
                for (l=0; l<26; l++)
                {
                    if (!strands[l]) continue;
                    working = strands[l];
                    g_chain = l+65;
                    working->set_pdb_chain(l+65);
                    working->save_pdb(pf);
                }
                working->end_pdb(pf);

                cout << "Wrote " << psz << "." << endl;

                fclose(pf);
                delete[] psz;

                if (words[2])
                {
                    if ( !strcmp("QUIT", words[2]) || !strcmp("EXIT", words[2]) || !strcmp("END", words[2]) )
                        return 0;
                    else raise_error((std::string)"Too many parameters given for SAVE: " + (std::string)words[2]);
                }
            }   // SAVE

            else if (!strcmp(words[0], "SCENV"))
            {
                l = 1;
                if (!words[l]) raise_error("Insufficient parameters given for SCENV.");
                int resno = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for SCENV.");

                std::string result = "";
                int bitmask = words[l+1] ? interpret_single_int(words[l++]) : 0xffffffff;

                // Find all residues close enough to interact side chains with resno.
                AminoAcid** nearby = working->get_residues_can_clash(resno);
                AminoAcid* residue = working->get_residue(resno);
                int ra = residue->get_atom_count();
                if (nearby)
                {
                    // For each side chain atom of nearby residues, each side chain atom of resno,
                    // if the atoms are close enough and the right kinds to interact,
                    // concatenate the letter and residue number to the output string.
                    for (i=0; nearby[i]; i++)
                    {
                        bool matched = false;
                        int na = nearby[i]->get_atom_count();
                        for (j=0; j<na; j++)
                        {
                            Atom* a = nearby[i]->get_atom(j);
                            if (a->is_backbone) continue;

                            for (k=0; k<ra; k++)
                            {
                                Atom* b = residue->get_atom(k);
                                if (b->is_backbone) continue;

                                // For debugging.
                                if (resno == 185 && nearby[i]->get_residue_no() == 193
                                    // && !strcmp(a->name, "1HE2")
                                    // && !strcmp(b->name, "OE2")
                                    )
                                {
                                    // cout << a->name << " " << b->name << endl;
                                }

                                float r = a->distance_to(b);
                                if (r < 7)
                                {
                                    // TODO:
                                    InteratomicForce** iff = InteratomicForce::get_applicable(a, b);
                                    if (iff) for (m=0; iff[m]; m++)
                                    {
                                        if (r < 1.333 * iff[m]->get_distance())
                                        {
                                            // TODO: Implement a types bitmask.
                                            switch (iff[m]->get_type())
                                            {
                                                case vdW:
                                                if (bitmask & 0x01) matched = true;
                                                break;

                                                case hbond:
                                                if (bitmask & 0x02) matched = true;
                                                break;

                                                case ionic:
                                                case mcoord:
                                                if (bitmask & 0x04) matched = true;
                                                break;

                                                case pi:
                                                case polarpi:
                                                if (bitmask & 0x08) matched = true;
                                                break;

                                                default:
                                                ;
                                            }
                                            
                                            if (matched) goto _exit_atomloop;
                                        }
                                    }
                                }
                            }
                        }

                        _exit_atomloop:
                        if (matched)
                        {
                            sprintf(buffer, "%c%d ", nearby[i]->get_letter(), nearby[i]->get_residue_no());
                            result += buffer;
                        }
                    }
                }
                
                Star sv;
                strcpy(buffer, words[l++]);
                strcpy(buffer1, result.c_str());
                sv.psz = buffer1;
                set_variable(buffer, sv);

                if (words[l]) raise_error("Too many parameters given for SCENV.");
            }   // SCENV

            else if (!strcmp(words[0], "SEARCH"))
            {
                l = 1;
                int sr, er, esr, sim;
                if (!words[l]) raise_error("Insufficient parameters given for SEARCH.");
                sr = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for SEARCH.");
                er = interpret_single_int(words[l++]);
                if (!words[l]) raise_error("Insufficient parameters given for SEARCH.");
                psz = interpret_single_string(words[l++]);
                esr = er; // - strlen(psz);

                int threshold = -1;
                int num_eq;

                if (!words[l]) raise_error("Insufficient parameters given for SEARCH.");
                // if (words[l+1]) threshold = interpret_single_int(words[l++]);
                if (!strcmp(words[l], "TH"))
                {
                    l++;
                    threshold = interpret_single_int(words[l++]);
                }

                k = working->search_sequence(sr, esr, psz, threshold, &sim);

                delete[] psz;

                n = find_var_index(words[l]);
                if (n<0) n = vars++;
                n &= _VARNUM_MASK;
                if (!words[l]) raise_error("Insufficient parameters given for SEARCH.");
                if (words[l+1] && words[l+2]) raise_error("Too many parameters given for SEARCH.");
                script_var[n].name = words[l];
                script_var[n].vt = type_from_name(words[l]);

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
                if (words[l])
                {
                    n = find_var_index(words[l]);
                    if (n<0) n = vars++;
                    n &= _VARNUM_MASK;
                    script_var[n].name = words[l];
                    script_var[n].vt = type_from_name(words[l]);

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

            else if (!strcmp(words[0], "STRAND"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for STRAND.");
                if (words[2])
                {
                    if (words[3]) raise_error("Too many parameters given for STRAND.");
                    chain = words[1][0] - 65;
                    char new_chain = words[2][0] - 65;
                    strands[new_chain] = strands[chain];
                    strands[chain] = nullptr;
                    working = strands[new_chain];
                    g_chain = new_chain+65;
                }
                else
                {
                    chain = words[1][0] - 65;
                    if (!strands[chain]) strands[chain] = new Protein("Another Prot");
                    working = strands[chain];
                    g_chain = chain+65;
                }
            }   // STRAND

            else if (!strcmp(words[0], "UNCHAIN"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for UNCHAIN.");
                chain = words[1][0] - 65;
                if (working == strands[chain]) raise_error("Cannot delete current working strand.");
                if (nullptr != strands[chain])
                {
                    delete strands[chain];
                    strands[chain] = nullptr;
                }
            }

            else if (!strcmp(words[0], "UPRIGHT"))
            {
                l = 1;
                if (words[l]) raise_error("Too many parameters given for UPRIGHT.");

                try
                {
                    working->upright();
                }
                catch (int ex)
                {
                    if (ex == 0xbad7312) raise_error("Cannot UPRIGHT protein without transmembrane regions named TMR{n}.");
                    else raise_error("Unknown error.");
                }
            }	// UPRIGHT

            else
            {
                raise_error( (std::string)"Unimplemented command: \"" + (std::string)words[0] + (std::string)"\"");
                return 0x51974c5;
            }
        }

    _pc_continue:
        delete[] owords;
        program_counter++;
    }

    return 0;
}










