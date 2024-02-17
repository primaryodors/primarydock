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
std::vector<int> call_stack;
int stack_pointer = 0;

Protein* strands[26];
Protein* working = nullptr;
char g_chain = 'A';
Molecule ligand("ligand");
char lig_chain = '\0';

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

#define debug_contact_energy 0

float contact_energy(Protein* a, Protein* b, bool repack = false, int a_from = 0, int a_to = 0, int b_from = 0, int b_to = 0)
{
    std::vector<AminoAcid*> vca = a->get_contact_residues(b, 2);
    float f=0;

    char cha = a->get_pdb_chain(), chb = b->get_pdb_chain();
    int i, j=0, n = vca.size();
    int aend = a->get_end_resno(), bend = b->get_end_resno();
    bool a_dirty[aend+1], b_dirty[bend+1];

    for (i=0; i<=aend; i++) a_dirty[i] = false;
    for (i=0; i<=bend; i++) b_dirty[i] = false;

    Molecule* mols[n+4];
    for (i=0; i<n; i++)
    {
        if (!vca[i]) continue;
        char c = vca[i]->get_pdb_chain();
        int resno = vca[i]->get_residue_no();

        if (c == cha)
        {
            if (a_from && (resno < a_from)) continue;
            if (a_to && (resno > a_to)) continue;
            if (a_dirty[resno]) continue;
            a_dirty[resno] = true;
        }
        else if (c == chb)
        {
            if (b_from && (resno < b_from)) continue;
            if (b_to && (resno > b_to)) continue;
            if (b_dirty[resno]) continue;
            b_dirty[resno] = true;
        }

        mols[j++] = (Molecule*)vca[i];
        #if debug_contact_energy
        cout << *vca[i] << endl;
        #endif
    }
    mols[j] = nullptr;

    if (repack) Molecule::conform_molecules(mols, 20);

    for (i=0; mols[i]; i++)
    {
        if (!mols[i]) continue;
        f -= mols[i]->get_intermol_binding(mols);
    }

    #if debug_contact_energy
    cout << endl;
    #endif

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

    if (c1 && *(c1-1) == '.') c1 = nullptr;
    if (c2 && *(c2-1) == '.') c2 = nullptr;

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
                    // TODO: Make this a stack allocation, instead of heap.
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
        if (script_var[n].value.psz) delete[] script_var[n].value.psz;
        script_var[n].value.psz = new char[strlen(vvalue.psz)+4];
        strcpy(script_var[n].value.psz, vvalue.psz);
    }
    else script_var[n].value = vvalue;

    return n;
}

int interpret_special_resno(const char* varname, char* strandid = nullptr)
{
    int rgno, member;
    Protein* p = working;

    if ( varname[0] >= '0' && varname[0] <= '9'
        && ( varname[1] == '.' || varname[1] == 'x' )
        && varname[2] >= '0' && varname[2] <= '9'
        )
    {
        rgno = varname[0] - '0';
        member = atoi(varname+2);
        if (strandid) *strandid = g_chain;
    }
    else if ( varname[0] >= '0' && varname[0] <= '9'
        && varname[1] == varname[0] + 1
        && ( varname[2] == '.' || varname[2] == 'x' )
        && varname[3] >= '0' && varname[3] <= '9'
        )
    {
        rgno = 10*(varname[0] - '0') + (varname[1] - '0');
        member = atoi(varname+3);
        if (strandid) *strandid = g_chain;
    }
    else if ( varname[0] >= 'A' && varname[0] <= 'Z'
        && varname[1] == '.'
        && varname[2] >= '0' && varname[2] <= '9'
        && ( varname[3] == '.' || varname[3] == 'x' )
        && varname[4] >= '0' && varname[4] <= '9'
        )
    {
        p = strands[varname[0] - 'A'];
        rgno = varname[2] - '0';
        member = atoi(varname+4);
        if (strandid) *strandid = varname[0];
    }
    else if ( varname[0] >= 'A' && varname[0] <= 'Z'
        && varname[1] == '.'
        && varname[2] >= '0' && varname[2] <= '9'
        && varname[3] == varname[2] + 1
        && ( varname[4] == '.' || varname[4] == 'x' )
        && varname[5] >= '0' && varname[5] <= '9'
        )
    {
        p = strands[varname[0] - 'A'];
        rgno = 10*(varname[2] - '0') + (varname[3] - '0');
        member = atoi(varname+5);
        if (strandid) *strandid = varname[0];
    }
    else return 0;

    int bw50 = p->get_bw50(rgno);
    if (!bw50) return 0;
    else return max(0, bw50 + member - 50);
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
    if (param[0] == '%')
    {
        int resno = interpret_special_resno(param+1);
        if (resno) return resno;
    }
    return round(interpret_single_float(param));
}

Point interpret_single_point(const char* param, Point old_value = Point(0,0,0))
{
    int n;
    Point pt = old_value;
    AminoAcid* aa;

    if (param[0] == '@')
    {
        char lchain = g_chain;
        int resno = interpret_special_resno(param+1, &lchain);
        if (resno)
        {
            AminoAcid* aa = strands[lchain - 'A']->get_residue(resno);
            if (!aa)
            {
                raise_error((std::string)"Residue " + std::to_string(resno) + (std::string)" is missing from strand.");
            }
            return aa->get_CA_location();
        }
    }

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
    if (!param) return nullptr;

    int n;
    char* buffer = new char[65536];
    for (n=0; n<65536; n++) buffer[n] = 0;
    char lchain = g_chain;
    int resno;

    switch (param[0])
    {
    case '%':
        sprintf(buffer, "%d", interpret_single_int(param));
        return buffer;

    case '&':
        sprintf(buffer, "%f", interpret_single_float(param));
        return buffer;

    case '@':
        n = find_var_index(param);
        if (n >= 0 &  _HAS_DOT)
        {
            sprintf(buffer, "%f", interpret_single_float(param));
            return buffer;
        }
        sprintf(buffer, "%s", interpret_single_point(param).printable().c_str());
        return buffer;

    case '$':
        resno = interpret_special_resno(param+1, &lchain);
        if (resno)
        {
            sprintf(buffer, "%c", strands[lchain - 'A']->get_residue(resno)->get_letter());
            return buffer;
        }

        n = find_var_index(param);
        if (n<0) return buffer;
        if (n &  _HAS_DOT)
        {
            param = strchr(param, '.');
            if (!param) return 0;
            param++;
            char lbuffer[1024];
            strcpy(lbuffer, param);
            char* dot2 = strchr(lbuffer, ',');
            n &= _VARNUM_MASK;
            if (!script_var[n].value.psz) return buffer;
            if (dot2)
            {
                *dot2 = 0;
                dot2++;
            }
            int i = interpret_single_int(lbuffer);
            if (!i) raise_error("Subscript out of range.");
            if (i > strlen(script_var[n].value.psz)) raise_error("Subscript out of range.");
            strcpy(buffer, script_var[n].value.psz+i-1);
            if (dot2) buffer[interpret_single_int(dot2)] = 0;
            return buffer;
        }
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
    Point pt3[4];
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
        bool include_ligand = false;
        Protein ptmp("RipMeAsunder");
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
                bool keeps = false, keepe = false;
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
                    if (!strcmp(words[l], "KEEP")) keeps = true;
                    sp = interpret_single_point(words[l++]);
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    aer = interpret_single_int(words[l++]);   
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    if (!strcmp(words[l], "KEEP")) keepe = true;
                    ep = interpret_single_point(words[l++]);
                }
                else
                {
                    asr = sr;
                    aer = er;                    
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    if (!strcmp(words[l], "KEEP")) keeps = true;
                    sp = interpret_single_point(words[l++]);
                    if (!words[l]) raise_error("Insufficient parameters given for ALIGN.");
                    if (!strcmp(words[l], "KEEP")) keepe = true;
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

                if (keeps) sp = sl;
                if (keepe)
                {
                    ep = Point(0,0,0);

                    for (i=0; i<eachend; i++)
                    {
                        AminoAcid* aa = working->get_residue(aer-i);
                        if (aa) ep = ep.add(aa->get_CA_location());
                    }

                    if (eachend > 1) ep.scale(ep.magnitude()/eachend);
                }

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

            else if (!strcmp(words[0], "ATOMTO"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for ATOMTO.");
                int resno = interpret_single_int(words[1]);
                if (!words[2]) raise_error("Insufficient parameters given for ATOMTO.");
                char* aname = interpret_single_string(words[2]);
                if (!words[3]) raise_error("Insufficient parameters given for ATOMTO.");
                Point target = interpret_single_point(words[3]);
                if (words[4]) raise_error("Too many parameters given for ATOMTO.");

                AminoAcid* aa = working->get_residue(resno);
                if (!aa) raise_error("Residue not found.");

                Atom* a = aa->get_atom(aname);
                if (!strcmp("EXTENT", aname)) a = aa->get_reach_atom();
                if (!a) raise_error("Atom not found.");

                aa->movability = MOV_FLEXONLY;
                aa->conform_atom_to_location(a->name, target);
            }	// ATOMTO

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

            else if (!strcmp(words[0], "BWCENTER"))
            {
                if (words[1]) raise_error("Too many parameters given for BWCENTER.");
                Point bw50_CA[10];
                for (l=1; l<=7; l++)
                {
                    int resno = working->get_bw50(l);
                    if (resno < 1) raise_error((std::string)"Strand does not have a " + to_string(l) + (std::string)".50 residue.");
                    bw50_CA[l] = working->get_residue(resno)->get_CA_location();
                    bw50_CA[l].weight = 1;
                }
                Point old_center = average_of_points(&bw50_CA[1], 7);
                Point new_center = working->get_region_center(1, working->get_end_resno()).subtract(old_center);
                working->move_piece(1, working->get_end_resno(), new_center);

                if (lig_chain == g_chain && ligand.get_atom_count())
                {
                    // cout << old_center << " -> " << new_center << endl;
                    ligand.move(Point(0,0,0).subtract(old_center));
                }
            }   // BWCENTER

            #if _dbg_point_avg
            else if (!strcmp(words[0], "DBGBWCEN"))
            {
                Point bw50_CA[10];
                for (l=1; l<=7; l++)
                {
                    int resno = working->get_bw50(l);
                    if (resno < 1) raise_error((std::string)"Strand does not have a " + to_string(l) + (std::string)".50 residue.");
                    bw50_CA[l] = working->get_residue(resno)->get_CA_location();
                    bw50_CA[l].weight = 1;
                    
                    AminoAcid* aa = working->get_residue(resno);
                    Atom* a = aa->get_atom("CA");
                    cout << "Residue " << resno << " CA is located at " << a->get_location()
                        << "; function returned " << bw50_CA[l] << endl;
                }
                Point old_center = average_of_points(&bw50_CA[1], 7);
                cout << old_center << endl;
            }
            #endif

            else if (!strcmp(words[0], "BWCOPY"))
            {
                if (!words[1] || !words[2]) raise_error("Insufficient parameters given for BWCOPY.");

                char c = words[1][0];
                if (c < 'A' || c > 'Z') raise_error("Invalid source strand given for BWCOPY.");
                Protein* source = strands[c-'A'];

                c = words[2][0];
                if (c < 'A' || c > 'Z') raise_error("Invalid destination strand given for BWCOPY.");
                Protein* dest = strands[c-'A'];

                for (l=1; l<=8; l++)
                {
                    std::string rgname = (std::string)"TMR" + std::to_string(l);
                    int sr = source->get_region_start(rgname), er = source->get_region_end(rgname), bw50 = source->get_bw50(l);
                    if (sr > 0 && er > 0) dest->set_region(rgname, sr, er);
                    if (bw50 > 0) dest->set_bw50(l, bw50);
                }
                for (l=1; l<=3; l++)
                {
                    std::string rgname = (std::string)"CYT" + std::to_string(l);
                    int bwhx = 10 + (l-1)*20 + (l*2);
                    int sr = source->get_region_start(rgname), er = source->get_region_end(rgname), bw50 = source->get_bw50(bwhx);
                    if (sr > 0 && er > 0) dest->set_region(rgname, sr, er);
                    if (bw50 > 0) dest->set_bw50(bwhx, bw50);

                    rgname = (std::string)"EXR" + std::to_string(l);
                    bwhx += 11;
                    sr = source->get_region_start(rgname), er = source->get_region_end(rgname), bw50 = source->get_bw50(bwhx);
                    if (sr > 0 && er > 0) dest->set_region(rgname, sr, er);
                    if (bw50 > 0) dest->set_bw50(bwhx, bw50);
                }
            }   // BWCOPY

            else if (!strcmp(words[0], "BWMOTIF"))
            {
                if (!words[1] || !words[2]) raise_error("Insufficient parameters given for BWMOTIF.");

                l = atoi(words[1]);
                char rgn[8];
                if (l < 8) sprintf(rgn, "TMR%d", l);
                else if (l >= 12 && l <= 78)
                {
                    int preced = floor(0.1*l);
                    int subseq = l % 10;
                    l = subseq >> 1;

                    if (preced & 1)
                    {
                        // Cytoplasmic.
                        sprintf(rgn, "CYT%d", l);
                    }
                    else
                    {
                        // Extracellular.
                        sprintf(rgn, "EXR%d", l);
                    }
                    l = atoi(words[1]);
                }
                else sprintf(rgn, "HXR%d", l);
                int sr = working->get_region_start(rgn), er = working->get_region_end(rgn);
                
                if (!sr || !er) raise_error((std::string)"Region " + std::to_string(atoi(words[1]))
                    + (std::string)"/" + (std::string)rgn
                    + (std::string)" not found in strand.");

                n = -1;
                for (i=0; words[2][i]; i++)
                {
                    if (words[2][i] >= 'A' && words[2][i] <= 'Z') n = i;
                    else words[2][i] &= 0x5f;
                }

                if (n < 0) raise_error("Motif must indicate BW50 residue by capitalization.");

                k = working->search_sequence(sr, er, words[2], 4, nullptr);
                if (!k) raise_error("Motif not found in region of strand.");
                // cout << l << ".50 = " << k << endl;

                working->set_bw50(l, k+n);
            }   // BWMOTIF

            else if (!strcmp(words[0], "CANMOVE"))
            {
                int sr1, er1, sr2, er2;

                if (!words[0] || !words[1] || !words[2] || !words[3] || !words[4] || !words[5] || !words[6])
                    raise_error("Insufficient parameters given for CANMOVE.");
                else if (words[7])
                    raise_error("Too many parameters given for CANMOVE.");

                sr1 = interpret_single_int(words[1]);
                er1 = interpret_single_int(words[2]);

                sr2 = interpret_single_int(words[4]);
                er2 = interpret_single_int(words[5]);

                if (!strcmp(words[3], "TO"))
                {
                    float f = 0, clash_limit = 0;
                    Point center1 = working->get_region_center(sr1, er1);
                    Point center2 = working->get_region_center(sr2, er2);
                    SCoord df = center2.subtract(center1);
                    float dfr = df.r - unconnected_residue_mindist;
                    df.r = 0.5;

                    Pose poses[er1+1];
                    for (l=sr1; l<=er1; l++)
                    {
                        AminoAcid* aa = working->get_residue(l);
                        if (aa)
                        {
                            poses[l].copy_state(aa);
                        }
                    }

                    clash_limit = clash_limit_per_aa;

                    for (l=0; l<200; l++)
                    {
                        working->move_piece(sr1, er1, df);
                        working->get_internal_clashes(sr1, er1, true);
                        float c = 0, c1;

                        int resno1, resno2;
                        AminoAcid *aa1, *aa2;
                        for (resno1 = sr1; resno1 <= er1; resno1++)
                        {
                            aa1 = working->get_residue(resno1);
                            if (!aa1) continue;
                            for (resno2 = sr2; resno2 <= er2; resno2++)
                            {
                                aa2 = working->get_residue(resno2);
                                c1 = aa1->get_intermol_clashes(aa2);
                                if (c1 > c) c = c1;
                            }
                        }
                        // cout << (f+df.r) << " " << c << endl;

                        if (c > clash_limit)
                        {
                            df.r *= -1;
                            working->move_piece(sr1, er1, df);
                            df.r *= -0.666;
                        }
                        else
                        {
                            f += df.r;
                        }

                        if (fabs(df.r) < 0.001) break;
                    }

                    for (l=sr1; l<=er1; l++)
                    {
                        AminoAcid* aa = working->get_residue(l);
                        if (aa) poses[l].restore_state(aa);
                    }

                    if (f > dfr) f = dfr;

                    Star s;
                    s.f = f;
                    set_variable(words[6], s);
                }
                else
                {
                    raise_error((std::string)"Unimplemented direction " + (std::string)words[3] + (std::string)".");
                }

            }   // CANMOVE

            else if (!strcmp(words[0], "CENTER"))
            {
                l = 1;
                Point newcen(0,0,0);
                if (words[l]) newcen = interpret_single_point(words[l++]);
                if (words[l]) raise_error("Too many parameters given for CENTER.");

                if (lig_chain == g_chain && ligand.get_atom_count())
                {
                    ligand.move(newcen.subtract(working->get_region_center(1, working->get_end_resno())));
                }

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

                if (ct > sr)
                {
                    // Line up the C and O of er to the expected prevaa C and O of ct.
                    a1 = era->get_atom("C");
                    a2 = era->get_atom("O");
                    cta->predict_previous_COCA(pt3);
                }
                else
                {
                    // Line up the N and HN (or substitute) of er to the expected nextaa N and HN of ct.
                    a1 = era->get_atom("N");
                    a2 = era->HN_or_substitute();
                     cta->predict_next_NHCA(pt3);
                }
                a3 = era->get_atom("CA");

                working->conform_backbone(sr, er, a1, pt3[0], a2, pt3[1], iters);
                working->backconnect(sr, er);

            _no_connect:
                goto _pc_continue;
            }	// CONNECT

            #define do_ctnrg_iters 1
            else if (!strcmp(words[0], "CTNRG"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for CTNRG.");
                for (n=1; words[n]; n++);       // count.
                n--;

                int asr = 1, aer = 9999, bsr = 1, ber = 9999;
                char *sa, *sb;
                SCoord dir(0,0,0);
                char* outvar = nullptr;
                float baseline;

                switch (n)
                {
                    case 1:
                    case 2:
                    raise_error("Insufficient parameters given for CTNRG.");

                    case 3:
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                    case 8:

                    l = 1;                  
                    sa = interpret_single_string(words[l++]);

                    if (n >= 7)
                    {
                        asr = interpret_single_int(words[l++]);
                        aer = interpret_single_int(words[l++]);
                    }

                    sb = interpret_single_string(words[l++]);

                    if (n > 4)
                    {
                        bsr = interpret_single_int(words[l++]);
                        ber = interpret_single_int(words[l++]);
                    }

                    outvar = words[l++];
                    if (words[l])
                    {
                        dir = interpret_single_point(words[l++]);
                    }
                    break;

                    default:
                    raise_error("Too many parameters given for CTNRG.");
                }

                strands[sa[0]-65]->set_pdb_chain(sa[0]);
                strands[sb[0]-65]->set_pdb_chain(sb[0]);
                sa[0] -= 65;
                sb[0] -= 65;

                if (dir.r)
                {
                    float former, latter;

                    baseline = contact_energy(strands[sa[0]], strands[sb[0]], true, 9999, 9999, bsr, ber);
                    baseline += contact_energy(strands[sa[0]], strands[sb[0]], true, asr, aer, 9999, 9999);
                    former = contact_energy(strands[sa[0]], strands[sb[0]], true, asr, aer, bsr, ber) - baseline;
                    dir.r = 0.1;

                    #if debug_contact_energy
                    cout << former << endl << flush;
                    #endif

                    for (i=0; i<200; i++)
                    {
                        strands[sb[0]]->move_piece(bsr, ber, dir);
                        #if !do_ctnrg_iters
                        break;
                        #else
                        latter = contact_energy(strands[sa[0]], strands[sb[0]], false, asr, aer, bsr, ber) - baseline;
                        #if _dbg_homology
                        cout << latter << " " << dir.r << endl;
                        #endif

                        if (latter > former)
                        {
                            strands[sb[0]]->undo();
                            dir.r *= 0.5;
                        }
                        else
                        {
                            former = latter;
                        }

                        if (latter < 100 || fabs(dir.r) < 0.01)
                        {
                            #if _dbg_homology
                            cout << "Got it in " << i << endl;
                            #endif
                            break;
                        }
                        #endif
                    }
                }

                Star s;
                s.f = contact_energy(strands[sa[0]], strands[sb[0]], true, asr, aer, bsr, ber) - baseline;

                set_variable(outvar, s);
            }   // CTNRG

            else if (!strcmp(words[0], "MOC"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for MOC.");
                int resno1 = interpret_single_int(words[1]);
                if (!words[2]) raise_error("Insufficient parameters given for MOC.");
                int resno2 = interpret_single_int(words[2]);
                if (!words[3]) raise_error("Insufficient parameters given for MOC.");
                char* outvar = words[3];

                Star s;
                s.n = 0;

                AminoAcid *aa1 = working->get_residue(resno1), *aa2 = working->get_residue(resno2);
                if (!aa1 || !aa2)
                {
                    set_variable(outvar, s);
                    continue;
                }
                SCoord optimize = aa1->motion_to_optimal_contact(aa2);

                /*Atom *atom1 = aa1->get_nearest_atom(aa2->get_CA_location()), *atom2 = aa2->get_nearest_atom(aa1->get_CA_location());
                atom1 = aa1->get_nearest_atom(atom2->get_location());
                atom2 = aa2->get_nearest_atom(atom1->get_location());
                float nearest_atom_distance = atom1->distance_to(atom2);
                // cout << *atom1 << " is " << nearest_atom_distance << "A from " << *atom2 << endl;

                if (optimize.r < fmax(nearest_atom_distance - 2.5, 0))
                {
                    if (!optimize.r) optimize = atom2->get_location().subtract(atom1->get_location());
                    else optimize.r = fmax(nearest_atom_distance - 2.5, 0);
                }*/

                switch(outvar[0])
                {
                    case '@':
                    s.ppt = new Point(optimize);
                    break;

                    case '&':
                    s.f = optimize.r;
                    break;

                    default:
                    raise_error("Ouput variable must be either Cartesian or float.");
                }
                set_variable(outvar, s);
            }   // MOC

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
                if (!words[3]) raise_error("Insufficient parameters for DOWNLOAD.");

                std::string url = "", retfmt = "", destfn = "";
                bool skip_if_exists = false;

                if (!strcmp(words[3], "ONCE")) raise_error("Insufficient parameters for DOWNLOAD.");
                else destfn = interpret_single_string(words[3]);

                if (words[4])
                {
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

                        delete dls;
                    }

                    if (!url.length()) raise_error("Download source not found in data file.");
                    j = url.find('%');
                    url.replace(j, 1, interpret_single_string(words[2]));

                    if (retfmt == "PDBDATA")
                    {
                        if (!download_file(url, destfn)) raise_error("Download failed.");
                    }
                    else if (retfmt.substr(0, 5) == "JSON:")
                    {
                        std::string key = (std::string)"\"" + retfmt.substr(5) + (std::string)"\":";
                        n = key.length();
                        if (!download_file(url, "tmp/.pepdl")) raise_error("JSON retrieval failed.");
                        FILE* fp = fopen("tmp/.pepdl", "rb");
                        if (!fp) raise_error("JSON retrieval failed.");
                        url = "";
                        while (!feof(fp))
                        {
                            fgets(buffer1, 1022, fp);
                            char* psz = strstr(buffer1, key.c_str());
                            if (psz)
                            {
                                psz += n+1;
                                char* quot = strchr(psz, '"');
                                if (quot)
                                {
                                    *quot = 0;
                                    url = psz;
                                    break;
                                }
                            }
                        }
                        fclose(fp);
                        remove("tmp/.pepdl");

                        if (!url.length()) raise_error("URL not found in JSON data.");
                        if (!download_file(url, destfn)) raise_error("Download failed.");
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

            else if (!strcmp(words[0], "GOSUB"))
            {
                call_stack.push_back(program_counter);
                stack_pointer = call_stack.size();
                goto _goto;
            }   // GOSUB

            else if (!strcmp(words[0], "GOTO"))
            {
                _goto:
                if (!words[1]) raise_error("Insufficient parameters given for GOTO.");
                sprintf(buffer1, "%s:", words[1]);
                for (n=0; n<script_lines.size(); n++)
                {
                    psz = new char[256];
                    strcpy(psz, script_lines[n].c_str());
                    if (!strcmp(psz, buffer1))
                    {
                        delete psz;
                        program_counter = n+1;
                        goto _found_goto_target;
                    }
                    delete psz;
                }
                raise_error( (std::string)"Label not found: \"" + (std::string)buffer1 + (std::string)"\"");
                return 0x51974c5;

            _found_goto_target:
                continue;
            }	// GOTO

            else if (!strcmp(words[0], "RETURN"))
            {
                if (!stack_pointer) raise_error("RETURN without GOSUB.");
                stack_pointer--;
                program_counter = call_stack[stack_pointer];
                call_stack.erase(call_stack.begin()+stack_pointer);
            }   // RETURN

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
                Protein *tpl, *ref;
                FILE* fp;

                bool delete_tpl = false, delete_ref = false;
                
                words[1] = interpret_single_string(words[1]);
                if (strlen(words[1]) == 1)
                {
                    tpl = strands[words[1][0] - 65];
                }
                else
                {
                    fp = fopen(words[1], "rb");
                    if (!fp) raise_error((std::string)"File not found: " + (std::string)words[1]);

                    tpl = new Protein("template");
                    delete_tpl = true;
                    tpl->load_pdb(fp);
                    fclose(fp);
                }

                if (words[2])
                {
                    words[2] = interpret_single_string(words[2]);
                    if (strlen(words[2]) == 1)
                    {
                        ref = strands[words[2][0] - 65];
                    }
                    else
                    {
                        fp = fopen(words[2], "rb");
                        if (!fp) raise_error((std::string)"File not found: " + (std::string)words[2]);

                        ref = new Protein("reference");
                        delete_ref = true;
                        ref->load_pdb(fp);
                        fclose(fp);
                    }
                    
                    working->homology_conform(tpl, ref);
                }
                else working->homology_conform(tpl, working);

                if (delete_tpl) delete tpl;
                if (delete_ref) delete ref;
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
                // if (lig_chain && ligand.get_atom_count()) ligand.hydrogenate();
            }   // HYDRO

            else if (!strcmp(words[0], "IF"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for IF.");
                if (!words[2]) raise_error("Insufficient parameters given for IF.");

                bool inverse_result = false;
                l = 2;
                if (!strcmp(words[1], "NOT"))
                {
                    inverse_result = true;
                    l++;
                }

                if (!strcmp(words[l-1], "EXISTS"))
                {
                    char* file_name = interpret_single_string(words[l]);
                    l--;
                    if (file_exists(file_name)) goto _evaluated_true;
                    else goto _evaluated_false;
                }

                // TODO: IF %var GOTO blah.
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
                if (inverse_result) goto _false_result;
            _true_result:
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
                if (inverse_result) goto _true_result;
            _false_result:
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

            else if (!strcmp(words[0], "INTC"))
            {
                l = 1;
                n = -1;
                if (!words[l])
                {
                    f = working->get_internal_clashes(1, working->get_end_resno());
                }
                else if (words[l][0] >= 'A' && words[l][0] <= 'Z')
                {
                    if (!strands[words[l][0]-'A']) raise_error("No strand with specified chain letter.");
                    else f = strands[words[l][0]-'A']->get_internal_clashes(1, strands[words[l][0]-'A']->get_end_resno());
                    l++;
                }
                else f = working->get_internal_clashes(1, working->get_end_resno());

                if (words[l] && (words[l][0] == '&'))
                {
                    n = find_var_index(words[l++]);
                    if (n<0)
                    {
                        n = vars++;
                        script_var[n].name = words[l-1];
                        script_var[n].vt = type_from_name(words[l-1]);
                    }
                }

                if (n < 0)
                {
                    bool b = true;
                    cout << f;
                    if (words[l])
                    {
                        if (words[l][0] == '~')
                        {
                            b = false;
                        }
                        else raise_error("Unknown argument.");
                    }
                    
                    if (b) cout << endl;
                }
                else script_var[n].value.f = f;
            }	// INTC

            else if (!strcmp(words[0], "WORST"))
            {
                l = 1;
                j = -1;
                if (words[l] && words[l][0] >= 'A' && words[l][0] <= 'Z')
                {
                    j = words[l++][0] - 'A';
                    if (!strands[j]) raise_error("Empty strand.");
                }

                if (words[l] && words[l][0] != '~')
                {
                    n = find_var_index(words[l]);
                    if (n<0)
                    {
                        n = vars++;
                        script_var[n].name = words[l];
                        script_var[n].vt = type_from_name(words[l]);
                    }

                    if (words[l][0] == '$')
                    {
                        script_var[n].value.psz = new char[65536];
                        strcpy(script_var[n].value.psz, (j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop1->get_name() : "null");
                        l++;
                    }
                    else if (words[l][0] == '%')
                    {
                        script_var[n].value.n = (j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop1->get_residue_no() : 0;
                        l++;
                    }
                }
                else cout << ((j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop1->get_name() : "null") << " ";

                if (words[l] && words[l][0] != '~')
                {
                    n = find_var_index(words[l]);
                    if (n<0)
                    {
                        n = vars++;
                        script_var[n].name = words[l];
                        script_var[n].vt = type_from_name(words[l]);
                    }

                    if (words[l][0] == '$')
                    {
                        script_var[n].value.psz = new char[65536];
                        strcpy(script_var[n].value.psz, (j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop2->get_name() : "null");
                        l++;
                    }
                    else if (words[l][0] == '%')
                    {
                        script_var[n].value.n = (j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop2->get_residue_no() : 0;
                        l++;
                    }
                }
                else
                {
                    bool b = false;
                    if (words[l] && words[l][0] == '~') b = true;
                    cout << ((j<0 ? working : strands[j])->stop1 ? (j<0 ? working : strands[j])->stop2->get_name() : "null");
                    if (!b) cout << endl;
                }
            }	// WORST

            else if (!strcmp(words[0], "MINC"))
            {
                l = 1;
                if (words[l]) raise_error("Too many parameters given for MINC.");
                working->get_internal_clashes(1, working->get_end_resno(), true, 50);
            }	// MINC

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
                        delete psz;
                        raise_error( (std::string)"Unimplemented operator " + (std::string)words[2] + (std::string)" for string assignment.");
                        return 0x51974c5;		// If you use your imagination, that says "syntax".
                    }
                    psz -= m;
                    delete psz;
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
                        f = interpret_single_float(words[5+l]);
                        if (!strcmp(words[4+l], "+")) *script_var[n].value.ppt = script_var[n].value.ppt->add(pt);
                        else if (!strcmp(words[4+l], "-")) *script_var[n].value.ppt = script_var[n].value.ppt->subtract(pt);
                        else if (!strcmp(words[4+l], "*")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() * f); // pt.magnitude());
                        else if (!strcmp(words[4+l], "/")) script_var[n].value.ppt->scale(script_var[n].value.ppt->magnitude() / f); // pt.magnitude());
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
                include_ligand = false;
                
				if (words[2]) chain = words[2][0];
                if (words[2] && words[3])
                {
                    new_chain = words[3][0];
                    if (words[4])
                    {
                        if (words[5]) raise_error("Too many parameters given for LOAD.");
                        else if (strcmp(words[4], "+")) raise_error("Too many parameters given for LOAD.");
                        else
                        {
                            ligand.delete_all_atoms();
                            include_ligand = true;
                            lig_chain = new_chain;
                        }
                    }

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
                l = working->load_pdb(pf, n, chain);
                if (!l) raise_error("No residues loaded.");
                working->set_name_from_pdb_name(words[1]);
                if (include_ligand) ligand.from_pdb(pf, true);

                fclose(pf);

            _prot_deets:

				char buffer[1024];
				char buffer1[1024];
				Star sv;

                sprintf(buffer, "$PDB%c", g_chain);
				sv.psz = new char[strlen(psz)+4];
                strcpy(sv.psz, psz);
				set_variable(buffer, sv);

                const char* pname = working->get_name().c_str();
                sprintf(buffer, "$PROTEIN%c", g_chain);
                sv.psz = new char[strlen(pname)+4];
                strcpy(sv.psz, pname);
				set_variable(buffer, sv);

                delete[] psz;

                int seqlen = working->get_seq_length();
                sprintf(buffer, "%cSEQLEN%c", '%', g_chain);
                sv.n = seqlen;
				set_variable(buffer, sv);

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
                    // cout << "ln " << program_counter << " set " << g_chain << ":" << words[3] << " to " << atoi(words[4]) << "-" << atoi(words[5]) << endl;

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
                bool uses_fancy_options = false;
                Atom** ba = nullptr;

                _yes_I_used_goto_for_this:
                if (!strcmp(words[l], "YO"))
                {
                    force_tyrosine_O = true;
                    uses_fancy_options = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(words[l], "YAr"))
                {
                    force_tyrosine_O = false;
                    uses_fancy_options = true;
                    l++;
                    goto _yes_I_used_goto_for_this;
                }
                else if (!strcmp(words[l], "Th8"))
                {
                    thiolate = true;
                    uses_fancy_options = true;
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
                        uses_fancy_options = true;
                        l++;
                        goto _another_goto;
                    }
                    else if (!strcmp(words[l], "YAr"))
                    {
                        local_O = false;
                        uses_fancy_options = true;
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

                if (ncr < 2) raise_error("MCOORD requires at least 2 coordinating atoms.");
                if (ncr < 3 && uses_fancy_options) raise_error("MCOORD with YO, YAr, or Th8 requires at least 3 coordinating atoms.");
                if (uses_fancy_options) working->coordinate_metal(ma, ncr, resnos, cratoms);
                else
                {
                    std::vector<MCoord> mtlcoords;
                    MCoord mc;
                    mc.Z = ma->get_Z();
                    mc.charge = elem_charge;
                    mc.mtl = ma;

                    for (l=0; l<ncr; l++)
                    {
                        ResiduePlaceholder rp;
                        rp.resno = resnos[l];
                        mc.coordres.push_back(rp);
                    }

                    mtlcoords.push_back(mc);

                    working->coordinate_metal(mtlcoords);
                }
            } // MCOORD

            else if (!strcmp(words[0], "MEASURE"))
            {
                for (l=1; words[l]; l++);       // Count params.

                if (l != 4 && l != 6) raise_error("Wrong number of parameters given for MEASURE.");

                bool atom_names = (l >= 5);
                int res1 = interpret_single_int(words[1]);
                std::string atom1 = atom_names ? interpret_single_string(words[2]) : "CA";
                int res2 = interpret_single_int(words[atom_names ? 3 : 2]);
                std::string atom2 = atom_names ? interpret_single_string(words[4]) : "CA";

                if (!res1 || !working->get_residue(res1)) raise_error("Residue A not found in protein.");
                if (!res2 || !working->get_residue(res2)) raise_error("Residue B not found in protein.");

                if (!working->get_atom(res1, atom1.c_str())) raise_error(std::to_string(res1) + (std::string)":" + atom1 + (std::string)" not found in protein.");
                if (!working->get_atom(res2, atom2.c_str())) raise_error(std::to_string(res2) + (std::string)":" + atom2 + (std::string)" not found in protein.");

                float r = working->get_atom_location(res1, atom1.c_str()).get_3d_distance(working->get_atom_location(res2, atom2.c_str()));

                Star s;
                s.f = r;
                set_variable(words[l-1], s);
            } // MEASURE

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

                // TODO: Option to include ligand in motion.
            }	// MOVE

            else if (!strcmp(words[0], "MOVEREL"))
            {
                l = 1;
                SCoord movamt(0,0,0);
                int sr, er;
                if (words[l]) sr = interpret_single_int(words[l++]);
                else raise_error("Not enough parameters given for MOVEREL.");
                if (words[l]) er = interpret_single_int(words[l++]);
                else raise_error("Not enough parameters given for MOVEREL.");
                if (words[l]) movamt = interpret_single_point(words[l++]);
                if (words[l]) raise_error("Too many parameters given for MOVEREL.");
                working->move_piece(sr, er, movamt);

                // TODO: Option to include ligand in motion.
            }	// MOVEREL

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
                int sr = working->get_start_resno(), er = working->get_end_resno(), piv = 0;

				if (!words[1]) raise_error("Insufficient parameters given for ROTATE.");
                SCoord axis = interpret_single_point(words[1]);
				if (!words[2]) raise_error("Insufficient parameters given for ROTATE.");
                float theta = interpret_single_float(words[2]) * fiftyseventh;
                if (words[3])
                {
                    piv = interpret_single_int(words[3]);
                    if (words[4] && !words[5])
                    {
                        piv = 0;
                        sr = interpret_single_int(words[3]);
                        er = interpret_single_int(words[4]);
                    }
                    else if (words[4] && words[5])
                    {
                        sr = interpret_single_int(words[4]);
                        er = interpret_single_int(words[5]);
                        if (words[6]) raise_error("Too many parameters given for ROTATE.");
                    }
                }

                LocatedVector lv = axis;
                AminoAcid* aapiv = nullptr;
                if (piv) aapiv = working->get_residue(piv);
                if (aapiv) lv.origin = aapiv->get_CA_location();
                else lv.origin = working->get_region_center(sr, er);

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
                    if (!working->get_seq_length()) continue;
                    g_chain = l+65;
                    working->set_pdb_chain(l+65);
                    working->save_pdb(pf, ligand.get_atom_count() ? &ligand : nullptr);
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
                                    InteratomicForce* iff[32];
                                    InteratomicForce::fetch_applicable(a, b, iff);
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
                esr = er;

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

                delete psz;

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

            else if (!strcmp(words[0], "SIDEREPL"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for SIDEREPL.");
                if (words[2] && words[3]) raise_error("Too many parameters given for SIDEREPL.");

                pf = fopen(interpret_single_string(words[1]), "rb");
                if (!pf) raise_error((std::string)"File not found" + (std::string)words[1] + (std::string)".");
                if (words[2]) words[2] = interpret_single_string(words[2]);
                ptmp.load_pdb(pf, 0, words[2] ? words[2][0] : 'A');

                working->replace_side_chains_from_other_protein(&ptmp);

                n = working->get_end_resno();
                for (l=1; l<=n; l++)
                {
                    f = working->get_internal_clashes(l, l);
                    if (f > clash_limit_per_aa) working->minimize_residue_clashes(l);
                }
            }   // SIDEREPL

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

            else if (!strcmp(words[0], "STRLEN"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for STRLEN.");
                if (!words[2]) raise_error("Insufficient parameters given for STRLEN.");
                if (words[3]) raise_error("Too many parameters given for STRLEN.");

                char* c = interpret_single_string(words[1]);
                if (!c) l = 0;
                else l = strlen(c);

                n = find_var_index(words[2]);
                if (n<0)
                {
                    n = vars++;
                    script_var[n].name = words[2];
                    script_var[n].vt = type_from_name(words[2]);
                }
                int flags = n & _VARFLAGS_MASK;
                n &= _VARNUM_MASK;

                script_var[n].value.n = l;
            }

            else if (!strcmp(words[0], "STRPOS"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for STRPOS.");
                if (!words[2]) raise_error("Insufficient parameters given for STRPOS.");
                if (!words[3]) raise_error("Insufficient parameters given for STRPOS.");
                if (words[4]) raise_error("Too many parameters given for STRPOS.");

                n = find_var_index(words[3]);
                if (n<0)
                {
                    n = vars++;
                    script_var[n].name = words[3];
                    script_var[n].vt = type_from_name(words[3]);
                }
                int flags = n & _VARFLAGS_MASK;
                n &= _VARNUM_MASK;

                char* haystack = interpret_single_string(words[1]);
                char* needle = interpret_single_string(words[2]);
                char* c = strstr(haystack, needle);
                if (!c) script_var[n].value.n = 0;
                else script_var[n].value.n = c - haystack + 1;
            }   // STRPOS

            else if (!strcmp(words[0], "UNCHAIN"))
            {
                if (!words[1]) raise_error("Insufficient parameters given for UNCHAIN.");
                chain = words[1][0] - 65;
                if (nullptr != strands[chain])
                {
                    delete strands[chain];
                    strands[chain] = nullptr;
                }
                if (g_chain == 65+chain) working = strands[chain] = new Protein("TheProt");
            }   // UNCHAIN

            else if (!strcmp(words[0], "UNLIG"))
            {
                if (words[1]) raise_error("Too many parameters given for UNLIG.");
                ligand.delete_all_atoms();
                lig_chain = '\0';
            }   // UNLIG

            else if (!strcmp(words[0], "UPRIGHT"))
            {
                if (words[1] && words[2]) raise_error("Too many parameters given for UPRIGHT.");

                try
                {
                    working->upright();
                    // cout << (Point)working->last_uprighted_xform << endl;

                    for (i=0; i<26; i++)
                    {
                        if (words[1] && !strchr(words[1], 65+i)) continue;

                        if (strands[i] && (strands[i] != working))
                        {
                            strands[i]->move_piece(1, 9999, working->last_uprighted_xform);
                            strands[i]->rotate_piece(1, 9999, working->last_uprighted_A.origin, working->last_uprighted_A.v, working->last_uprighted_A.a);
                            strands[i]->rotate_piece(1, 9999, working->last_uprighted_B.origin, working->last_uprighted_B.v, working->last_uprighted_B.a);
                        }

                        if (lig_chain == 65+i && ligand.get_atom_count())
                        {
                            // cout << "Upright with ligand: " << ligand.get_barycenter() << " x " << (Point)working->last_uprighted_xform << endl;
                            ligand.move(working->last_uprighted_xform);
                            LocatedVector lv = working->last_uprighted_A.v;
                            lv.origin = working->last_uprighted_A.origin;
                            ligand.rotate(lv, working->last_uprighted_A.a);
                            lv = working->last_uprighted_B.v;
                            lv.origin = working->last_uprighted_B.origin;
                            ligand.rotate(lv, working->last_uprighted_B.a);
                        }
                    }
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










