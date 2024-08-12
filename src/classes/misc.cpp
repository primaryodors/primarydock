
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include "point.h"

using namespace std;

float _INTERA_R_CUTOFF = _DEFAULT_INTERA_R_CUTOFF;
const char* Greek = "ABGDEZHQIKLMNJOPRSTYFXCW";
std::ofstream *debug = nullptr;
bool last_iter = false;
float pre_ligand_multimol_radius = default_pre_ligand_multimol_radius;
float pre_ligand_flex_radius = default_pre_ligand_flex_radius;

char asterisk[5] = "*";

int in_array(void* needle, void** haystack)
{
    int i;

    for (i=0; haystack[i]; i++)
    {
        if (haystack[i] == needle) return i;
    }

    return -1;
}

int in_array(int needle, int* haystack)
{
    int i;

    for (i=0; haystack[i]; i++)
    {
        if (haystack[i] == needle) return i;
    }

    return -1;
}

int in_array(Star needle, Star* haystack)
{
    int i;

    for (i=0; haystack[i].n; i++)
    {
        if (haystack[i].n == needle.n) return i;
    }

    return -1;
}

char** chop_spaced_words(char* line, char separator)
{
    if (!line[0]) return nullptr;

    char** retval = new char*[100];
    int i, j=0;

    if (separator == ' ' && line[0] == '\t') line[0] = separator;
    if (line[0] != separator) retval[j++] = line;
    for (i=1; line[i] && (line[i] != '\n'); i++)
    {
        if (line[i] == '"')
        {
            retval[j++] = &line[i++];
            for (; line[i] && line[i] != '"'; i++);		// Get to next quote character.
            bool last = false;
            if (!line[i+1]) last = true;
            if (line[i] == '"') line[++i] = '\0';
            if (last) break;
            continue;
        }

        if (separator == ' ' && line[i] == '\t') line[i] = separator;
        if (line[i-1] == separator && line[i] != separator) retval[j++] = line+i;
        else if (line[i-1] == 0 && line[i] != separator) retval[j++] = line+i;
        else if (line[i-1] != separator && line[i] == separator) line[i] = 0;
        if (line[i] == '\n') line[i] = 0;
        if (j >= 99) break;
    }
    if (line[i] == '\n') line[i] = 0;
    retval[j] = 0;

    return retval;
}


// http://www.zedwood.com/article/cpp-str_pad
// Code by zedwood.com, CC-BY-SA 3.0.
string str_pad(const string &str, int pad_length, string pad_string, STR_PAD pad_type)
{
    int i,j,x;
    int str_size = str.size();
    int pad_size = pad_string.size();
    if (pad_length<=str_size || pad_size<1)
        return str;

    string o;
    o.reserve(pad_length);//allocate enough memory only once

    if (pad_type==STR_PAD_RIGHT)
    {
        for(i=0,x=str_size; i<x; i++)
            o.push_back( str[i] );
        for(i=str_size; i<pad_length; )
            for(j=0; j<pad_size && i<pad_length; j++,i++)
                o.push_back( pad_string[j] );
    }
    else if (pad_type==STR_PAD_LEFT)
    {
        int a1= pad_length-str_size;
        for(i=0; i<a1; )
            for(j=0; j<pad_size && i<a1; j++,i++)
                o.push_back( pad_string[j] );
        for(i=0,x=str_size; i<x; i++)
            o.push_back( str[i] );
    }
    else if (pad_type==STR_PAD_BOTH)
    {
        int a1= (pad_length-str_size)/2;
        int a2= pad_length-str_size-a1;
        for(i=0; i<a1; )
            for(j=0; j<pad_size && i<a1; j++,i++)
                o.push_back( pad_string[j] );
        for(i=0,x=str_size; i<x; i++)
            o.push_back( str[i] );
        for(i=0; i<a2; )
            for(j=0; j<pad_size && i<a2; j++,i++)
                o.push_back( pad_string[j] );
    }
    return o;
}

// End zedwood.com code.

Star* array_unique(Star* inp)
{
    int i, k=0;

    Star buffer[16384];
    for (i=0; inp[i].n; i++)
    {
        if (!in_array(inp[i], buffer))
            buffer[k++].n = inp[i].n;
    }
    buffer[k].n = 0;

    Star* retval = new Star[k+4];
    for (i=0; i<k; i++) retval[i] = buffer[i];
    retval[k].n = 0;

    return retval;
}

int greek_from_aname(const char* aname)
{
    int i, j;
    if (!aname) return -1;
    for (i=0; aname[i]; i++);	// Get length.
    i--;
    for (; i>=0 && aname[i] < 'A'; i--);	// Walk back until find a letter.
    if (i<0) return i;
    for (j=0; Greek[j]; j++)
        if (Greek[j] == aname[i]) return j;
    return -1;
}

int randsgn()
{
    return (rand()&1) ? 1 : -1;
}

float frand(float lmin, float lmax)
{
    int r = rand();
    float f = (float)r / RAND_MAX;
    f *= (lmax-lmin);
    return f+lmin;
}

// http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
float Pearson_correlation(float* xarr, float* yarr, int length)
{
    if (!xarr || !yarr || !length) return 0;

    int i;
    float xavg=0, yavg=0;

    for (i=0; i<length; i++)
    {
        // cout << i << ".";
        xavg += xarr[i];
        yavg += yarr[i];
    }
    xavg /= length;
    yavg /= length;

    float ssxx = 0.0, ssyy = 0.0, ssxy = 0.0;

    for (i=0; i<length; i++)
    {
        float x = xarr[i];
        float y = yarr[i];

        ssxx += pow(x - xavg, 2);
        ssyy += pow(y - yavg, 2);
        ssxy += (x - xavg)*(y - yavg);
    }

    if (ssxx == 0 || ssyy == 0) return -1;		// insufficient data.

    float r = ssxy/pow(ssxx*ssyy,0.5);
    return r;
}

std::string cardinality_printable(float card)
{
    std::string retval = ".";
    if (card > 1 && card < 2) retval = ":";
    else switch ((int)card)
        {
        case 0:
            break;

        case 1:
            retval = "-";
            break;

        case 2:
            retval = "=";
            break;

        case 3:
            retval = "#";
            break;

        case 4:
            retval = "$";
            break;

        default:
            retval = "!?";
        }

    return retval;
}

const float pH_minus_pKa[22] = {-1.9, -0.9, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.6, 2.1, 2.6, 3.1};
const float f_protonated[22] = {0.988, 0.888, 0.715, 0.666, 0.613, 0.557, 0.5, 0.443, 0.387, 0.334, 0.285, 0.240, 0.201, 0.166, 0.137, 0.112, 0.091, 0.074, 0.024, 0.008, 0.002, 0.001};
float protonation(float pKa)
{
    float d = pH - pKa;
    int i;

    if (d < pH_minus_pKa[0]) return 1;
    else if (d > pH_minus_pKa[21]) return 0;

    float f;
    for (i=0; i<21; i++)
    {
        f = (d - pH_minus_pKa[i]) / (pH_minus_pKa[i+1] - pH_minus_pKa[i]);
        if (f <= 1)
        {
            float e = 1.0 - f;
            return e*f_protonated[i] + f*f_protonated[i+1];
        }
    }

    throw 0xbadbca;
}

float larger(float v1, float v2)
{
    float v1a = fabs(v1), v2a = fabs(v2);
    int sign = (v2a > v1a) ? sgn(v2) : sgn(v1);
    return fmax(v1a, v2a) * sign;
}

bool file_exists(std::string fname)
{
    struct stat s;
    if (stat(fname.c_str(), &s) == 0) return true;
    else return false;
}


void colorrgb(int r, int g, int b)
{
    r = max(0, min(255, r));
    g = max(0, min(255, g));
    b = max(0, min(255, b));

    cout << "\x1b[38;2;" << r << ";" << g << ";" << b << "m";
}

void colorize(float f)
{
    float red, green, blue;

    if (f >= 0)
    {
        f = sqrt(f/5);
        blue = 128 + 128 * f;
        green = fmax(48, (f-1) * 255);
        red = fmax(64, (f-2) * 255);
    }
    else
    {
        f = sqrt(-f)*3;
        f = fmax(0,fmin(128,f*16));
        red = 128+f;
        blue = 128-f;
        green = 0.333 * red + 0.666 * blue;
    }

    colorrgb(red, green, blue);
}

void colorless()
{
    cout << "\x1b[0m";
}











