
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <vector>

#include "point.h"

using namespace std;

float _INTERA_R_CUTOFF = _DEFAULT_INTERA_R_CUTOFF;
const char* Greek = "ABGDEZHQIKLMNJOPRSTYFXCW";
std::ofstream *debug = nullptr;
bool last_iter = false;
bool differential_dock = false;
float pre_ligand_multimol_radius = default_pre_ligand_multimol_radius;
float pre_ligand_flex_radius = default_pre_ligand_flex_radius;

char asterisk[5] = "*";

#if active_persistence
int active_persistence_resno[active_persistence_limit];
#endif

#if active_persistence_noflex
bool allow_ligand_flex = true;
#endif

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

int in_array(int needle, std::vector<int> haystack)
{
    int i;

    for (i=0; i < haystack.size(); i++)
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

char** chop_spaced_fields(char* line, char separator)
{
    char** retval = new char*[100];
    int i, j=0;

    if (!line[0]) return nullptr;

    if (separator == ' ' && line[0] == '\t') line[0] = separator;
    if (line[0] != separator) retval[j++] = line;
    for (i=1; line[i] && (line[i] != '\n'); i++)
    {
        if (line[i] == '"')
        {
            retval[j++] = &line[i++];
            for (; line[i] && line[i] != '"'; i++);		// Get to next quote character.
            if (line[i] == '"') line[++i] = '\0';
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

#if active_persistence
float residue_binding_multiplier(int resno)
{
    int i;
    for (i=0; i<active_persistence_limit; i++)
    {
        if (!active_persistence_resno[i]) return 1;
        if (active_persistence_resno[i] == resno) return active_persistence_ratio;
    }

    return 1;
}
#endif


// Modified from here: https://gist.github.com/TheRayTracer/2644387
size_t levenshtein_distance(const char* s, size_t n, const char* t, size_t m)
{
   ++n; ++m;
   size_t* d = new size_t[n * m];

   memset(d, 0, sizeof(size_t) * n * m);

   for (size_t i = 1, im = 0; i < m; ++i, ++im)
   {
      for (size_t j = 1, jn = 0; j < n; ++j, ++jn)
      {
         if (s[jn] == t[im])
         {
            d[(i * n) + j] = d[((i - 1) * n) + (j - 1)];
         }
         else
         {
            d[(i * n) + j] =
                min(d[(i - 1) * n + j] + 1, /* A deletion. */
                min(d[i * n + (j - 1)] + 1, /* An insertion. */
                d[(i - 1) * n + (j - 1)] + 1)); /* A substitution. */
         }
      }
   }

#ifdef debug_levenshtein
   for (size_t i = 0; i < m; ++i)
   {
      for (size_t j = 0; j < n; ++j)
      {
         cout << d[(i * n) + j] << " ";
      }
      cout << endl;
   }
#endif

   size_t r = d[n * m - 1];
   delete [] d;
   return r;
}

size_t levenshtein_distance(std::string s1, std::string s2)
{
    return levenshtein_distance(s1.c_str(), s1.length(), s2.c_str(), s2.length());
}









