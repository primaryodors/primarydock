
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <strings.h>
#include "protein.h"

#define _DBG_REACHLIG true

using namespace std;

float *g_rgnxform_r = nullptr, *g_rgnxform_theta = nullptr, *g_rgnxform_y = nullptr;
float *g_rgnrot_alpha = nullptr, *g_rgnrot_w = nullptr, *g_rgnrot_u = nullptr;

Protein::Protein()
{
    aaptrmin.n = aaptrmax.n = 0;

    residues = nullptr;
    sequence = nullptr;
    ca = nullptr;
    res_reach = nullptr;
    metals = nullptr;

    int i;
    for (i=0; i<79; i++) Ballesteros_Weinstein[i] = 0;
}

Protein::Protein(const char* lname)
{
    name = lname;
    aaptrmin.n = aaptrmax.n = 0;

    residues = nullptr;
    sequence = nullptr;
    ca = nullptr;
    res_reach = nullptr;
    metals = nullptr;

    int i;
    for (i=0; i<79; i++) Ballesteros_Weinstein[i] = 0;
}

Protein::~Protein()
{
    connections.clear();

    delete[] remarks;
    remarksz = 0;
}

bool Protein::add_residue(const int resno, const char aaletter)
{
    int i;

    if (!residues)
    {
        arrlimit = resno+1024;
        residues = new AminoAcid*[arrlimit];
        sequence = new char[arrlimit];
        ca = new Atom*[arrlimit];
        res_reach = new float[arrlimit];

        for (i=0; i<arrlimit; i++)
        {
            residues[i] = NULL;
            sequence[i] = 0;
            ca[i] = NULL;
            res_reach[i] = 0;
        }
    }
    else if (resno >= arrlimit)
    {
        AminoAcid** oldres = residues;
        char* oldseq = sequence;
        Atom** oldca = ca;
        float* oldreach = res_reach;
        arrlimit = resno+1024;
        residues = new AminoAcid*[arrlimit];
        sequence = new char[arrlimit];
        ca = new Atom*[arrlimit];
        res_reach = new float[arrlimit];

        for (i=0; i<arrlimit; i++)
        {
            residues[i] = NULL;
            sequence[i] = 0;
            ca[i] = NULL;
            res_reach[i] = 0;
        }

        for (i=0; oldres[i]; i++)
        {
            residues[i] = oldres[i];
            sequence[i] = oldseq[i];
            ca[i] = oldca[i];
            res_reach[i] = oldreach[i];
        }
        residues[i] = 0;
        sequence[i] = 0;
    }

    i = resno-1;

    if (i)
    {
        Point pts[4];
        residues[i-1]->predict_next_NHCA(pts);
        residues[i] = new AminoAcid(aaletter, residues[i-1], false);
        residues[i]->ensure_pi_atoms_coplanar();

        residues[i]->establish_internal_clash_baseline();

        Atom* a = residues[i]->get_atom("N");
        if (a)
        {
            float r = pts[0].get_3d_distance(a->get_location());
            if (fabs(r) >= 0.01) cout << "Warning: Residue " << resno << " N location is off by " << r << "A." << endl;

            Atom* b = residues[i]->get_atom("HN");
            if (!b) b = residues[i]->get_atom("H");
            if (b)
            {
                Point HN = b->get_location().subtract(a->get_location());
                HN.scale(1);
                Point pt = pts[1].subtract(pts[0]);
                pt.scale(1);
                r = pt.get_3d_distance(HN);
                if (fabs(r) >= 0.01) cout << "Warning: Residue " << resno << " HN location is off by " << r << "A." << endl;
            }
        }
        else cout << "Warning: Residue " << resno << " has no N atom." << endl << flush;
    }
    else
    {
        residues[i] = new AminoAcid(aaletter, 0, false);
        residues[i]->ensure_pi_atoms_coplanar();
    }

    sequence[i] = aaletter;
    ca[i] = residues[i]->get_atom("CA");
    res_reach[i] = residues[i]->get_aa_definition()->reach;
    residues[i+1] = 0;
    sequence[i+1] = 0;

    if (!aaptrmin.n || residues[i] < aaptrmin.paa) aaptrmin.paa = residues[i];
    if (!aaptrmax.n || residues[i] > aaptrmax.paa) aaptrmax.paa = residues[i];

    return true;
}

void BallesterosWeinstein::from_string(const char* inpstr)
{
    helix_no = atoi(inpstr);
    const char* dot = strchr(inpstr, '.');
    if (dot) member_no = atoi(dot+1);
}

AminoAcid* Protein::get_residue(int resno)
{
    if (!resno) return 0;
    if (!residues) return 0;

    int i;
    for (i=0; residues[i]; i++)
    {
        if (residues[i]->get_residue_no() == resno) return residues[i];
    }

    return NULL;
}

AminoAcid* Protein::get_residue(BallesterosWeinstein bw)
{
    return get_residue_bw(bw.helix_no, bw.member_no);
}

AminoAcid* Protein::get_residue_bw(const char* bwno)
{
    char buffer[16];
    strcpy(buffer, bwno);
    char* dot = strchr(buffer, '.');
    if (!dot) return get_residue(atoi(bwno));
    *dot = 0;
    return get_residue_bw(atoi(buffer), atoi(dot+1));
}

AminoAcid* Protein::get_residue_bw(int hxno, int bwno)
{
    if (!Ballesteros_Weinstein[3])
    {
        cout << "Call to get_residue_bw() on a protein with no BW numbers set." << endl;
        throw -1;
    }

    int bw50 = get_bw50(hxno);
    if (bw50 < 1)
    {
        cout << "BW number not found: " << hxno << ".50." << endl;
        throw -1;
    }

    int resno = bw50 - 50 + bwno;
    return get_residue(resno);
}

Atom* Protein::get_atom(int resno, const char* aname)
{
    AminoAcid* aa = get_residue(resno);

    if (!aa) return NULL;
    return aa->get_atom(aname);
}

Point Protein::get_atom_location(int resno, const char* aname)
{
    AminoAcid* aa = get_residue(resno);

    if (!aa)
    {
        Point pt;
        return pt;
    }
    return aa->get_atom_location(aname);
}

bool Protein::add_sequence(const char* lsequence)
{
    if (!lsequence) return false;

    int i;
    for (i=0; lsequence[i]; i++)
    {
        add_residue(i+1, lsequence[i]);
    }

    for (i=0; lsequence[i]; i++)
    {
        if (get_atom(i+1, "CB"))		// Don't check glycine.
        {
            float r = get_atom_location(i+1, "CA").get_3d_distance(get_atom_location(i+1, "CB"));
            if (fabs(r-1.54) > 0.5) cout << "Error: " << i+1 << lsequence[i] << ":CA-CB = " << r << endl;
        }
    }

    int seql = get_seq_length();
    Molecule* aas[seql+4];
    for (i=1; i<=seql; i++)
    {
        aas[i-1] = get_residue(i);
        get_residue(i)->ensure_pi_atoms_coplanar();
        aas[i-1]->clear_all_bond_caches();
        aas[i] = 0;
    }
    Molecule::conform_molecules(aas, 25);
    for (i=1; i<=seql; i++)
    {
        get_residue(i)->ensure_pi_atoms_coplanar();
    }

    set_clashables();
    allocate_undo_poses();

    return true;
}

char Protein::set_pdb_chain(char c)
{
    if (c >= 'A' && c <= 'Z')
    {
        pdbchain = c;
    }

    int i;
    for (i=0; residues[i]; i++)
    {
        residues[i]->set_pdb_chain(pdbchain);
    }

    return pdbchain;
}

void Protein::save_pdb(FILE* os, Molecule* lig)
{
    int i, offset=0;

    if (remarksz)
    {
        for (i=0; i<remarksz && remarks[i]; i++)
        {
            fprintf(os, "%s", remarks[i]);
        }
    }

    if (regions_from == rgn_manual && regions)
    {
        for (i=0; i<PROT_MAX_RGN; i++)
        {
            if (!regions[i].start || regions[i].end <= regions[i].start) break;
            fprintf(os, "REMARK 650 HELIX %s %d %d\n", regions[i].name.c_str(), regions[i].start, regions[i].end);
        }
    }

    if (!residues) return;
    for (i=0; residues[i]; i++)
    {
        residues[i]->set_pdb_chain(pdbchain);
        residues[i]->save_pdb(os, offset);
        offset += residues[i]->get_atom_count();
    }
    if (m_mcoord)
    {
        for (i=0; m_mcoord[i]; i++)
        {
            cout << "Saving " << m_mcoord[i]->metal->name << endl;
            m_mcoord[i]->metal->save_pdb_line(os, ++offset);
        }
    }

    if (lig)
    {
        int ac = lig->get_atom_count();
        for (i=0; i<ac; i++)
        {
            Atom* a = lig->get_atom(i);
            if (a) a->save_pdb_line(os, ++offset);
        }
    }

    last_saved_atom_number = offset;

    // Example CONECT syntax:
    // CONECT  487 1056
    if (connections.size())
    {
        for (i=0; i<connections.size(); i++)
        {
            if (!connections[i]->atom || !connections[i]->btom) continue;

            fprintf(os, "CONECT ");
            int a, b;
            a = connections[i]->atom->pdbidx;
            b = connections[i]->btom->pdbidx;

            if (a < 1000) fprintf(os, " ");
            if (a <  100) fprintf(os, " ");
            if (a <   10) fprintf(os, " ");
            fprintf(os, "%d ", a);

            if (b < 1000) fprintf(os, " ");
            if (b <  100) fprintf(os, " ");
            if (b <   10) fprintf(os, " ");
            fprintf(os, "%d \n", b);
        }
    }

    fprintf(os, "\nTER\n");
}

void Protein::end_pdb(FILE* os)
{
    fprintf(os, "END\n");
}

float Protein::get_internal_clashes(int sr, int er, bool repack, int repack_iters)
{
    if (!residues) return 0;
    if (repack) save_undo_state();
    int i, j, l, m;
    float result = 0, r;
    Point clashttl(0,0,0);
    float clash_worst = 0;
    for (i=0; residues[i]; i++)
    {
        int resno = residues[i]->get_residue_no();
        if (sr && resno < sr) continue;
        if (er && resno > er) continue;

        if (repack)
        {
            AminoAcid** laa = get_residues_can_clash(resno);
            if (laa)
            {
                int n;
                for (n=0; laa[n]; n++);         // Get count.
                Molecule* interactors[n+4];
                Molecule* backdrop[n+4];
                MovabilityType wasmov[n+4];

                l = m = 0;
                interactors[l++] = residues[i];
                #if _dbg_repack
                std::string dbgresstr;
                #endif
                for (j=0; j<n; j++)
                {
                    if (residues[i] == laa[j]) continue;
                    // if (residues[i]->get_intermol_clashes(laa[j]) >= 1)
                    if (fabs(residues[i]->get_intermol_binding(laa[j])) >= 1)
                    {
                        #if _dbg_repack
                        dbgresstr += (std::string)" " + (std::string)laa[j]->get_name();
                        #endif
                        interactors[l] = laa[j];
                        wasmov[l] = laa[j]->movability;
                        if (laa[j]->movability != MOV_FLXDESEL) laa[j]->movability = MOV_FLEXONLY;
                        l++;
                    }
                    else
                    {
                        Atom* laaa = laa[j]->get_nearest_atom(residues[i]->get_CA_location());
                        if (laaa->distance_to(residues[i]->get_atom("CA")) < residues[i]->get_reach())
                        {
                            #if _dbg_repack
                            dbgresstr += (std::string)" [" + (std::string)laa[j]->get_name() + (std::string)"]";
                            #endif
                            backdrop[m++] = laa[j];
                        }
                    }
                }

                interactors[l] = nullptr;
                backdrop[m] = nullptr;
                if (l > 1)
                {
                    #if _dbg_repack
                    cout << "Repacking " << residues[i]->get_name() << " with" << dbgresstr << "..." << endl;
                    #endif

                    Molecule::conform_molecules(interactors, backdrop, repack_iters);
                }

                for (l=0; interactors[l]; l++)
                    interactors[l]->movability = wasmov[l];
            }
        }

        for (j=i; residues[j]; j++)
        {
            int resno2 = residues[j]->get_residue_no();
            if (abs(resno - resno2) < 5) continue;
            if (j==i) result += residues[i]->get_internal_clashes();
            else
            {
                if (repack)
                {
                    r = residues[i]->get_CA_location().get_3d_distance(residues[j]->get_CA_location());
                    float rr = residues[i]->get_reach() + residues[j]->get_reach();
                    if (r > rr) continue;
                }

                Atom* ia = residues[i]->get_nearest_atom(residues[j]->get_CA_location());
                Atom* ja = residues[j]->get_nearest_atom(ia->get_location());
                r = ia->distance_to(ja);
                if (r > ia->get_vdW_radius() + ja->get_vdW_radius()) continue;

                if (abs(resno2-resno) <= 1) continue;
                if (resno2 >= sr && resno2 <= er) continue;

                float f = residues[i]->get_intermol_clashes(residues[j]);
                result += f;
                SCoord dpos = residues[j]->get_CA_location().subtract(residues[i]->get_CA_location());
                dpos.r = f;
                clashttl = clashttl.add(dpos);
                float limit = unconnected_residue_mindist + sqrt(residues[i]->get_reach()) + sqrt(residues[j]->get_reach());
                f = fmax(limit - dpos.r, 0);
                if (f > clash_worst)
                {
                    clash_worst = f;
                    stop1 = residues[i];
                    stop2 = residues[j];
                }
            }
        }
    }

    last_int_clash_dir = clashttl;
    last_int_clash_dir.r = clash_worst;

    return result;
}

float Protein::get_rel_int_clashes()
{
    return get_internal_clashes(0, 0, false) - initial_int_clashes;
}

float Protein::get_internal_binding()
{
    if (!residues) return 0;
    int i, j;
    float result = 0;
    for (i=0; residues[i]; i++)
    {
        for (j=i; residues[j]; j++)
        {
            /*if (j==i) result += residues[i]->get_internal_binding();
            else*/ result += residues[i]->get_intermol_binding(residues[j]);
        }
    }
    result += initial_int_clashes;              // Compensate for AminoAcid::get_intermol_binding() which factors in clashes.
    return result;
}

float Protein::get_intermol_clashes(Molecule* ligand)
{
    AminoAcid* laminos[100];
    Point size(0,0,0);
    int cres = get_residues_can_clash_ligand(laminos, ligand, ligand->get_barycenter(), size, nullptr);
    if (!cres) return 0;
    int i;
    float result = 0;
    for (i=0; i<cres; i++)
    {
        result += laminos[i]->Molecule::get_intermol_clashes(ligand);
    }
    return result;
}

float Protein::get_intermol_binding(Molecule* ligand)
{
    AminoAcid** laminos = new AminoAcid*[SPHREACH_MAX];
    Point size(0,0,0);
    int cres = get_residues_can_clash_ligand(laminos, ligand, ligand->get_barycenter(), size, nullptr);
    if (!cres) return 0;
    int i;
    float result = 0;
    for (i=0; i<cres; i++)
    {
        result += laminos[i]->Molecule::get_intermol_binding(ligand);
    }
    delete laminos;
    return result;
}

void Protein::find_residue_initial_bindings()
{
    if (!residues) return;

    aabridges.clear();
    
    int i, j, k;
    for (i=0; residues[i]; i++)
    {
        AminoAcid** aa = residues; // get_residues_can_clash(residues[i]->get_residue_no());

        #if _ALLOW_PROTONATE_PNICTOGENS
        // If the current residue has one or more negatively charged non-backbone atoms,
        // and a nearby residue has one or more neutral non-backbone pnictogen atoms not part of an amide,
        // put a positive charge on the pnictogen(s).
        if (residues[i]->get_charge() < -0.5)
        {
            for (j=0; aa[j]; j++)
            {
                Atom* a = residues[i]->get_nearest_atom(aa[j]->get_barycenter());
                Atom* b = aa[j]->get_nearest_atom(a->get_location());
                float r = a->distance_to(b);

                if (r < 4)
                {
                    int ac = aa[j]->get_atom_count();
                    for (k=0; k<ac; k++)
                    {
                        Atom* a = aa[j]->get_atom(k);
                        if ( !a->is_backbone && a->get_family() == PNICTOGEN && !a->get_charge() )
                        {
                            Atom* C = a->is_bonded_to(TETREL, 1);
                            if (C && C->is_bonded_to(CHALCOGEN, 2)) continue;
                            else a->increment_charge(0.75);
                        }
                    }
                }
            }
        }
        #endif

        float ib = 0, maxb = 0;
        AABridge aab;
        aab.aa1 = residues[i];
        for (j=0; aa[j]; j++)
        {
            if (aa[j] == residues[i]) continue;
            float f = residues[i]->get_intermol_binding(aa[j]);
            ib += f;
            if (j > i && f > maxb)
            {
                aab.aa2 = aa[j];
                maxb = f;
            }
        }

        if (ib >= 5)
        {
            if (residues[i]->movability != MOV_FORCEFLEX) residues[i]->movability = min(residues[i]->movability, MOV_PINNED);
            if (maxb >= 5) aabridges.push_back(aab);
        }

        #if _debug_locks
        if (residues[i]->get_residue_no() == _dbg_lock_res)
        {
            for (j=0; aa[j]; j++)
            {
                cout << _dbg_lock_res << " - " << aa[j]->get_residue_no() << ": " << residues[i]->get_intermol_binding(aa[j]) << endl;
            }
        }
        #endif

        // delete[] aa;
    }
}

void Protein::set_name_from_pdb_name(const char* pdb_name)
{
    char buffer[1024], *begin = buffer;
    strcpy(buffer, pdb_name);

    char* slash = strchr(buffer, '/');
    while (slash && slash[1])
    {
        begin = slash+1;
        slash = strchr(begin, '/');
    }

    char* dot = strchr(begin, '.');
    if (dot) *dot = 0;

    name = begin;
}

int Protein::load_pdb(FILE* is, int rno, char chain)
{
    AminoAcid* restmp[65536];
    char buffer[1024];
    Atom* a;

    if (residues) delete[] residues;
    if (sequence) delete[] sequence;
    if (ca) delete[] ca;
    // if (res_reach) delete res_reach;         // This was causing a segfault.
    if (metals) delete[] metals;

    origpdb_residues.clear();
    connections.clear();

    Atom* pdba[65536];

    AminoAcid useless('#');		// Feed it nonsense just so it has to load the data file.
    AminoAcid* prevaa = nullptr;

    int i, j, rescount=0;

    for (i=0; i<65536; i++) pdba[i] = nullptr;

    bool got_atoms = false;
    while (!feof(is))
    {
        try
        {
            int told = ftell(is);
            fgets(buffer, 1003, is);

            if (got_atoms &&
                buffer[0] == 'T' &&
                buffer[1] == 'E' &&
                buffer[2] == 'R'
               )
            {
                break;
            }

            if (buffer[0] == 'A' &&
                buffer[1] == 'N' &&
                buffer[2] == 'I' &&
                buffer[3] == 'S'
               )
                continue;

            else if (buffer[0] == 'A' &&
                     buffer[1] == 'T' &&
                     buffer[2] == 'O' &&
                     buffer[3] == 'M' &&
                     (buffer[22] == ' ' || (buffer[22] >= '0' && buffer[22] <= '9')) &&
                     (buffer[23] == ' ' || (buffer[23] >= '0' && buffer[23] <= '9')) &&
                     (buffer[24] == ' ' || (buffer[24] >= '0' && buffer[24] <= '9')) &&
                     (buffer[25] == ' ' || (buffer[25] >= '0' && buffer[25] <= '9'))
               )
            {
                if (buffer[21] != ' ' && buffer[21] != chain) continue;

                buffer[16] = ' ';
                fseek(is, told, SEEK_SET);

                char tmp3let[5];
                for (i=0; i<3; i++)
                    tmp3let[i] = buffer[17+i];
                tmp3let[3] = 0;
                tmp3let[4] = 0;
                got_atoms = true;

                for (i=0; i<256; i++)
                {
                    if (aa_defs[i].name[0] && !strcasecmp(aa_defs[i]._3let, tmp3let))
                    {
                        AminoAcid* aa = new AminoAcid(is, prevaa, rno);

                        int n = aa->get_atom_count();
                        for (j=0; j<n; j++)
                        {
                            Atom* la = aa->get_atom(j);
                            if (la && la->pdbidx) pdba[la->pdbidx] = la;
                        }

                        restmp[rescount++] = aa;
                        restmp[rescount] = NULL;
                        prevaa = aa;
                        Pose aap, filler;
                        aap.copy_state(aa);
                        for (; origpdb_residues.size() < aa->get_residue_no(); origpdb_residues.push_back(filler));
                        origpdb_residues.push_back(aap);
                        goto _found_AA;
                    }
                }

                // If no AA, load a metal.
                if (!metcount % 16)
                {
                    Atom** mtmp = new Atom*[metcount+20];
                    if (metals)
                    {
                        for (i=0; metals[i]; i++)
                            mtmp[i] = metals[i];
                        mtmp[i] = NULL;
                        delete metals;
                    }
                    metals = mtmp;
                }

                a = new Atom(is);
                metals[metcount++] = a;
                metals[metcount] = NULL;
                // cout << "M" << metcount << " ";

                for (i=0; i<26; i++)
                {
                    if (aa_defs[i]._1let && !strcasecmp(aa_defs[i]._3let, a->aa3let ))
                    {
                        a->aaletter = aa_defs[i]._1let;
                        break;
                    }
                }
            }
            else if (buffer[0] == 'R' &&
                     buffer[1] == 'E' &&
                     buffer[2] == 'M' &&
                     buffer[3] == 'A' &&
                     buffer[4] == 'R' &&
                     buffer[5] == 'K'
                )
            {
                add_remark(buffer);
                // cout << "Found remark " << buffer;

                if (buffer[7] == '6' && buffer[8] < '!')
                {
                    char** words = chop_spaced_words(buffer);
                    if (words[2] && !words[3])
                        name = words[2];
                }
            }
            else if (buffer[0] == 'C'
                  && buffer[1] == 'O'
                  && buffer[2] == 'N'
                  && buffer[3] == 'E'
                  && buffer[4] == 'C'
                  && buffer[5] == 'T'
                    )
            {
                buffer[6] = buffer[11] = buffer[16] = 0;
                int a1 = atoi(&buffer[7]);
                int a2 = atoi(&buffer[12]);
                if (pdba[a1] && pdba[a2])
                {
                    pdba[a1]->bond_to(pdba[a2], 1);
                    Bond* b = pdba[a1]->get_bond_between(pdba[a2]);
                    if (b) connections.push_back(b);
                }
            }
            

        _found_AA:
            ;
        }
        catch (unsigned int ex)
        {
            // cout << "Exception " << ex << endl;
            switch (ex)
            {
            case ATOM_NOT_OF_AMINO_ACID:
                if (!metcount % 16)
                {
                    Atom** mtmp = new Atom*[metcount+20];
                    if (metals)
                    {
                        for (i=0; metals[i]; i++)
                            mtmp[i] = metals[i];
                        mtmp[i] = NULL;
                        delete metals;
                    }
                    metals = mtmp;
                }

                a = new Atom(is);
                metals[metcount++] = a;
                metals[metcount] = NULL;

                for (i=0; i<26; i++)
                {
                    if (aa_defs[i]._1let && !strcasecmp(aa_defs[i]._3let, a->aa3let ))
                    {
                        a->aaletter = aa_defs[i]._1let;
                        break;
                    }
                }
                break;

            case NOT_ATOM_RECORD:
                fgets(buffer, 1003, is);
                cout << buffer << endl;
                if (buffer[0] == 'R' &&
                        buffer[1] == 'E' &&
                        buffer[2] == 'M' &&
                        buffer[3] == 'A' &&
                        buffer[4] == 'R' &&
                        buffer[5] == 'K'
                   )
                    add_remark(buffer);
                break;

            default:
                throw 0xbadca22;
            }
        }
    }

    int arrlimit = rescount+1;
    residues 	= new AminoAcid*[arrlimit];
    sequence 	= new char[arrlimit];
    ca       	= new Atom*[arrlimit];
    res_reach	= new float[arrlimit];

    for (i=0; i<arrlimit; i++)
    {
        residues[i] = NULL;
        sequence[i] = 0;
        ca[i] = NULL;
        res_reach[i] = 0;
    }

    for (i=0; i<rescount; i++)
    {
        try
        {
            residues[i] = restmp[i];
            residues[i]->clear_cache();
            residues[i]->establish_internal_clash_baseline();

            Atom *atom = residues[i]->get_atom("N"), *btom;
            AminoAcid* prev = get_residue(residues[i]->get_residue_no()-1);
            if (prev)
            {
                btom = prev->get_atom("C");
                if (atom && btom) atom->bond_to(btom, 1.5);
            }

            if (!aaptrmin.n || residues[i] < aaptrmin.paa) aaptrmin.paa = residues[i];
            if (!aaptrmax.n || residues[i] > aaptrmax.paa) aaptrmax.paa = residues[i];

            AADef* raa   = restmp[i]->get_aa_definition();
            if (!raa) cout << "Warning: Residue " << (i+1) << " has no AADef." << endl;
            sequence[i]  = raa ? raa->_1let : '?';
            ca[i]		 = restmp[i]->get_atom("CA");
            res_reach[i] = raa ? raa->reach : 2.5;
        }
        catch (int e)
        {
            cout << "Residue " << (i+1) << " threw an error." << endl;
            throw e;
        }
    }
    // cout << "Read residue " << *residues[rescount] << endl;
    residues[rescount] = 0;

    res_can_clash = nullptr;
    set_clashables();

    int l;
    std::vector<std::string> rem_hx = get_remarks("650 HELIX");
    for (l=0; l<rem_hx.size(); l++)
    {
        char buffer[1024];
        char buffer1[1024];
        strcpy(buffer, rem_hx[l].c_str());
        char** words = chop_spaced_words(buffer);

        set_region(words[3], atoi(words[4]), atoi(words[5]));
        regions_from = rgn_pdb;
        delete[] words;
    }

    std::vector<std::string> rem_st = get_remarks("800 SITE");
    for (l=0; l<rem_st.size(); l++)
    {
        char buffer[1024];
        char buffer1[1024];
        strcpy(buffer, rem_st[l].c_str());
        char** words = chop_spaced_words(buffer);

        if (!words[0] || !words[1] || !words[2] || !words[3] || !words[4]) continue;

        int f4 = atoi(words[4]);
        if (!strcmp(words[3], "BW"))
        {
            Ballesteros_Weinstein[f4] = atoi(words[5]);
        }
        delete[] words;
    }

    initial_int_clashes = get_internal_clashes();
    allocate_undo_poses();
    set_conditional_basicities();

    return rescount;
}

void Protein::allocate_undo_poses()
{
    int i, l;
    if (undo_poses) delete[] undo_poses;

    l = get_end_resno();
    undo_poses = new Pose*[l+4];
    undo_poses[0] = nullptr;
    for (i=1; i<l+4; i++)
    {
        AminoAcid* aa = get_residue(i);
        undo_poses[i] = aa ? new Pose(aa) : nullptr;
    }
}

void Protein::save_undo_state()
{
    if (mass_undoable) return;
    if (!undo_poses) allocate_undo_poses();

    int i, n;
    n = get_end_resno();

    for (i=1; i<=n; i++)
    {
        if (undo_poses[i])
        {
            AminoAcid* aa = get_residue(i);
            if (aa) undo_poses[i]->copy_state(aa);
        }
    }
}

void Protein::undo()
{
    mass_undoable = false;
    int i, n;
    n = get_end_resno();

    for (i=1; i<=n; i++)
    {
        if (undo_poses[i])
        {
            AminoAcid* aa = get_residue(i);
            if (aa) undo_poses[i]->restore_state(aa);
        }
    }
}

void Protein::revert_to_pdb()
{
    if (!residues) return;
    int i, n;
    n = origpdb_residues.size();
    if (!n) return;

    for (i=1; i<n; i++)
    {
        AminoAcid* aa = get_residue(i);
        if (aa) origpdb_residues[i].restore_state(aa);
    }
    allocate_undo_poses();
}

int Protein::get_bw50(int helixno)
{
    if (helixno < 1 || helixno > 78) return -1;
    return Ballesteros_Weinstein[helixno];
}

void Protein::set_bw50(int hxno, int resno)
{
    if (hxno < 1 || hxno > 78) return;
    // if (Ballesteros_Weinstein[hxno]) return;
    Ballesteros_Weinstein[hxno] = resno;
}

int Protein::get_seq_length()
{
    if (!sequence) return 0;
    int i;
    for (i=0; sequence[i]; i++);	// Obtain the count.
    return i;
}

std::string Protein::get_sequence()
{
    if (!sequence) return 0;
    return std::string(sequence);
}

int Protein::get_start_resno()
{
    if (!residues || !residues[0]) return 0;
    else return residues[0]->get_residue_no();
}

int Protein::get_end_resno()
{
    int retval = 0;
    if (!residues || !residues[0]) return retval;
    int i;
    for (i=0; residues[i]; i++) retval = residues[i]->get_residue_no();
    return retval;
}

void Protein::add_remark(const char* remark)
{
    if (!remarks) remarks = new char*[65536];

    remarks[remarksz] = new char[strlen(remark)+2];
    strcpy(remarks[remarksz++], remark);

    char buffer[1024];
    strcpy(buffer, remark);
    char** words = chop_spaced_words(buffer);

    if (words[1]
        && !strcmp(words[1], "800")
        && words[2]
        && !strcmp(words[2], "SITE")
        && !strcmp(words[3], "BW")
        )
    {
        int f4 = atoi(words[4]);
        Ballesteros_Weinstein[f4] = atoi(words[5]);
    }
}

std::vector<std::string> Protein::get_remarks(std::string search_for)
{
    std::vector<string> retval;
    int i;
    for (i=0; i<remarksz; i++)
    {
        if (!search_for.length() || strstr(remarks[i], search_for.c_str()))
        {
            char buffer[1024];
            strcpy(buffer, remarks[i]);
            retval.push_back((std::string)buffer);
        }
    }

    return retval;
}

void Protein::add_remark(std::string new_remark)
{
    add_remark(new_remark.c_str());

    // TODO: Sort remarks by number, becarefuling to preserve the sequence of same numbered remarks.
}

void Protein::set_clashables(int resno, bool recursed)
{
    int i, j, k;

    // cout << "Setting clashables." << endl;


    int maxres = get_end_resno();
    if (!res_can_clash)
    {
        res_can_clash = new AminoAcid**[maxres+8];
        for (i=0; i<=maxres; i++) res_can_clash[i] = nullptr;
    }

    int sr = get_start_resno(), er = get_end_resno();
    int sr1 = sr, er1 = er;
    if (resno > 0) sr = er = resno;
    for (i=sr; i<=er; i++)
    {
        AminoAcid* resi = get_residue(i);
        if (!resi) continue;

        if (debug) *debug << endl << "Testing residue " << resi->get_residue_no() << endl;
        AminoAcid* temp[maxres+1];
        k=0;
        for (j=sr1; j<=er1; j++)
        {
            if (j == i) continue;
            AminoAcid* resj = get_residue(j);
            if (!resj) continue;
            if (resi->can_reach(resj))
            {
                temp[k++] = resj;
                if (debug) *debug << *resj << " can reach " << *resi << endl;
            }
        }

        if (!res_can_clash[i])
            res_can_clash[i] = new AminoAcid*[maxres+8];
        for (j=0; j<k; j++)
        {
            res_can_clash[i][j] = temp[j];
            if (resno > 0 && !recursed)
            {
                set_clashables(temp[j]->get_residue_no(), true);
            }
            /*cout << residues[i]->get_aa_definition()->_3let << residues[i]->get_residue_no()
            	 << " can clash with "
            	 << res_can_clash[i][j]->get_aa_definition()->_3let << res_can_clash[i][j]->get_residue_no()
            	 << endl;*/
        }
        res_can_clash[i][k] = 0;
    }

    res_can_clash[maxres+1] = 0;

    /*for (i=0; i<=maxres; i++)
    {
    	cout << i << ": ";
    	for (j=0; res_can_clash[i][j]; j++)
    		cout << *res_can_clash[i][j] << " ";
    	cout << endl;
    }*/
}

void Protein::set_conditional_basicities()
{
    if (!residues) return;

    int i, j, l, n;
    for (i=0; residues[i]; i++)
    {
        if (!residues[i]->conditionally_basic()) continue;

        std::vector<AminoAcid*> nearby = get_residues_near(residues[i]->get_CA_location(), 10, true);
        n = nearby.size();
        if (!n) continue;

        Molecule* mols[256];
        j = 0;
        for (l=0; l<n; l++)
        {
            if (nearby[l]->get_charge() > -0.25) continue;
            mols[j++] = (Molecule*)nearby[l];
        }
        mols[j] = nullptr;

        residues[i]->set_conditional_basicity(mols);
    }
}

std::vector<AminoAcid*> Protein::get_residues_near(Point pt, float maxr, bool facing)
{
    std::vector<AminoAcid*> retval;

    float cb_tolerance_angle = 44 * fiftyseventh;
    float tolerance_sine = sin(cb_tolerance_angle);

    if (!residues) return retval;
    int i, j;

    #if _DBG_TUMBLE_SPHERES
    cout << "Protein::get_residues_near(" << pt << ", " << maxr << ")" << endl;
    #endif

    for (i=0; residues[i]; i++)
    {
        float r = pt.get_3d_distance(residues[i]->get_atom_location("CA"));

        if (facing && residues[i]->get_atom("CB"))
        {
            float r1 = pt.get_3d_distance(residues[i]->get_atom_location("CB"));
            float r2 = residues[i]->get_atom_location("CA").get_3d_distance(residues[i]->get_atom_location("CB"));
            float tolerance = r2 * tolerance_sine;
            if (r1 > r+tolerance) continue;
        }

        if (r <= maxr)
        {
            retval.push_back(residues[i]);
            #if _DBG_TUMBLE_SPHERES
            cout << residues[i]->get_3letter() << residues[i]->get_residue_no() << " ";
            #endif
        }
    }
    #if _DBG_TUMBLE_SPHERES
    cout << endl << endl;
    #endif

    return retval;
}

std::vector<AminoAcid*> Protein::get_contact_residues(Protein* op, float cd)
{
    std::vector<AminoAcid*> retval;

    if (!residues) return retval;
    if (!op->residues) return retval;

    int m = get_end_resno(), n = op->get_end_resno();
    bool dirty[n + 4];
    
    int i, j;
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=1; i<m; i++)
    {
        AminoAcid* a = get_residue(i);
        if (!a) continue;

        bool adirty = false;

        for (j=1; j<n; j++)
        {
            AminoAcid* b = op->get_residue(j);

            if (b)
            {
                if (dirty[j]) continue;

                float f = a->get_reach() + b->get_reach() + cd;
                float r = a->get_CA_location().get_3d_distance(b->get_CA_location());

                if (r < f)
                {
                    if (!adirty) retval.push_back(a);
                    retval.push_back(b);
                    adirty = dirty[j] = true;
                }
            }
        }
    }

    return retval;
}

AminoAcid** Protein::get_residues_can_clash(int resno)
{
    if (!residues) return 0;
    if (!res_can_clash) set_clashables();

    int i, j;
    for (i=0; residues[i]; i++)
    {
        if (residues[i]->get_residue_no() == resno)
        {
            if (!res_can_clash[i] || !res_can_clash[i][0]) set_clashables();
            residues[i]->mclashables = (Molecule**)res_can_clash[i];
            /*cout << i << ": " << flush;
            if (res_can_clash[i] && res_can_clash[i][0])
            {
            	cout << *res_can_clash[i][0] << endl;
            	for (j=0; res_can_clash[i][j]; j++)
            		cout << j << " " << flush;
            	cout << endl;
            }*/
            return res_can_clash[i];
        }
    }

    return 0;
}

std::vector<AminoAcid*> Protein::get_residues_can_clash(int sr, int er)
{
    std::vector<AminoAcid*> result;
    int i, j;
    bool already[get_seq_length()+4];
    for (i=0; residues[i]; i++) already[i] = false;

    for (j=sr; j<=er; j++)
    {
        AminoAcid* aa = get_residue(j);
        for (i=0; residues[i]; i++)
        {
            if (already[i]) continue;

            int resno = residues[i]->get_residue_no();
            if (resno >= sr && resno <= er) continue;

            if (aa->can_reach(residues[i]))
            {
                result.push_back(residues[i]);
                already[i] = true;
            }
        }
    }

    return result;
}

int Protein::get_residues_can_clash_ligand(AminoAcid** reaches_spheroid,
        Molecule* ligand,
        const Point nodecen,
        const Point size,
        const int* addl_resno
        )
{
    int i, j, sphres = 0;
    int seql = get_end_resno();
    bool resno_already[8192];
    for (i=0; i<8192; i++) resno_already[i] = false;

    for (i=0; i<SPHREACH_MAX; i++) reaches_spheroid[i] = NULL;

    for (i=1; i<seql && residues[i]; i++)
    {
        AminoAcid* aa = residues[i-1];
        if (!aa) continue;

        int resno = aa->get_residue_no();

        if (addl_resno)
        {
            for (j=0; addl_resno[j]; j++)
            {
                if (addl_resno[j] == resno)
                {
                    if (!resno_already[resno])
                    {
                        reaches_spheroid[sphres++] = aa;
                        resno_already[resno] = true;
                        #if _DBG_REACHLIG
                        if (debug)
                        {
                            Star s;
                            s.paa = aa;
                            *debug << std::hex << s.n << std::dec << " " << flush;
                        }
                        #endif
                    }
                    continue;
                }
            }
        }

        Atom* ca = aa->get_atom("CA");
        if (!ca) continue;
        Atom* cb = aa->get_atom("CB");

        Point pt = ca->get_location();
        Atom* a = ligand->get_nearest_atom(pt);
        Point pt2;
        if (a) pt2 = a->get_location();
        else   pt2 = ligand->get_barycenter();

        if (cb)
        {
            Point pt1 = pt.subtract(&pt2);
            float angle = find_3d_angle(cb->get_location(), pt2, ca->get_location());
            if (pt1.magnitude() < 3 || angle < _can_clash_angle)
            {
                if (pt1.magnitude() < (aa->get_reach()+2))
                {
                    if (!resno_already[resno])
                    {
                        reaches_spheroid[sphres++] = aa;
                        resno_already[resno] = true;
                        #if _DBG_REACHLIG
                        if (debug)
                        {
                            Star s;
                            s.paa = aa;
                            *debug << std::hex << s.n << std::dec << " " << flush;
                        }
                        #endif
                    }
                    continue;
                }
            }

            angle = find_3d_angle(cb->get_location(), nodecen, ca->get_location());
            if (angle > _can_clash_angle) continue;
        }

        pt = pt.subtract(&nodecen);
        Point pt1 = pt;
        float dist = pt.magnitude();
        pt1.scale(fmax(dist - aa->get_reach(), 0));

        pt1.x /= size.x;
        pt1.y /= size.y;
        pt1.z /= size.z;

        SCoord dir(&pt1);

        if (dir.r <= 1)
        {
            if (!resno_already[resno])
            {
                reaches_spheroid[sphres++] = aa;
                resno_already[resno] = true;
                #if _DBG_REACHLIG
                if (debug)
                {
                    Star s;
                    s.paa = aa;
                    *debug << std::hex << s.n << std::dec << " " << flush;
                }
                #endif
            }
        }

        aa->reset_conformer_momenta();
    }

    reaches_spheroid[sphres] = NULL;
    #if _DBG_REACHLIG
    if (debug) *debug << endl << flush;
    #endif

    return sphres;
}

bool Protein::aa_ptr_in_range(AminoAcid* aaptr)
{
    if (!aaptr) return false;
    if (aaptr < aaptrmin.paa || aaptr > aaptrmax.paa) return false;
    else return true;
}

Molecule* Protein::metals_as_molecule()
{
    Molecule* met=NULL;
    if (metals) met = new Molecule("(metals)", metals);

    InteratomicForce f;

    // Associate coordinating residue atoms with the metal ions so that the ions can have geometry.
    // Otherwise we end up with an exception later on in the InteratomicForce::total_binding() function.
    if (mcoord_resnos && mcoord_resnos[0])
    {
        int i, j, k, l, m, n;
        for (i=0; metals[i]; i++)
        {
            Point mloc = metals[i]->get_location();
            k = 0;
            for (j=0; mcoord_resnos[j]; j++)
            {
                AminoAcid* caa = get_residue(mcoord_resnos[j]);		// caa = coordinating amino acid.
                if (!caa) continue;
                caa->movability = MOV_NONE;
                Atom* mca = caa->get_nearest_atom(mloc, mcoord);
                if (!mca) continue;

                float r = mloc.get_3d_distance(mca->get_location());

                if (r < 2 * InteratomicForce::coordinate_bond_radius(metals[i], mca, mcoord))
                {
                    Bond* b = metals[i]->get_bond_by_idx(k++);
                    if (!b) break;
                    b->btom = mca;
                    b->cardinality = 0.5;
                    // cout << metals[i]->name << " coordinates to " << *caa << ":" << mca->name << endl;
                }
            }
        }
    }

    return met;
}

int Protein::search_sequence(const int sr, const int esr, const char* psz, const int threshold, int* psim)
{
    int i, j, k = 0, l, num_eq;
    float m = 0, n = 0, sim;
    char aac;
    AminoAcid* aa;

    for (i=sr; i<esr; i++)
    {
        m = num_eq = 0;
        char lc = 0;
        l = 0;
        for (j=0; psz[j]; j++)
        {
            if (psz[j] == '^') j++;
            char c = psz[j];
            aa = get_residue(i+l);
            if (!aa) continue;
            aac = aa->get_letter();

            if (lc == '^' && aa->get_residue_no() != get_start_resno()) goto _wrong_place;
            if (psz[j+1] == '$' && aa->get_residue_no() != get_end_resno()) goto _wrong_place;

            if (c == 'X')
            {
                m += 1;
                num_eq++;
            }
            else
            {
                if (c == aac)
                {
                    num_eq++;
                    sim = 1;
                }
                else sim = aa->similarity_to(c);
                m += sim;
            }

            lc = c;
            l++;
        }

        if (m > n && num_eq >= threshold)
        {
            k = i;
            n = m;
        }

        _wrong_place:
        ;
    }
    sim = n;
    if (psim) *psim = sim;

    return k;
}

void Protein::rotate_backbone(int resno, bb_rot_dir dir, float angle)
{
    save_undo_state();

    AminoAcid* bendy = get_residue(resno);
    if (!bendy) return;
    LocatedVector lv = bendy->rotate_backbone(dir, angle);
    bendy->ensure_pi_atoms_coplanar();

    if (lv.r)
    {
        int i, inc;
        AminoAcid* movable;
        inc = (dir == CA_desc || dir == C_desc) ? -1 : 1;

        for (i=resno+inc; movable = get_residue(i); i+=inc)
        {
            // cout << "Rotating " << i << endl;
            movable->rotate(lv, angle);
            movable->ensure_pi_atoms_coplanar();
            set_clashables(i);
        }
    }
}

void Protein::rotate_backbone_partial(int startres, int endres, bb_rot_dir dir, float angle)
{
    save_undo_state();
    if (startres == endres) return;
    int inc = (dir == CA_desc || dir == C_desc) ? -1 : 1;
    if (sgn(endres - startres) != sgn(inc))
    {
        cout << "ERROR: direction mismatch " << startres << "->" << endres
             << " but direction is " << inc << endl;
        return;
    }

    AminoAcid* bendy = get_residue(startres);
    if (!bendy) return;
    LocatedVector lv = bendy->rotate_backbone(dir, angle);
    bendy->ensure_pi_atoms_coplanar();

    if (lv.r)
    {
        int i;
        AminoAcid* movable;

        for (i=startres+inc; movable = get_residue(i); i+=inc)
        {
            movable->rotate(lv, angle);
            movable->ensure_pi_atoms_coplanar();
            if (i == endres) break;
        }
    }

    set_clashables();
}

void Protein::conform_backbone(int startres, int endres, int iters, bool backbone_atoms_only)
{
    Point pt;
    conform_backbone(startres, endres, NULL, pt, NULL, pt, iters, backbone_atoms_only);
}

void Protein::conform_backbone(int startres, int endres, Atom* a, Point target, int iters)
{
    Point pt;
    conform_backbone(startres, endres, a, target, NULL, pt, iters, false);
}

void Protein::conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, int iters, bool backbone_atoms_only)
{
    conform_backbone(startres, endres, a1, target1, a2, target2, nullptr, Point(), iters, backbone_atoms_only);
}

#define DBGCONF 0

void Protein::conform_backbone(int startres, int endres,
                               Atom* a1, Point target1,
                               Atom* a2, Point target2,
                               Atom* a3, Point target3,
                               int iters, bool backbone_atoms_only
                              )
{
    cout << "conform_backbone() has never worked and probably never will." << endl;
    throw -1;

    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;

    int inc = sgn(endres-startres);
    int res, i, j, iter;
    bb_rot_dir dir1 = (inc>0) ? N_asc : CA_desc,
               dir2 = (inc>0) ? CA_asc : C_desc;

    int am = abs(endres-startres), minres = (inc>0) ? startres : endres;
    float momenta1o[am+4], momenta2o[am+4], momenta1e[am+4], momenta2e[am+4];
    int eando_res[am+4];
    float eando_mult[am+4];
    float r = 0, lastr = 99999;
    int iters_since_improvement = 0;

    for (res = startres; res; res += inc)
    {
        int residx = res-minres;
        momenta1o[residx] = randsgn()*_fullrot_steprad;
        momenta2o[residx] = randsgn()*_fullrot_steprad;
        momenta1e[residx] = randsgn()*_fullrot_steprad;
        momenta2e[residx] = randsgn()*_fullrot_steprad;
        eando_res[residx] = min(res + (rand() % 5) + 1, endres);
        if (eando_res[residx] == res) eando_res[residx] = 0;
        eando_mult[residx] = 1;
        if (res == endres) break;
    }

    set_clashables();
    float tolerance = 1.2, alignfactor = 100, reversal = -0.81, enhance = 1.5;
    int ignore_clashes_until = iters*0.666;
    for (iter=0; iter<iters; iter++)
    {
        #if DBGCONF
        // cout << "Iteration " << iter << endl;
        cout << " " << iter << flush;
        #endif

        for (res = startres; res != endres; res += inc)
        {
            int residx = res-minres;

            // Get the preexisting nearby residues and inter-residue binding/clash value.
            // These will likely have changed since last iteration.
            float bind=0, bind1=0, angle;

            for (i=res; i != endres; i += inc)
            {
                AminoAcid* aa = get_residue(i);
                if (!aa) continue;
                AminoAcid** rcc = get_residues_can_clash(i);
                if (!rcc) cout << "No clashables." << endl;
                if (a1 && (iter >= ignore_clashes_until)) bind -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
                else bind += aa->get_intermol_binding(rcc, backbone_atoms_only);
            }
            if (a1 && iter>10)
            {
                Point pt = a1->get_location();
                bind += alignfactor/(pt.get_3d_distance(target1)+0.001);
            }
            if (a2 && iter>10)
            {
                Point pt = a2->get_location();
                bind += alignfactor/(pt.get_3d_distance(target2)+0.001);
            }

            if (reinterpret_cast<long>(get_residue(res)) < 0x1000) cout << "Warning missing residue " << res << endl << flush;
            else if (strcmp(get_residue(res)->get_3letter(), "Pro"))		// TODO: Don't hard code this to proline, but check bond flexibility.
            {
                // Rotate the first bond a random amount. TODO: use angular momenta.
                angle = (iter & 1) ? momenta1o[residx] : momenta1e[residx]; // frand(-_fullrot_steprad, _fullrot_steprad);
                rotate_backbone_partial(res, endres, dir1, angle);
                if ((iter & 1) && eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir2, -angle*eando_mult[residx]);

                // Have bindings/clashes improved?
                for (i=res; i != endres; i += inc)
                {
                    AminoAcid* aa = get_residue(i);
                    AminoAcid** rcc = get_residues_can_clash(i);
                    if (a1 && (iter >= ignore_clashes_until)) bind1 -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
                    else bind1 += aa->get_intermol_binding(rcc, backbone_atoms_only);
                }
                if (a1)
                {
                    Point pt = a1->get_location();
                    bind1 += alignfactor/(pt.get_3d_distance(target1)+0.001);
                }
                if (a2)
                {
                    Point pt = a2->get_location();
                    bind1 += alignfactor/(pt.get_3d_distance(target2)+0.001);
                }
                if (a3)
                {
                    Point pt = a3->get_location();
                    bind1 += alignfactor/(pt.get_3d_distance(target3)+0.001);
                }

                // If no, put it back.
                // if (res == startres) cout << bind << " v. " << bind1 << endl;
                if (bind1 < tolerance*bind || (a1 && iters_since_improvement > 10 && frand(0,1)<0.25))
                {
                    rotate_backbone_partial(res, endres, dir1, -angle);
                    if (eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir2, angle*eando_mult[residx]);
                    if (iter & 1) momenta1o[residx] *= reversal;
                    else momenta1e[residx] *= reversal;
                }
                else
                {
                    if (bind1 < bind) bind = bind1;
                    if (iter & 1) momenta1o[residx] *= enhance;
                    else momenta1e[residx] *= enhance;
                }
            }

            // Rotate the second bond.
            angle = (iter & 1) ? momenta2o[residx] : momenta2e[residx]; // frand(-_fullrot_steprad, _fullrot_steprad);
            rotate_backbone_partial(res, endres, dir2, angle);
            if ((iter & 1) && eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir1, -angle*eando_mult[residx]);

            // Improvement?
            bind1 = 0;
            for (i=res; i != endres; i += inc)
            {
                AminoAcid* aa = get_residue(i);
                if (!aa) continue;
                AminoAcid** rcc = get_residues_can_clash(i);
                if (a1 && (iter >= ignore_clashes_until)) bind1 -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
                else bind1 += aa->get_intermol_binding(rcc, backbone_atoms_only);
            }
            if (a1)
            {
                Point pt = a1->get_location();
                bind1 += alignfactor/(pt.get_3d_distance(target1)+0.001);
            }
            if (a2)
            {
                Point pt = a2->get_location();
                bind1 += alignfactor/(pt.get_3d_distance(target2)+0.001);
            }
            if (a3)
            {
                Point pt = a3->get_location();
                bind1 += alignfactor/(pt.get_3d_distance(target3)+0.001);
            }

            // If no, put it back.
            // if (res == startres) cout << bind << " vs. " << bind1 << endl;
            if (bind1 < tolerance*bind || (a1 && iters_since_improvement > 10 && frand(0,1)<0.25))
            {
                rotate_backbone_partial(res, endres, dir2, -angle);
                if ((iter & 1) && eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir1, angle*eando_mult[residx]);
                if (iter & 1) momenta2o[residx] *= reversal;
                else momenta2e[residx] *= reversal;
            }
            else
            {
                if (iter & 1) momenta2o[residx] *= enhance;
                else momenta2e[residx] *= enhance;
            }

            alignfactor *= 1.003;
            tolerance = ((tolerance-1)*0.97)+1;
        }

        r = 0;
        if (a1)
        {
            Point pt = a1->get_location();
            r += pt.get_3d_distance(target1);
        }
        if (a2)
        {
            Point pt = a2->get_location();
            r += pt.get_3d_distance(target2);
        }
        if (a3)
        {
            Point pt = a3->get_location();
            r += pt.get_3d_distance(target3);
        }

        if (r < 0.999*lastr) iters_since_improvement = 0;
        else iters_since_improvement++;
        lastr = r;

        #if DBGCONF
        if (r) cout << "." << r << flush;
        #endif
    }
    #if DBGCONF
    cout << endl;
    #endif
    if (r > 2.5) cout << "Warning - protein strand alignment anomaly outside of tolerance." << endl << "# Anomaly is " << r << " Angstroms." << endl;

    mass_undoable = wmu;
}

#define DBG_BCKCONN 0
#define _INCREMENTAL_BKCONN 1
void Protein::backconnect(int startres, int endres)
{
    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;

    int i;
    int inc = sgn(endres-startres);

    #if _INCREMENTAL_BKCONN
    int iter;
    for (iter=24; iter>=0; iter--)
    {
    #endif

        // Glue the last residue onto the target.
        // Then adjust its inner bonds so the other end points as closely to the previous residue as possible.
        // Then do the same for the previous residue, and the one before, etc,
        // all the way back to the starting residue.
        // Give a warning if the starting residue has an anomaly > 0.1A.
        AminoAcid *next, *curr, *prev;
        int pointer = endres;
        float movfactor = 1, decrement, anomaly = 0;

        decrement = 1.0 / fabs(endres - startres);		// if exterior is the opposite of interior, what's the opposite of increment?

        next = get_residue(endres+inc);
        curr = get_residue(pointer);
        prev = get_residue(pointer-inc);

        #if DBG_BCKCONN
        cout << "backconnect( " << startres << ", " << endres << ")" << endl;
        #endif
        while (next && curr)
        {
            #if DBG_BCKCONN
            cout << pointer << ":";
            #endif

            Point pts[4];
            (inc > 0) ? next->predict_previous_COCA(pts) : next->predict_next_NHCA(pts);
            Point ptsc[4];

            #if _INCREMENTAL_BKCONN
            ptsc[0] = (inc > 0) ? curr->get_atom_location("C") : curr->get_atom_location("N");
            ptsc[1] = (inc > 0) ? curr->get_atom_location("O") : curr->HN_or_substitute_location();
            ptsc[2] = curr->get_atom_location("CA");

            Point _4avg[3];
            for (i=0; i<3; i++)
            {
                _4avg[0] = pts[i];
                _4avg[1] = ptsc[i];
                _4avg[0].weight = movfactor;
                _4avg[1].weight = 1.0 - movfactor;
                pts[i] = average_of_points(_4avg, 2);
            }
            #endif

            curr->attach_to_prediction(pts, inc > 0);
            // break;

            #if DBG_BCKCONN
            cout << "g";
            #endif

            if (prev)
            {
                MovabilityType fmov = curr->movability;
                curr->movability = MOV_ALL;

                #if DBG_BCKCONN
                cout << "a";
                #endif

                pts[4];
                (inc < 0) ? prev->predict_previous_COCA(pts) : prev->predict_next_NHCA(pts);
                Point target_heavy = pts[0];
                Point target_pole = pts[1];

                #if DBG_BCKCONN
                cout << "b";
                #endif

                float theta, step = fiftyseventh*1.0, r, btheta = 0, bestr;
                for (theta=0; theta < M_PI*2; theta += step)
                {
                    r = target_heavy.get_3d_distance( (inc > 0) ? curr->get_atom_location("N") : curr->get_atom_location("C") );
                    // r -= target_pole.get_3d_distance( (inc > 0) ? curr->HN_or_substitute_location() : curr->get_atom_location("O") );
                    if (inc > 0) r -= target_pole.get_3d_distance(curr->HN_or_substitute_location());

                    if (!theta || (r < bestr))
                    {
                        bestr = r;
                        btheta = theta;
                    }

                    curr->rotate_backbone( (inc > 0) ? CA_desc : N_asc, step );
                }
                curr->rotate_backbone( (inc > 0) ? CA_desc : N_asc, btheta );

                #if DBG_BCKCONN
                cout << "c";
                #endif

                btheta=0;
                for (theta=0; theta < M_PI*2; theta += step)
                {
                    r = target_heavy.get_3d_distance( (inc > 0) ? curr->get_atom_location("N") : curr->get_atom_location("C") );
                    // r -= target_pole.get_3d_distance( (inc > 0) ? curr->HN_or_substitute_location() : curr->get_atom_location("O") );
                    if (inc < 0) r -= target_pole.get_3d_distance(curr->get_atom_location("O"));

                    if (!theta || (r < bestr))
                    {
                        bestr = r;
                        btheta = theta;
                    }

                    curr->rotate_backbone( (inc > 0) ? C_desc : CA_asc, step );
                }
                curr->rotate_backbone( (inc > 0) ? C_desc : CA_asc, btheta );
                anomaly = bestr;

                #if DBG_BCKCONN
                cout << "d";
                #endif

                curr->movability = fmov;
            }

            if (pointer == startres)
            {
                #if DBG_BCKCONN
                cout << ". ";
                #endif
                break;
            }

            pointer -= inc;
            next = curr;
            curr = prev;
            prev = get_residue(pointer-inc);
            movfactor -= decrement;

            #if DBG_BCKCONN
            cout << ", ";
            #endif
        };

        #if DBG_BCKCONN
        cout << endl;
        #endif

        /*if (anomaly > 0.1) cout << "Warning! conform_backbone( " << startres << ", " << endres << " ) anomaly out of range." << endl
        						<< "# " << (startres-inc) << " anomaly: " << anomaly << endl;*/
        #if _INCREMENTAL_BKCONN
    }
        #endif

    mass_undoable = wmu;
}

void Protein::make_helix(int startres, int endres, float phi, float psi)
{
    make_helix(startres, endres, endres, phi, psi);
}

#define DBG_ASUNDER_HELICES 0
void Protein::make_helix(int startres, int endres, int stopat, float phi, float psi)
{
    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;

    int inc = sgn(endres-startres);
    if (!inc) inc = sgn(stopat-startres);
    else if (inc != sgn(stopat-startres)) return;
    if (inc > 0 && stopat < endres) endres = stopat;
    if (inc < 0 && stopat > endres) endres = stopat;
    int res, i, j, iter;
    int phis[365], psis[365];
    bb_rot_dir dir1 = (inc>0) ? N_asc : CA_desc,
               dir2 = (inc>0) ? CA_asc : C_desc;

    #if DBG_ASUNDER_HELICES
    cout << "make_helix( " << startres << ", " << endres << ", " << stopat << ", " << phi << ", " << psi << " )" << endl;
    #endif
    for (res = startres; inc; res += inc)
    {
        AminoAcid* aa = get_residue(res);

        LocRotation lr = aa->enforce_peptide_bond();
        aa->ensure_pi_atoms_coplanar();

        if (lr.v.r)
        {
            AminoAcid* movable;

            if (res != endres)
                for (i=res+inc; movable = get_residue(i); i+=inc)
                {
                    LocatedVector lv = lr.get_lv();
                    movable->rotate(lv, lr.a);
                    movable->ensure_pi_atoms_coplanar();
                    #if DBG_ASUNDER_HELICES
                    cout << i << " ";
                    #endif
                    if (i == stopat) break;
                }
            #if DBG_ASUNDER_HELICES
            cout << endl;
            #endif
        }

        LocRotation* lr2 = aa->flatten();

        for (j=0; j<5; j++)
        {
            // cout << "Rotating " << *aa << " " << lr2[j].a*fiftyseven << " degrees." << endl;
            if (lr2[j].v.r && lr2[j].a)
            {
                AminoAcid* movable;

                if (res != endres) for (i=res+1; movable = get_residue(i); i+=inc)
                    {
                        LocatedVector lv = lr2[j].get_lv();
                        movable->rotate(lv, lr2[j].a);
                        movable->ensure_pi_atoms_coplanar();
                        #if DBG_ASUNDER_HELICES
                        cout << i << " ";
                        #endif
                        if (i == stopat) break;
                    }
                #if DBG_ASUNDER_HELICES
                cout << endl;
                #endif

                int round_theta = (int)(lr2[j].a*fiftyseven+0.5);
                while (round_theta < 0) round_theta += 360;
                if (round_theta <= 360)
                {
                    if (j == 2) phis[round_theta]++;
                    if (j == 3) psis[round_theta]++;
                }
            }
        }
        delete[] lr2;

        // cout << "Rotating " << *aa << " phi " << (phi*fiftyseven) << " degrees." << endl;
        lr = aa->rotate_backbone_abs(dir1, phi);
        aa->ensure_pi_atoms_coplanar();

        if (lr.v.r)
        {
            AminoAcid* movable;

            if (res != endres) for (i=res+inc; movable = get_residue(i); i+=inc)
            {
                // cout << i << " ";
                LocatedVector lv = lr.get_lv();
                movable->rotate(lv, lr.a);
                movable->ensure_pi_atoms_coplanar();
                #if DBG_ASUNDER_HELICES
                cout << i << " ";
                #endif
                if (i == stopat) break;
            }
            #if DBG_ASUNDER_HELICES
            cout << endl;
            #endif
        }
        // cout << endl;

        // cout << "Rotating " << *aa << " psi " << (psi*fiftyseven) << " degrees." << endl;
        lr = aa->rotate_backbone_abs(dir2, psi);
        aa->ensure_pi_atoms_coplanar();

        if (lr.v.r)
        {
            AminoAcid* movable;

            if (res != endres) for (i=res+inc; movable = get_residue(i); i+=inc)
                {
                    // cout << i << " ";
                    LocatedVector lv = lr.get_lv();
                    movable->rotate(lv, lr.a);
                    movable->ensure_pi_atoms_coplanar();
                    #if DBG_ASUNDER_HELICES
                    cout << i << " ";
                    #endif
                    if (i == stopat) break;
                }
            #if DBG_ASUNDER_HELICES
            cout << endl;
            #endif
        }
        // cout << endl;

        if (res == endres) break;
    }

    // Hang onto this line, might want it later.
    // if (phi > 0.1 || psi > 0.1) conform_backbone(startres, endres, 20, true);

    #if 0
    // This is for finding values of phi and psi for helices.
    for (j=0; j<=360; j++)
    {
        if (phis[j])
        {
            cout << "φ=" << j << ": ";
            for (i=0; i<phis[j]; i++) cout << "*";
            cout << endl;
        }
    }
    for (j=0; j<=360; j++)
    {
        if (psis[j])
        {
            cout << "ψ=" << j << ": ";
            for (i=0; i<psis[j]; i++) cout << "*";
            cout << endl;
        }
    }
    #endif

    set_clashables();

    int seql = get_seq_length();
    Molecule* aas[seql+4];
    for (i=startres; i<=endres; i++)
    {
        aas[i-startres] = get_residue(i);
        aas[i-startres]->movability = MOV_FLEXONLY;
        aas[i-startres+1] = 0;
    }
    aas[endres-startres+1] = 0;
    Molecule::conform_molecules(aas, 25);

    mass_undoable = wmu;
}

void Protein::delete_residue(int resno)
{
    if (!resno) return;
    // if (resno > get_seq_length()) return;
    if (!residues) return;

    int i, j;
    for (i=0; residues[i]; i++)
    {
        if (residues[i]->get_residue_no() == resno)
        {
            for (j=i+1; residues[j-1]; j++)
            {
                residues[j-1] = residues[j];
                sequence[j-1] = sequence[j];
            }
            // TODO: sequence.
            ca = 0;
            res_reach = 0;

            return;
        }
    }
}

void Protein::delete_residues(int startres, int endres)
{
    int i;
    for (i=startres; i<=endres; i++) delete_residue(i);
}

void Protein::renumber_residues(int startres, int endres, int new_startres)
{
    int i, j;
    j = new_startres - startres;
    for (i=0; residues[i]; i++) 
    {
        residues[i]->renumber(residues[i]->get_residue_no() + j);
    }
    allocate_undo_poses();
}

void Protein::delete_sidechain(int resno)
{
    if (!resno) return;
    // if (resno > get_seq_length()) return;
    if (!residues) return;

    AminoAcid* aa = get_residue(resno);
    aa->delete_sidechain();
}

void Protein::delete_sidechains(int startres, int endres)
{
    if (!residues) return;
    int i;
    for (i=0; residues[i]; i++)
    {
        int res = residues[i]->get_residue_no();
        if (res >= startres && res <= endres) residues[i]->delete_sidechain();
    }
}

Protein* gmprot;
Point gmtgt;

void ext_mtl_coord_cnf_cb(int iter)
{
    gmprot->mtl_coord_cnf_cb(iter);
}

int Protein::get_metals_count()
{
    if (!metals) return 0;

    int i;
    for (i=0; metals[i]; i++);      // Get count.
    return i;
}

MetalCoord* Protein::coordinate_metal(Atom* metal, int residues, int* resnos, std::vector<string> res_anames)
{
    int i, j=0, k=0, l, n;
    if (!m_mcoord)
    {
        m_mcoord = new MetalCoord*[2];
        m_mcoord[0] = new MetalCoord();
        m_mcoord[1] = NULL;
    }
    else
    {
        for (j=0; m_mcoord[j]; j++);	// Get count.
        MetalCoord** nmc = new MetalCoord*[j+2];
        for (i=0; i<j; i++) nmc[i] = m_mcoord[i];
        nmc[j] = new MetalCoord();
        nmc[j+1] = NULL;
        delete[] m_mcoord;
        m_mcoord = nmc;
    }

    if (!metals)
    {
        metals = new Atom*[2];
        metals[0] = metal;
        metals[1] = NULL;
    }
    else
    {
        for (k=0; metals[k]; k++);	// Get count.
        Atom** nma = new Atom*[k+2];
        for (i=0; i<k; i++) nma[i] = metals[i];
        nma[k] = metal;
        nma[k+1] = NULL;
        delete[] metals;
        metals = nma;
    }

    m_mcoord[j]->metal = metal;
    m_mcoord[j]->coord_res = new AminoAcid*[residues+2];
    m_mcoord[j]->coord_atoms = new Atom*[residues+2];

    int maxres = 0, minres = 0;
    for (i=0; i<residues; i++)
    {
        if (resnos[i] > maxres) maxres = resnos[i];
        if (!minres || resnos[i] < minres) minres = resnos[i];

        char buffer[256];
        sprintf(buffer, "REMARK 800 SITE MCOORD %d %s\n", resnos[i], metal->name);
        add_remark(buffer);

        m_mcoord[j]->coord_res[i] = get_residue(resnos[i]);
        if (!m_mcoord[j]->coord_res[i])
        {
            cout << "Attempt to bind metal to residue " << resnos[i] << " not found in protein!" << endl;
            throw 0xbad12e5d;
        }
        m_mcoord[j]->coord_res[i]->m_mcoord = m_mcoord[j];
        if (i < res_anames.size())
        {
            m_mcoord[j]->coord_atoms[i] = m_mcoord[j]->coord_res[i]->get_atom(res_anames[i].c_str());
            if (!m_mcoord[j]->coord_atoms[i])
            {
                cout << "Attempt to bind metal to " << resnos[i] << ":" << res_anames[i] << " not found in protein!" << endl;
                throw 0xbada70b;
            }
        }
    }
    m_mcoord[j]->coord_res[residues] = NULL;
    m_mcoord[j]->coord_atoms[residues] = NULL;

    add_remark("REMARK 800\n");

    // Get the plane of the coordinating atoms, then get the normal.
    Point ptarr[3] =
    {
        m_mcoord[j]->coord_atoms[0]->get_location(),
        m_mcoord[j]->coord_atoms[1]->get_location(),
        m_mcoord[j]->coord_atoms[2]->get_location()
    };
    Point coordcen = average_of_points(ptarr, 3);
    SCoord normal = compute_normal(ptarr[0], ptarr[1], ptarr[2]);
    normal.r = 3;
    Point pnormal = coordcen.add(normal), pantinormal = coordcen.subtract(normal);
    int nc = 0, anc = 0;

    // Iterate from minres to maxres, counting how many CA atoms are on each side of the plane and are within a threshold distance.
    for (i = minres; i <= maxres; i++)
    {
        Atom* la = ca[i];
        if (!la) continue;
        Point lpt = la->get_location();
        float nr = lpt.get_3d_distance(pnormal);
        float anr = lpt.get_3d_distance(pantinormal);

        // Don't sweat residues on the opposite side of the binding pocket if the coord atoms are on different helices.
        if (nr > 7 || anr > 7) continue;

        // If equidistant, don't count it.
        if (nr < anr) nc++;
        if (nr > anr) anc++;
    }

    // Choose the side of the plane with the fewest CA atoms and set gmtgt about 1A in that direction of plane center.
    // If the two sides are the same, set gmtgt to coordcen.
    // cout << "nc " << nc << " | anc " << anc << endl;
    if (nc < anc)
    {
        normal.r = 7;
        gmtgt = pnormal.add(normal);
        metal->move(gmtgt);
        gmtgt = pnormal;

    }
    else if (nc > anc)
    {
        normal.r = 7;
        gmtgt = pantinormal.subtract(normal);
        metal->move(gmtgt);
        gmtgt = pantinormal;
    }
    else gmtgt = coordcen;

    // Create an array of Molecules containing the coordinating residues and surrounding residues, with the metal as its own Molecule.
    Atom* ma[2];
    ma[0] = metal;
    ma[1] = NULL;
    Molecule m("Metal", ma);
    m.movability = MOV_ALL;
    Molecule* lmols[maxres-minres+8];
    Molecule* lbkg[maxres-minres+8];
    int lmolc=0, lbkgc=0;
    lmols[lmolc++] = &m;
    lbkg[lbkgc++] = &m;
    for (n=0; n<residues; n++)
    {
        if (n>0 && resnos[n] == resnos[n-1]) continue;
        AminoAcid* aa = get_residue(i);
        if (aa)
        {
            lmols[lmolc++] = aa;
        }
    }
    for (i=minres; i<=maxres; i++)
    {
        AminoAcid* aa = get_residue(i);
        if (aa)
        {
            lbkg[lbkgc++] = aa;
            if (i >= minres && i <= maxres)
                aa->m_mcoord = m_mcoord[j];		// Sets the coordinating residues, and all in between residues, to backbone immovable.
        }
    }
    lmols[lmolc] = nullptr;
    lbkg[lbkgc] = nullptr;

    // Make sure to set the coordinating residues so that their side chains are flexible.
    for (i=0; m_mcoord[j]->coord_res[i]; i++)
        m_mcoord[j]->coord_res[i]->movability = MOV_FLEXONLY;    

    // Move the metal to the new center of all coordinating residues' CA.
    Point pt4avg[residues+2];
    l=0;
    for (n=0; n<residues; n++)
    {
        Point respt = get_atom_location(resnos[n], "CA");

        if (n>0 && resnos[n] == resnos[n-1]) continue;
        else pt4avg[l++] = respt;
    }
    Point ptmtl = average_of_points(pt4avg, l);
    metal->move(ptmtl);

    // Flex the side chains to all be close to the alpha center.
    int iter;
    for (iter=0; iter<50; iter++)
    {
        for (i=0; i<residues; i++)
        {
            AminoAcid* aa = get_residue(resnos[i]);
            if (aa)
            {
                Bond** bb = aa->get_rotatable_bonds();
                if (bb)
                {
                    float rad = 0, step = 10*fiftyseventh;
                    float bestrad, bestr, r;

                    for (l=0; bb[l]; l++)
                    {
                        bestrad = 0;
                        bestr = 999999;
                        for (; rad < M_PI*2; rad += step)
                        {
                            bb[l]->rotate(step);
                            r = 0;
                            float clashes = 0;
                            for (n=0; n<residues; n++)
                            {
                                // if (resnos[n] == resnos[i]) continue;
                                // r += get_atom_location(resnos[n], res_anames[n].c_str()).get_3d_distance(get_atom_location(resnos[i], res_anames[i].c_str()));
                                r += fabs(get_atom_location(resnos[n], res_anames[n].c_str()).get_3d_distance(ptmtl) - 2);      // VERY rough approximation.
                                clashes += get_residue(resnos[n])->get_intermol_clashes(lbkg);
                            }   // for n

                            /*cout << iter << " " << *aa << ":"
                            	 << bb[l]->atom->name << "-" << bb[l]->btom->name
                            	 << " " << rad*fiftyseven << "deg, r=" << r
                            	 << ", clash=" << clashes << endl;*/

                            r += 0.05*clashes;

                            if (r < bestr)
                            {
                                bestrad = rad;
                                bestr = r;
                            }
                        }

                        if (bestrad) bb[l]->rotate(bestrad);

                    }		// for (l=0; bb[l]; l++)
                }		// if (bb)
            }		// if (aa)
        }		// for (i=0; i<residues; i++)
    }		// for iter

    // Multimol conform the array.
    gmprot = this;
    Molecule::conform_molecules(lmols, lbkg, 50); // &ext_mtl_coord_cnf_cb);
    // metal->move(ptmtl);

    // Flex the side chains to all be close to the metal.
    ptmtl = metal->get_location();
    for (iter=0; iter<50; iter++)
    {
        for (i=0; i<residues; i++)
        {
            AminoAcid* aa = get_residue(resnos[i]);
            if (aa)
            {
                Bond** bb = aa->get_rotatable_bonds();
                if (bb)
                {
                    float rad = 0, step = 10*fiftyseventh;
                    float bestrad, bestr, r;

                    for (l=0; bb[l]; l++)
                    {
                        bestrad = 0;
                        bestr = 999999;
                        for (; rad < M_PI*2; rad += step)
                        {
                            bb[l]->rotate(step);
                            r = 0;
                            float clashes = 0;
                            for (n=0; n<residues; n++)
                            {
                                r += fabs(get_atom_location(resnos[n], res_anames[n].c_str()).get_3d_distance(ptmtl) - 2);      // VERY rough approximation.
                                clashes += get_residue(resnos[n])->get_intermol_clashes(lbkg);
                            }   // for n

                            /*cout << iter << " " << *aa << ":"
                            	 << bb[l]->atom->name << "-" << bb[l]->btom->name
                            	 << " " << rad*fiftyseven << "deg, r=" << r
                            	 << ", clash=" << clashes << endl;*/

                            r += 0.05*clashes;

                            if (r < bestr)
                            {
                                bestrad = rad;
                                bestr = r;
                            }
                        }

                        if (bestrad) bb[l]->rotate(bestrad);

                    }		// for (l=0; bb[l]; l++)
                }		// if (bb)
            }		// if (aa)
        }		// for (i=0; i<residues; i++)
    }		// for iter

    // Move the metal to the new center of all coordinating atoms.
    #if 0
    l=0;
    for (n=0; n<residues; n++)
    {
        Point respt = get_atom_location(resnos[n], res_anames[n].c_str());

        if (n>0 && resnos[n] == resnos[n-1])
        {
            pt4avg[l-1].x = (pt4avg[l-1].x + respt.x)/2;
            pt4avg[l-1].y = (pt4avg[l-1].y + respt.y)/2;
            pt4avg[l-1].z = (pt4avg[l-1].z + respt.z)/2;
        }
        else
            pt4avg[l++] = respt;
    }
    ptmtl = average_of_points(pt4avg, l);
    #endif
    metal->move(ptmtl);

    // Set the coordinating residues' sidechains to immovable.
    for (i=0; m_mcoord[j]->coord_res[i]; i++)
        m_mcoord[j]->coord_res[i]->movability = MOV_NONE;

    m_mcoord[j]->locked = true;
    return m_mcoord[j];
}

void Protein::mtl_coord_cnf_cb(int iter)
{
    int i;
    for (i=0; m_mcoord[i]; i++)
    {
        // SCoord delta(m_mcoord[i]->coord_atom_avg_loc().subtract(m_mcoord[i]->metal->get_location()));
        SCoord delta(gmtgt.subtract(m_mcoord[i]->metal->get_location()));
        delta.r *= 0.1;
        m_mcoord[i]->metal->move_rel(&delta);
    }
}

float Protein::get_helix_orientation(int startres, int endres)
{
    int i, j;

    AminoAcid* aa;
    Atom* a;
    int rescount = (endres - startres)+1;
    int acount = 3*rescount;
    if (acount < 11)
    {
        cout << "Helix too short for determining orientation. " << startres << "-" << endres << endl;
        return 0;
    }
    Point pt, ptarr[acount+8];
    for (i=0; i<=rescount; i++)
    {
        /*Star s;
        s.pprot = this;
        cout << i << ":" << hex << s.n << dec << " " << flush;*/
        int resno = i+startres;
        // cout << resno << " ";
        aa = get_residue(resno);
        if (aa) a = aa->get_atom("N");
        if (a) pt = a->get_location();
        ptarr[i*3] = pt;

        if (aa) a = aa->get_atom("CA");
        if (a) pt = a->get_location();
        ptarr[i*3+1] = pt;

        if (aa) a = aa->get_atom("C");
        if (a) pt = a->get_location();
        ptarr[i*3+2] = pt;
    }

    // Take a running average of 10-atom blocks, then measure the average radians from vertical for imaginary lines between consecutive averages.
    int bcount = acount-10;
    Point blkavg[bcount+8];
    float retval = 0;
    j=0;

    for (i=0; i<bcount; i++)
    {
        blkavg[i] = average_of_points(&ptarr[i], 10);
        if (i>0)
        {
            SCoord v(blkavg[i].subtract(blkavg[i-1]));
            retval += v.theta;
            j++;
        }
    }

    return retval/j;
}

float Protein::orient_helix(int startres, int endres, int stopat, float angle, int iters)
{
    AminoAcid* aa = get_residue(startres-1);
    float n_am = 0.1, ca_am = 0.1;
    int iter;
    float ha;

    ha = get_helix_orientation(startres, endres);
    for (iter = 0; iter < iters; iter++)
    {
        rotate_backbone_partial(startres, stopat, N_asc, n_am);
        float nha = get_helix_orientation(startres, endres);

        if (fabs(nha-angle) <= fabs(ha-angle))
        {
            ha = nha;
            n_am *= 1.1;
            cout << "+";
        }
        else
        {
            rotate_backbone_partial(startres, stopat, N_asc, -n_am);
            n_am *= -0.75;
            cout << "x";
        }

        rotate_backbone_partial(startres, stopat, CA_asc, ca_am);
        nha = get_helix_orientation(startres, endres);

        if (fabs(nha-angle) <= fabs(ha-angle))
        {
            ha = nha;
            ca_am *= 1.1;
            cout << "+";
        }
        else
        {
            rotate_backbone_partial(startres, stopat, CA_asc, -ca_am);
            ca_am *= -0.75;
            cout << "x";
        }
    }
    cout << " ";

    return ha;
}

SCoord Protein::get_region_axis(int startres, int endres)
{
    int rglen = endres-startres;
    if (rglen < 4) throw 0xbadc0de;     // TODO
    Point N[4], C[4];
    int i;
    for (i=0; i<4; i++)
    {
        AminoAcid* aa = get_residue(startres + i);
        if (!aa) throw 0xdeadac1d;      // TODO
        N[i] = aa->get_CA_location();

        aa = get_residue(endres - i);
        if (!aa) throw 0xdeadac1d;      // TODO
        C[i] = aa->get_CA_location();
    }

    Point nterm = average_of_points(N, 4), cterm = average_of_points(C, 4);

    return (SCoord)(cterm.subtract(nterm));
}

void Protein::set_region(std::string rgname, int start, int end)
{
    int i;
    for (i=0; i<PROT_MAX_RGN; i++) if (!regions[i].start || !strcmp(regions[i].name.c_str(), rgname.c_str())) break;
    if (i >= PROT_MAX_RGN) return;		// Nope.

    regions[i].name = rgname;
    regions[i].start = start;
    regions[i].end = end;
    regions_from = rgn_manual;
}

Region Protein::get_region(const std::string rgname)
{
    int i;
    for (i=0; i<PROT_MAX_RGN; i++) if (regions[i].name == rgname) return regions[i];
    return Region();
}


int Protein::get_region_start(const std::string name)
{
    Region rgn = get_region(name);
    return rgn.start;
}

int Protein::get_region_end(const std::string name)
{
    Region rgn = get_region(name);
    return rgn.end;
}

Point Protein::get_region_center(int startres, int endres)
{
    int rglen = endres-startres;
    Point range[rglen+4];

    int i;
    int samples = 0;
    for (i=0; i<rglen; i++)
    {
        // This is slow but that's okay.
        AminoAcid* aa = get_residue(startres+i);
        if (!aa) continue;
        range[samples] = aa->get_barycenter();
        samples++;
    }

    return average_of_points(range, samples);
}

void Protein::move_piece(int start_res, int end_res, Point new_center)
{
    Point old_center = get_region_center(start_res, end_res);
    SCoord move_amt = new_center.subtract(old_center);

    move_piece(start_res, end_res, move_amt);
}

void Protein::move_piece(int start_res, int end_res, SCoord move_amt)
{
    save_undo_state();

    int i;
    for (i=start_res; i<=end_res; i++)
    {
        AminoAcid* aa = get_residue(i);
        if (!aa) continue;
        MovabilityType mov = aa->movability;
        aa->movability = MOV_ALL;
        aa->aamove(move_amt);
        aa->movability = mov;
        set_clashables(i);
    }
}

LocRotation Protein::rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res)
{
    Point pivot = pivot_res ? get_residue(pivot_res)->get_CA_location() : get_region_center(start_res, end_res);
    Point align = get_residue(align_res)->get_CA_location();
    Rotation rot = align_points_3d(&align, &align_target, &pivot);
    return rotate_piece(start_res, end_res, rot, pivot_res);
}

LocRotation Protein::rotate_piece(int start_res, int end_res, Rotation rot, int pivot_res)
{
    AminoAcid* aa = pivot_res ? get_residue(pivot_res) : nullptr;
    Point pivot = aa ? aa->get_CA_location() : get_region_center(start_res, end_res);

    return rotate_piece(start_res, end_res, pivot, rot.v, rot.a);
}

LocRotation Protein::rotate_piece(int start_res, int end_res, Point pivot, SCoord axis, float theta)
{
    save_undo_state();

    LocatedVector lv(axis);
    lv.origin = pivot;
    int i;
    for (i=start_res; i<=end_res; i++)
    {
        AminoAcid* aa = get_residue(i);
        if (!aa) continue;
        MovabilityType mov = aa->movability;
        aa->movability = MOV_ALL;
        aa->rotate(lv, theta);
        aa->ensure_pi_atoms_coplanar();
        aa->movability = mov;
        set_clashables(i);
    }

    LocRotation retval(lv);
    retval.a = theta;
    return retval;
}



Point Protein::find_loneliest_point(Point cen, Point sz)
{
    if (!residues) return cen;

    float x, y, z, xp, yp, zp, xa, ya, za, r, bestr = 0, step = 0.25;
    int i;
    Point retval = cen;

    sz.x /= 2;
    sz.y /= 2;
    sz.z /= 2;

    /*if (fabs(sz.x) > 4) sz.x = 4;
    if (fabs(sz.y) > 4) sz.y = 4;
    if (fabs(sz.z) > 4) sz.z = 4;*/

    for (x = -sz.x; x <= sz.x; x += step)
    {
        xp = x / sz.x;
        xp *= xp;
        xa = cen.x + x;
        for (y = -sz.y; y <= sz.y; y += step)
        {
            yp = y / sz.y;
            yp *= yp;
            ya = cen.y + y;
            for (z = -sz.z; z <= sz.z; z += step)
            {
                zp = z / sz.z;
                zp *= zp;
                za = cen.z + z;
                r = sqrt(xp+yp+zp);
                if (r > 1) continue;

                Point maybe(xa, ya, za);
                float minr = Avogadro;

                for (i=0; residues[i]; i++)
                {
                    Atom* a = residues[i]->get_nearest_atom(maybe);
                    if (a)
                    {
                        r = a->get_location().get_3d_distance(maybe);
                        if (r < minr) minr = r;
                    }
                }

                if (minr < 1e9 && minr > bestr)
                {
                    retval = maybe;
                    bestr = minr;
                }
            }
        }
    }

    return retval;
}

Point Protein::estimate_pocket_size(std::vector<AminoAcid*> ba)
{
    int i, n = ba.size();
    if (!n) return Point();
    float cx, cy, cz;

    cx = cy = cz = 0;
    for (i=0; i<n; i++)
    {
        Point pt = ba[i]->get_atom_location("CA");
        cx += pt.x;
        cy += pt.y;
        cz += pt.z;
    }

    Point center(cx/n, cy/n, cz/n);

    float sx, sy, sz, wx, wy, wz;
    sx = sy = sz = wx = wy = wz = 0;
    for (i=0; i<n; i++)
    {
        Point pt = ba[i]->get_atom_location("CA").subtract(center);
        float mag = pt.magnitude();
        mag -= 0.666 * ba[i]->get_reach();
        pt.scale(mag);
        float lwx = (fabs(pt.x) / sqrt(pt.y*pt.y + pt.z*pt.z)) / mag;
        float lwy = (fabs(pt.y) / sqrt(pt.x*pt.x + pt.z*pt.z)) / mag;
        float lwz = (fabs(pt.z) / sqrt(pt.x*pt.x + pt.y*pt.y)) / mag;

        sx += lwx * fabs(pt.x);
        sy += lwy * fabs(pt.y);
        sz += lwz * fabs(pt.z);

        wx += lwx;
        wy += lwy;
        wz += lwz;
    }

    Point size(sx/wx, sy/wy, sz/wz);

    return size;
}

Molecule** Protein::all_residues_as_molecules()
{
    if (!residues) return nullptr;
    Molecule** retval = new Molecule*[get_seq_length()+4];

    int i;
    for (i=0; residues[i]; i++) retval[i] = reinterpret_cast<Molecule*>(residues[i]);
    retval[i] = nullptr;

    return retval;
}

Molecule** Protein::all_residues_as_molecules_except(Molecule** mm)
{
    if (!residues) return nullptr;
    Molecule** retval = new Molecule*[get_seq_length()+4];

    int i, j, k;
    k=0;
    for (i=0; residues[i]; i++)
    {
        if (mm)
        {
            for (j=0; mm[j]; j++)
            {
                if (mm[j] == residues[i]) goto _exclude;
            }
        }
        retval[k++] = reinterpret_cast<Molecule*>(residues[i]);
        _exclude:
        ;
    }
    retval[k] = nullptr;

    return retval;
}

bool Protein::disulfide_bond(int resno1, int resno2)
{
    AminoAcid* res1, *res2;

    save_undo_state();

    res1 = get_residue(resno1);
    res2 = get_residue(resno2);

    if (!res1 || !res2) return false;

    int i, j, k, l;
    Atom *S1 = nullptr, *S2 = nullptr, *H1 = nullptr, *H2 = nullptr;
    bool result = false;

    j = res1->get_atom_count();
    l = res2->get_atom_count();
    for (i=0; i<j; i++)
    {
        S1 = res1->get_atom(i);
        if (!S1) continue;
        if (S1->get_Z() == 16)
        {
            H1 = S1->is_bonded_to("H");
            if (H1)
            {
                for (k=0; k<l; k++)
                {
                    S2 = res2->get_atom(k);
                    if (!S2) continue;
                    if (S2->get_Z() == 16)
                    {
                        H2 = S2->is_bonded_to("H");
                        if (H2)
                        {
                            float r = S1->get_location().get_3d_distance(S2->get_location());
                            if (fabs(r - 2.07) < 0.25)
                            {
                                res1->delete_atom(H1);
                                res2->delete_atom(H2);
                                if (S1->bond_to(S2, 1))
                                {
                                    connections.push_back(S1->get_bond_between(S2));
                                    result = true;
                                }
                                goto _next_S1;
                            }
                        }
                    }
                }
            }
        }
        _next_S1:
        ;
    }

    return result;
}

void Protein::upright()
{
    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;
    Point oldcen = get_region_center(1, 9999).negate();
    last_uprighted_xform = oldcen;
    move_piece(1, 9999, Point(0,0,0));

    Point extracellular[256], cytoplasmic[256];
    int i, j, exr_n=0, cyt_n=0;

    for (i=1; i<=7; i++)
    {
        int sr = get_region_start((std::string)"TMR" + std::to_string(i));
        if (!sr) continue;
        int er = get_region_end((std::string)"TMR" + std::to_string(i));

        for (j=0; j<4; j++)
        {
            if (i & 1)			// TMR1, TMR3, TMR5, TMR7 begin on the extracellular side and descend.
            {
                extracellular[exr_n++] = get_atom_location(sr+j, "CA");
                cytoplasmic[cyt_n++] = get_atom_location(er-j, "CA");
            }
            else				// TMR2, TMR4, TMR6 ascend from the cytoplasmic side.
            {
                cytoplasmic[cyt_n++] = get_atom_location(sr+j, "CA");
                extracellular[exr_n++] = get_atom_location(er-j, "CA");
            }
        }
    }

    if (!exr_n || !cyt_n) throw 0xbad7312;

    Point exrdir = average_of_points(extracellular, exr_n);
    Point cytdir = average_of_points(cytoplasmic, cyt_n);

    Rotation rot = align_points_3d(&exrdir, new Point(0,1e6,0), &cytdir);

    rotate_piece(1, 9999, rot, 0);
    last_uprighted_A.v = rot.v;
    last_uprighted_A.a = rot.a;
    last_uprighted_A.origin = get_region_center(1, 9999);

    // Rotate to place TMR4 in the +Z direction relative to TMR1.
    int sr = get_region_start("TMR4");
    if (sr)
    {
        int er = get_region_end("TMR4");

        Point tmr1[64], tmr4[64];
        int tmr1_n=0, tmr4_n=0;

        for (i=sr; i<=er; i++)
        {
            tmr4[tmr4_n++] = get_atom_location(i, "CA");
        }

        sr = get_region_start("TMR1");
        er = get_region_end("TMR1");

        for (i=sr; i<=er; i++)
        {
            tmr1[tmr1_n++] = get_atom_location(i, "CA");
        }

        Point tmr1dir = average_of_points(tmr1, tmr1_n);
        Point tmr4dir = average_of_points(tmr4, tmr4_n);


        tmr1dir.y = tmr4dir.y = 0;

        rot = align_points_3d(&tmr4dir, new Point(0,0,1e9), &tmr1dir);

        rotate_piece(1, 9999, rot, 0);
        last_uprighted_B.v = rot.v;
        last_uprighted_B.a = rot.a;
        last_uprighted_B.origin = get_region_center(1, 9999);
    }

    mass_undoable = wmu;
}

void Protein::homology_conform(Protein* target, Protein* reference)
{
    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;

    if (!reference) reference = this;

    // Check that both proteins have TM helices and BW numbers set. If not, error out.
    if (!get_region_start("TMR6") || !target->get_region_start("TMR6")) throw 0xbadbeb7;
    if (!Ballesteros_Weinstein[3] || !target->Ballesteros_Weinstein[3]) throw 0xbadbeb7;
    if (!reference->get_region_start("TMR6")) throw 0xbadbeb7;
    if (!reference->Ballesteros_Weinstein[3]) throw 0xbadbeb7;

    upright();
    target->upright();

    // Get the average location delta for all CA atoms in the TM helices. Match them by BW number.
    // Include only BW numbers that are inside a TMR for both proteins.
    int hxno, resno1, resno2;
    char buffer[256];
    Point xform_delta(0,0,0), center(0,0,0);
    int count = 0;

    std::vector<int> resnos1, resnos2;

    for (hxno = 1; hxno <= 7; hxno++)
    {
        sprintf(buffer, "TMR%d", hxno);
        int rgend1 = get_region_end(buffer);
        int rgstart2 = target->get_region_start(buffer);
        int rgend2 = target->get_region_end(buffer);
        int bw50a = get_bw50(hxno), bw50b = target->get_bw50(hxno);
        for (resno1 = get_region_start(buffer); resno1 <= rgend1; resno1++)
        {
            int i = resno1 - bw50a;
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                Point caloc = get_atom_location(resno1, "CA");
                Point ptdiff = caloc.subtract(target->get_atom_location(resno2, "CA"));
                xform_delta = xform_delta.add(ptdiff);
                center = center.add(caloc);
                count++;

                resnos1.push_back(resno1);
                resnos2.push_back(resno2);
            }
        }
    }
    if (count)
    {
        xform_delta.scale(xform_delta.magnitude() / count);
        center.scale(center.magnitude() / count);
    }

    // Transform the target to bring its TM center to coincide with that of the current protein.
    SCoord move_amt = xform_delta;
    target->move_piece(1, 9999, move_amt);

    // Get the average necessary rotation, about the +Y axis centered on the TM center, to match
    // the TM CA atoms as closely as possible.
    int i, l;
    float theta = 0;
    Point axis(0,1,0);
    count = 0;
    for (i=0; i<resnos1.size(); i++)
    {
        resno1 = resnos1[i];
        resno2 = resnos2[i];
        Point pt1 = get_atom_location(resno1, "CA"), pt2 = target->get_atom_location(resno2, "CA");
        pt1.y = pt2.y = 0;
        Rotation rot = align_points_3d(pt2, pt1, center);
        Point rotv = rot.v;
        if (rotv.y < 0) theta -= rot.a;
        else theta += rot.a;
        count++;
    }

    if (count) theta /= count;

    // Perform the rotation.
    target->rotate_piece(1, 9999, center, axis, theta);
    #if _dbg_homology
    cout << "Rotated about Y axis " << theta*fiftyseven << "deg." << endl;
    #endif

    if (reference != this)
    {
        count = 0;
        xform_delta = center = Point(0,0,0);
        resnos1.clear();
        resnos2.clear();

        for (hxno = 1; hxno <= 7; hxno++)
        {
            sprintf(buffer, "TMR%d", hxno);
            int rgend1 = get_region_end(buffer);
            int rgstart2 = reference->get_region_start(buffer);
            int rgend2 = reference->get_region_end(buffer);
            int bw50a = get_bw50(hxno), bw50b = reference->get_bw50(hxno);
            for (resno1 = get_region_start(buffer); resno1 <= rgend1; resno1++)
            {
                int i = resno1 - bw50a;
                resno2 = bw50b + i;
                if (resno2 >= rgstart2 && resno2 <= rgend2)
                {
                    Point caloc = get_atom_location(resno1, "CA");
                    Point ptdiff = caloc.subtract(reference->get_atom_location(resno2, "CA"));
                    xform_delta = xform_delta.add(ptdiff);
                    center = center.add(caloc);
                    count++;

                    resnos1.push_back(resno1);
                    resnos2.push_back(resno2);
                }
            }
        }
        if (count)
        {
            xform_delta.scale(xform_delta.magnitude() / count);
            center.scale(center.magnitude() / count);
        }

        move_amt = xform_delta;
        reference->move_piece(1, 9999, move_amt);

        // Get the average necessary rotation, about the +Y axis centered on the TM center, to match
        // the TM CA atoms as closely as possible.
        int i;
        float theta = 0;
        Point axis(0,1,0);
        count = 0;
        for (i=0; i<resnos1.size(); i++)
        {
            resno1 = resnos1[i];
            resno2 = resnos2[i];
            Point pt1 = get_atom_location(resno1, "CA"), pt2 = reference->get_atom_location(resno2, "CA");
            pt1.y = pt2.y = 0;
            Rotation rot = align_points_3d(pt2, pt1, center);
            Point rotv = rot.v;
            if (rotv.y < 0) theta -= rot.a;
            else theta += rot.a;
            count++;
        }

        if (count) theta /= count;

        // Perform the rotation.
        reference->rotate_piece(1, 9999, center, axis, theta);
    }

    // Find the rotations and transformations for each TM region to bring its CA atoms as close as
    // possible to those of the target.
    for (hxno = 1; hxno <= 7; hxno++)
    {
        sprintf(buffer, "TMR%d", hxno);
        int rgstart0 = get_region_start(buffer);
        int rgend0 = get_region_end(buffer);
        int rgstart1 = reference->get_region_start(buffer);
        int rgend1 = reference->get_region_end(buffer);
        int rgstart2 = target->get_region_start(buffer);
        int rgend2 = target->get_region_end(buffer);
        #if _dbg_homology
        cout << "Region " << hxno << " target " << rgstart2 << "-" << rgend2 << ", reference " << rgstart1 << "-" << rgend1 << endl;
        #endif

        int bw50a = reference->get_bw50(hxno), bw50b = target->get_bw50(hxno);

        #if homology_phi_psi_rotations
        // Do the phi and psi rotations.
        float theta_carry = 0;
        for (resno1 = rgstart1; resno1 <= rgend1; resno1++)
        {
            i = resno1 - bw50a;
            int resno0 = get_bw50(hxno) + i;
            AminoAcid* aa = get_residue(resno0);
            if (!aa) continue;

            i = resno1 - bw50a;
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                AminoAcid* aa1 = reference->get_residue(resno1);
                if (!aa1) continue;
                AminoAcid* aa2 = target->get_residue(resno2);
                if (!aa2) continue;

                float phi1 = aa1->get_phi(), phi2 = aa2->get_phi();
                if (phi1 >= M_PI) phi1 -= M_PI*2;
                if (phi2 >= M_PI) phi2 -= M_PI*2;
                float theta = phi2 - phi1;

                bond_rotation_fail_reason phi_fail = aa->rotate_phi(theta + theta_carry);
                if (phi_fail != bf_none) theta_carry += theta;
                else theta_carry = 0;

                #if _dbg_homology
                cout << resno0 << " (t " << resno2 << ", r " << resno1 << ") phi " << (theta*fiftyseven) << " (" << (phi2*fiftyseven) << " - " << (phi1*fiftyseven) << ") deg. " << phi_fail << "." << endl;
                #endif

                AminoAcid* aaprev = aa->get_prev();
                if (aaprev)
                {
                    Point* nhca = aaprev->predict_next_NHCA();
                    aa->attach_to_prediction(nhca);
                }
 
                float psi1 = aa1->get_psi(), psi2 = aa2->get_psi();
                if (psi1 >= M_PI) psi1 -= M_PI*2;
                if (psi2 >= M_PI) psi2 -= M_PI*2;
                theta = psi2 - psi1;

                bond_rotation_fail_reason psi_fail = aa->rotate_psi(theta + theta_carry);
                if (psi_fail != bf_none) theta_carry += theta;
                else theta_carry = 0;

                #if _dbg_homology
                cout << resno0 << " (t " << resno2 << ", r " << resno1 << ") psi " << (theta*fiftyseven) << " (" << (psi2*fiftyseven) << " - " << (psi1*fiftyseven) << ") deg. " << psi_fail << "." << endl;
                #endif

                /*AminoAcid* aanext = aa->get_next();
                if (!aanext) continue;
                Point* nhca = aa->predict_next_NHCA();
                aanext->attach_to_prediction(nhca);*/
            }
        }
        #endif

        // Compute the TM region transformation.
        Point rcen1(0,0,0);
        Point rcen2(0,0,0);
        count = 0;
        for (resno1 = rgstart1; resno1 <= rgend1; resno1++)
        {
            i = resno1 - bw50a;
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                rcen1 = rcen1.add(reference->get_atom_location(resno1, "CA"));
                rcen2 = rcen2.add(target->get_atom_location(resno2, "CA"));
                count++;
            }
        }

        if (count)
        {
            rcen1.scale(rcen1.magnitude()/count);
            rcen2.scale(rcen2.magnitude()/count);
        }

        // Perform the TM region transformation.
        Point transformation = rcen2.subtract(rcen1);
        move_piece(rgstart0, rgend0, (SCoord)transformation);
        #if _dbg_homology
        cout << "Region " << hxno << " (" << rgstart0 << "-" << rgend0 << ")" << " transformed by " << transformation.magnitude() << "A" << endl;
        #endif

        // Compute the TM region rotation.
        Point axis(0,0,0);
        Point rcen = get_region_center(rgstart2, rgend2);
        float theta = 0;
        count = 0;
        for (resno1 = rgstart1; resno1 <= rgend1; resno1++)
        {
            i = resno1 - bw50a;
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                Point caloc1 = reference->get_atom_location(resno1, "CA"),
                      caloc2 = target->get_atom_location(resno2, "CA");
                Rotation rot = align_points_3d(caloc1, caloc2, rcen);

                axis = axis.add(rot.v);
                theta += rot.a;
                count++;
            }
        }

        if (count)
        {
            axis.scale(axis.magnitude()/count);
            theta/=count;
        }

        // Perform the TM region rotation.
        rotate_piece(rgstart0, rgend0, rcen, axis, theta);
        #if _dbg_homology
        cout << "Region " << hxno << " rotated " << theta*fiftyseven << "deg." << endl;
        #endif


        #if homology_long_axis_rotations
        // Compute the TM region long axes.
        Point rgn_begin(0,0,0), rgn_end(0,0,0);
        l = 0;
        bool starting = true;
        for (resno1 = rgstart1; resno1 <= rgend1; resno1++)
        {
            AminoAcid* aa = get_residue(resno1);
            if (!aa) continue;

            i = resno1 - bw50a;
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                Point caloc2 = target->get_atom_location(resno2, "CA");
                if (starting) rgn_begin = rgn_begin.add(caloc2);
                else rgn_end = rgn_end.add(caloc2);

                l++;
                if (l >= 4)
                {
                    if (starting)
                    {
                        rgn_begin.scale(rgn_begin.magnitude()/l);
                        resno1 = rgend1 - 4;
                        l = 0;
                        starting = false;
                    }
                    else
                    {
                        rgn_end.scale(rgn_end.magnitude()/l);
                        break;
                    }
                }
            }
        }

        if (starting || !l) continue;

        // Compute the TM region long-axis rotation.
        axis = rgn_end.subtract(rgn_begin);
        theta = 0;
        count = 0;
        for (resno1 = rgstart0; resno1 <= rgend0; resno1++)
        {
            i = resno1 - get_bw50(hxno);
            resno2 = bw50b + i;
            if (resno2 >= rgstart2 && resno2 <= rgend2)
            {
                Point caloc1 = get_atom_location(resno1, "CA"),
                      caloc2 = target->get_atom_location(resno2, "CA");

                theta += find_angle_along_vector(caloc1, caloc2, rgn_begin, axis);
                count++;
            }
        }

        if (count)
        {
            theta/=count;
        }

        // Perform the TM region long-axis rotation.
        rotate_piece(rgstart0, rgend0, rgn_begin, axis, theta);
        #if _dbg_homology
        cout << "Region " << hxno << " long-axis rotated " << theta*fiftyseven << "deg." << endl;
        #endif
        #endif

        // if (helix) backconnect(rgend1+1, helix-1);
    }

    #if _dbg_homology
    cout << "Minimizing clashes..." << endl;
    #endif
    get_internal_clashes(1, 9999, true, 25);

    #if homology_region_optimization
    for (l=0; l<20; l++)
    {
        #if _dbg_homology
        cout << endl << flush;
        #endif

        float total_motion = 0;
        for (hxno = 1; hxno <= 7; hxno++)
        {
            sprintf(buffer, "TMR%d", hxno);
            int rgstart = get_region_start(buffer);
            int rgend = get_region_end(buffer);

            SCoord transform;
            Rotation rotation;

            region_optimal_positioning(rgstart, rgend, &transform, &rotation, nullptr);

            total_motion += transform.r + rotation.a*fiftyseven;

            #if _dbg_homology
            Point ptx = transform;
            cout << buffer << " should move " << ptx << " and rotate " << rotation << endl << flush;
            #endif

            transform.r *= 0.2;
            rotation.a *= 0.333;

            Point rgncenter = get_region_center(rgstart, rgend);

            move_piece(rgstart, rgend, transform);

            // rotate_piece(rgstart, rgend, rgncenter, rotation.v, rotation.a);
        }

        if (total_motion < 0.1) break;
    }
    #else
    for (l=0; l<50; l++)
    {
        int nummoved = 0;
        for (hxno = 1; hxno <= 7; hxno++)
        {
            sprintf(buffer, "TMR%d", hxno);
            int rgstart0 = get_region_start(buffer);
            int rgend0 = get_region_end(buffer);
            float threshold = clash_limit_per_aa * (rgend0 - rgstart0);
            float f = get_internal_clashes(rgstart0, rgend0, true, 10);

            #if _dbg_homology
            cout << "Helix " << hxno << " is clashy " << f << " of threshold " << threshold << endl;
            #endif

            if (f < threshold) continue;
            if (hxno == 3) continue;

            SCoord clashmov = last_int_clash_dir;
            clashmov.r = -0.003 * (f - threshold);
            #if _dbg_homology
            cout << "Motion is " << clashmov.r << " A." << endl;
            #endif
            if (fabs(clashmov.r) < 0.01) continue;

            #if _dbg_homology
            cout << "Moving " << rgstart0 << "-" << rgend0 << " by " << clashmov << endl;
            #endif

            move_piece(rgstart0, rgend0, clashmov);
            nummoved++;
        }

        if (!nummoved) break;
    }
    #endif

    mass_undoable = wmu;
}

void Protein::bridge(int resno1, int resno2)
{
    save_undo_state();
    AminoAcid *aa1 = get_residue(resno1), *aa2 = get_residue(resno2);
    if (!aa1 || !aa2) return;

    aa1->movability = aa2->movability = MOV_FLEXONLY;

    _INTERA_R_CUTOFF = aa1->get_CA_location().get_3d_distance(aa2->get_CA_location())
        + aa1->get_reach() + aa2->get_reach() + _DEFAULT_INTERA_R_CUTOFF;
    
    float r = aa1->get_CA_location().get_3d_distance(aa2->get_CA_location()) - aa1->get_reach() - aa2->get_reach();
    if (r > 2.5)
    {
        aa1->conform_atom_to_location(aa1->get_reach_atom()->name, aa2->get_reach_atom()->get_location());
        aa2->conform_atom_to_location(aa2->get_reach_atom()->name, aa1->get_reach_atom()->get_location());
        aa1->conform_atom_to_location(aa1->get_reach_atom()->name, aa2->get_reach_atom()->get_location());
        aa2->conform_atom_to_location(aa2->get_reach_atom()->name, aa1->get_reach_atom()->get_location());
    }

    Molecule* mols[3];
    mols[0] = aa1;
    mols[1] = aa2;
    mols[2] = nullptr;

    Molecule::conform_molecules(mols, 15);

    Molecule** mols2;

    mols2 = (Molecule**)get_residues_can_clash(resno1);
    Molecule::conform_molecules(mols, mols2, 10);

    mols2 = (Molecule**)get_residues_can_clash(resno2);
    Molecule::conform_molecules(mols, mols2, 10);

    Molecule::conform_molecules(mols, 15);

    aa1->movability = MOV_PINNED;
    aa2->movability = MOV_PINNED;

    _INTERA_R_CUTOFF = _DEFAULT_INTERA_R_CUTOFF;
}

void Protein::soft_iteration(std::vector<Region> l_soft_rgns, Molecule* ligand)
{
    save_undo_state();
    bool wmu = mass_undoable;
    mass_undoable = true;

    int l;
    float prebind;

    int sz = l_soft_rgns.size();
    if (sz)
    {
        for (l=0; l<sz; l++)
        {
            if (!l) prebind = (ligand ? get_intermol_binding(ligand)*soft_ligand_importance : 0) + get_internal_binding()*_kJmol_cuA;         // /'kʒmɑɫ.kju.ə/
            
            #if _dbg_soft
            cout << iter << ": from " << prebind;
            #endif

            int tweak = rand() % 6;
            float amount = nanf("unbiased");

            if (isnan(amount) || !amount) amount = frand(-1, 1);
            else
            {
                if (amount > 0) amount = frand(-amount*soft_bias_overlap, amount);
                else amount = frand(amount, -amount*soft_bias_overlap);
            }

            Point bary = get_region_center(1, 9999);
            Point ptrgn = get_region_center(l_soft_rgns[l].start, l_soft_rgns[l].end);
            SCoord r = ptrgn.subtract(bary);
            r.theta = 0;
            r.r = amount;
            Point pr1 = bary.add(r);
            Point pr2 = pr1;
            pr2.y += 20;
            SCoord normal = compute_normal(bary, pr1, pr2);
            normal.r = amount;
            SCoord alpha = get_region_axis(l_soft_rgns[l].start, l_soft_rgns[l].end);
            alpha.r = amount;
            Point rgncen = get_region_center(l_soft_rgns[l].start, l_soft_rgns[l].end);

            switch (tweak)
            {
                case 0:
                move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, r);
                break;

                case 1:
                move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, normal);
                break;

                case 2:
                move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, alpha);
                break;

                case 3:
                rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, alpha, amount/10);
                break;

                case 4:
                rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, r, amount/50);
                break;

                case 5:
                rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, normal, amount/30);
                break;

                default:
                ;
            }

            float postbind = (ligand ? get_intermol_binding(ligand)*soft_ligand_importance : 0) + get_internal_binding()*_kJmol_cuA;
            #if _dbg_soft
            cout << " to " << postbind;
            #endif

            switch (tweak)
            {
                case 0:
                if (ligand && postbind > prebind) g_rgnxform_r[l] += amount;
                else
                {
                    r.r = -r.r;
                    move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, r);
                    #if _dbg_soft
                    cout << " reverting r.";
                    #endif
                }
                break;

                case 1:
                if (ligand && postbind > prebind) g_rgnxform_theta[l] += amount;
                else
                {
                    normal.r = -normal.r;
                    move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, normal);
                    #if _dbg_soft
                    cout << " reverting normal.";
                    #endif
                }
                break;

                case 2:
                if (ligand && postbind > prebind) g_rgnxform_y[l] += amount;
                else
                {
                    alpha.r = -alpha.r;
                    move_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, alpha);
                    #if _dbg_soft
                    cout << " reverting alpha.";
                    #endif
                }
                break;

                case 3:
                if (ligand && postbind > prebind) g_rgnrot_alpha[l] += amount/10;
                else
                {
                    rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, alpha, -amount/10);
                    #if _dbg_soft
                    cout << " reverting alpha rot.";
                    #endif
                }
                break;

                case 4:
                if (ligand && postbind > prebind) g_rgnrot_w[l] += amount/50;
                else
                {
                    rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, r, -amount/50);
                    #if _dbg_soft
                    cout << " reverting r rot.";
                    #endif
                }
                break;

                case 5:
                if (ligand && postbind > prebind) g_rgnrot_u[l] += amount/30;
                else
                {
                    rotate_piece(l_soft_rgns[l].start, l_soft_rgns[l].end, rgncen, normal, -amount/30);
                    #if _dbg_soft
                    cout << " reverting normal rot.";
                    #endif
                }
                break;

                default:
                ;
            }

            if (postbind > prebind) prebind = postbind;

            #if _dbg_soft
            cout << endl;
            #endif
        }
    }

    mass_undoable = wmu;
}

void Protein::minimize_residue_clashes(int resno)
{
    AminoAcid** caa = get_residues_can_clash(resno);
    if (!caa) return;

    int i;
    for (i=0; caa[i]; i++);     // count.

    Molecule* mols[i+4];
    MovabilityType mt[i+4];
    AminoAcid* aa = get_residue(resno);
    mols[0] = (Molecule*)aa;
    if (!aa) return;
    mt[0] = aa->movability;
    aa->movability = MOV_FLEXONLY;

    for (i=0; caa[i]; i++)
    {
        mt[i+1] = caa[i]->movability;
        mols[i+1] = (Molecule*)caa[i];
        caa[i]->movability = MOV_NONE;
    }
    mols[i+1] = nullptr;

    Molecule::conform_molecules(mols, 20);

    for (i=0; caa[i]; i++)
    {
        caa[i]->movability = mt[i+1];
    }
    aa = get_residue(resno);
    aa->movability = mt[0];
}

float Protein::binding_to_nearby_residues(int resno)
{
    AminoAcid** caa = get_residues_can_clash(resno);
    if (!caa) return 0;

    AminoAcid* aa = get_residue(resno);
    if (!aa) return 0;

    return aa->get_intermol_binding(caa);
}

#define dbg_region_can_move 0

float Protein::region_can_move(int startres, int endres, SCoord direction, bool repack, int isr, int ier)
{
    int n = get_end_resno();
    Pose revert_to[n+1];
    std::string debug_msg;

    stop1 = stop2 = nullptr;
    stop1a = stop2a = nullptr;

    int i, j, l;

    for (i=1; i<=n; i++)
    {
        AminoAcid *aa = get_residue(i);
        if (aa) revert_to[i].copy_state(aa);
    }

    float result = 0, increment = 1, clash;

    for (l=0; l<200; l++)
    {
        direction.r = increment;
        move_piece(startres, endres, direction);
        if (repack) get_internal_clashes(startres, endres, true);

        clash = 0;
        AminoAcid *aa1, *aa2;
        for (i=startres; i<=endres; i++)
        {
            aa1 = get_residue(i);
            if (!aa1) continue;

            for (j=1; j<=n; j++)
            {
                if (j >= startres-7 && j <= endres+7) continue;
                if (j >= isr && j <= ier) continue;
                aa2 = get_residue(j);
                if (!aa2->can_reach(aa1)) continue;
                float c = aa1->get_intermol_clashes(aa2);
                if (c > clash)
                {
                    clash = c;
                    debug_msg = (std::string)aa1->get_3letter() + std::to_string(aa1->get_residue_no())
                        + (std::string)" clashes by " + std::to_string(c)
                        + (std::string)" with "
                        + (std::string)aa2->get_3letter() + std::to_string(aa2->get_residue_no());
                    stop1 = aa1;
                    stop2 = aa2;
                    stop1a = aa1->clash1;
                    stop2a = aa1->clash2;
                }
            }
        }

        if (clash > clash_limit_per_aa)
        {
            #if dbg_region_can_move
            cout << "At " << (result+increment) << ", " << debug_msg << "." << endl;
            #endif

            direction.r *= -1;
            move_piece(startres, endres, direction);
            increment *= 0.81;
        }
        else result += increment;

        if (fabs(increment) < 0.01) break;
    }

    for (i=1; i<=n; i++)
    {
        AminoAcid *aa = get_residue(i);
        if (aa) revert_to[i].restore_state(aa);
    }

    return result;
}

float Protein::region_can_rotate(int startres, int endres, LocatedVector axis, bool repack, float eca, int isr, int ier)
{
    int n = get_end_resno();
    Pose revert_to[n+1];
    std::string debug_msg;

    stop1 = stop2 = nullptr;
    stop1a = stop2a = nullptr;

    int i, j, l;

    for (i=1; i<=n; i++)
    {
        AminoAcid *aa = get_residue(i);
        if (aa) revert_to[i].copy_state(aa);
    }

    float result = 0, increment = fiftyseventh*5, clash=0, initclash=0;

    for (l=0; l<200; l++)
    {
        rotate_piece(startres, endres, axis.origin, axis, increment);

        clash = 0;
        AminoAcid *aa1, *aa2;
        for (i=startres; i<=endres; i++)
        {
            aa1 = get_residue(i);
            if (!aa1) continue;

            for (j=1; j<=n; j++)
            {
                if (j >= startres-7 && j <= endres+7) continue;
                if (j >= isr && j <= ier) continue;
                aa2 = get_residue(j);
                if (!aa2->can_reach(aa1)) continue;

                float c = aa1->get_intermol_clashes(aa2);

                if (repack && c > clash)
                {
                    Molecule* mols[3];
                    mols[0] = (Molecule*)aa1;
                    mols[1] = (Molecule*)aa2;
                    mols[2] = nullptr;

                    Molecule::conform_molecules(mols, 20);

                    c = aa1->get_intermol_clashes(aa2);
                }

                if (c > clash)
                {
                    clash = c;
                    debug_msg = (std::string)aa1->get_3letter() + std::to_string(aa1->get_residue_no())
                        + (std::string)" clashes by " + std::to_string(c)
                        + (std::string)" with "
                        + (std::string)aa2->get_3letter() + std::to_string(aa2->get_residue_no());
                    stop1 = aa1;
                    stop2 = aa2;
                    stop1a = aa1->clash1;
                    stop2a = aa1->clash2;
                }
            }
        }

        if (!l) initclash = clash;
        if (clash > clash_limit_per_aa+eca && clash > 1.1*initclash)
        {
            // cout << "At " << ((result+increment)*fiftyseven) << "deg, " << debug_msg << "." << endl;
            rotate_piece(startres, endres, axis.origin, axis, -increment);
            increment *= 0.81;
        }
        else
        {
            result += increment;
            initclash = clash;
        }

        if (fabs(increment) < 0.01) break;
    }

    for (i=1; i<=n; i++)
    {
        AminoAcid *aa = get_residue(i);
        if (aa) revert_to[i].restore_state(aa);
    }

    return result;
}

void Protein::region_optimal_positioning(int sr, int er, SCoord* x, Rotation* r, Protein** p)
{
    if (!x || !r)
    {
        cout << "Protein::region_optimal_positioning() output parameters cannot be null." << endl;
        throw -1;
    }

    int i, j, n;
    for (j=0; p && p[j]; j++);          // Count other strands.

    Protein* strands[j+2];

    strands[0] = this;
    for (i=0; i<j; i++) strands[i+1] = p[i];
    strands[i+1] = nullptr;
    AminoAcid* reaching[1024];

    int resno, divisor = 0;
    x->phi = x->theta = x->r = 0;

    Point sum_pt(0,0,0), sum_aln(0,0,0);
    Point center = get_region_center(sr, er);

    for (resno = sr; resno <= er; resno++)
    {
        Star aa;
        aa.paa = get_residue(resno);
        if (!aa.n) continue;
        Point CA = aa.paa->get_CA_location();
        bool residue_counted_yet = false;

        for (i=0; strands[i]; i++)
        {
            n = strands[i]->get_residues_can_clash_ligand(reaching, aa.pmol, aa.paa->get_CA_location(), Point(7,7,7), nullptr);
            for (j=0; j<n; j++)
            {
                /*float e = -aa.pmol->get_intermol_binding(reaching[j]);
                if (e < clash_limit_per_aa) continue;*/

                SCoord motion = aa.pmol->motion_to_optimal_contact(reaching[j]);
                if (fabs(motion.r) < 0.1) continue;

                Point pt_temp = *x;
                Point pt_add = motion;
                Point pt_diff = aa.paa->get_CA_location().subtract(reaching[j]->get_CA_location());

                if (pt_add.x * sgn(pt_diff.x) > pt_temp.x * sgn(pt_diff.x)) pt_temp.x = pt_add.x;
                if (pt_add.y * sgn(pt_diff.y) > pt_temp.y * sgn(pt_diff.y)) pt_temp.y = pt_add.y;
                if (pt_add.z * sgn(pt_diff.z) > pt_temp.z * sgn(pt_diff.z)) pt_temp.z = pt_add.z;

                *x = pt_temp;

                if (!residue_counted_yet)
                {
                    divisor++;
                    residue_counted_yet = true;
                }

                Point CA_new = CA.add(motion);

                if (CA.y > center.y)
                {
                    sum_pt = sum_pt.add(CA.subtract(center));
                    sum_aln = sum_aln.add(CA_new.subtract(center));
                }
                else
                {
                    sum_pt = sum_pt.add(center.subtract(CA));
                    sum_aln = sum_aln.add(center.subtract(CA_new));
                }
            }
        }
    }

    // if (divisor) x->r /= divisor;
    *r = align_points_3d(sum_pt, sum_aln, Point(0,0,0));
}

BallesterosWeinstein Protein::get_bw_from_resno(int resno)
{
    int i, b=0, w=0;

    i=1;
    do
    {
        int j = this->get_bw50(i);
        if (j<0) continue;

        int x = 50 + resno-j;
        if (!w || abs(x-50) < abs(w-50))
        {
            w = x;
            b = i;
        }

        if (i > 10) i %= 10;
        else i = i*10 + i+1;
    } while (i != 8);

    return BallesterosWeinstein(b, w);
}

int Protein::replace_side_chains_from_other_protein(Protein* other)
{
    int i, j, l, n = other->get_end_resno(), num_applied = 0;

    for (i=1; i<=n; i++)
    {
        AminoAcid* source = other->get_residue(i);
        if (!source) continue;
        BallesterosWeinstein bw = other->get_bw_from_resno(i);
        if (!bw.helix_no) continue;
        if (get_bw50(bw.helix_no) < 1) continue;
        AminoAcid* dest = get_residue(bw);
        if (!dest) continue;
        if (source->get_aa_definition() == dest->get_aa_definition()) continue;

        // cout << *source << " -> " << *dest << endl;

        Atom *sN = source->get_atom("N"),
            *sHN = source->HN_or_substitute(),
            *sCA = source->get_atom("CA"),
            *sC  = source->get_atom("C"),
            *sO  = source->get_atom("O"),
            *dN  = dest->get_atom("N"),
            *dHN = dest->HN_or_substitute(),
            *dCA = dest->get_atom("CA"),
            *dC  = dest->get_atom("C"),
            *dO  = dest->get_atom("O");
        
        if (!sN || !sHN || !sCA || !sC || !sO || !dN || !dHN || !dCA || !dC || !dO) continue;
        
        // TODO: Write a way to copy an amino acid, including all atoms and bonds. For now:
        // Move the source CA to coincide with the dest CA.
        SCoord motion = dCA->get_location().subtract(sCA->get_location());
        source->movability = MOV_ALL;
        source->move(motion, true);

        // Rotate the source so that the N coincides with the dest N.
        LocatedVector lv;
        Rotation rot = align_points_3d(sN->get_location(), dN->get_location(), sCA->get_location());
        lv = rot.v;
        lv.origin = sCA->get_location();
        source->rotate(lv, rot.a);

        // Rotate the source about the N-CA axis so that the C coincides with the dest C.
        lv = (SCoord)sCA->get_location().subtract(sN->get_location());
        lv.origin = sCA->get_location();
        float theta = find_angle_along_vector(sC->get_location(), dC->get_location(), sCA->get_location(), lv);
        source->rotate(lv, theta);

        // Rotate the CA-N bond so that the HN coincides with the dest HN.
        lv = (SCoord)sN->get_location().subtract(sCA->get_location());
        theta = find_angle_along_vector(sHN->get_location(), dHN->get_location(), sN->get_location(), lv);
        Bond* b = sCA->get_bond_between(sN);
        if (!b) cout << "WARNING - " << *source << ":" << sCA->name << " is not bonded to " << sN->name << endl;
        b->rotate(theta);

        // Rotate the CA-C bond so that the O coincides with the dest O.
        lv = (SCoord)sC->get_location().subtract(sCA->get_location());
        theta = find_angle_along_vector(sO->get_location(), dO->get_location(), sC->get_location(), lv);
        b = sCA->get_bond_between(sC);
        if (!b) cout << "WARNING - " << *source << ":" << sCA->name << " is not bonded to " << sC->name << endl;
        b->rotate(theta);

        // Replace the dest residue with the source residue and delete the source residue from the other protein.
        // Repoint the source residue's prevaa and nextaa, as well as those of its neighbors.
        for (j=0; residues[j]; j++)
        {
            if (residues[j] == dest)
            {
                residues[j] = source;
                source->set_prev(dest->get_prev());
                source->set_next(dest->get_next());

                break;
            }
        }

        for (j=0; other->residues[j]; j++)
        {
            if (other->residues[j] == source)
            {
                for (l=j+1; other->residues[l-1]; l++)
                {
                    other->residues[l-1] = other->residues[l];
                }
                break;
            }
        }

        num_applied++;
    }

    return num_applied;
}
