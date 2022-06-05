
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "protein.h"

#define _DBG_REACHLIG true

using namespace std;

Protein::Protein(const char* lname)
{
    name = lname;
    aaptrmin.n = aaptrmax.n = 0;
}

bool Protein::add_residue(const int resno, const char aaletter)
{
    int i;

    if (!residues)
    {
        residues = new AminoAcid*[resno+256];
        sequence = new char[resno+256];
        ca = new Atom*[resno+256];
        res_reach = new float[resno+256];
    }
    else if (resno % 256)
    {
        AminoAcid** oldres = residues;
        char* oldseq = sequence;
        residues = new AminoAcid*[resno+261];
        sequence = new char[resno+261];
        ca = new Atom*[resno+261];
        res_reach = new float[resno+261];
        for (i=0; oldres[i]; i++)
        {
            residues[i] = oldres[i];
            sequence[i] = oldseq[i];
        }
        residues[i] = 0;
        sequence[i] = 0;
    }

    i = resno-1;
    residues[i] = new AminoAcid(aaletter, i ? residues[i-1] : 0);
    sequence[i] = aaletter;
    ca[i] = residues[i]->get_atom("CA");
    res_reach[i] = residues[i]->get_aa_definition()->reach;
    residues[i+1] = 0;
    sequence[i+1] = 0;

    // TODO: Create a peptide bond to the previous residue.

    if (!aaptrmin.n || residues[i] < aaptrmin.paa) aaptrmin.paa = residues[i];
    if (!aaptrmax.n || residues[i] > aaptrmax.paa) aaptrmax.paa = residues[i];

    return true;
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

Atom* Protein::get_atom(int resno, const char* aname)
{	AminoAcid* aa = get_residue(resno);
	
	if (!aa) return NULL;
	return aa->get_atom(aname);
}

Point Protein::get_atom_location(int resno, const char* aname)
{	AminoAcid* aa = get_residue(resno);
	
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
        float r = get_atom_location(i+1, "CA").get_3d_distance(get_atom_location(i+1, "CB"));
        if (fabs(r-1.54) > 0.2) cout << i+1 << lsequence[i] << ":CA-CB = " << r << endl;
    }
    
    int seql = get_seq_length();
    Molecule* aas[seql+4];
    for (i=1; i<=seql; i++)
    {
    	aas[i-1] = get_residue(i);
    	aas[i] = 0;
    }
    Molecule::multimol_conform(aas, 25);
    
    set_clashables();

    return true;
}

void Protein::save_pdb(FILE* os)
{
    if (!residues) return;
    int i, offset=0;
    for (i=0; residues[i]; i++)
    {
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

    fprintf(os, "\nTER\n");
}

void Protein::end_pdb(FILE* os)
{
    fprintf(os, "END\n");
}

int Protein::load_pdb(FILE* is)
{
    AminoAcid* restmp[65536];

    int i, rescount=0;

    while (!feof(is))
    {
        try
        {
            AminoAcid* aa = new AminoAcid(is);
            restmp[rescount++] = aa;
        }
        catch (int ex)
        {
            // cout << "Exception " << ex << endl;
            if (ex == ATOM_NOT_OF_AMINO_ACID)
            {
                if (!metcount % 16)
                {
                    Atom** mtmp = new Atom*[metcount+20];
                    if (metals)
                    {
                        for (i=0; metals[i]; i++)
                            mtmp[i] = metals[i];
                        delete metals;
                    }
                    metals = mtmp;
                }

                Atom* a = new Atom(is);
                metals[metcount++] = a;

                for (i=0; i<26; i++)
                {
                    if (aa_defs[i]._1let && !strcmp(aa_defs[i]._3let, a->aa3let ))
                    {
                        a->aaletter = aa_defs[i]._1let;
                        break;
                    }
                }

            }
            else throw 0xbadca22;
        }
    }

    residues 	= new AminoAcid*[rescount+1];
    sequence 	= new char[rescount+1];
    ca       	= new Atom*[rescount+1];
    res_reach	= new float[rescount+1];

    for (i=0; i<rescount; i++)
    {
        try
        {
            residues[i] = restmp[i];

            Atom *atom = residues[i]->get_atom("N"), *btom;
            AminoAcid* prev = get_residue(residues[i]->get_residue_no()-1);
            if (prev)
            {
                btom = prev->get_atom("C");
                if (atom && btom) atom->bond_to(btom, 1.5);
            }

            if (!aaptrmin.n || residues[i] < aaptrmin.paa) aaptrmin.paa = residues[i];
            if (!aaptrmax.n || residues[i] > aaptrmax.paa) aaptrmax.paa = residues[i];

            AADef* raa = restmp[i]->get_aa_definition();
            if (!raa) cout << "Warning: Residue " << (i+1) << " has no AADef." << endl;
            sequence[i] = raa ? raa->_1let : '?';
            ca[i]		= restmp[i]->get_atom("CA");
            res_reach[i]= raa ? raa->reach : 2.5;
        }
        catch (int e)
        {
            cout << "Residue " << (i+1) << " threw an error." << endl;
            throw e;
        }
    }
    // cout << "Read residue " << *residues[rescount] << endl;
    residues[rescount] = 0;

    set_clashables();

    return rescount;
}

int  Protein::get_seq_length()
{
    if (!sequence) return 0;
    int i;
    for (i=0; sequence[i]; i++)
        ;
    return i;
}

int  Protein::get_start_resno()
{
    if (!residues) return 0;
    else return residues[0]->get_residue_no();
}


void Protein::set_clashables()
{
    int i, j, k;

    // cout << "Setting clashables." << endl;

    if (res_can_clash)
    {
        delete[] res_can_clash;
    }

    int seqlen = get_seq_length();
    res_can_clash = new AminoAcid**[seqlen+1];

    // cout << "seqlen is " << seqlen << endl;

    for (i=0; i<seqlen; i++)
    {
        if (debug) *debug << endl << "Testing residue " << residues[i]->get_residue_no() << endl;
        AminoAcid* temp[seqlen+1];
        k=0;
        for (j=0; j<seqlen; j++)
        {
            if (j == i) continue;
            if (residues[i]->can_reach(residues[j]))
            {
                temp[k++] = residues[j];
                if (debug) *debug << *residues[j] << " can reach " << *residues[i] << endl;
            }
        }
        res_can_clash[i] = new AminoAcid*[k+1];
        for (j=0; j<k; j++)
        {
            res_can_clash[i][j] = temp[j];
            /*cout << residues[i]->get_aa_definition()->_3let << residues[i]->get_residue_no()
            	 << " can clash with "
            	 << res_can_clash[i][j]->get_aa_definition()->_3let << res_can_clash[i][j]->get_residue_no()
            	 << endl;*/
        }
        res_can_clash[i][k] = 0;
    }

    res_can_clash[seqlen] = 0;
    
    /*for (i=0; i<seqlen; i++)
    {
    	cout << i << ": ";
    	for (j=0; res_can_clash[i][j]; j++)
    		cout << *res_can_clash[i][j] << " ";
    	cout << endl;
    }*/
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

int Protein::get_residues_can_clash_ligand(AminoAcid** reaches_spheroid,
        const Molecule* ligand,
        const Point nodecen,
        const Point size,
        const int* mcoord_resno
                                            )
{
    int i, j, sphres = 0;
    int seql = get_seq_length();
    bool resno_already[8192];

    for (i=0; i<SPHREACH_MAX; i++) reaches_spheroid[i] = NULL;

    for (i=1; i<=seql; i++)
    {
        AminoAcid* aa = residues[i-1];
        if (!aa) continue;

        int resno = aa->get_residue_no();

        if (mcoord_resno)
            for (j=0; mcoord_resno[j]; j++)
            {
                if (mcoord_resno[j] == resno)
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
            float angle = find_3d_angle(cb->get_location(), pt2, ca->get_location());
            if (angle < M_PI/1.5)
            {
                Point pt1 = pt;
                pt1 = pt1.subtract(&pt2);
                if (pt1.magnitude() <= aa->get_reach()*1.25)
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
            if (angle > M_PI/1.5) continue;
        }

        pt = pt.subtract(&nodecen);
        Point pt1 = pt;
        float dist = pt.magnitude();
        pt1.scale(fmax(dist - aa->get_reach(), 0));

        pt1.x /= size.x;
        pt1.y /= size.y;
        pt1.z /= size.z;

        SCoord dir(&pt1);

        if (dir.r <= 1.25)
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

void Protein::rotate_backbone(int resno, bb_rot_dir dir, float angle)
{
    AminoAcid* bendy = get_residue(resno);
    if (!bendy) return;
    LocatedVector lv = bendy->rotate_backbone(dir, angle);

    if (lv.r)
    {
        int i, inc;
        AminoAcid* movable;
        inc = (dir == CA_desc || dir == C_desc) ? -1 : 1;

        for (i=resno+inc; movable = get_residue(i); i+=inc)
        {
            // cout << "Rotating " << i << endl;
            movable->rotate(lv, angle);
        }
    }
}

void Protein::rotate_backbone_partial(int startres, int endres, bb_rot_dir dir, float angle)
{
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

    if (lv.r)
    {
        int i;
        AminoAcid* movable;

        for (i=startres+inc; movable = get_residue(i); i+=inc)
        {
            movable->rotate(lv, angle);
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

void Protein::conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, int iters)
{
    conform_backbone(startres, endres, a1, target1, a2, target2, iters, false);
}

void Protein::conform_backbone(int startres, int endres,
                               Atom* a1, Point target1,
                               Atom* a2, Point target2,
                               int iters, bool backbone_atoms_only
                              )
{
    int inc = sgn(endres-startres);
    int res, i, j, iter;
    bb_rot_dir dir1 = (inc>0) ? N_asc : CA_desc,
               dir2 = (inc>0) ? CA_asc : C_desc;
    
    int am = abs(endres-startres), minres = (inc>0) ? startres : endres;
    float momenta1[am+4], momenta2[am+4];
    int eando_res[am+4];
    float eando_mult[am+4];
    
    for (res = startres; res <= endres; res += inc)
    {
    	int residx = res-minres;
    	momenta1[residx] = randsgn()*_fullrot_steprad;
    	momenta2[residx] = randsgn()*_fullrot_steprad;
    	eando_res[residx] = min(res + (rand() % 5) + 1, endres);
    	if (eando_res[residx] == res) eando_res[residx] = 0;
    	eando_mult[residx] = 1;
    }

	set_clashables();
    float tolerance = 1.2, alignfactor = 100;
    for (iter=0; iter<iters; iter++)
    {
        cout << "Iteration " << iter << endl;
        for (res = startres; res != endres; res += inc)
        {
        	int residx = res-minres;
        	
            // Get the preexisting nearby residues and inter-residue binding/clash value.
            // These will likely have changed since last iteration.
            float bind=0, bind1=0, angle;

            for (i=res; i != endres; i += inc)
            {
                AminoAcid* aa = get_residue(i);
                AminoAcid** rcc = get_residues_can_clash(i);
                if (!rcc) cout << "No clashables." << endl;
                if (a1)		bind -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
                else		bind += aa->get_intermol_binding(rcc, backbone_atoms_only);
            }
            if (a1)
            {
                Point pt = a1->get_location();
                bind += alignfactor/pt.get_3d_distance(target1);
            }
            if (a2)
            {
                Point pt = a2->get_location();
                bind += alignfactor/pt.get_3d_distance(target2);
            }

            if (strcmp(get_residue(res)->get_3letter(), "PRO"))		// TODO: Don't hard code this to proline, but check bond flexibility.
            {
                // Rotate the first bond a random amount. TODO: use angular momenta.
                angle = momenta1[residx]; // frand(-_fullrot_steprad, _fullrot_steprad);
                rotate_backbone_partial(res, endres, dir1, angle);
                if (eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir2, -angle*eando_mult[residx]);

                // Have bindings/clashes improved?
                for (i=res; i != endres; i += inc)
                {
                    AminoAcid* aa = get_residue(i);
                    AminoAcid** rcc = get_residues_can_clash(i);
                    if (a1)		bind1 -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
                	else		bind1 += aa->get_intermol_binding(rcc, backbone_atoms_only);
                }
                if (a1)
                {
                    Point pt = a1->get_location();
                    bind1 += alignfactor/pt.get_3d_distance(target1);
                }
                if (a2)
                {
                    Point pt = a2->get_location();
                    bind1 += alignfactor/pt.get_3d_distance(target2);
                }

                // If no, put it back.
                if (res == startres) cout << bind << " v. " << bind1 << endl;
                if (bind1 < tolerance*bind)
                {
                    rotate_backbone_partial(res, endres, dir1, -angle);
                    if (eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir2, angle*eando_mult[residx]);
                    momenta1[residx] *= -0.666;
                }
                else
                {
                    if (bind1 < bind) bind = bind1;
                    momenta1[residx] *= 1.05;
                }
            }

            // Rotate the second bond.
            angle = momenta2[residx]; // frand(-_fullrot_steprad, _fullrot_steprad);
            rotate_backbone_partial(res, endres, dir2, angle);
            if (eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir1, -angle*eando_mult[residx]);

            // Improvement?
            bind1 = 0;
            for (i=res; i != endres; i += inc)
            {
                AminoAcid* aa = get_residue(i);
                AminoAcid** rcc = get_residues_can_clash(i);
                if (a1)		bind1 -= aa->get_intermol_clashes(AminoAcid::aas_to_mols(rcc));
            	else		bind1 += aa->get_intermol_binding(rcc, backbone_atoms_only);
            }
            if (a1)
            {
                Point pt = a1->get_location();
                bind1 += alignfactor/pt.get_3d_distance(target1);
            }
            if (a2)
            {
                Point pt = a2->get_location();
                bind1 += alignfactor/pt.get_3d_distance(target2);
            }

            // If no, put it back.
            if (res == startres) cout << bind << " vs. " << bind1 << endl;
            if (bind1 < tolerance*bind)
            {
                rotate_backbone_partial(res, endres, dir2, -angle);
                if (eando_res[residx]) rotate_backbone_partial(eando_res[residx], endres, dir1, angle*eando_mult[residx]);
                momenta2[residx] *= -0.666;
            }
            else
            {
                momenta2[residx] *= 1.05;
            }

            alignfactor *= 1.003;
            tolerance = ((tolerance-1)*0.97)+1;
        }
    }
}

void Protein::make_helix(int startres, int endres, float phi, float psi)
{
    make_helix(startres, endres, endres, phi, psi);
}

void Protein::make_helix(int startres, int endres, int stopat, float phi, float psi)
{
    int inc = sgn(endres-startres);
    if (inc != sgn(stopat-startres)) return;
    if (stopat < endres) endres = stopat;
    int res, i, j, iter;
    int phis[365], psis[365];
    bb_rot_dir dir1 = (inc>0) ? N_asc : CA_desc,
               dir2 = (inc>0) ? CA_asc : C_desc;

    for (res = startres; inc; res += inc)
    {
        AminoAcid* aa = get_residue(res);

        LocRotation* lr2 = aa->flatten();

        for (j=0; j<5; j++)
        {
            if (lr2[j].v.r && lr2[j].a)
            {
                AminoAcid* movable;

                for (i=res+1; movable = get_residue(i); i+=inc)
                {
                    LocatedVector lv = lr2[j].get_lv();
                    movable->rotate(lv, lr2[j].a);
                    if (i == stopat) break;
                }

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

        LocRotation lr = aa->rotate_backbone_abs(dir1, phi);

        if (lr.v.r)
        {
            AminoAcid* movable;

            for (i=res+inc; movable = get_residue(i); i+=inc)
            {
                LocatedVector lv = lr.get_lv();
                movable->rotate(lv, lr.a);
                if (i == stopat) break;
            }
        }

        lr = aa->rotate_backbone_abs(dir2, psi);

        if (lr.v.r)
        {
            AminoAcid* movable;

            for (i=res+inc; movable = get_residue(i); i+=inc)
            {
                LocatedVector lv = lr.get_lv();
                movable->rotate(lv, lr.a);
                if (i == stopat) break;
            }
        }

        if (res >= endres) break;
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
    }
    aas[endres-startres+1] = 0;
    Molecule::multimol_conform(aas, 25);
}

void Protein::delete_residue(int resno)
{
    if (!resno) return;
    if (resno > get_seq_length()) return;
    if (!residues) return;

    int i, j;
    for (i=0; residues[i]; i++)
    {
        if (residues[i]->get_residue_no() == resno)
        {
            for (j=i+1; residues[j-1]; j++)
            {
                residues[j-1] = residues[j];
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

void Protein::delete_sidechain(int resno)
{
    if (!resno) return;
    if (resno > get_seq_length()) return;
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

MetalCoord* Protein::coordinate_metal(Atom* metal, int residues, int* resnos, char** res_anames)
{
    int i, j=0, k=0;
    if (!m_mcoord)
    {
        m_mcoord = new MetalCoord*[2];
        m_mcoord[0] = new MetalCoord();
        m_mcoord[1] = NULL;
    }
    else
    {
        for (j=0; m_mcoord[j]; j++);
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
        for (k=0; metals[k]; k++);
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
        m_mcoord[j]->coord_res[i] = get_residue(resnos[i]);
        if (!m_mcoord[j]->coord_res[i])
        {
            cout << "Attempt to bind metal to residue " << resnos[i] << " not found in protein!" << endl;
            throw 0xbad12e5d;
        }
        m_mcoord[j]->coord_res[i]->m_mcoord = m_mcoord[j];
        m_mcoord[j]->coord_atoms[i] = m_mcoord[j]->coord_res[i]->get_atom(res_anames[i]);
        if (!m_mcoord[j]->coord_atoms[i])
        {
            cout << "Attempt to bind metal to " << resnos[i] << ":" << res_anames[i] << " not found in protein!" << endl;
            throw 0xbada70b;
        }
    }
    m_mcoord[j]->coord_res[residues] = NULL;
    m_mcoord[j]->coord_atoms[residues] = NULL;

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
    cout << "nc " << nc << " | anc " << anc << endl;
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
    int lmolc=0;
    lmols[lmolc++] = &m;
    for (i=minres; i<=maxres; i++)
    {
        AminoAcid* aa = get_residue(i);
        if (aa)
        {
            lmols[lmolc++] = aa;
            if (i >= minres && i <= maxres)
                aa->m_mcoord = m_mcoord[j];		// Sets the coordinating residues, and all in between residues, to backbone immovable.
        }
    }
    lmols[lmolc] = NULL;

    // Make sure to set the coordinating residues so that their side chains are flexible.
    for (i=0; m_mcoord[j]->coord_res[i]; i++)
        m_mcoord[j]->coord_res[i]->movability = MOV_FLEXONLY;

    // Multimol conform the array.
    gmprot = this;
    Molecule::multimol_conform(lmols, 250, &ext_mtl_coord_cnf_cb);

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

void Protein::set_region(std::string rgname, int start, int end)
{
    int i;
    for (i=0; i<PROT_MAX_RGN; i++) if (!regions[i].start) break;
    if (i >= PROT_MAX_RGN) return;		// Nope.

    regions[i].name = rgname;
    regions[i].start = start;
    regions[i].end = end;
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
	for (i=0; i<rglen; i++)
	{
		// This is slow but that's okay.
		range[i] = get_residue(startres+i)->get_barycenter();
	}
	
	return average_of_points(range, rglen);
}

void Protein::move_piece(int start_res, int end_res, Point new_center)
{
	Point old_center = get_region_center(start_res, end_res);
	SCoord move_amt = new_center.subtract(old_center);
	
	int i;
	for (i=start_res; i<=end_res; i++)
	{
		AminoAcid* aa = get_residue(i);
		aa->move(move_amt);
	}
}

void Protein::rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res)
{
	Point pivot = pivot_res ? get_residue(pivot_res)->get_barycenter() : get_region_center(start_res, end_res);
	Point align = get_residue(align_res)->get_barycenter();
	Rotation rot = align_points_3d(&align, &align_target, &pivot);
	
	LocatedVector lv(rot.v);
	lv.origin = pivot;
	int i;
	for (i=start_res; i<=end_res; i++)
	{
		AminoAcid* aa = get_residue(i);
		aa->rotate(lv, rot.a);
	}
}















