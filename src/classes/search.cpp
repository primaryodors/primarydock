
#include "search.h"

Point loneliest;
Point size(10,10,10);
std::vector<int> exclusion;
AtomGroup ligand_groups[3];
ResidueGroup sc_groups[3];

std::vector<std::shared_ptr<AtomGroup>> agc;
std::vector<AminoAcid*> cs_res;
std::vector<intera_type> cs_bt;
std::vector<AtomGroup*> cs_lag;
int cs_idx = -1;

void Search::do_tumble_spheres(Protein* protein, Molecule* ligand, Point l_pocket_cen)
{
    int i, j, l, n;
    float lig_min_int_clsh = ligand->get_internal_clashes();

    // Begin tumble sphere behavior.
    std::vector<AminoAcid*> tsphres = protein->get_residues_near(l_pocket_cen, size.magnitude()+6);
    int tsphsz = tsphres.size();
    float outer_sphere[tsphsz+4], inner_sphere[tsphsz+4];

    for (i=0; i<tsphsz+4; i++) outer_sphere[i] = inner_sphere[i] = 0;

    Point pocketsize = protein->estimate_pocket_size(tsphres);
    Point ligbbox = ligand->get_bounding_box();

    for (i=0; !ligbbox.fits_inside(pocketsize) && i<100; i++)
    {
        ligand->crumple(fiftyseventh*30);
    }

    for (i=0; i<tsphsz; i++)
    {
        #if use_exclusions
        if (exclusion.size()
                &&
                std::find(exclusion.begin(), exclusion.end(), tsphres[i]->get_residue_no())!=exclusion.end()
        )
        {
            tsphres.erase(tsphres.begin()+i);
            tsphsz--;
            continue;
        }
        #endif

        // TODO: Algorithmically determine more accurate values based on interaction type, etc.
        outer_sphere[i] = tsphres[i]->get_reach() + 2.5;
        inner_sphere[i] = tsphres[i]->get_reach() / 3 + 1;
    }

    const SCoord xaxis = Point(1,0,0), yaxis = Point(0,1,0), zaxis = Point(0,0,1);
    float loneliness=0, blone=0, xrad, yrad, zrad, lrad, step, bestxr, bestyr, bestzr, score, worth, weight, bestscore;
    const int ac = ligand->get_atom_count();
    Pose besp(ligand);
    #if _DBG_TUMBLE_SPHERES
    std::string tsdbg = "", tsdbgb = "";
    #endif

    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) ligand->rotate(zaxis, square);
    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) ligand->rotate(yaxis, square);
    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) ligand->rotate(zaxis, square);
    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) ligand->rotate(xaxis, square);
    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) ligand->rotate(yaxis, square);
    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) ligand->rotate(xaxis, square);

    step = fiftyseventh*30;
    bestscore = -Avogadro;
    float lonely_step = 1.0 / loneliest.get_3d_distance(l_pocket_cen);
    #if _DBG_LONELINESS
    cout << "Loneliest point " << loneliest << " is " << loneliest.get_3d_distance(l_pocket_cen) << "A from pocketcen " << l_pocket_cen << "." << endl;
    cout << "Pocket size is " << pocketsize << " vs. ligand bounding box " << ligbbox << endl;
    #endif
    if (isnan(lonely_step) || lonely_step < 0.1) lonely_step = 0.1;

    #if pocketcen_is_loneliest
    if (1)
    {
        ligand->recenter(l_pocket_cen);
    #else
    for (loneliness=0; loneliness <= 1; loneliness += lonely_step)
    {
        float centeredness = 1.0 - loneliness;
        Point tmpcen(loneliest.x * loneliness + l_pocket_cen.x * centeredness,
                     loneliest.y * loneliness + l_pocket_cen.y * centeredness,
                     loneliest.z * loneliness + l_pocket_cen.z * centeredness
                    );
        ligand->recenter(tmpcen);
    #endif

        #if _DBG_LONELINESS && !pocketcen_is_loneliest
        cout << "Ligand is " << loneliness << " lonely centered at " << tmpcen << "." << endl;
        #endif

        for (xrad=0; xrad <= M_PI*2; xrad += step)
        {
            for (yrad=0; yrad <= M_PI*2; yrad += step)
            {
                for (zrad=0; zrad <= M_PI*2; zrad += step)
                {
                    ligbbox = ligand->get_bounding_box();

                    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) continue;
                    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) continue;
                    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) continue;
                    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) continue;
                    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) continue;
                    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) continue;

                    Bond** rb = ligand->get_rotatable_bonds();

                    if (!rb) n = 0;
                    else for (n=0; rb[n]; n++);		// Get count.

                    l = 0;
                    lrad = 0;
                _xyzl_loop:
                    if (ligand->get_internal_clashes() >= lig_min_int_clsh*5+5) goto _xyzl_skip_loop;

                    score = 0;
                    #if _DBG_TUMBLE_SPHERES
                    tsdbg = "";
                    // cout << ligand->get_internal_clashes() << " vs. " << lig_min_int_clsh << endl;
                    #endif
                    for (i=0; i<ac; i++)
                    {
                        Atom* a = ligand->get_atom(i);
                        intera_type it = vdW;

                        for (j=0; j<tsphsz; j++)
                        {
                            worth = 0.4;
                            if (a->get_charge() && tsphres[j]->get_charge()
                                    &&
                                    sgn(a->get_charge()) == -sgn(tsphres[j]->get_charge())
                            )
                            {
                                it = ionic;
                                worth = 100;
                            }
                            else if (a->get_charge() || a->is_polar())
                            {
                                it = hbond;
                                worth = 40;
                            }
                            else if (a->is_pi())
                            {
                                it = pi;
                                worth = 7;
                            }

                            #if active_persistence
                            worth *= residue_binding_multiplier(tsphres[j]->get_residue_no());
                            #endif

                            if (tsphres[j]->capable_of_inter(it))
                            {
                                float r = a->get_location().get_3d_distance(tsphres[j]->get_atom_location("CA"));
                                if (r <= outer_sphere[j])
                                {
                                    if (r > inner_sphere[j])
                                    {
                                        weight = 1;

                                        if (tsphres[j]->priority)
                                        {
                                            weight = ts_priority_coefficient;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                                        }

                                        #if !tumble_spheres_include_vdW
                                        if ((worth*weight) < 1) continue;
                                        #endif

                                        score += worth * weight;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("+ ")
                                                +  std::string(a->name) + std::string(" reaches ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                    else
                                    {
                                        score -= 200;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("- ")
                                                +  std::string(a->name) + std::string(" clashes ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                }
                            }
                        }
                    }

                    #if !pocketcen_is_loneliest
                    if (score > 0) score *= 1.0 + 0.1 * centeredness;
                    #endif

                    if (score > bestscore)
                    {
                        besp.copy_state(ligand);
                        blone = loneliness;
                        bestxr = xrad;
                        bestyr = yrad;
                        bestzr = zrad;
                        bestscore = score;

                        #if _DBG_TUMBLE_SPHERES
                        tsdbgb = tsdbg;

                        cout << "Tumble score " << score << " for ligand box " << ligand->get_bounding_box() << endl;


                        #if output_tumble_debug_docs
                        int u, v, w;
                        char protfttl[1000];
                        strcpy(protfttl, protfname);

                        char** lwords = chop_spaced_words(protfttl, '/');

                        for (u=0; lwords[u]; u++);
                        u--;

                        char fname[1000];
                        sprintf(fname, "output/tumble_%s_%d_%d_%d_%f.dock",
                                lwords[u],
                                (int)(xrad*fiftyseven),
                                (int)(yrad*fiftyseven),
                                (int)(zrad*fiftyseven),
                                score);
                        cout << fname << endl;
                        std::ofstream tspdbdat(fname, std::ofstream::out);

                        tspdbdat << "PDB file: " << protfname << endl;
                        tspdbdat << "Pose: 1\nNode: 0\nPDBDAT:\n";

                        int lac = ligand->get_atom_count();
                        for (u=0; u<lac; u++) ligand->get_atom(u)->stream_pdb_line(tspdbdat, 9000+u);

                        int pseql = protein->get_seq_length();
                        v = 1;
                        for (u = 1; u < pseql; u++)
                        {
                            AminoAcid* dbgaa = protein->get_residue(u);
                            if (dbgaa)
                            {
                                int aaac = dbgaa->get_atom_count();
                                for (w=0; w<aaac; w++)
                                {
                                    Atom* dbga = dbgaa->get_atom(w);
                                    if (!strcmp(dbga->name, "CA") || !strcmp(dbga->name, "CB")) dbga->stream_pdb_line(tspdbdat, v++);
                                }
                            }
                        }
                        tspdbdat << "END" << endl;
                        tspdbdat.close();
                        #endif
                        #endif
                    }

                _xyzl_skip_loop:

                    if (rb && rb[l])
                    {
                        rb[l]->rotate(step);

                        lrad += step;
                        if (lrad >= M_PI*2)
                        {
                            l++;
                            if (l < n) goto _xyzl_loop;
                        }
                        else goto _xyzl_loop;
                    }

                    ligand->rotate(zaxis, step);
                }		// zrad.
                ligand->rotate(yaxis, step);
            }			// yrad.
            ligand->rotate(xaxis, step);
        }				// xrad.

        #if !pocketcen_is_loneliest
        if (bestscore >= (ligand->get_atom_count()*13)) break;
        #endif

        #if _DBG_LONELINESS
        cout << "Best score: " << bestscore << endl;
        #endif
    }					// loneliness.

    #if _DBG_TUMBLE_SPHERES
    cout << "Tumble sphere best score " << bestscore << " for "
        << "x" << bestxr*fiftyseven << "deg, "
        << "y" << bestyr*fiftyseven << "deg, "
        << "z" << bestzr*fiftyseven << "deg."
        << " (" << blone << " lonely)."
        << endl;
    cout << tsdbgb << endl;
    #endif

    // Load the best ligand conformer.
    besp.restore_state(ligand);

    // Minimize ligand clashes.
    #if prerot_sidechains_from_ligand
    for (i=0; i<tsphsz; i++)
    {
        Bond** tsphb = tsphres[i]->get_rotatable_bonds();
        if (tsphb)
        {
            for (j=0; tsphb[j]; j++)
            {
                float rad=0, bestrad=0, clash, bestclash=6.25e24;
                for (rad=0; rad < M_PI*2; rad += step)
                {
                    clash = tsphres[i]->get_intermol_clashes(ligand);

                    if (clash < bestclash)
                    {
                        bestrad = rad;
                        bestclash = clash;
                    }

                    tsphb[j]->rotate(step);
                }

                tsphb[j]->rotate(bestrad);
            }
            // delete[] tsphb;
        }
    }
    #endif
    // End tumble sphere behavior.
}

void Search::do_best_binding(Protein* protein, Molecule* ligand, Point l_pocket_cen, AminoAcid** reaches_spheroid)
{
    #if _dbg_groupsel
    cout << "Candidate binding residues: ";
    for (i=0; i<sphres; i++) cout << *reaches_spheroid[i] << " ";
    cout << endl;
    #endif

    std::vector<std::shared_ptr<AtomGroup>> agc = AtomGroup::get_potential_ligand_groups(ligand, mtlcoords.size() > 0);
    std::vector<std::shared_ptr<ResidueGroup>> scg = ResidueGroup::get_potential_side_chain_groups(reaches_spheroid, l_pocket_cen);
    global_pairs = GroupPair::pair_groups(agc, scg, l_pocket_cen);

    if (global_pairs.size() > 2)
    {
        // If the 2nd group is closer to the 1st group than the 3rd group is, swap the 2nd and 3rd groups.
        Point grpcen1 = global_pairs[0]->ag->get_center(),
            grpcen2 = global_pairs[1]->ag->get_center(),
            grpcen3 = global_pairs[2]->ag->get_center();

        float r12 = grpcen1.get_3d_distance(grpcen2),
            r13 = grpcen1.get_3d_distance(grpcen3);
        
        // cout << r12 << " " << r13 << endl;

        if (r12 < r13)
        {
            std::shared_ptr<GroupPair> tmpg = global_pairs[2];
            global_pairs[2] = global_pairs[1];
            global_pairs[1] = tmpg;
        }
    }

    ligand->recenter(l_pocket_cen);

    GroupPair::align_groups(ligand, global_pairs, false, 1);

    #if _dbg_groupsel
    cout << endl;
    #endif

    int l;
    Molecule* lmols[256];
    for (l=0; reaches_spheroid[l]; l++)
    {
        lmols[l] = (Molecule*)reaches_spheroid[l];
    }

    #if bb_enable_residue_disqualifications
    if (ligand->get_intermol_clashes(lmols) >= bb_disqualification_energy)
    {
        #if _dbg_groupsel
        cout << "Primary residue group is disqualified." << endl;
        #endif

        global_pairs[0]->disqualify();
        std::vector<std::shared_ptr<ResidueGroup>> scg = ResidueGroup::get_potential_side_chain_groups(reaches_spheroid, l_pocket_cen);
        global_pairs = GroupPair::pair_groups(agc, scg, l_pocket_cen);
        GroupPair::align_groups(ligand, global_pairs, false, 1);
    }
    #endif
}

void Search::prepare_constrained_search(Protein* protein, Molecule* ligand, Point l_pocket_cen)
{
    // Enumerate all binding pocket residues into an array.
    AminoAcid* baa[256];
    int nba = protein->get_residues_can_clash_ligand(baa, ligand, l_pocket_cen, size, nullptr);
    agc = AtomGroup::get_potential_ligand_groups(ligand, mtlcoords.size() > 0);

    const int num_allowed_type = 5;
    const intera_type allowed_types[num_allowed_type] = {mcoord, ionic, hbond, pi, vdW};

    // For each pocket residue,
    int i, j, l, n;
    for (i=0; i<nba; i++)
    {
        bool res_has_nonvdw = false;
        bool res_has_ionic = false;

        // For each of the binding types: mcoord, ionic, hbond, pi, vdW:
        for (j=0; j<num_allowed_type; j++)
        {
            // If other binding types have already been found for this residue, skip vdW.
            if (allowed_types[j] == vdW && res_has_nonvdw) continue;
            else if (res_has_ionic) continue;
        
            // Can the residue's side chain form this type of bond with the ligand?
            bool can_bind = false;
            float lc, rc;
            switch (allowed_types[j])
            {
                case mcoord:
                if (baa[i]->coordmtl && ligand->count_atoms_by_element("S"))
                    can_bind = res_has_nonvdw = true;
                break;

                case ionic:
                rc = baa[i]->get_charge();
                lc = ligand->get_charge();
                if (rc && lc && sgn(rc) == -sgn(lc)) can_bind = res_has_nonvdw = res_has_ionic = true;
                break;

                case hbond:
                if (res_has_ionic) continue;
                if (baa[i]->has_hbond_donors() && ligand->has_hbond_acceptors()) can_bind = res_has_nonvdw = true;
                else if (baa[i]->has_hbond_acceptors() && ligand->has_hbond_donors()) can_bind = res_has_nonvdw = true;
                break;

                case pi:
                if (baa[i]->has_pi_atoms() && ligand->has_pi_atoms()) can_bind = res_has_nonvdw = true;
                break;

                case vdW:
                default:
                if (baa[i]->hydrophilicity() < 1.5*hydrophilicity_cutoff && ligand->hydrophilicity() < 1.5*hydrophilicity_cutoff) can_bind = true;
            }
            
            if (can_bind)
            {
                // If so, record the residue, the binding type, and the ligand atom group with the strongest potential.
                AtomGroup* ag = nullptr;
                float agbb = 0;
                n = agc.size();
                for (l=0; l<n; l++)
                {
                    float mpb = agc[l]->max_potential_binding(allowed_types[j]);
                    if (mpb > agbb)
                    {
                        agbb = mpb;
                        ag = agc[l].get();
                    }
                }
                if (!ag) cout << "BAD CSAG BINDING: " << baa[i]->get_name() << " " << allowed_types[j] << endl;
                if (!ag) continue;

                cs_res.push_back(baa[i]);
                cs_bt.push_back(allowed_types[j]);
                cs_lag.push_back(ag);
            }
        }
    }
}

void Search::do_constrained_search(Protein* protein, Molecule* ligand)
{
    int i, j, l, n;

    // Choose a residue-type-group combination, randomly but weighted by binding energy of binding type.
    n = cs_res.size();
    while (true)
    {
        for (j=0; j<n; j++)
        {
            float b;
            switch (cs_bt[j])
            {
                case mcoord: b = 200; break;
                case ionic: b = 60; break;
                case hbond: b = 25; break;
                case pi: b = 12; break;
                case vdW: default: b = 4;
            }

            float r = fmax(2.8, cs_res[j]->get_CA_location().get_3d_distance(loneliest) - cs_res[j]->get_reach());
            float w = pow(b/500, cs_bondweight_exponent) / pow(r, 3) * 10;
            if (cs_res[j]->priority) w *= 3.7;

            if (frand(0,1) < w) goto chose_residue;
        }
    }
    chose_residue:
    if (j<n) cs_idx = j;

    // Place the ligand so that the atom group is centered in the binding pocket.
    Point agp = cs_lag[j]->get_center();
    SCoord mov = loneliest.subtract(agp);
    ligand->move(mov);
    
    // Move the ligand so that the atom group is at the optimal distance to the residue.
    Point resna = cs_res[j]->get_nearest_atom(loneliest)->get_location();
    Point agcen = cs_lag[j]->get_center();
    mov = resna.subtract(agcen);
    mov.r = cs_lag[j]->distance_to(resna) - 2;              // One size fits all for now. It's the iterations' job to correct this distance.
    
    // Rotate the ligand about the residue so that its barycenter aligns with the "loneliest" point.
    Rotation rot = align_points_3d(ligand->get_barycenter(), loneliest, agcen);
    LocatedVector lv = rot.v;
    lv.origin = agcen;
    ligand->rotate(lv, rot.a);
    
    // Conform the side chain and ligand to each other, ignoring other residues.
    Molecule* mm[3];
    mm[0] = cs_res[j];
    mm[1] = ligand;
    mm[2] = nullptr;
    Molecule::conform_molecules(mm, 30);

    // return;
    
    // Perform a monaxial 360° rotation about the residue and the imaginary line between ligand barycenter and residue,
    // and look for the rotamer with the smallest clash total.
    lv = (SCoord)resna.subtract(ligand->get_barycenter());
    lv.origin = agcen;
    Pose best(ligand);
    best.copy_state(ligand);
    float least_clash;
    float theta = 0;
    for (; theta < M_PI*2; theta += cs_360_step)
    {
        AminoAcid* cc[256];
        protein->get_residues_can_clash_ligand(cc, ligand, ligand->get_barycenter(), size, nullptr);
        float f = ligand->get_intermol_clashes(reinterpret_cast<Molecule**>(cc));
        if (!theta || f < least_clash)
        {
            least_clash = f;
            best.copy_state(ligand);
        }

        ligand->rotate(lv, cs_360_step);
    }
            
    best.restore_state(ligand);
}
