
#include "group.h"
#include "scoring.h"

float init_total_binding_by_type[_INTER_TYPES_LIMIT];
float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

float* initial_binding;
#if compute_vdw_repulsion
float* initial_vdWrepl;
#endif

DockResult::DockResult()
{
    ;
}

DockResult::DockResult(Protein* protein, Molecule* ligand, Point size, int* addl_resno, int drcount, Molecule** waters)
{
    int end1 = SPHREACH_MAX+4;
    AminoAcid* reaches_spheroid[end1];
    int sphres = protein->get_residues_can_clash_ligand(reaches_spheroid, ligand, ligand->get_barycenter(), size, addl_resno);
    // cout << "sphres " << sphres << endl;
    Molecule* met = protein->metals_as_molecule();
    int i, j, k, l, n;

    char metrics[end1][20];
    float lmkJmol[end1];
    float limkJmol[end1];
    #if compute_vdw_repulsion
    float lmvdWrepl[end1];
    float limvdWrepl[end1];
    #endif
    int metcount = 0;
    float btot = 0;
    float pstot = 0;
    char* lmba1n[end1];
    char* lmba2n[end1];
    char* lmca1n[end1];
    char* lmca2n[end1];
    float lmc[end1];

    ligpos.copy_state(ligand);

    for (i=0; i<end1; i++)
    {
        #if compute_vdw_repulsion
        lmkJmol[i] = limkJmol[i] = lmvdWrepl[i] = limvdWrepl[i] = lmc[i] = 0;
        #else
        lmkJmol[i] = limkJmol[i] = lmc[i] = 0;
        #endif
    }

    worst_energy = worst_nrg_aa = 0;
    worst_clash_1 = worst_clash_2 = nullptr;

    // if (debug) *debug << "Pose " << pose << " pathnode " << nodeno /*<< " clashes " << clash*/ << endl;

    ligand->clear_atom_binding_energies();

    float final_binding[end1];
    #if compute_vdw_repulsion
    float final_vdWrepl[end1];
    for (i=0; i<end1; i++) final_binding[i] = final_vdWrepl[i] = 0;
    #else
    for (i=0; i<end1; i++) final_binding[i] = 0;
    #endif

    std::vector<AminoAcid*> allres = protein->get_residues_near(ligand->get_barycenter(), 100000, false);
    int qpr = allres.size();
    Molecule* postaa[qpr];
    postaa[0] = ligand;
    for (i=0; i<qpr; i++)
    {
        postaa[i+1] = reinterpret_cast<Molecule*>(allres[i]);
    }

    for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = total_binding_by_type[i];

    #if _peratom_audit
    interaudit.clear();
    interauditing = true;
    #endif

    for (i=0; i<_INTER_TYPES_LIMIT; i++)
    {
        init_total_binding_by_type[i] = 0;
        fin_total_binding_by_type[i] = 0;
        total_binding_by_type[i] = 0;
    }

    int prot_seq_len = protein->get_end_resno();
    #if compute_clashdirs
    residue_clash = new float[prot_seq_len+8];
    res_clash_dir = new SCoord[prot_seq_len+8];

    for (i=0; i<=prot_seq_len; i++)
    {
        residue_clash[i] = 0;
        res_clash_dir[i] = SCoord(0,0,0);
    }
    #endif

    for (i=0; i<sphres; i++)
    {
        if (!reaches_spheroid[i]) continue;
        if (!protein->aa_ptr_in_range(reaches_spheroid[i])) continue;
        reaches_spheroid[i]->clear_atom_binding_energies();
        int resno = reaches_spheroid[i]->get_residue_no();

        if ((reaches_spheroid[i]->has_hbond_acceptors() && ligand->has_hbond_donors())
            || (reaches_spheroid[i]->has_hbond_donors() && ligand->has_hbond_acceptors())
           )
        {
            int aai, li, aan = reaches_spheroid[i]->get_atom_count(), ln = ligand->get_atom_count();
            for (li=0; li<ln; li++)
            {
                Atom* a = ligand->get_atom(li);
                if (!a) continue;
                int aZ = a->get_Z();

                if (aZ > 10 && a->get_family() == CHALCOGEN)
                {
                    if (met && !a->get_charge() && met->get_atom_count())
                    {
                        Atom* mtl = met->get_nearest_atom(a->get_location());
                        if (mtl && mtl->is_metal() && mtl->distance_to(a) < _INTERA_R_CUTOFF)
                        {
                            Atom* H = a->is_bonded_to("H");
                            if (H)
                            {
                                ligand->delete_atom(H);
                                a->increment_charge(-1);
                            }
                        }
                    }
                }

                float apol = a->is_polar();
                if (fabs(apol) < hydrophilicity_cutoff) continue;

                if (sgn(apol) > 0 && aZ > 1) continue;
                for (aai=0; aai<aan; aai++)
                {
                    Atom* b = reaches_spheroid[i]->get_atom(aai);
                    if (!b) continue;
                    float bpol = b->is_polar();
                    if (fabs(bpol) < hydrophilicity_cutoff) continue;
                    int bZ = b->get_Z();
                    if (sgn(bpol) > 0 && bZ > 1) continue;
                    if (sgn(apol) == sgn(bpol)) continue;
                    float r = a->distance_to(b);
                    if (r > 6) continue;

                    Atom* H = (aZ==1) ? a : b;
                    Atom* target = (aZ==1) ? b : a;
                    Atom* targetH = target->is_bonded_to("H");
                    Bond* bond = H->get_bond_by_idx(0);
                    if (!bond) continue;
                    Atom* heavy = bond->atom2;
                    if (!heavy) continue;
                    if (heavy->get_family() != CHALCOGEN) continue;
                    if (heavy->is_pi()) continue;

                    bond = heavy->get_bond_by_idx(0);
                    if (!bond || !bond->atom2 || !bond->can_rotate) continue;

                    SCoord axis = heavy->get_location().subtract(bond->atom2->get_location());
                    float theta = find_angle_along_vector(H->get_location(), target->get_location(), heavy->get_location(), axis);
                    if (targetH && heavy->distance_to(targetH) < target->distance_to(H)) theta += M_PI;

                    Point newHloc1 = rotate3D(H->get_location(), heavy->get_location(), axis, theta);
                    Point newHloc2 = rotate3D(H->get_location(), heavy->get_location(), axis, -theta);

                    if (newHloc1.get_3d_distance(target->get_location()) < newHloc2.get_3d_distance(target->get_location()))
                        H->move(newHloc1);
                    else H->move(newHloc2);
                }
            }
        }

        mc_bpotential = 0;
        Interaction lb = ligand->get_intermol_binding(reaches_spheroid[i], false);
        float clash = lb.repulsive;
        if (ligand->clash_worst > worst_energy)
        {
            worst_energy = ligand->clash_worst;
            worst_clash_1 = ligand->clash1;
            worst_clash_2 = ligand->clash2;
        }
        if (/* lb.summed() < 0 && */ clash > worst_nrg_aa)
        {
            worst_nrg_aa = clash;
            #if _dbg_worst_energy
            cout << reaches_spheroid[i]->get_name() << " binding strength " << lb << " clashes " << clash << " updates worst_nrg_aa to " << worst_nrg_aa << endl;
            #endif
        }

        #if compute_clashdirs
        if (lb > 0 && ligand->clash1 && ligand->clash2)
        {
            residue_clash[resno] += lb;
            SCoord clashdir = ligand->clash2->get_location().subtract(ligand->clash1->get_location());
            clashdir.r = ligand->clash1->get_vdW_radius() + ligand->clash2->get_vdW_radius() - ligand->clash1->distance_to(ligand->clash2);
            res_clash_dir[resno] = res_clash_dir[resno].add(clashdir);
        }
        #endif

        #if include_residue_eclipses
        lb -= fmax(reaches_spheroid[i]->total_eclipses() - reaches_spheroid[i]->initial_eclipses, 0);
        #endif

        #if _dbg_51e2_ionic
        if (resno == 262)
        {
            cout << endl << resno << " charge " << reaches_spheroid[i]->get_charge()
                << " vs. ligand charge " << ligand->get_charge()
                << ": " << lb << endl << endl;
        }
        #endif

        if (lb.summed() > 500) lb = 0;
        lmkJmol[metcount] = lb.summed();
        lmc[metcount] = -mc_bpotential / missed_connection.r;

        lmba1n[metcount] = ligand->best_intera ? ligand->best_intera->name : nullptr;
        lmba2n[metcount] = ligand->best_other_intera ? ligand->best_other_intera->name : nullptr;
        lmca1n[metcount] = ligand->clash1 ? ligand->clash1->name : nullptr;
        lmca2n[metcount] = ligand->clash2 ? ligand->clash2->name : nullptr;

        BallesterosWeinstein bw = protein->get_bw_from_resno(resno);
        if (bw.helix_no)
            sprintf(metrics[metcount], "%s%d(%d.%d)", reaches_spheroid[i]->get_3letter(), resno, bw.helix_no, bw.member_no);
        else
            sprintf(metrics[metcount], "%s%d", reaches_spheroid[i]->get_3letter(), resno);
        // cout << metrics[metcount] << ": " << lb << " . ";

        #if compute_vdw_repulsion
        lmvdWrepl[metcount] = 0;
        lmvdWrepl[metcount] += ligand->get_vdW_repulsion(reaches_spheroid[i]);
        /*for (j=0; j<sphres; j++)
        {
            if (j == i) continue;
            mvdWrepl[metcount] += reaches_spheroid[i]->get_vdW_repulsion(reaches_spheroid[j]);
        }*/
        limvdWrepl[metcount] = 0;
        #endif
        limkJmol[metcount] = 0;
        
        metcount++;
        btot += lb.summed();
        // cout << *(reaches_spheroid[i]) << " adds " << lb << " to btot, making " << btot << endl;

        float lf = ligand->get_intermol_polar_sat(reaches_spheroid[i]);
        pstot += lf;

        #if _dbg_polsat
        cout << *(reaches_spheroid[i]) << " adds " << lf << " to pstot, making " << pstot << endl;
        #endif
    }
    // cout << btot << endl;

    int mcn;
    if (mcn = mtlcoords.size())         // Assignment, not comparison.
    {
        strcpy(metrics[metcount], "Metals");
        float lmb = 0;

        for (i=0; i<mcn; i++)
        {
            if (!mtlcoords[i].mtl) continue;
            Molecule lm("MTL");
            lm.add_existing_atom(mtlcoords[i].mtl);
            float f = ligand->get_intermol_binding(&lm).summed();
            btot += f;
            lmb += f;
        }

        lmkJmol[metcount] = lmb;
    }

    #if _peratom_audit
    cout << endl << "Interatomic Audit:" << endl;
    cout << "Total energy: " << -btot << endl;
    int ian = interaudit.size(), iai;
    for (iai=0; iai<ian; iai++) cout << interaudit[iai] << endl;
    cout << endl << endl;
    interauditing = false;
    #endif

    if (btot > 100*ligand->get_atom_count()) btot = 0;

    kJmol          = btot;
    ikJmol          = 0;
    polsat           = pstot;
    metric            = new char*[metcount+4];
    this->mkJmol       = new float[metcount];
    this->imkJmol       = new float[metcount];
    #if compute_vdw_repulsion
    this->mvdWrepl       = new float[metcount];
    this->imvdWrepl       = new float[metcount];
    #endif
    this->mb_atom1_name     = new const char*[metcount];
    this->mb_atom2_name      = new const char*[metcount];
    this->mc_atom1_name       = new const char*[metcount];
    this->mc_atom2_name        = new const char*[metcount];
    #if compute_missed_connections
    this->missed_connections = new float[metcount];
    #endif
    ligand_self = ligand->get_intermol_binding(ligand).summed() - ligand->get_base_clashes() - ligand->total_eclipses();
    A100 = protein->A100();
    kJmol += ligand_self;
    #if _dbg_internal_energy
    cout << "Ligand internal = " << ligand_self << endl;
    #endif
    protclash = protein->get_rel_int_clashes();

    int itn;
    for (itn=0; itn<_INTER_TYPES_LIMIT; itn++)
    {
        i = itn;
        bytype[i] = total_binding_by_type[i];
        ibytype[i] = init_total_binding_by_type[i];
        ikJmol += init_total_binding_by_type[i];
        kJmol += fin_total_binding_by_type[i];
        // cout << drcount << "|" << i << " ";
    }

    // Populate the array.
    for (i=0; i<metcount; i++)
    {
        metric[i] = new char[max(8,(int)strlen(metrics[i])+4)];
        strcpy(metric[i], metrics[i]);
        mkJmol[i] = lmkJmol[i];
        imkJmol[i] = limkJmol[i];
        #if compute_vdw_repulsion
        mvdWrepl[i] = lmvdWrepl[i];
        imvdWrepl[i] = limvdWrepl[i];
        #endif
        mb_atom1_name[i] = lmba1n[i];
        mb_atom2_name[i] = lmba2n[i];
        mc_atom1_name[i] = lmca1n[i];
        mc_atom2_name[i] = lmca2n[i];
        #if compute_missed_connections
        missed_connections[i] = lmc[i];
        #endif
        // cout << "*" << metric[i] << ": " << mkJmol[i] << endl;
    }

    // Terminate with an empty string and a null pointer.
    metric[i] = new char[2];
    metric[i][0] = 0;
    metric[i+1] = 0;

    std::ostringstream lpdbdat;

    // Prepare a partial PDB of the ligand atoms and all involved residue sidechains.
    n = ligand->get_atom_count();
    int offset = n;
    for (l=0; l<n; l++)
    {
        Atom* a = ligand->get_atom(l);
        if (!a) continue;
        a->residue = pose;
        a->stream_pdb_line(lpdbdat, 9000+l, true);
    }

    if (waters)
    {
        for (k=0; waters[k]; k++)
        {
            for (l=0; l<3; l++)
            {
                Atom* a = waters[k]->get_atom(l);
                if (!a) continue;
                a->residue = pose;
                a->stream_pdb_line(lpdbdat, 9000+offset+l+3*k, true);
            }
        }
    }

    int en = protein->get_end_resno();
    int resno;
    for (resno = protein->get_start_resno(); resno <= en; resno++)
    {
        AminoAcid* laa = protein->get_residue(resno);
        if (!laa) continue;
        if (!laa->been_flexed)
        {
            if (laa->distance_to(ligand) > 5) continue;
            for (k=0; reaches_spheroid[k]; k++)
            {
                if (!protein->aa_ptr_in_range(reaches_spheroid[k])) continue;
                if (reaches_spheroid[k] == laa) goto _afterall;
            }
            continue;
        }
        _afterall:
        n = laa->get_atom_count();
        for (l=0; l<n; l++)
        {
            laa->get_atom(l)->stream_pdb_line(
                lpdbdat,
                laa->atno_offset+l
            );
        }
    }

    if (mtlcoords.size())
    {
        for (l=0; l<mtlcoords.size(); l++)
        {
            mtlcoords[l].mtl->stream_pdb_line(
                lpdbdat,
                9900+l
            );
        }
    }

    this->pdbdat = lpdbdat.str();
}

std::ostream& operator<<(std::ostream& output, const DockResult& dr)
{
    int l;

    if (dr.isomer.length())
    {
        output << "Isomer: " << dr.isomer << endl << endl;
    }

    if (dr.out_per_res_e || dr.out_per_btyp_e)
    {
        output << "# Binding energies:" << endl;
    }

    if (dr.out_per_res_e)
    {
        output << "BENERG:" << endl;
        for (	l=0;

                dr.mkJmol
                && dr.metric
                && dr.metric[l]
                && dr.metric[l][0]
                ;

                l++
            )
        {
            if (fabs(dr.mkJmol[l]) < dr.out_itemized_e_cutoff) continue;

            if (dr.do_output_colors) colorize(dr.mkJmol[l]);
            output << dr.metric[l] << ": " << -dr.mkJmol[l]*dr.energy_mult;
            if (dr.display_binding_atoms) output << " |";
            if (dr.display_binding_atoms) output << " " << (dr.mb_atom1_name[l] ? dr.mb_atom1_name[l] : "-");
            if (dr.display_binding_atoms) output << " " << (dr.mb_atom2_name[l] ? dr.mb_atom2_name[l] : "-");
            if (dr.display_clash_atoms) output << " |";
            if (dr.display_clash_atoms) output << " " << (dr.mc_atom1_name[l] ? dr.mc_atom1_name[l] : "-");
            if (dr.display_clash_atoms && dr.mc_atom1_name[l] && dr.mc_atom2_name[l]) output << " clash";
            if (dr.display_clash_atoms) output << " " << (dr.mc_atom2_name[l] ? dr.mc_atom2_name[l] : "-");
            output << endl;
            if (dr.do_output_colors) colorless();
        }
        output << endl;
        // output << "Worst energy: " << dr.worst_energy << endl;
    }

    if (dr.out_per_btyp_e) for (l=0; l<_INTER_TYPES_LIMIT; l++)
    {
        if (fabs(dr.bytype[l]) < dr.out_itemized_e_cutoff) continue;

        char lbtyp[64];
        switch (l+covalent)
        {
            case covalent:
                continue; /*strcpy(lbtyp, "Total covalent: ");		break;*/
            case ionic:
                strcpy(lbtyp, "Total ionic: ");
                break;
            case hbond:
                strcpy(lbtyp, "Total H-bond: ");
                break;
            case pi:
                strcpy(lbtyp, "Total pi stack: ");
                break;
            case polarpi:
                strcpy(lbtyp, "Total polar-pi and cation-pi: ");
                break;
            case mcoord:
                strcpy(lbtyp, "Total metal coordination: ");
                break;
            case vdW:
                strcpy(lbtyp, "Total van der Waals: ");
                break;
            default:
                goto _btyp_unassigned;
        }

        if (dr.do_output_colors) colorize(dr.bytype[l]);
        output << lbtyp << -dr.bytype[l]*dr.energy_mult << endl;
        if (dr.do_output_colors) colorless();
    }
    output << endl;

_btyp_unassigned:

    if (dr.do_output_colors) colorize(dr.kJmol);
    output << "Total: " << -dr.kJmol*dr.energy_mult << endl;
    output << "Worst atom clash: " << dr.worst_energy*dr.energy_mult << endl;
    if (dr.do_output_colors) colorless();

    #if include_eclipses
    if (dr.out_lig_int_e) output << "Ligand internal energy: " << -dr.ligand_self*dr.energy_mult << endl << endl;
    #endif

    if (dr.miscdata.size())
    {
        output << dr.miscdata << endl;
    }

    output << "A100 score: " << dr.A100 << endl;
    output << endl;

    if (dr.out_lig_pol_sat)
    {
        output << "Ligand polar satisfaction: " << dr.polsat << endl;
        output << endl;
    }

    if (dr.out_prox) output << "Proximity: " << dr.proximity << endl << endl;

    if (dr.out_pro_clash) output << "Protein clashes: " << dr.protclash << endl << endl;

    #if compute_missed_connections
    if (dr.out_mc)
    {
        output << "# Missed Connections" << endl << "MC:" << endl;
        output << endl;
        
        for (	l=0;

                dr.metric
                && dr.metric[l]
                && dr.metric[l][0];

                l++
            )
        {
            if (fabs(dr.missed_connections[l]) < dr.out_itemized_e_cutoff) continue;
            if (dr.metric[l]) output << dr.metric[l] << ": " << dr.missed_connections[l]*dr.energy_mult << endl;
        }
        output << endl;
    }
    #endif

    #if compute_vdw_repulsion
    if (dr.out_vdw_repuls)
    {
        output << "# van der Waals repulsion" << endl << "vdWRPL:" << endl;
        output << endl;
        
        for (	l=0;

                dr.metric
                && dr.metric[l]
                && dr.metric[l][0];

                l++
            )
        {
            if (fabs(dr.mvdWrepl[l]) < dr.out_itemized_e_cutoff) continue;
            if (dr.metric[l]) output << dr.metric[l] << ": " << dr.mvdWrepl[l]*dr.energy_mult << endl;
        }
        output << endl;
    }
    #endif

    if (dr.include_pdb_data)
    {
        if (!dr.pdbdat.length())
        {
            output << "WARNING: Missing PDB data." << endl;
        }
        else
        {
            output << "# PDB Data" << endl << "PDBDAT:" << endl;

            output << dr.pdbdat << endl;

            output << "TER" << endl << "END" << endl << endl << endl;
        }
    }

    return output;
}
