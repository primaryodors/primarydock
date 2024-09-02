
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

DockResult::DockResult(Protein* protein, Molecule* ligand, Point size, int* addl_resno, int drcount)
{
    int end1 = SPHREACH_MAX+4;
    AminoAcid* reaches_spheroid[end1];
    int sphres = protein->get_residues_can_clash_ligand(reaches_spheroid, ligand, ligand->get_barycenter(), size, addl_resno);
    // cout << "sphres " << sphres << endl;
    Molecule* met = protein->metals_as_molecule();
    int i, j;

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
    char* lma1n[end1];
    char* lma2n[end1];
    float lmc[end1];

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

    std::vector<AminoAcid*> allres = protein->get_residues_near(ligand->get_barycenter(), 100000);
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

        mc_bpotential = 0;
        float lb = ligand->get_intermol_binding(reaches_spheroid[i], false);
        float clash = ligand->get_intermol_clashes(reaches_spheroid[i]);
        if (ligand->clash_worst > worst_energy)
        {
            worst_energy = ligand->clash_worst;
            worst_clash_1 = ligand->clash1;
            worst_clash_2 = ligand->clash2;
        }
        if (lb < 0 && clash > worst_nrg_aa)
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

        if (lb > 500) lb = 0;
        lmkJmol[metcount] = lb;
        lmc[metcount] = -mc_bpotential / missed_connection.r;

        lma1n[metcount] = ligand->clash1 ? ligand->clash1->name : nullptr;
        lma2n[metcount] = ligand->clash2 ? ligand->clash2->name : nullptr;

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
        btot += lb;
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
            float f = ligand->get_intermol_binding(&lm);
            btot += f;
            lmb += f;
        }

        lmkJmol[metcount] = lmb;
        metcount++;
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
    this->m_atom1_name     = new const char*[metcount];
    this->m_atom2_name      = new const char*[metcount];
    #if compute_missed_connections
    this->missed_connections = new float[metcount];
    #endif
    #if include_residue_eclipses
    ligand_self = ligand->get_intermol_binding(ligand) - ligand->total_eclipses();
    #endif
    A100 = protein->A100();
    lig_occlusion = ligand->occlusion;
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
    // cout << endl;
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Filled btypes." << endl;
    #endif

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
        m_atom1_name[i] = lma1n[i];
        m_atom2_name[i] = lma2n[i];
        #if compute_missed_connections
        missed_connections[i] = lmc[i];
        #endif
        // cout << "*" << metric[i] << ": " << mkJmol[i] << endl;
    }

    // Terminate with an empty string and a null pointer.
    metric[i] = new char[2];
    metric[i][0] = 0;
    metric[i+1] = 0;
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "More metrics or something idfk." << endl;
    #endif

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
            if (dr.display_clash_atom1) output << " " << (dr.m_atom1_name[l] ? dr.m_atom1_name[l] : "-");
            if (dr.display_clash_atom2) output << " " << (dr.m_atom2_name[l] ? dr.m_atom2_name[l] : "-");
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
    output << "Total: " << -dr.kJmol*dr.energy_mult << endl << endl;
    if (dr.do_output_colors) colorless();

    #if include_eclipses
    if (dr.out_lig_int_e) output << "Ligand internal energy: " << -dr.ligand_self*dr.energy_mult << endl << endl;
    #endif

    if (dr.miscdata.size())
    {
        output << dr.miscdata << endl;
    }

    output << "A100 score: " << dr.A100 << endl;
    output << "Ligand occlusion: " << dr.lig_occlusion << endl;
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
