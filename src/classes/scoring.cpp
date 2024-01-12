
#include "group.h"
#include "scoring.h"

float init_total_binding_by_type[_INTER_TYPES_LIMIT];
float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

float* initial_binding;
float* initial_vdWrepl;

DockResult::DockResult()
{
    ;
}

DockResult::DockResult(Protein* protein, Molecule* ligand, Point size, int* addl_resno, int drcount, bool differential_dock)
{
    int end1 = protein->get_end_resno()+1;
    AminoAcid* reaches_spheroid[end1];
    int sphres = protein->get_residues_can_clash_ligand(reaches_spheroid, ligand, ligand->get_barycenter(), size, addl_resno);
    // cout << "sphres " << sphres << endl;
    float maxclash = 0;
    Molecule* met = protein->metals_as_molecule();
    int i, j;

    char metrics[end1][10];
    float lmkJmol[end1];
    float limkJmol[end1];
    float lmvdWrepl[end1];
    float limvdWrepl[end1];
    int metcount = 0;
    float btot = 0;
    float pstot = 0;
    char* lma1n[end1];
    char* lma2n[end1];

    for (i=0; i<end1; i++)
    {
        lmkJmol[i] = limkJmol[i] = lmvdWrepl[i] = limvdWrepl[i] = 0;
    }

    worst_energy = 0;

    // if (debug) *debug << "Pose " << pose << " pathnode " << nodeno /*<< " clashes " << clash*/ << endl;

    ligand->clear_atom_binding_energies();

    if (met)
    {
        met->clear_atom_binding_energies();
        float lb = ligand->get_intermol_binding(met);
        strcpy(metrics[metcount], "Metals");
        lmkJmol[metcount] = lb;
        limkJmol[metcount] = 0;								// TODO

        lmvdWrepl[metcount] = 0;
        lmvdWrepl[metcount] += ligand->get_vdW_repulsion(met);		// TODO: Include repulsions with non-mcoord side chains.

        limvdWrepl[metcount] = 0;							// TODO

        metcount++;

        btot += lb;
        // cout << "Metal adds " << lb << " to btot, making " << btot << endl;
    }

    float final_binding[end1];
    float final_vdWrepl[end1];
    for (i=0; i<end1; i++) final_binding[i] = final_vdWrepl[i] = 0;

    std::vector<AminoAcid*> allres = protein->get_residues_near(ligand->get_barycenter(), 100000);
    int qpr = allres.size();
    Molecule* postaa[qpr];
    postaa[0] = ligand;
    for (i=0; i<qpr; i++)
    {
        postaa[i+1] = reinterpret_cast<Molecule*>(allres[i]);
    }

    if (differential_dock)
    {
        for (i=0; i<qpr+1; i++)
        {
            int resno = i ? (allres[i-1]->get_residue_no()) : 0;
            #if _DBG_TOOLARGE_DIFFNUMS
            std::string ibdbg = to_string(resno) + (std::string)" ibdbg:\n";
            #endif

            #if use_trip_switch
            bool is_trip_i = false;
            for (n=0; n<tripswitch_clashables.size(); n++)
                if (tripswitch_clashables[n] == resno)
                {
                    is_trip_i = true;
                    break;
                }
            #endif

            for (j=0; j<qpr+1; j++)
            {
                if (j == i) continue;
                int jres = j ? allres[j-1]->get_residue_no() : 0;

                #if use_trip_switch
                bool is_trip_j = false;
                if (is_trip_i && j)
                    for (n=0; n<tripswitch_clashables.size(); n++)
                        if (tripswitch_clashables[n] == jres)
                        {
                            is_trip_j = true;
                            break;
                        }
                #endif

                float f = postaa[i]->get_intermol_binding(postaa[j], j==0);
                #if use_trip_switch
                if (f < 0 && is_trip_j)
                {
                    tripclash -= f;
                    f = 0;
                }
                #endif
                final_binding[resno] += f;

                #if _DBG_TOOLARGE_DIFFNUMS
                if (f) ibdbg += to_string(postaa[j]->get_atom(0)->residue) + (std::string)" " + to_string(f) + (std::string)"\n";
                #endif

                final_vdWrepl[resno] += postaa[i]->get_vdW_repulsion(postaa[j]);
            }

            #if _DBG_TOOLARGE_DIFFNUMS
            if (fabs(final_binding[resno]) >= 200) cout << ibdbg << endl;
            #endif
        }
    }
    else
    {      
        #if use_trip_switch          
        for (i=0; i<qpr; i++)
        {
            int resno = allres[i]->get_residue_no();

            bool is_trip_i = false;
            for (n=0; n<tripswitch_clashables.size(); n++)
                if (tripswitch_clashables[n] == resno)
                {
                    is_trip_i = true;
                    break;
                }

            for (j=0; j<qpr; j++)
            {
                if (j == i) continue;
                int jres = allres[j]->get_residue_no();

                bool is_trip_j = false;
                if (is_trip_i)
                    for (n=0; n<tripswitch_clashables.size(); n++)
                        if (tripswitch_clashables[n] == jres)
                        {
                            is_trip_j = true;
                            break;
                        }
                if (!is_trip_j) continue;

                float f = postaa[i]->get_intermol_binding(postaa[j], j==0);
                if (f < 0)
                {
                    tripclash -= f;
                    f = 0;
                }
            }
        }
        #endif
    }

    for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = total_binding_by_type[i];

    #if active_persistence
    float res_kJmol[end1];
    for (i=0; i<end1; i++) res_kJmol[i] = 0;
    #endif

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
    residue_clash = new float[prot_seq_len+8];
    res_clash_dir = new SCoord[prot_seq_len+8];

    for (i=0; i<=prot_seq_len; i++)
    {
        residue_clash[i] = 0;
        res_clash_dir[i] = SCoord(0,0,0);
    }

    for (i=0; i<sphres; i++)
    {
        if (!reaches_spheroid[i]) continue;
        if (!protein->aa_ptr_in_range(reaches_spheroid[i])) continue;
        reaches_spheroid[i]->clear_atom_binding_energies();
        int resno = reaches_spheroid[i]->get_residue_no();

        float lb = ligand->get_intermol_binding(reaches_spheroid[i], false);
        if (lb < -maxclash) maxclash = -lb;

        if (lb > 0 && ligand->clash1 && ligand->clash2)
        {
            residue_clash[resno] += lb;
            SCoord clashdir = ligand->clash2->get_location().subtract(ligand->clash1->get_location());
            clashdir.r = ligand->clash1->get_vdW_radius() + ligand->clash2->get_vdW_radius() - ligand->clash1->distance_to(ligand->clash2);
            res_clash_dir[resno] = res_clash_dir[resno].add(clashdir);
        }

        lb -= reaches_spheroid[i]->total_eclipses() - reaches_spheroid[i]->initial_eclipses;

        #if _dbg_51e2_ionic
        if (resno == 262)
        {
            cout << endl << resno << " charge " << reaches_spheroid[i]->get_charge()
                << " vs. ligand charge " << ligand->get_charge()
                << ": " << lb << endl << endl;
        }
        #endif

        if (differential_dock)
        {
            lmkJmol[metcount] = final_binding[resno] + lb;
        }
        else
        {
            if (lb > 500) lb = 0;
            lmkJmol[metcount] = lb;

            if (-lb > worst_energy) worst_energy = -lb;
        }

        lma1n[metcount] = ligand->clash1 ? ligand->clash1->name : nullptr;
        lma2n[metcount] = ligand->clash2 ? ligand->clash2->name : nullptr;

        #if active_persistence
        res_kJmol[resno] = lb;
        #endif

        sprintf(metrics[metcount], "%s%d", reaches_spheroid[i]->get_3letter(), resno);
        // cout << metrics[metcount] << ": " << lb << " . ";

        if (differential_dock)
        {
            limkJmol[metcount] = initial_binding[resno];
            lmvdWrepl[metcount] = final_vdWrepl[resno];
            limvdWrepl[metcount] = initial_vdWrepl[resno];
        }
        else
        {
            lmvdWrepl[metcount] = 0;
            lmvdWrepl[metcount] += ligand->get_vdW_repulsion(reaches_spheroid[i]);
            /*for (j=0; j<sphres; j++)
            {
                if (j == i) continue;
                mvdWrepl[metcount] += reaches_spheroid[i]->get_vdW_repulsion(reaches_spheroid[j]);
            }*/
            limkJmol[metcount] = 0;
            limvdWrepl[metcount] = 0;
        }
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
        for (i=0; i<mcn; i++)
        {
            if (!mtlcoords[i].mtl) continue;
            Molecule lm("MTL");
            lm.add_existing_atom(mtlcoords[i].mtl);
            float f = ligand->get_intermol_binding(&lm);
            btot += f;
        }
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
    if (differential_dock && (maxclash > individual_clash_limit)) btot = -Avogadro;

    kJmol = (differential_dock && (maxclash > individual_clash_limit)) ? -Avogadro : btot;
    ikJmol = 0;
    polsat  = pstot;
    metric   = new char*[metcount+4];
    this->mkJmol    = new float[metcount];
    this->imkJmol    = new float[metcount];
    this->mvdWrepl    = new float[metcount];
    this->imvdWrepl    = new float[metcount];
    this->m_atom1_name  = new const char*[metcount];
    this->m_atom2_name  = new const char*[metcount];
    ligand_self = ligand->get_intermol_binding(ligand) - ligand->total_eclipses();
    kJmol += ligand_self;
    #if _dbg_internal_energy
    cout << "Ligand internal = " << ligand_self << endl;
    #endif
    #if use_trip_switch
    tripswitch  = tripclash;
    #endif
    protclash = protein->get_rel_int_clashes();

    int itn;
    for (itn=0; itn<_INTER_TYPES_LIMIT; itn++)
    {
        i = itn;
        bytype[i] = differential_dock ? fin_total_binding_by_type[i] : total_binding_by_type[i];
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
        mvdWrepl[i] = lmvdWrepl[i];
        imvdWrepl[i] = limvdWrepl[i];
        m_atom1_name[i] = lma1n[i];
        m_atom2_name[i] = lma2n[i];
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

    if (differential_dock)
    {
        output << "# Binding energies: delta = with ligand minus without ligand." << endl;
    }
    else
    {
        output << "# Binding energies:" << endl;
    }

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
        if (fabs(dr.mkJmol[l]) < 0.001) continue;

        if (differential_dock)
        {
            if (dr.metric[l]) output << dr.metric[l]
                        << ": " << -(dr.mkJmol[l] - dr.imkJmol[l])*dr.energy_mult
                        << " = " << -dr.mkJmol[l]*dr.energy_mult
                        << " minus " << -dr.imkJmol[l]*dr.energy_mult
                        << endl;
        }
        else
        {
            if (dr.do_output_colors) colorize(dr.mkJmol[l]);
            output << dr.metric[l] << ": " << -dr.mkJmol[l]*dr.energy_mult;
            if (dr.display_clash_atom1) output << " " << (dr.m_atom1_name[l] ? dr.m_atom1_name[l] : "-");
            if (dr.display_clash_atom2) output << " " << (dr.m_atom2_name[l] ? dr.m_atom2_name[l] : "-");
            output << endl;
            if (dr.do_output_colors) colorless();
        }
    }
    output << endl;

    for (l=0; l<_INTER_TYPES_LIMIT; l++)
    {
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

        if (differential_dock)
        {
            output << lbtyp << -(dr.bytype[l] - dr.ibytype[l])*dr.energy_mult
                << " = " << -dr.bytype[l]*dr.energy_mult
                << " minus " << -dr.ibytype[l]*dr.energy_mult
                << endl;
        }
        else
        {
            if (dr.do_output_colors) colorize(dr.bytype[l]);
            output << lbtyp << -dr.bytype[l]*dr.energy_mult << endl;
            if (dr.do_output_colors) colorless();
        }
    }
    output << endl;

_btyp_unassigned:

    if (differential_dock)
    {
        output << "Total: " << -(dr.kJmol - dr.ikJmol)*dr.energy_mult
            << " = " << -dr.kJmol*dr.energy_mult
            << " minus " << -dr.ikJmol*dr.energy_mult
            << endl << endl;
    }
    else
    {
        if (dr.do_output_colors) colorize(dr.kJmol);
        output << "Total: " << -dr.kJmol*dr.energy_mult << endl << endl;
        if (dr.do_output_colors) colorless();
    }

    output << "Ligand internal energy: " << -dr.ligand_self*dr.energy_mult << endl << endl;

    if (dr.miscdata.size())
    {
        output << endl << dr.miscdata << endl;
    }

    if (dr.softrock.size())
    {
        output << "Active Helix Soft Rotations:" << endl
                << dr.softrock << endl;
    }

    output << "Ligand polar satisfaction: " << dr.polsat << endl;
    output << endl;

    output << "Proximity: " << dr.proximity << endl << endl;

    output << "Protein clashes: " << dr.protclash << endl << endl;

    #if use_trip_switch
    if (tripswitch_clashables.size())
    {
        output << "Trip switch: " << dr.tripswitch << endl << endl;
    }
    #endif

    output << "# van der Waals repulsion" << endl << "vdWRPL:" << endl;
    output << endl;
    
    for (	l=0;

            dr.metric
            && dr.metric[l]
            && dr.metric[l][0];

            l++
        )
    {
        if (fabs(dr.mvdWrepl[l]) < 0.001) continue;

        if (differential_dock)
        {
            if (dr.metric[l]) output << dr.metric[l]
                << ": " << (dr.mvdWrepl[l] - dr.imvdWrepl[l])*dr.energy_mult
                << " = " << dr.mvdWrepl[l]*dr.energy_mult
                << " minus " << dr.imvdWrepl[l]*dr.energy_mult
                << endl;
        }
        else
        {
            if (dr.metric[l]) output << dr.metric[l] << ": " << dr.mvdWrepl[l]*dr.energy_mult << endl;
        }
    }
    output << endl;

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
