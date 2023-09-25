
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

    for (i=0; i<end1; i++)
    {
        lmkJmol[i] = limkJmol[i] = lmvdWrepl[i] = limvdWrepl[i] = 0;
    }

    for (i=0; i<_INTER_TYPES_LIMIT; i++) init_total_binding_by_type[i] = 0;
    for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = 0;
    for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

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
        for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

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
        for (i=0; i<qpr; i++)
        {
            int resno = allres[i]->get_residue_no();

            #if use_trip_switch
            bool is_trip_i = false;
            for (n=0; n<tripswitch_clashables.size(); n++)
                if (tripswitch_clashables[n] == resno)
                {
                    is_trip_i = true;
                    break;
                }
            #endif

            for (j=0; j<qpr; j++)
            {
                if (j == i) continue;
                int jres = allres[j]->get_residue_no();

                #if use_trip_switch
                bool is_trip_j = false;
                if (is_trip_i)
                    for (n=0; n<tripswitch_clashables.size(); n++)
                        if (tripswitch_clashables[n] == jres)
                        {
                            is_trip_j = true;
                            break;
                        }
                if (!is_trip_j) continue;
                #endif

                float f = postaa[i]->get_intermol_binding(postaa[j], j==0);
                #if use_trip_switch
                if (f < 0)
                {
                    tripclash -= f;
                    f = 0;
                }
                #endif
            }
        }
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

    for (i=0; i<sphres; i++)
    {
        if (!reaches_spheroid[i]) continue;
        if (!protein->aa_ptr_in_range(reaches_spheroid[i])) continue;
        reaches_spheroid[i]->clear_atom_binding_energies();
        int resno = reaches_spheroid[i]->get_residue_no();

        float lb = ligand->get_intermol_binding(reaches_spheroid[i], false);
        if (lb < -maxclash) maxclash -= lb;

        if (differential_dock)
        {
            mkJmol[metcount] = final_binding[resno] + lb;
        }
        else
        {
            if (lb > 500) lb = 0;
            mkJmol[metcount] = lb;
        }

        #if active_persistence
        res_kJmol[resno] = lb;
        #endif

        sprintf(metrics[metcount], "%s%d", reaches_spheroid[i]->get_3letter(), resno);
        // cout << metrics[metcount] << ": " << lb << " . ";

        if (differential_dock)
        {
            imkJmol[metcount] = initial_binding[resno];
            mvdWrepl[metcount] = final_vdWrepl[resno];
            imvdWrepl[metcount] = initial_vdWrepl[resno];
        }
        else
        {
            mvdWrepl[metcount] = 0;
            mvdWrepl[metcount] += ligand->get_vdW_repulsion(reaches_spheroid[i]);
            /*for (j=0; j<sphres; j++)
            {
                if (j == i) continue;
                mvdWrepl[metcount] += reaches_spheroid[i]->get_vdW_repulsion(reaches_spheroid[j]);
            }*/
            imkJmol[metcount] = 0;
            imvdWrepl[metcount] = 0;
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

    if (btot > 60*ligand->get_atom_count()) btot = 0;
    if (differential_dock && (maxclash > individual_clash_limit)) btot = -Avogadro;

    kJmol = (differential_dock && (maxclash > individual_clash_limit)) ? -Avogadro : btot;
    ikJmol = 0;
    polsat  = pstot;
    metric   = new char*[metcount+4];
    this->mkJmol    = new float[metcount];
    this->imkJmol    = new float[metcount];
    this->mvdWrepl    = new float[metcount];
    this->imvdWrepl    = new float[metcount];
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
