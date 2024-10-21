
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    char tstname[1024] = "sdf/floralozone.sdf";
    if (argc > 1) strcpy(tstname, argv[1]);
    Molecule m(tstname);
    int nloaded;
    float density = (argc > 2) ? atoi(argv[2]) : 20;

    char buffer[65536];
    FILE* pf = fopen(tstname, "rb");
    if (!pf)
    {
        // cout << "Failed to open " << tstname << " for reading." << endl;
        // return -1;
        m.from_smiles(tstname);
    }
    else
    {
        int fyrw = fread(buffer, 1, 65535, pf);
        fclose(pf);
        cout << "Read data.\n";
        nloaded = m.from_sdf(buffer);
    }
    cout << "Loaded " << nloaded << " atoms of " << tstname << " into molecule.\n";

    const Point* vertices = m.obtain_vdW_surface(density);
    Atom** va = m.get_vdW_vertex_atoms();
    if (!vertices)
    {
        cout << "No verts returned." << endl;
        return -1;
    }

    int i;
    for (i=0; vertices[i].x || vertices[i].y || vertices[i].z; i++)
    {
        cout << va[i]->name << " " << vertices[i] << endl;
    }
    int atomZ[i+4];

    pf = fopen("tmp/surface", "w");
    fprintf(pf, "var shape = new NGL.Shape( \"shape\" );\n");
    fprintf(pf, "var sphereBuffer = new NGL.SphereBuffer(\n{\n");
    fprintf(pf, "    position: new Float32Array( [ ");
    for (i=0; vertices[i].x || vertices[i].y || vertices[i].z; i++)
    {
        if (i) fprintf(pf, ", ");
        fprintf(pf, "%f, %f, %f", vertices[i].x, vertices[i].y, vertices[i].z);
        atomZ[i] = va[i]->get_Z();
    }
    fprintf(pf, " ] ),\n");
    fprintf(pf, "    color: new Float32Array( [ ");
    for (i=0; vertices[i].x || vertices[i].y || vertices[i].z; i++)
    {
        if (i) fprintf(pf, ", ");
        switch (atomZ[i])
        {
            case 1: 
            if (fabs(va[i]->is_polar()) > hydrophilicity_cutoff) fprintf(pf, "0.7, 0.9, 1");
            else fprintf(pf, "0.9, 0.95, 1");
            break;
            case 5: fprintf(pf, "0.5, 0.4, 0.3"); break;
            case 6: fprintf(pf, "0.6, 0.6, 0.6"); break;
            case 7: fprintf(pf, "0, 0, 1"); break;
            case 8: fprintf(pf, "0, 1, 1"); break;
            case 9: fprintf(pf, "0.7, 1, 0.1"); break;
            case 14: fprintf(pf, "0.1, 0.1, 0.5"); break;
            case 15: fprintf(pf, "1, 0.5, 0.6"); break;
            case 16: fprintf(pf, "1, 0.8, 0.1"); break;
            case 17: fprintf(pf, "0.1, 1, 0.1"); break;
            default: fprintf(pf, "1, 0.1, 1");
        }
    }
    fprintf(pf, " ] ),\n");
    fprintf(pf, "    radius: new Float32Array( [ ");
    for (i=0; vertices[i].x || vertices[i].y || vertices[i].z; i++)
    {
        if (i) fprintf(pf, ", ");
        fprintf(pf, "0.03");
    }
    fprintf(pf, " ] ),\n");
    fprintf(pf, "} );\n");
    fprintf(pf, "shape.addBuffer( sphereBuffer );\n");
    fprintf(pf, "var shapeComp = stage.addComponentFromObject( shape );\n");
    fprintf(pf, "shapeComp.addRepresentation( \"buffer\" );\n");
    fprintf(pf, "shapeComp.autoView();\n");
    fprintf(pf, "\n");
    fclose(pf);

    return 0;
}

