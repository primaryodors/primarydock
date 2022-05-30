
#include <cstring>
#include <stdio.h>
#include "molecule.h"

#ifndef _AMINACID
#define _AMINACID

enum bb_rot_dir
{
 	N_asc,
 	CA_asc,
 	CA_desc,
 	C_desc
};

class AABondDef
{
public:
 	char aname[7];
 	char bname[7];
 	float cardinality;
 	float acharge;
 	bool can_rotate;
};

class AADef
{
public:
 	char _1let = '\0';
 	char _3let[4] = {};
 	char name[20] = {};
 	char smiles[50] = {};
 	float reach = 2.5;
 	AABondDef** aabonds=0;
 	bool proline_like = false;
};

class MetalCoord
{
public:
 	Atom* metal;
 	Atom** coord_atoms;
 	AminoAcid** coord_res;
 	bool locked = false;

 	Point coord_atom_avg_loc();
};

class AminoAcid : public MolVirtual, public Molecule
{
public:
 	// Constructors.
 	AminoAcid(FILE* instream);
 	AminoAcid(const char letter, AminoAcid* prev_res=0);

 	// Getters and setters.
 	AminoAcid* get_prev() const 				{ 	 	return prev_aa; 						}
 	AminoAcid* get_next() const 				{ 	 	return next_aa; 						}
 	int get_residue_no() const 					{ 	 	return residue_no; 						}
 	char get_letter() const 					{ 	 	return aadef ? aadef->_1let : 0; 		}
 	char* get_3letter() const 					{ 	 	return aadef ? aadef->_3let : 0; 		}
 	AADef* get_aa_definition() const 			{ 	 	return aadef; 							}
 	float get_reach() const 					{ 	 	return aadef ? aadef->reach : 2.5; 		}
 	bool can_reach(AminoAcid* other) const;
 	void set_region(const char* regname) 		{ 	 	strcpy(region, regname); 				}
 	Atom* previous_residue_C();
 	Atom* next_residue_N();

 	// Serialization.
 	int from_pdb(FILE* instream);							// returns number of atoms loaded.
 	void save_pdb(FILE* outstream, int atomno_offset=0);

 	// Spatial functions.
 	void move(Vector move_amt) 					{ 	 	return; 	}
 	void move(Point move_amt) 					{ 	 	return; 	}
 	void recenter(Point new_location) 			{ 	 	return; 	}
 	void rotate(Vector* vector, float theta) 	{ 	 	return; 	}
 	void rotate(LocatedVector vector, float theta);
 	LocatedVector rotate_backbone(bb_rot_dir direction, float angle);	// Return the origin and direction of the rotation axis.
 	LocRotation rotate_backbone_abs(bb_rot_dir direction, float angle);
 	LocRotation* flatten();		// Ensure the peptide bond is coplanar and the CA lies in the same plane. Return LocRotation[2].

 	// Bond functions.
 	bool disulfide_bond(const AminoAcid* bond_to);
 	Bond** get_rotatable_bonds();
 	bool capable_of_inter(intera_type inter);

 	// Intermol functions.
 	float get_intermol_binding(AminoAcid** neighbors, bool backbone_atoms_only = false);

 	// Misc.
 	void delete_sidechain();
 	static Molecule** aas_to_mols(AminoAcid** aas);

 	// Public properties.
 	int strand;
 	int atno_offset=0;
 	MetalCoord* m_mcoord=0;

protected:
 	void load_aa_defs();
 	void copy_loaded_to_object(char letter, int tbdctr, AABondDef** tmpbdefs, bool proline_like);

 	int residue_no=0;
 	char region[25] {};
 	AADef* aadef=0;
 	AminoAcid *prev_aa=0, *next_aa=0;
};

extern AADef aa_defs[26];		// Indexed by letter so 0=A ala, 2=C cys, 3=D asp, etc.
extern char* override_aminos_dat;

std::ostream& operator<<(std::ostream& os, const AminoAcid& aa);

#endif

