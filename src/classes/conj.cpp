
#include "conj.h"


Conjugation::Conjugation()
{
    atoms = new Atom*[1024];
    atoms[0] = nullptr;
}

Conjugation::Conjugation(Atom* from_atom)
{
    atoms = new Atom*[1024];
    atoms[0] = nullptr;
    add_atom(from_atom);
}

Conjugation::~Conjugation()
{
    delete[] atoms;
}

void Conjugation::add_atom(Atom* add)
{
    add_atom(add, add, add);
}

void Conjugation::add_atom(Atom* a, Atom* prev, Atom* orig)
{
    int i, j;

    if (!a->is_pi()) return;

    for (j=0; atoms[j]; j++)
    {
        if (atoms[j] == a) return;
    }

    atoms[j] = a;
    atoms[j+1] = nullptr;
    net_charge_known = false;
    // if (a->conjugation) delete[] a->conjugation;
    a->conjugation = std::shared_ptr<Conjugation>(this);

    Bond* b[16];
    for (i=0; i<16; i++) b[i] = nullptr;
    a->fetch_bonds(b);

    for (i=0; b[i]; i++)
    {
        if (!b[i]->btom) continue;
        if (!b[i]->btom->is_pi()) continue;
        if (b[i]->btom == prev) continue;
        if (b[i]->btom == orig) continue;

        add_atom(b[i]->btom, a, orig);              // RECURSION!
    }
}

float Conjugation::get_net_charge()
{
    if (!net_charge_known)
    {
        int i;
        net_charge = 0;
        for (i=0; atoms[i]; i++)
        {
            net_charge += atoms[i]->get_orig_charge();
        }

        net_charge_known = true;
    }

    return net_charge;
}

Point Conjugation::get_barycenter()
{
    if (!atoms || !atoms[0]) return Point(0,0,0);

    Point result, pt;
    float divisor = 0;
    int i;

    for (i=0; atoms[i]; i++)
    {
        pt = atoms[i]->get_location();
        float atwt = atoms[i]->get_atomic_weight();
        pt.scale(atwt);
        result = result.add(pt);
        divisor += atwt;
    }

    if (divisor) result.multiply(1.0 / divisor);

    return result;
}

Atom* Conjugation::get_nearest_atom(Point pt)
{
    if (!atoms || !atoms[0]) return nullptr;

    float r, best_r;
    Atom* best_a;
    int i;

    for (i=0; atoms[i]; i++)
    {
        r = atoms[i]->get_location().get_3d_distance(pt);
        if (!i || r < best_r)
        {
            best_a = atoms[i];
            best_r = r;
        }
    }

    return best_a;
}

Atom* Conjugation::get_nearest_atom(Conjugation* conj)
{
    /*if (mutual == conj)
    {
        mutual->mutual = nullptr;
        mutual = nullptr;
        return mutual_nearest;
    }*/

    Point bcen1 = get_barycenter(), bcen2 = conj->get_barycenter();
    Atom* a = get_nearest_atom(bcen2), *b = conj->get_nearest_atom(bcen1);

    a = get_nearest_atom(b->get_location());
    b = get_nearest_atom(a->get_location());
    a = get_nearest_atom(b->get_location());
    b = get_nearest_atom(a->get_location());
    a = get_nearest_atom(b->get_location());
    b = get_nearest_atom(a->get_location());

    conj->mutual = this;
    conj->mutual_nearest = b;

    return a;
}
