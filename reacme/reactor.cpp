#include "reactor.hpp"

Reactor::Reactor(string r, string p, int confs)
{
  numconfs = confs;
  rmol = new OpenBabel::OBMol;
  pmol = new OpenBabel::OBMol;
  rsmi = r;
  psmi = p;
}

Reactor::~Reactor()
{
  delete rmol;
  delete pmol;
}

bool Reactor::Init()
{
  OpenBabel::OBConversion conv;
  if(conv.SetInFormat("smi"))
  {
    if(conv.ReadString(rmol, rsmi) && conv.ReadString(pmol, psmi))
    {
      // generate structure
      if(!Structurize())
      {
        cout << "number of atoms differ between reactant (N=" << rmol->NumAtoms() << ") and product (N=" << pmol->NumAtoms() << ")" << endl;
        return false;
      }
      // obtain class maps
      if(!ObtainClassMaps())
      {
        cout << "atom class count differ between reactant (N=" << rmap.size() << ") and product (N=" << pmap.size() << ")" << endl;
        return false;
      }
      cout << "reactant, # atoms=" << rmol->NumAtoms() << ", # confs=" << rmol->NumConformers() << endl;
      cout << "product,  # atoms=" << pmol->NumAtoms() << ", # confs=" << pmol->NumConformers() << endl;
      return true;
    }
  }
  return false;
}

bool Reactor::Structurize()
{
  SetupMol(rmol);
  SetupMol(pmol);
  return rmol->NumAtoms() == pmol->NumAtoms();
}

void Reactor::SetupMol(OpenBabel::OBMol* mol)
{
  OpenBabel::OBForceField *pFF = OpenBabel::OBForceField::FindForceField("MMFF94");
  if(pFF)
  {
    mol->AddHydrogens();
    builder.Build(*mol);
    if(pFF->Setup(*mol))
    {
      pFF->SteepestDescent(250);
      pFF->WeightedRotorSearch(numconfs, 250);
      pFF->GetConformers(*mol);
    }
  }
}

bool Reactor::ObtainClassMaps()
{
  ObtainClassMapForMol(rmol,&rmap);
  ObtainClassMapForMol(pmol,&pmap);
  return rmap.size() == pmap.size();
}

void Reactor::ObtainClassMapForMol(OpenBabel::OBMol* mol, map<int,int>* mymap)
{
  int idx = -1;
  int cls = -1;
  if(mol->HasData("Atom Class"))
  {
    OpenBabel::OBAtomClassData* adc = (OpenBabel::OBAtomClassData*)mol->GetData("Atom Class");
    FOR_ATOMS_OF_MOL(atom,mol)
    {
      idx = atom->GetIdx();
      cls = adc->GetClass(idx);
      if(cls >=0) (*mymap)[cls] = idx;
    }
  }
}

