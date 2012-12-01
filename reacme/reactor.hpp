#ifndef REACTOR_H
#define REACTOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <openbabel/mol.h>
#include <openbabel/builder.h>
#include <openbabel/atomclass.h>
#include <openbabel/generic.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>

using namespace std;

class Reactor
{
  private:
    int numconfs;
    string rsmi;
    string psmi;
    map<int,int> rmap;
    map<int,int> pmap;

    OpenBabel::OBMol* rmol;
    OpenBabel::OBMol* pmol;
    OpenBabel::OBBuilder builder;

    void SetupMol(OpenBabel::OBMol*);
    bool Structurize();
    bool ObtainClassMaps();
    void ObtainClassMapForMol(OpenBabel::OBMol*, map<int,int>*);

  public:
    Reactor(string, string, int);
   ~Reactor();
    bool Init();
};

#endif
