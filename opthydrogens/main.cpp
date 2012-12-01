#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/forcefield.h>


using namespace std;
using namespace OpenBabel;

int main(int argc, char **argv) {
  // some default parameters and properties
  char *program_name = argv[0];
  int n_steps = 2500;
  double convergence=1e-5;
  string option;
  string basename, filename_in, extension = "";

  if (argc < 2) {
    cout << "usage: " << program_name << " [options] filename" << endl;
    cout << "    [options]:" << endl;
    cout << "        -n    " << " number of steps in minimization (default=" << n_steps << ")" << endl;
    cout << endl;
    cout << "        -c    " << " convergence criteria (default=" << convergence << ")" << endl;
    cout << endl;
    return 0;
  }

  // parse parameters, heavily inspired by obminimize
  int ifile = 1;
  for (int i = 0; i < argc; i++) {
    option = argv[i];
    if ((option == "-n") && (argc > (i+1))) {
      n_steps = atoi(argv[i+1]);
      ifile += 2;
    }
    if ((option == "-c") && (argc > (i+1))) {
      convergence = atof(argv[i+1]);
      ifile += 2;
    }
  }
  // store the filename
  basename = filename_in = argv[ifile];
  size_t extPos = filename_in.rfind('.');
  if( extPos != string::npos ) {
    basename = filename_in.substr(0, extPos);
    extension = filename_in.substr(extPos+1, string::npos);
  }

  // get the formats, always output the same format as we get in
  OBConversion conv;
  OBFormat *format = conv.FindFormat(extension);

  if( !format ) {
    cerr << program_name << ": cannot read input/output format '" << extension << "'!" << endl;
    exit(-1);
  }

  // attempt to open input file to see if exists
  ifstream ifs;

  ifs.open(filename_in.c_str());
  if(!ifs) {
    cerr << program_name << ": cannot read input file '" << filename_in.c_str() << "'!" << endl;
  }
  

  if(conv.SetInAndOutFormats(format, format)) {
    OBMol mol;
    if(conv.Read(&mol, &ifs)) {
      // setup force field
      OBForceField* pFF = OBForceField::FindForceField("MMFF94");
      if( pFF ) {
        // set some properties for the force field
        pFF->SetLogFile(&cerr);
        pFF->SetLogLevel(OBFF_LOGLVL_LOW);

        // fix all atoms which are not hydrogens
        OBFFConstraints constraints;
        FOR_ATOMS_OF_MOL(current_atom, mol) {
          int atom_idx = current_atom->GetIdx();
          if( !current_atom->IsHydrogen() )
            constraints.AddAtomConstraint( atom_idx );
        }

        // carry out the minimization if the constraints were successfull
        if( pFF->Setup( mol, constraints) ) {
          pFF->ConjugateGradients(n_steps, convergence);
          pFF->GetCoordinates( mol );
        }
        conv.Write(&mol, &cout );
      } else
        cout << program_name << ": could not locate force field!";
    } else
      cout << program_name << ": error while reading file";
  }
  return 0;
}
