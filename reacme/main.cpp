#include <string>
#include "reactor.hpp"
using namespace std;

struct reacme_settings
{
  int n_steps;
  int n_conformers;
};

void defaults(reacme_settings* settings)
{
  settings->n_steps = 250;
  settings->n_conformers = 5;
}

int main(int argc, char* argv[])
{
  reacme_settings settings;
  defaults(&settings);

  string program_name = argv[0];
  if (argc < 2)
  {
    cout << "usage: " << program_name << " [options] <REACTANT> <PRODUCT>" << endl;
    cout << "    [options]:" << endl;
    cout << "        -n    " << " number of steps in minimization (default=" << settings.n_steps << ")" << endl;
    cout << endl;
    cout << "        -c    " << " conformer count (default=" << settings.n_conformers << ")" << endl;
    cout << endl;
    return 0;
  }

  // parse parameters, heavily inspired by obminimize
  int ireactant = 1;
  string option;
  for (int i = 0; i < argc-1; i++) {
    option = argv[i];
    if ((option == "-n") && (argc > (i+1))) {
      settings.n_steps = atoi(argv[i+1]);
      ireactant += 2;
    }
    if ((option == "-c") && (argc > (i+1))) {
      settings.n_conformers = atoi(argv[i+1]);
      ireactant += 2;
    }
  }
  // the last two occurences of input arguments is the
  // reactant and the product
  if(ireactant +2 != argc)
  {
    cout << "reactants or products not specified. aborting." << endl;
    return 2;
  }

  string s_reactant = argv[ireactant];
  string s_product  = argv[ireactant+1];

  Reactor r(s_reactant, s_product, settings.n_conformers);
  if(r.Init())
  {
    cout << "initialized successfully." << endl;
  }
  return 0;
}
