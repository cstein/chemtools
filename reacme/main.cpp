#include "reactor.hpp"
using namespace std;

int main(int argc, char* argv[])
{
  char* program_name = argv[0];
  if(argc != 3)
  {
    cout << "usage: " << program_name << " <REACTANT> <PRODUCT>" << endl;
    return 0;
  }
  Reactor r(argv[1], argv[2], 5);
  if(r.Init())
  {
    cout << "initialized successfully." << endl;
  }
  return 0;
}
