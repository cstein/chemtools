# Chemtools for OpenBabel
The Chemtools for chemistry based on mainly the [Open Babel](http://openbabel.org/wiki/Main_Page) [API](http://openbabel.org/api/)

## Tools

* opthydrogens optimizes hydrogens of a molecule structure leaving the rest frozen. Uses the MMFF94 force field to minimize.

* reacme takes reactant and product "SMILES":http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html strings and constructs. Please note that this is still a work in progress.

##t Setup

each subdirectory has a simple Makefile. You must change the OB_BASE variable to point to your Open Babel installation directory followed by

  @make
