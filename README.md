# SmilesCombiner
Simple script to automatically combine SMILES structures (https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system). It relies on pysmiles (pip install pysmiles) and networkx (installed as part of pysmiles). Written for Python 3.9; undefined behavior on other versions but will probably work.

It looks for two files in the same directory as it, titled "structures.txt" and "subgroups.txt", which each have one SMILES per line. It performs the substitutions and outputs an "output.csv" file like the example in the repository.

The main structures can have any number of atoms whose elements are R. These must be bonded to only one other atom, and they will be replaced. If no R is found, then instead of combining with subgroups, one line will be output, specifying that there is no replacement. If there are multiple R atoms, all will be replaced by the same subgroup; there is currently no way to replace some with a subgroup and other with another.

The subgroups must have exactly one atom whose element is \*. It must be bonded to exactly one other element, and it will be replaced.

Thus, whatever was bonded to the R and whatever was bonded to the * will be bonded to each other in the output molecules, since R and * will be deleted.

In the example, 3 structures are combined with 2 substitution groups for a total of 5 molecules (because one of the structures has no R spots). The output smiles aren't pretty, but they work.
