from pysmiles import read_smiles, write_smiles
import networkx as nx
from collections import namedtuple


# Written by D S H, but all the heavy lifting is done by external libraries


def networkx_graph_from_smiles(smiles_string, should_add_all_h_around_structure):
    return read_smiles(smiles_string, explicit_hydrogen=should_add_all_h_around_structure)


def get_mol(smile_string):
    manually_specified_h = False
    networkx_mol = networkx_graph_from_smiles(smile_string, not manually_specified_h)
    return networkx_mol


def substitute(main_structure, symbol, replacement):
    Replace_Site = namedtuple('Replace_Site', ['atom', 'bonder'])
    U = main_structure.copy()
    replace_sites = list()
    nodes_to_delete = list()

    for idx in U.nodes:
        atom_dict = U.nodes[idx]
        if atom_dict['element'] == symbol:
            bonders = list(U.adj[idx])
            assert len(bonders) == 1
            bonder = bonders[0]
            replace_sites.append(Replace_Site(idx, bonder))

    for replace_site in replace_sites:
        U = nx.disjoint_union(U, replacement)
        # A stands for Asterisk
        A_atom = None
        A_bonder = None
        for idx in U.nodes:
            atom_dict = U.nodes[idx]
            if 'element' not in atom_dict:
                A_atom = idx
                bonders = list(U.adj[idx])
                assert len(bonders) == 1
                A_bonder = bonders[0]
                break
        assert A_atom is not None
        U.add_edge(replace_site.bonder, A_bonder)
        U.remove_edge(A_atom, A_bonder)
        U.remove_edge(replace_site.atom, replace_site.bonder)
        nodes_to_delete.append(A_atom)
        U.nodes[A_atom]['element'] = '*'
        nodes_to_delete.append(replace_site.atom)

    for node in nodes_to_delete:
        U.remove_node(node)

    return U


def make_smiles(networkx_mol):
    sm_str = write_smiles(networkx_mol)
    return sm_str


if __name__ == '__main__':
    output = open("output.csv", 'w')
    with open("structures.txt") as structureFile:
        structures_text = structureFile.readlines()
    with open("subgroups.txt") as subgroupFile:
        subgroups_text = subgroupFile.readlines()
    subgroups = [get_mol(s) for s in subgroups_text]
    num_structure = 0
    num_subgroup = 0
    for smiles in structures_text:
        num_structure += 1
        structure = get_mol(smiles)
        if 'R' in smiles:
            num_subgroup = 0
            for subgroup in subgroups:
                num_subgroup += 1
                result_mol = substitute(structure, 'R', subgroup)
                result_str = make_smiles(result_mol)

                check_mol = get_mol(result_str)
                assert check_mol.order() == result_mol.order()
                assert check_mol.size() == result_mol.size()
                output.write(f"base structure {num_structure}, sub group {num_subgroup}, {result_str}\n")
        else:
            output.write(f"base structure {num_structure}, no replacement, {smiles.strip()}\n")
    print(f"Finished {num_structure} main structures and {num_subgroup} sub groups")
