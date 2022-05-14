import sys
import csv
from pathlib import Path
from build_graphs import load_structures, parse_structures, write_graphs
from pglib.ms2_tool.Helper_Funcs import Helper_Funcs
from pglib.ms2_tool.Bio_Graph import Bio_Graph
from tempfile import TemporaryDirectory

MASSES = "Data/Constants/masses_table.csv"
MODS = "Data/Constants/mods_table.csv"


# FIXME: This is almost certainly in the wrong place...
def pretty_print_nodes(nodes):
    return "-".join(sorted(nodes, key=lambda n: int(n[1:])))


if __name__ == "__main__":
    ms1_file = sys.argv[1]
    out_dir = sys.argv[2]
    mols = load_structures(ms1_file)
    with TemporaryDirectory() as tmp:
        for mol in parse_structures(mols):
            write_graphs(mol, tmp)
        p = Path(tmp)
        for node_file in p.glob("* NL.csv"):
            parent = node_file.parent
            structure_name = node_file.name.removesuffix(" NL.csv")
            edge_file = parent / (structure_name + " EL.csv")
            hf = Helper_Funcs(node_file, edge_file)
            masses = hf.generate_dict(MASSES)
            mods = hf.generate_dict(MODS)
            bg = Bio_Graph(hf.nodes_df(), hf.edges_df(), masses, mods)
            graph = bg.construct_graph()
            fragments = bg.fragmentation(graph)
            nfrags, cfrags, ifrags = [
                sorted(f, key=lambda i: len(fragments[i].nodes))
                for f in bg.sort_fragments(fragments)
            ]

            tagged_frags = [
                *[(i, "N-Terminal", fragments[i].nodes) for i in nfrags],
                *[(i, "C-Terminal", fragments[i].nodes) for i in cfrags],
                *[(i, "Internal", fragments[i].nodes) for i in ifrags],
            ]
            with open(Path(out_dir) / f"{structure_name} Fragments.csv", "w") as f:
                ff = csv.writer(f)
                ff.writerow(["Type", "Ion", "Mass", "Parts"])
                for id, type, frag in tagged_frags:
                    mass = bg.monoisotopic_mass_calculator(fragments, [id])[0][0]
                    ff.writerow(
                        [type, mass + mods["Proton"], mass, pretty_print_nodes(frag)]
                    )
