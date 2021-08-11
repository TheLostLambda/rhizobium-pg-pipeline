import os.path as path
import sys
import csv
import re


class Node:
    # Some helpful mappings
    chain_abbr = {"G": "GlcNAc", "M": "MurNAc_alditol"}
    stem_abbr = {"J": "mDAP"}
    # A shared, class-wide counter
    id_counter = 0

    @classmethod
    def chain(cls, residue):
        "Create a new glycan chain residue"
        return cls.__new__(residue, Node.chain_abbr)

    @classmethod
    def stem(cls, residue):
        "Create a new peptide stem residue"
        return cls.__new__(residue, Node.stem_abbr)

    @classmethod
    def __new__(cls, residue, abbr):
        # Create an instance of `Node` using the `__new__` of its parent
        self = super().__new__(cls)
        # Resolve residue abbreviations using the specified mapping
        self.residue = abbr.get(residue) or residue
        # Generate a unique ID from the residue abbreviation and a counter
        self.id = residue + str(Node.id_counter)
        # By default there is no modification (encoded as "zero")
        self.mod = "zero"
        # Increment the counter so that no two IDs are the same
        Node.id_counter += 1
        return self

    def __repr__(self):
        return self.residue


class Edge:
    def __init__(self, id_a, id_b, type=2):
        "Create an Edge from two Node IDs and a linktype"
        self.id_a = id_a
        self.id_b = id_b
        self.type = type

    def __repr__(self):
        return f"{self.id_a}->{self.id_b}"


class Monomer:
    """
    Rules:
      - A monomer always contains both a glycan chain and peptide stem,
        separated by a hyphen (`-`)
      - The stem may contain at most one lateral chain, offset by brackets `()`
      - The lateral chain is bound to the side-chain of the preceding residue
    """

    lat_re = re.compile(r"\((\w*)\)")

    def __init__(self, ms1_str):
        # Save the MS1 string for later
        self.structure = ms1_str
        # Split the MS1 string into the glycan chain and peptide stem
        chain, stem = ms1_str.split("-")
        # Create chain nodes
        self.chain = [Node.chain(r) for r in chain]
        # Split the stem on the lateral chain (if present)
        stem = Monomer.lat_re.split(stem)
        # Was a lateral chain present?
        if len(stem) > 1:
            # If so, extract, resolve, and remove it
            self.lat = [Node.stem(r) for r in stem.pop(1)]
        else:
            # Otherwise, initialise the attribute with `[]`
            self.lat = []
        # Create stem nodes from the remaining, non-empty parts
        self.stem = [[Node.stem(r) for r in p] for p in stem if p]
        # The terminal stem residue should have a "Hydroxy" modification
        self.stem[-1][-1].mod = "Hydroxy"
        # The terminal lateral residue should have a "Hydrogen" modification
        if self.lat:
            self.lat[-1].mod = "Hydrogen"

    def __repr__(self):
        return self.structure

    def nodes(self):
        # Flatten the `self.stem` nested list
        flat_stem = sum(self.stem, [])
        # Return the chain, stem, and lateral chain nodes + structure name
        return {self.structure: [*self.chain, *flat_stem, *self.lat]}

    def edges(self):
        # Generates unidirectional bonds for a linear portion of the molecule
        def bond_segment(nodes):
            # Stop looping one node before the end of the segment
            span = range(0, len(nodes) - 1)
            # A 2-node sliding window is used to build edges
            return [Edge(nodes[i].id, nodes[i + 1].id) for i in span]

        # Flatten the `self.stem` nested list
        flat_stem = sum(self.stem, [])
        # Generate edges for the linear chain->stem segment
        edges = bond_segment(self.chain + flat_stem)
        # Check if a lateral chain is present
        if self.lat:
            # If so, generate a segment using the lateral chain and the stem
            # node that it's connected to
            lat = [self.stem[0][-1], *self.lat]
            lat.reverse()
            # Generate edges for that segment, but with bonds running in the
            # opposite direction
            edges += bond_segment(lat)
        # Return the completed edge-list (with the structure name)
        return {self.structure: edges}


class Dimer:
    """
    Rules:
      - Dimers consist of monomers separated by either `~` or `=`
      - Dimers with `~` are connected via a glycosidic bond
      - Dimers with `=` are connected via a peptide bond (3-3 or 3-4)
    """

    sep_re = re.compile(r"([~=])")

    def __init__(self, ms1_str):
        # Save the MS1 string for later
        self.structure = ms1_str
        # Split the MS1 monomers
        monomers = Dimer.sep_re.split(ms1_str)
        # Remove and record the dimer separator
        self.glycosidic = monomers.pop(1) == "~"
        # If the monomers joined by a glycosidic bond
        if self.glycosidic:
            # Then convert the rest of the string to monomers
            self.monomers = [Monomer(m) for m in monomers]
            # And oxidise / re-bond the MurNAc of the first monomer
            self.monomers[0].chain[-1].residue = "MurNAc"
            self.monomers[0].chain[-1].mod = "negH"
            self.monomers[1].chain[0].mod = "negHOxy"
        else:
            # Otherwise, split into 3-3 and 3-4 bonded monomers
            self.monomers = {
                3: [Monomer(m) for m in monomers],
                4: [Monomer(m) for m in monomers],
            }
            # Both result in a negH mod at position 3 of the first monomer
            self.monomers[3][0].stem[0][2].mod = "negH"
            self.monomers[4][0].stem[0][2].mod = "negH"
            # 3-3 bonds modify position 3 of the second monomer as well
            self.monomers[3][1].stem[0][2].mod = "negHOxy"
            # 3-4 bonds remove a mod from the second monomer stem terminus
            self.monomers[4][1].stem[-1][-1].mod = "zero"

    def __repr__(self):
        return self.structure

    def nodes(self):
        # If the bond is glycosidic
        if self.glycosidic:
            # Then just splice both monomers' node lists together
            return {self.structure: Dimer.splice(self.monomers, "nodes")}
        else:
            # Otherwise, return two different node lists for 3-3 and 3-4 bonds
            key_leader = self.structure + " (3-"
            return {
                key_leader + "3)": Dimer.splice(self.monomers[3], "nodes"),
                key_leader + "4)": Dimer.splice(self.monomers[4], "nodes"),
            }

    def edges(self):
        # If the bond is glycosidic
        if self.glycosidic:
            # Then splice both monomers' edge lists together
            edges = Dimer.splice(self.monomers, "edges")
            # And add the glycosidic bond
            mono_a, mono_b = self.monomers
            edges.append(Edge(mono_a.chain[-1].id, mono_b.chain[0].id, 2))
            # Finally, return the full, named edge-list
            return {self.structure: edges}
        else:
            # Otherwise, create two different edge lists for 3-3 and 3-4 bonds
            edges_3 = Dimer.splice(self.monomers[3], "edges")
            edges_4 = Dimer.splice(self.monomers[4], "edges")
            # Add a 3-3 peptide bond
            mono_3a, mono_3b = self.monomers[3]
            bond_33 = Edge(mono_3b.stem[0][2].id, mono_3a.stem[0][2].id, 2)
            edges_3.append(bond_33)
            # Add a 3-4 peptide bond
            mono_4a, mono_4b = self.monomers[4]
            bond_34 = Edge(mono_4b.stem[0][2].id, mono_4a.stem[-1][-1].id, 2)
            edges_4.append(bond_34)
            # Return the named edge lists
            return {
                self.structure + " (3-3)": edges_3,
                # FIXME: 3-4 bonds should only be generated for dimers where
                # one of the constituent monomers actually has 4 stem residues
                self.structure + " (3-4)": edges_4,
            }

    @staticmethod
    def splice(monomers, getter):
        # Extract the actual node list from the monomers
        node_lists = [getattr(m, getter)()[m.structure] for m in monomers]
        # Return the flattened result
        return sum(node_lists, [])


def write_graphs(molecule, out_dir):
    for name, nodes in molecule.nodes().items():
        with open(path.join(out_dir, f"{name} NL.csv"), "w") as nlf:
            nl = csv.writer(nlf)
            nl.writerow(["node", "molecule_ID", "colour", "mods"])
            nl.writerows([[n.id, n.residue, "cyan", n.mod] for n in nodes])
    for name, edges in molecule.edges().items():
        with open(path.join(out_dir, f"{name} EL.csv"), "w") as elf:
            el = csv.writer(elf)
            el.writerow(["node1", "node2", "chain_colour_ID", "linktype"])
            el.writerows([[e.id_a, e.id_b, "violet", e.type] for e in edges])


def load_structures(ms1_file):
    # Load the MS1 structures from CSV
    with open(ms1_file) as f:
        rows = csv.reader(f)
        ms1 = [r[0].split("|")[0] for r in rows][1:]
    # Return generated graphs for every MS1 structure
    return [Dimer(s) if Dimer.sep_re.search(s) else Monomer(s) for s in ms1]


if __name__ == "__main__":
    # Extract the MS1 input filepath and graph output path from the arguments
    ms1_file = sys.argv[1]
    out_dir = sys.argv[2]
    # Load the MS1 structures from CSV
    mols = load_structures(ms1_file)
    # Write graphs for every MS1 structure
    for mol in mols:
        write_graphs(mol, out_dir)
