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

    # TODO: Is this really the most elegant solution?
    @classmethod
    def reset_ids(cls):
        "Reset residue numbering"
        cls.id_counter = 0

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

    def __init__(self, ms1_str, reset_ids=True):
        # Save the MS1 string for later
        self.structure = ms1_str
        # Split the MS1 string into the glycan chain and peptide stem
        chain, stem = ms1_str.split("-")
        # Conditionally reset node numbering
        if reset_ids:
            Node.reset_ids()
        # Create chain nodes
        self.chain = [Node.chain(r) for r in chain]
        # The first chain residue should have a "Hydrogen" modification
        self.chain[0].mod = "Hydrogen"
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
        # The terminal stem residue should have a "Hydroxyl" modification
        self.stem[-1][-1].mod = "Hydroxy"
        # The terminal lateral residue should have a "Hydrogen" modification
        # and the branch-point should have a "negH"
        if self.lat:
            self.stem[0][2].mod = "negH"
            self.lat[-1].mod = "Hydrogen"

    def __repr__(self):
        return self.structure

    def nodes(self):
        # Return the chain, stem, and lateral chain nodes + structure name
        return {self.structure: [*self.chain, *self.flat_stem(), *self.lat]}

    def edges(self):
        # Generates unidirectional bonds for a linear portion of the molecule
        def bond_segment(nodes):
            # Stop looping one node before the end of the segment
            span = range(0, len(nodes) - 1)
            # A 2-node sliding window is used to build edges
            return [Edge(nodes[i].id, nodes[i + 1].id) for i in span]

        # Flatten the `self.stem` nested list
        # Generate edges for the linear chain->stem segment
        edges = bond_segment(self.chain + self.flat_stem())
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

    def flat_stem(self):
        # Flatten the `self.stem` nested list
        return sum(self.stem, [])


class Dimer:
    """
    Rules:
      - Dimers consist of monomers separated by either `~` or `=`
      - Dimers with `~` are connected via a glycosidic bond
      - Dimers with `=` are connected via a peptide bond (3-3 or 3-4)
    """

    sep_re = re.compile(r"([~=])")
    crosslinks = ["3-3", "3-4"]

    def __init__(self, ms1_str, flipped=False):
        # Split the MS1 monomers
        monomers = Dimer.sep_re.split(ms1_str)
        # Optionally swap the acceptor and donor monomers
        if flipped:
            monomers.reverse()
        # Save the MS1 string for later
        self.structure = "".join(monomers)
        # Remove and record the dimer separator
        self.glycosidic = monomers.pop(1) == "~"
        # If the monomers joined by a glycosidic bond
        if self.glycosidic:
            # Reset node numbering
            Node.reset_ids()
            # Then convert the rest of the string to monomers
            self.monomers = [Monomer(m, False) for m in monomers]
            # And oxidise / re-bond the MurNAc of the first monomer
            self.monomers[0].chain[-1].residue = "MurNAc"
            self.monomers[0].chain[-1].mod = "negH"
            self.monomers[1].chain[0].mod = "negHOxy"
        else:
            # Otherwise, split into crosslinked sets of monomers
            self.monomers = {
                type: monos
                for type in Dimer.crosslinks
                if (monos := self.crosslink(type, monomers))
            }

    def __repr__(self):
        return self.structure

    def crosslink(self, type, monomers):
        # Since the name `3-3` can only exist as a string, we need to check
        # that the provided string `type` is a valid mode of crosslinking
        assert type in Dimer.crosslinks
        # Reset node numbering
        Node.reset_ids()
        # Parse the monomer strings into Monomer objects
        acceptor, donor = [Monomer(m, False) for m in monomers]
        # Get monomer stem lengths
        a_len, d_len = (len(acceptor.flat_stem()), len(donor.flat_stem()))
        # Check for illegal crosslinking combinations
        if (
            a_len < 3
            or (type == "3-3" and d_len != 3)
            or (type == "3-4" and (d_len != 4 or donor.stem[-1][-1].residue != "A"))
        ):
            return None
        # In both types, a hydrogen is lost from residue 3 of the acceptor stem
        if acceptor.stem[0][2].mod == "Hydroxy":
            acceptor.stem[0][2].mod = "Oxygen"
        else:
            acceptor.stem[0][2].mod = "negH"
        # 3-3 bonds lose an OH on the donor stem
        if type == "3-3":
            if donor.stem[0][2].mod == "Hydroxy":
                donor.stem[0][2].mod = "zero"
            else:
                donor.stem[0][2].mod = "negHOxy"
        # 3-4 bonds remove a hydroxy from the donor terminus
        # TODO: Take a closer look at this!
        if type == "3-4":
            donor.stem[-1][-1].mod = "zero"
        # Return the crosslinked monomers
        return [acceptor, donor]

    def nodes(self):
        # If the bond is glycosidic
        if self.glycosidic:
            # Then just splice both monomers' node lists together
            return {self.structure: Dimer.splice(self.monomers, "nodes")}
        else:
            # Otherwise, return different node lists for the crosslinking types
            return {
                f"{self.structure} ({type})": Dimer.splice(monos, "nodes")
                for type, monos in self.monomers.items()
            }

    def edges(self):
        # If the bond is glycosidic
        if self.glycosidic:
            # Then splice both monomers' edge lists together
            edges = Dimer.splice(self.monomers, "edges")
            # And add the glycosidic bond
            acceptor, donor = self.monomers
            edges.append(Edge(acceptor.chain[-1].id, donor.chain[0].id, 2))
            # Finally, return the full, named edge-list
            return {self.structure: edges}
        else:
            # Otherwise, create two different edge lists for 3-3 and 3-4 bonds
            edges = {
                type: Dimer.splice(monos, "edges")
                for type, monos in self.monomers.items()
            }
            # Add a 3-3 peptide bond if possible
            if "3-3" in edges:
                acceptor, donor = self.monomers["3-3"]
                bond = Edge(donor.stem[0][2].id, acceptor.stem[0][2].id, 2)
                edges["3-3"].append(bond)
            # Add a 3-4 peptide bond if possible
            if "3-4" in edges:
                acceptor, donor = self.monomers["3-4"]
                bond = Edge(donor.stem[-1][-1].id, acceptor.stem[0][2].id, 2)
                edges["3-4"].append(bond)
            # Return the named edge lists
            return {
                f"{self.structure} ({type})": edges for type, edges in edges.items()
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
    return ms1


def parse_structures(ms1):
    # Return generated graphs for every MS1 structure
    structs = [
        [Dimer(s), Dimer(s, flipped=True)] if Dimer.sep_re.search(s) else [Monomer(s)]
        for s in ms1
    ]
    # Flatten the nested lists
    return sum(structs, [])


if __name__ == "__main__":
    # Extract the MS1 input filepath and graph output path from the arguments
    ms1_file = sys.argv[1]
    out_dir = sys.argv[2]
    # Load the MS1 structures from CSV
    mols = load_structures(ms1_file)
    # Write graphs for every MS1 structure
    for mol in mols:
        write_graphs(mol, out_dir)
