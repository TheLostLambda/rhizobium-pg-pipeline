import networkx as nx
from networkx.classes.function import set_node_attributes
import pandas as pd
import numpy as np
from datetime import datetime
import decimal as dec
import hashlib as hl
import itertools as itool
import pylab as plt
import copy as cp


class Bio_Graph:

    sugars = ['GlcNAc', 'GlcN', 'MurNAc', 'MurN', 'MurNAc_alditol']

    def __init__(self, nodes_from_file: pd.DataFrame, edges_from_file: pd.DataFrame, mass_dict: dict, mods_dict: dict):
        self.nodes_from_file = nodes_from_file
        self.edges_from_file = edges_from_file
        self.mass_dict = mass_dict
        self.mods_dict = mods_dict

    # Generate edge order dict
    def construct_edge_order(self):

        # Read edge order defined in edge dataframe (read from csv file)
        edge_order_df = self.edges_from_file[['node1', 'node2']]
        # dictionary for containing edge order information in the form of a tuple (u,v)
        edge_order_dict = {}
        # Populate edge_order_dict with edge order information (key:edge_order) where edge order is a tuple of the form (node1,node2) as (u,v)
        for row in edge_order_df.itertuples(index=False):
            u = getattr(row, 'node1')
            v = getattr(row, 'node2')
            # Map both underorder values of node tuples (u,v and v,u) to ordered tuple (u,v)
            edge_order_dict[u, v] = (u, v)
            edge_order_dict[v, u] = (u, v)

        self.edge_order = edge_order_dict

        return edge_order_dict

    # construct graph from pandas dataframes
    def construct_graph(self):

        # create master graph from edge list
        master_graph = nx.from_pandas_edgelist(
            self.edges_from_file, 'node1', 'node2', edge_attr=True)

        # add node attributes from node list
        node_attributes = self.nodes_from_file.set_index(
            'node').to_dict('index').items()
        master_graph.add_nodes_from(node_attributes)
        nodeskeys = dict.fromkeys(master_graph, [])
        nx.set_node_attributes(master_graph, nodeskeys, "flag")
        for id in master_graph.nodes:
            mod = self.mods_dict[master_graph.nodes[id]["mods"]]
            master_graph.nodes[id]["mods"] = mod
        # generate edge order list (used for fragmentation methods)
        self.construct_edge_order()

        return master_graph

    # Get ordered edge tuple
    def get_ordered_edge(self, node1, node2):
        key_tuple = (node1, node2)
        ordered_edge = self.edge_order[key_tuple]
        return ordered_edge

    # draw graph in window (Debug function)
    def draw_graph(self, graph):

        # Get colours values for nodes from node attribute : 'colour'
        ncolours = [n[1]['colour'] for n in graph.nodes(data=True)]

        # Get Colour values for edges from edge attribute : 'chain_colour_ID'
        ecolours = [e[2]['chain_colour_ID'] for e in graph.edges(data=True)]

        plt.figure('1')
        nx.draw_kamada_kawai(graph, with_labels=True, font_weight='bold', node_color=ncolours,
                             edge_color=ecolours,
                             node_size=500)
        plt.show()

    # Generate unique ID values for graphs
    def graph_hash(self, graph):
        # List all node IDs as string values
        node_ids = str(sorted(graph.nodes))
        # Set hash algorithm
        h = hl.blake2b()
        # Generate hash value for graph (utf8 encoding)
        h.update(node_ids.encode('utf8'))
        # Digest hash to create unique 6 letter ID
        graph_hash_id = h.hexdigest()

        return graph_hash_id

    # calculate monoisotopic mass of graph (Debug function)
    def calc_check(self, graph):

        node_residue_mass = []
        mod_mass = []
        # get node ID values
        node_ids = [n[1]['molecule_ID'] for n in graph.nodes(data=True)]
        # get mod ID values
        mod_ids = [m[1]['mods'] for m in graph.nodes(data=True)]
        print(mod_ids)
        # List masses for each node ID
        for ID in node_ids:
            node_residue_mass.append(self.mass_dict[ID])
        # List masses for each mod ID
        for ID in mod_ids:
            mod_mass.append(self.mods_dict[ID])
        # Sum node masses
        nodes_total = sum(np.array(node_residue_mass))
        # Sum mod masses
        mods_total = sum(np.array(mod_mass))
        # Sum mod masses & node masses to calculate neutral monoisotopic mass of molecule
        output_mass = nodes_total + mods_total

        return output_mass

    # Calculate monoisoptopic mass of graphs and return unique masses along with all corresponding graphs
    def monoisotopic_mass_calculator(self, graph_fragments, graph_IDs):

        # List of node mass values
        node_residue_mass = []

        # List of graph masses (will contain duplicates)
        graph_masses = []

        # Generate residue mass of each graph
        for f_graph in graph_IDs:
            # Create graph intance of f_graph
            g = nx.Graph(graph_fragments[f_graph])
            # List of node IDs in f_graph
            node_molecule_ID = [n[1]['molecule_ID']
                                for n in g.nodes(data=True)]
            # List of node mods in f_graph
            node_mods = [n[1]['mods'] for n in g.nodes(data=True)]
            # Generate list of node masses
            for ID in node_molecule_ID:
                node_residue_mass.append(self.mass_dict[ID])

            nodes_total = sum(np.array(node_residue_mass))
            mods_total = sum(np.array(node_mods))
            # Sum node and mod mass totals to generate neutral monoisotopic mass of f_graph
            total_mono_mass = nodes_total + mods_total

            graph_masses.append((total_mono_mass, f_graph))

            node_residue_mass.clear()

        return graph_masses

    # Cuts bonds in given graph
    def bond_cutter(self, intact_graph):
        # dict for holding generated fragments
        fragment_graphs = {}
        # Gets edges from intact graph (intact graph represents parent graph)
        for (u, v) in intact_graph.edges:
            # Generate copy of graph to apply cut operation to (prevents intact graph object used in line above from changing during execution)
            intact_copy = cp.deepcopy(intact_graph)

            # N1 is always N-terminal node
            (n1, n2) = self.get_ordered_edge(u, v)

            intact_copy = self.set_fragmentation_flags(intact_copy, n1, n2)

            # Remove edge (break bond)
            intact_copy.remove_edge(n1, n2)

            node1 = intact_copy.nodes[n1]
            node2 = intact_copy.nodes[n2]

            node1["mods"] += self.mods_dict["negProton"]
            node2["mods"] += self.mods_dict["Hydrogen"]

            # FIXME: Check the residues being separated and update any modifications. e.g
            # add the Hydrogen to a MurNAc terminal or add a Hydrogen to B ions, N-terminal
            # residues)
            # Generates subgraphs based on the intact_copy graph
            fragments = [cp.deepcopy(intact_copy.subgraph(comps))
                         for comps in nx.connected_components(intact_copy)]
            for fragment in fragments:
                graph_ID = self.graph_hash(fragment)
                # Add fragment to fragment list if the hash value is not already present in list.
                if graph_ID not in fragment_graphs:
                    fragment_graphs.update({graph_ID: fragment})

        return fragment_graphs

    # Set fragmentation flags (helps determine which ion type to calculate later)
    def set_fragmentation_flags(self, graph, node1, node2):

        if graph.edges[node1, node2]['linktype'] == 2:
            if graph.nodes[node1]['flag'] == []:
                nx.set_node_attributes(graph, {node1: ['N']}, 'flag')
            else:
                flags_list = graph.nodes[node1]['flag']
                flags_list.append('N')
                nx.set_node_attributes(graph, {node1: flags_list}, 'flag')

            if graph.nodes[node2]['flag'] == []:
                nx.set_node_attributes(graph, {node2: ['C']}, 'flag')
            else:
                flags_list = graph.nodes[node2]['flag']
                flags_list.append('C')
                nx.set_node_attributes(graph, {node1: flags_list}, 'flag')

        elif graph.edges[node1, node2]['linktype'] == 1:
            if graph.nodes[node1]['flag'] == []:
                nx.set_node_attributes(graph, {node1: ['C']}, 'flag')
            else:
                flags_list = graph.nodes[node1]['flag']
                flags_list.append('C')
                nx.set_node_attributes(graph, {node1: flags_list}, 'flag')
            if graph.nodes[node2]['flag'] == []:
                nx.set_node_attributes(graph, {node2: ['N']}, 'flag')
            else:
                flags_list = graph.nodes[node1]['flag']
                flags_list.append('N')
                nx.set_node_attributes(graph, {node1: flags_list}, 'flag')

        return graph

    # Generates multi cut fragments from parent graph
    def fragmentation(self, molecule):
        # Define a local helper function for recursion
        def recurse(molecule, seen=set()):
            # Add the current graph as a fragment to the dictionary
            all_fragments = {self.graph_hash(molecule): molecule}
            # If the molecule contains 1 or fewer nodes, return the current graph
            if len(molecule.nodes) <= 1:
                return all_fragments
            # Otherwise, recursively fragment the molecule
            else:
                # Generate a dictionary of single-cut products
                cut_products = {id: graph for id, graph in self.bond_cutter(
                    molecule).items() if id not in seen}
                seen |= cut_products.keys()
                # Loop through the cut products, fragmenting them and adding the
                # generated products to `all_fragments`
                for graph in cut_products.values():
                    all_fragments = recurse(
                        nx.Graph(graph), seen) | all_fragments
                return all_fragments

        # Recursively compute all fragments
        result = recurse(molecule)
        # Remove the unfragmented parent molecule from the results
        del result[self.graph_hash(molecule)]
        return result

    # Sorts fragments into N/C terminal sides or internal fragment
    def sort_fragments(self, fragment_graphs):

        # flag total counter
        ftotal = 0
        # list for sorted graphs
        nfrag = []
        cfrag = []
        ifrag = []
        # Counts total flags in graph
        for graph_ID in fragment_graphs:
            graph = nx.Graph(fragment_graphs[graph_ID])
            flagged = graph.nodes(data='flag')
            for n in flagged:
                if n[1]:
                    ftotal = ftotal + len(n[1])
        # Sorts graphs with only 1 flags into nfrag or cfrag lists all other fragments are stored in ifrag
            if ftotal <= 1:
                for n in flagged:
                    if n[1] == ['N']:
                        nfrag.append(graph_ID)
                    elif n[1] == ['C']:
                        cfrag.append(graph_ID)
            else:
                ifrag.append(graph_ID)
            ftotal = 0

        return nfrag, cfrag, ifrag

    # Generate m/z masses
    def generate_mass_to_charge_masses(self, fgraphs, n_terminal_fragments: list, c_terminal_fragments: list,
                                       internal_fragments: list, enabled_ions: list, charge_limit: int = 1):
        proton_mass = dec.Decimal('1.0073')
        mzlist_graph = []
        b_mz = []
        y_mz = []
        i_mz = []

        proton_mass = dec.Decimal('1.0073')

        # FIXME: This is obviously a bit repetitive, but I'm scared to refactor
        # critical code when I'm half-asleep. In reality, I should write tests.
        if 'y' in enabled_ions:
            for (mono_mass, graph) in c_terminal_fragments:
                for charge_state in range(1, charge_limit + 1):
                    y_mz.append(
                        (mono_mass + (proton_mass * charge_state)) / charge_state)
                mzlist_graph.append((y_mz, "y_ion", graph))
                y_mz = []

        if 'b' in enabled_ions:
            for (mono_mass, graph) in n_terminal_fragments:
                for charge_state in range(1, charge_limit + 1):
                    b_mz.append(
                        (mono_mass + (proton_mass * charge_state)) / charge_state)
                mzlist_graph.append((b_mz, "b_ion", graph))
                b_mz = []

        if 'i' in enabled_ions:
            for (mono_mass, graph) in internal_fragments:
                for charge_state in range(1, charge_limit + 1):
                    i_mz.append(
                        (mono_mass + (proton_mass * charge_state)) / charge_state)
                mzlist_graph.append((i_mz, "i_ion", graph))
                i_mz = []

        mzlist_graph_df = pd.DataFrame(
            mzlist_graph, columns=['mz', 'Ion_type', 'Graph_ID'])
        return mzlist_graph_df

    def group_by_mass(self, graph_masses: list):
        # Index position of matched graphs
        matched_graph_indexes = []

        # Container for unique values already added to graph masses (no duplicate values)
        seen = set()

        # Index position of graph(s) matched to unique mass value
        graph_index_group = []

        # Generate list of unque graph mass values from graph masses list
        for idx, graph_mass in enumerate(graph_masses):
            if graph_mass[1] not in seen:
                seen.add(graph_mass[1])

        # Group graphs according to their masses (graphs with the same mass are grouped together)
        for s_idx, unique_mass in enumerate(seen):
            for g_index, graph_mass in enumerate(graph_masses):
                if graph_mass[1] == unique_mass:
                    graph_index_group.append(graph_mass[0])
            # list of graphs grouped by mass - each item is a tuple containing (mass, list of graphs)
            matched_graph_indexes.append((unique_mass, graph_index_group))
            graph_index_group = []
