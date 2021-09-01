import os
from pglib.ms2_tool import Charge_Mono_Caller as CMC
from pglib.ms2_tool import Byspec_Reader as BR
import networkx as nx
import decimal as dec
from pglib.ms2_tool import Helper_Funcs as hfunc
from pglib.ms2_tool import Bio_Graph as BG
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

# Change to relevant graph files and data
data_file = Path("Data/Inputs/20210618_RhiLeg_ndslt_TY_1.raw.byspec2")
graph_folder = Path("Data/Outputs/Graphs/")
output_dir = Path("Data/Outputs/MS2/")

# Do not change
mass_table = "Data/Constants/masses_table.csv"
mod_table = "Data/Constants/mods_table.csv"
selected_ions = ['y', 'b', 'i']


def calculate_ppm_tolerance(mass, ppm_tol):
    return (mass*ppm_tol) / 1000000


def ppm_error(obs_mass, theo_mass):
    return (1-(dec.Decimal(obs_mass)/theo_mass))*1000000


def autosearch(graph_name, user_set_charge=2, intact_ppm_tol='10', frag_ppm='20'):
    nl = Path(graph_folder) / (graph_name + " NL.csv")
    el = Path(graph_folder) / (graph_name + " EL.csv")
    hf = hfunc.Helper_Funcs(nl, el)
    nodes_from_file = hf.nodes_df()
    edges_from_file = hf.edges_df()
    mass_dict = hf.generate_dict(dict_table_filepath=mass_table)
    mod_dict = hf.generate_dict(dict_table_filepath=mod_table)

    bio_graph = BG.Bio_Graph(
        nodes_from_file, edges_from_file, mass_dict, mod_dict)

    molecules = {}
    molecule_IDs = []
    frag_structure = []
    matched_output = []

    NN_SIMPLE_CHARGE_WEIGHTS = "Data/Constants/Models/simple_weights.nn"
    NN_MONO_WEIGHTS = "Data/Constants/Models/"
    charge_mono_caller = CMC.Charge_Mono_Caller(
        NN_SIMPLE_CHARGE_WEIGHTS, NN_MONO_WEIGHTS)
    byspec_reader = BR.Byspec_Reader(data_file)
    scan_mz_charges = byspec_reader.get_scan_mz_charge()
    i_ppm = dec.Decimal(intact_ppm_tol)
    f_ppm = dec.Decimal(frag_ppm)

    master_graph = bio_graph.construct_graph()

    for components in nx.connected_components(master_graph):
        molecule = nx.subgraph(master_graph, components)
        molecule_hash = bio_graph.graph_hash(molecule)
        molecule_IDs.append(molecule_hash)
        molecules.update({molecule_hash: molecule})

    molecule_momo_mass = bio_graph.monoisotopic_mass_calculator(
        molecules, molecule_IDs)

    for (mass, graph_ID) in molecule_momo_mass:
        print(f"{graph_name}: {mass}")
        scans_to_search = []
        upper_mass_lim = mass + calculate_ppm_tolerance(mass, i_ppm)
        lower_mass_lim = mass - calculate_ppm_tolerance(mass, i_ppm)

        for scan_mz_charge_tuple in scan_mz_charges:
            scan = byspec_reader.get_scan_by_scan_number(
                scan_mz_charge_tuple[0])
            try:
                caller_result = charge_mono_caller.process(
                    scan, scan_mz_charge_tuple[1])
                if caller_result['monoisotopic_mass'] > lower_mass_lim:
                    if caller_result['monoisotopic_mass'] < upper_mass_lim:
                        print('Valid scan added')
                        scans_to_search.append(scan_mz_charge_tuple[3])

            except IndexError:
                print('-' * 20)
                print('parent ion not found in scan')
                print('scan: ',  scan_mz_charge_tuple[0])
                print('mz: ', scan_mz_charge_tuple[1])
                print('charge: ', scan_mz_charge_tuple[2])
                print('#' * 20)

        if not scans_to_search:
            print('scan_to_search is empty')

        output_path = output_dir / graph_name
        os.mkdir(output_path)

        graph = nx.Graph(molecules[graph_ID])
        fragments = bio_graph.fragmentation(graph)
        total_theo_frags = len(fragments)
        n_frag, c_frag, i_frag = bio_graph.sort_fragments(fragments)
        nlist = bio_graph.monoisotopic_mass_calculator(fragments, n_frag)
        clist = bio_graph.monoisotopic_mass_calculator(fragments, c_frag)
        ilist = bio_graph.monoisotopic_mass_calculator(fragments, i_frag)
        frag_ions_df = bio_graph.generate_mass_to_charge_masses(
            fragments, nlist, clist, ilist, selected_ions, user_set_charge)

        all_obs_frags = {tuple(sorted(f.nodes)): 0 for f in fragments.values()}

        for scan_number in scans_to_search:
            scan = byspec_reader.get_scan_by_scan_number(scan_number)
            frags_in_scan = set()
            for obs_mz, count in scan:
                for row in frag_ions_df.values:
                    ions = row[0]
                    ion_type = row[1]
                    graph_key = row[2]
                    for idx, ion in enumerate(ions):
                        upper_frag_lim = ion + \
                            calculate_ppm_tolerance(ion, f_ppm)
                        lower_frag_lim = ion - \
                            calculate_ppm_tolerance(ion, f_ppm)
                        if lower_frag_lim < obs_mz < upper_frag_lim:
                            ppm_diff = ppm_error(obs_mz, ion)
                            charge = idx + 1
                            fragment = nx.Graph(fragments[graph_key])
                            frags_in_scan.add(tuple(sorted(fragment.nodes)))
                            fragment_nodes = str(fragment.nodes)
                            frag_structure.append(fragment_nodes)
                            matched_output.append(
                                (obs_mz, ion, charge, count, ppm_diff, ion_type, frag_structure))
                            frag_structure = []

            for frag in frags_in_scan:
                all_obs_frags[frag] += 1

            matched_output_df = pd.DataFrame(matched_output, columns=[
                                             'Observered Ion', 'Theoretical Ion', 'Charge', 'count', 'PPM Error', 'Ion Type', 'Structure'])
            matched_num = len(matched_output_df)
            score = round(((matched_num/total_theo_frags)*100))
            matched_output_df.to_csv(
                output_path / f'{scan_number} ({score}%).csv')
            matched_output.clear()
        df = pd.DataFrame(all_obs_frags.items(), columns=['Fragment', 'Count'])
        df.to_csv(
            output_path / f'Observed Fragments ({len([c for c in all_obs_frags.values() if c > 0])} of {total_theo_frags}).csv', index=False)


if __name__ == "__main__":
    with Pool() as p:
        structures = [nf.name.removesuffix(" NL.csv")
                      for nf in graph_folder.glob("* NL.csv")]
        p.map(autosearch, structures)
