from pglib.ms2_tool import Charge_Mono_Caller as CMC
from pglib.ms2_tool import Byspec_Reader as BR
import networkx as nx
import decimal as dec
from pglib.ms2_tool import Helper_Funcs as hfunc
from pglib.ms2_tool import Bio_Graph as BG
import pandas as pd
from pathlib import Path
from multiprocessing import Pool
from build_graphs import load_structures, parse_structures, write_graphs
import re
from datetime import datetime

# Change to relevant graph files and data
data_file = Path("Data/Inputs/20210618_RhiLeg_ndslt_TY_1.raw.byspec2")
ms1_file = Path("Data/Outputs/MS1.csv")
selected_structures = Path("Data/Outputs/selected_structures.csv")
graph_folder = Path("Data/Outputs/Graphs/")
output_dir = Path("Data/Outputs/MS2/")

# Do not change
mass_table = "Data/Constants/masses_table.csv"
mod_table = "Data/Constants/mods_table.csv"
selected_ions = ['y', 'b', 'i']

# Break this file up into sub-files in pglib/ms2_tool


def split_structures(st):
    return re.findall(r"([A-Z-=~]+)", st)


def filter_structures(df, structures):
    return df[df["inferredStructure"].apply(
        lambda ss: not set(structures).isdisjoint(split_structures(ss)))]


def all_structures(df):
    return {struct
            for structs in df["inferredStructure"].apply(split_structures)
            for struct in structs}


def calculate_ppm_tolerance(mass, ppm_tol):
    return (mass*ppm_tol) / 1000000


def ppm_error(obs_mass, theo_mass):
    return (1-(dec.Decimal(obs_mass)/theo_mass))*1000000


def resolve_row(ms1_row, set_charge=1, intact_ppm_tol='10', frag_ppm='20'):
    start = ms1_row.scanStart
    end = ms1_row.scanEnd
    structures = [nf.name.removesuffix(" NL.csv")
                  for s in all_structures(ms1_row.to_frame().T)
                  for nf in graph_folder.glob(f"{s} *NL.csv")]
    coverage_percents = {}
    for structure in structures:
        nl = Path(graph_folder) / (structure + " NL.csv")
        el = Path(graph_folder) / (structure + " EL.csv")
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
        scan_mz_charges = byspec_reader.filter_children_with_parent_in_range(
            start, end).get_scan_mz_charge()
        i_ppm = dec.Decimal(intact_ppm_tol)
        f_ppm = dec.Decimal(frag_ppm)

        master_graph = bio_graph.construct_graph()

        for components in nx.connected_components(master_graph):
            molecule = nx.subgraph(master_graph, components)
            molecule_hash = bio_graph.graph_hash(molecule)
            molecule_IDs.append(molecule_hash)
            molecules.update({molecule_hash: molecule})

        mass, graph_ID = bio_graph.monoisotopic_mass_calculator(
            molecules, molecule_IDs)[0]

        print(f"{structure}: {mass}")
        scans_to_search = []
        min_mass = mass - calculate_ppm_tolerance(mass, i_ppm)
        max_mass = mass + calculate_ppm_tolerance(mass, i_ppm)

        for i, scan_mz_charge_tuple in enumerate(scan_mz_charges, start=1):
            print(f"Scanning {i}/{len(scan_mz_charges)}...", end="\r")
            scan = byspec_reader.get_scan_by_scan_number(
                scan_mz_charge_tuple[0])
            try:
                caller_result = charge_mono_caller.process(
                    scan, scan_mz_charge_tuple[1])
                if min_mass < caller_result['monoisotopic_mass'] < max_mass:
                    # Padding hack!
                    print('Valid scan added' + ' ' * 10)
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
            continue

        output_path = output_dir / f"{start}-{end}" / structure
        output_path.mkdir(parents=True, exist_ok=True)

        graph = nx.Graph(molecules[graph_ID])
        fragments = bio_graph.fragmentation(graph)
        total_theo_frags = len(fragments)
        n_frag, c_frag, i_frag = bio_graph.sort_fragments(fragments)
        nlist = bio_graph.monoisotopic_mass_calculator(fragments, n_frag)
        clist = bio_graph.monoisotopic_mass_calculator(fragments, c_frag)
        ilist = bio_graph.monoisotopic_mass_calculator(fragments, i_frag)
        frag_ions_df = bio_graph.generate_mass_to_charge_masses(
            fragments, nlist, clist, ilist, selected_ions, set_charge)

        all_obs_frags = [[tuple(sorted(f.nodes)), bio_graph.monoisotopic_mass_calculator(
            fragments, [id])[0][0], 0] for id, f in fragments.items()]
        all_obs_frags = pd.DataFrame(all_obs_frags, columns=[
                                     'Structure', 'Mass', 'Count'])
        all_obs_frags = all_obs_frags.set_index('Structure')

        print("Writing results..." + ' ' * 10)
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
                all_obs_frags.at[frag, 'Count'] += 1

            matched_output_df = pd.DataFrame(matched_output, columns=[
                                             'Observered Ion', 'Theoretical Ion', 'Charge', 'count', 'PPM Error', 'Ion Type', 'Structure'])
            matched_num = len(matched_output_df)
            score = round(((matched_num/total_theo_frags)*100))
            matched_output_df.to_csv(
                output_path / f'{scan_number} ({score}%).csv')
            matched_output.clear()
        total_observed = len(all_obs_frags[all_obs_frags['Count'] > 0])
        coverage_percents[structure] = round(
            100 * sum(all_obs_frags['Count']) / (total_theo_frags * len(scans_to_search)), 3)
        all_obs_frags.sort_values(by=['Mass']).to_csv(
            output_path / f'Observed Fragments ({total_observed} of {total_theo_frags}).csv')
    ms1_row["coveragePercents"] = str(coverage_percents)
    return ms1_row


if __name__ == "__main__":
    # FIXME: Add the ability to automatically run MS2 on 100% of structures!
    mols = load_structures(selected_structures)
    ms1_df = pd.read_csv(ms1_file)
    filtered_ms1 = filter_structures(ms1_df, mols)
    for struct in parse_structures(all_structures(filtered_ms1)):
        write_graphs(struct, graph_folder)
    output_dir /= datetime.now().strftime("%F (%R)")
    #print(*[x for x in map(resolve_row,
    #      [r for _, r in filtered_ms1.iterrows()])], sep='\n')
    #quit()
    with Pool() as p:
        rows = p.map(resolve_row, [r for _, r in filtered_ms1.iterrows()])
    df = pd.concat(rows, axis=1)[1:].T
    df.to_csv(output_dir / "MS2.csv", index=False)
