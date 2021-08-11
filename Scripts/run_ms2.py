import sys
sys.path.append("Libraries")

import pandas as pd
from ms2_tool import Bio_Graph as BG
from ms2_tool import Helper_Funcs as hfunc
import decimal as dec
import networkx as nx
from ms2_tool import Byspec_Reader as BR
from ms2_tool import Charge_Mono_Caller as CMC
import os

# Change to relevant graph files and data
data_file = "Data/Inputs/20210618_RhiLeg_ndslt_TY_1.raw.byspec2"
TestNL = "Data/Outputs/Graphs/GM-AEJ=GM-AEJA (3-4) NL.csv"
TestEL = "Data/Outputs/Graphs/GM-AEJ=GM-AEJA (3-4) EL.csv"

# Do not change
mass_table = "Data/Constants/masses_table.csv"
mod_table = "Data/Constants/mods_table.csv"


hf = hfunc.Helper_Funcs(TestNL, TestEL)
nodes_from_file = hf.nodes_df()
edges_from_file = hf.edges_df()
mass_dict = hf.generate_dict(dict_table_filepath=mass_table)
mod_dict = hf.generate_dict(dict_table_filepath=mod_table)


# test_graph = BG.Bio_Graph(nodes_from_file,edges_from_file,mass_dict,mod_dict)
# mg1 = test_graph.construct_graph()
# test_graph.draw_graph(mg1)
#
# output = test_graph.fragmentation(mg1,cut_limit=3)
#
# l1,l2,l3 = test_graph.sort_fragments(output)
#
# nlist = test_graph.monoisotopic_mass_calculator(graph_fragments=output,graph_IDs=l1)
# clist = test_graph.monoisotopic_mass_calculator(graph_fragments=output,graph_IDs=l2)
# ilist = test_graph.monoisotopic_mass_calculator(graph_fragments=output,graph_IDs=l3)
# enabled = ['y']
# mzlist_graph = test_graph.generate_mass_to_charge_masses(nlist,clist,ilist,enabled_ions=enabled,charge_limit=2)
#
# print(mzlist_graph)
#

bio_graph = BG.Bio_Graph(nodes_from_file, edges_from_file, mass_dict, mod_dict)


def calculate_ppm_tolerance(mass, ppm_tol):
    return (mass*ppm_tol) / 1000000


def ppm_error(obs_mass, theo_mass):
    return (1-(dec.Decimal(obs_mass)/theo_mass))*1000000


def autosearch(selected_ions, user_set_cuts=1, user_set_charge=1, intact_ppm_tol: str = '10', frag_ppm: str = '20'):

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
        print(mass[0])
        scans_to_search = []
        upper_mass_lim = mass[0] + calculate_ppm_tolerance(mass[0], i_ppm)
        lower_mass_lim = mass[0] - calculate_ppm_tolerance(mass[0], i_ppm)

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

            except:
                print('-' * 20)
                print('parent ion not found in scan')
                print('scan: ',  scan_mz_charge_tuple[0])
                print('mz: ', scan_mz_charge_tuple[1])
                print('charge: ', scan_mz_charge_tuple[2])
                print('#' * 20)

        if not scans_to_search:
            print('scan_to_search is empty')

        for scan_number in scans_to_search:
            scan = byspec_reader.get_scan_by_scan_number(scan_number)
            graph = nx.Graph(molecules[graph_ID])
            fragments = bio_graph.fragmentation(graph, cut_limit=user_set_cuts)
            total_theo_frags = len(fragments)
            n_frag, c_frag, i_frag = bio_graph.sort_fragments(fragments)
            nlist = bio_graph.monoisotopic_mass_calculator(fragments, n_frag)
            clist = bio_graph.monoisotopic_mass_calculator(fragments, c_frag)
            ilist = bio_graph.monoisotopic_mass_calculator(fragments, i_frag)
            frag_ions_df = bio_graph.generate_mass_to_charge_masses(
                nlist, clist, ilist, selected_ions, user_set_charge)
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
                            fragment_nodes = str(fragment.nodes)
                            frag_structure.append(fragment_nodes)
                            matched_output.append(
                                (obs_mz, ion, charge, count, ppm_diff, ion_type, frag_structure))
                            frag_structure = []

            matched_output_df = pd.DataFrame(matched_output, columns=[
                                             'Observered Ion', 'Theoretical Ion', 'Charge', 'count', 'PPM Error', 'Ion Type', 'Structure'])
            matched_num = len(matched_output_df)
            score = round(((matched_num/total_theo_frags)*100))
            desired_width = 3840
            pd.set_option('display.width', desired_width)
            print(matched_output_df)
            output_dir_path = "Data/Outputs/MS2/"
            output_folder = "/" + str(scan_number) + ' Score ' + str(score)
            os.mkdir(output_dir_path + output_folder)
            matched_output_df.to_csv(output_dir_path + output_folder +
                                     "/scan" + str(scan_number) + "score" + str(score) + ".csv")

            matched_output.clear()

    return None


enabled_ions = ['y', 'b', 'i']
autosearch(selected_ions=enabled_ions, user_set_cuts=3,
           user_set_charge=2, intact_ppm_tol=10)
