import pandas as pd
from decimal import *

class Helper_Funcs:

    def __init__(self, nodes_filepath: str, edges_filepath: str):
        self.nodes_filepath = nodes_filepath
        self.edges_filepath = edges_filepath

    def nodes_df(self):
        nodes_df = pd.read_csv(self.nodes_filepath)
        return nodes_df

    def edges_df(self):
        edges_df = pd.read_csv(self.edges_filepath)
        return edges_df

    def generate_dict(self, dict_table_filepath):
        # read csv files
        dict_df = pd.read_csv(dict_table_filepath)
        # Set Value column to string (.apply decimal will convert entire float value to decimal type otherwise)
        dict_df['Value'] = dict_df['Value'].apply(str)
        # Set Value column to Decimal (all values have exactly the same decimal precision as values in input csv)
        dict_df['Value'] = dict_df['Value'].apply(Decimal)
        # Set Key column as index
        dict_df = dict_df.set_index('Key')
        # Tranpose dataframe so index value is header for each column (to_dict takes column headers as key value)
        transposed_df = dict_df.transpose()
        # Create dictionary
        output_dict = transposed_df.to_dict('list')
        return output_dict


