import pandas as pd
import struct
import sqlite3
import decimal as dec
import copy


pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)


class Byspec_Reader:

    def __init__(self, file_path: str):

        # variables
        self.file_path = file_path

        # These Data structures need to be saved for use in next class
        self.scans_details = None
        self.scans = None

        # function calls
        self.read_byspec2_file()

    def read_byspec2_file(self):

        with sqlite3.connect(self.file_path) as conn:

            sql = "SELECT * FROM Spectra ORDER BY Id ASC"
            self.scans_details = pd.read_sql(sql, conn)

            self.scans_details['ObservedMz'] = self.scans_details['ObservedMz'].astype(
                str)

            sql = "SELECT * FROM Peaks"
            self.scans = pd.read_sql(sql, conn)

            # print(self.scans_details.dtypes)

            # print(self.scans_details.columns)
            # print(self.scans_details[self.scans_details.MSLevel == 2])

    def convert_binary_arrays_to_scan_points(self, peaks_mz, peak_intensity) -> list:

        d_amount = ('d' * int(len(peaks_mz) / 8))
        i_amount = ('f' * int(len(peak_intensity) / 4))
        mz_decoded = struct.unpack(d_amount, peaks_mz)
        intensity_decoded = struct.unpack(i_amount, peak_intensity)

        return list(zip(mz_decoded, intensity_decoded))

    def get_scan_by_scan_number(self, scan_number: int) -> list:

        scan = self.scans.loc[self.scans.Id == scan_number]
        scan_converted \
            = self.convert_binary_arrays_to_scan_points(scan.PeaksMz.values[0], scan.PeaksIntensity.values[0])
        return(scan_converted)

    def get_parent_scan_by_scan_number(self, scan_number: int) -> int:
        return int(self.scans_details["ParentScanNumber"][self.scans_details["ScanNumber"] == scan_number].iloc[0])

    def filter_children_with_parent_in_range(self, start, end):
        self_copy = copy.deepcopy(self)
        scan_nums = self.scans_details["ScanNumber"]
        self_copy.scans_details = self.scans_details[(start <= scan_nums) & (scan_nums <= end)]
        return self_copy

    def get_scan_by_observed_mz(self, mz_number: str) -> list:

        for index, row in self.scans_details.iterrows():
            if row['ObservedMz'] == mz_number:
                scan_number = row['ScanNumber']
                charge_state = row['ChargeList']
                fragmentation_type = row['FragmentationType']
                print(scan_number)
                print(charge_state)
                print(fragmentation_type)
                scan = self.scans.loc[self.scans.Id == scan_number]
                scan_converted \
                    = self.convert_binary_arrays_to_scan_points(scan.PeaksMz.values[0], scan.PeaksIntensity.values[0])
                print(scan_converted)

    def get_scan_mz_charge(self):
        dec.getcontext().prec = 8

        scan_mz_charge = []
        ms2_scans = self.scans_details[self.scans_details.MSLevel == 2]

        for index, row in ms2_scans.iterrows():
            parent_scan_number = int(row['ParentScanNumber'])
            ms2_scan = int(row['ScanNumber'])
            observedMz = dec.Decimal(row['ObservedMz'])
            charge_state = int(row['ChargeList'])
            scan_mz_charge.append(
                [parent_scan_number, observedMz, charge_state, ms2_scan])

        return(scan_mz_charge)

# Usage
# byspec = Byspec_Reader(r"C:\Users\Hyperion\Documents\GitHub\ms2_graph_tool\OT_190122_APatel_Efaecalis_EnpA_10mAU.raw.byspec2")
# byspec.get_scan_mz_charge()
# byspec.get_scan_by_scan_number(668)
# byspec.get_scan_by_observed_mz("471.713958740234")
