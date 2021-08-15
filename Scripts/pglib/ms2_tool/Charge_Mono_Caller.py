import numpy as np
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from keras.models import load_model
import os.path
# from matplotlib import pyplot as plt

NEUTRON = 1.003
PROTON = 1.0078
MAX_CHARGE_STATE = 10
COMB_FILTER_TEETH_CHARGE = 4
COMB_FILTER_TEETH_MONO = 10
APEX_CHARGE_CLUSTERING = 4
GRANULARITY = 0.002
VECTOR_GRANULARITY = int(1 / GRANULARITY)
MZ_MATCH_INDEX = int(APEX_CHARGE_CLUSTERING * VECTOR_GRANULARITY)
VECTOR_ARRAY_LENGTH = int((MZ_MATCH_INDEX * 2) + 1)
ERROR_RANGE = 0.02
ERROR_RANGE_HASHED = int(ERROR_RANGE * VECTOR_GRANULARITY)
NN_SIMPLE_CHARGE_WEIGHTS = r"C:\Users\ankur\Documents\Code\PMI module\FOUR_APEX_OLD\Models\simple_weights.nn"
NN_MONO_WEIGHTS = r"C:/Users/ankur/Documents/Code/PMI module/FOUR_APEX_OLD/Models/"
MZ_END = 2000
APEX_CHARGE_CLUSTERING_HASHED = int(APEX_CHARGE_CLUSTERING * VECTOR_GRANULARITY)


def hash_mz(mz: float, granularity: float = VECTOR_GRANULARITY) -> int:
    return int(round(mz * granularity))


class Charge_Determinator:

    def __init__(self, neural_network_model):
        self.model = load_model(neural_network_model)
        self.successive_comb_filters = self.build_successive_comb_filter()

    def determine_charge(self, scan_slice):
        scan_slice = self.charge_preprocessor(scan_slice, self.successive_comb_filters)
        charge = self.model.predict(np.array([scan_slice]))[0]
        return np.where(charge == charge.max())[0][0] + 1, charge.max()

    def build_successive_comb_filter(self) -> np.array:

        succesive_comb_filters = []

        for charge in range(1, (MAX_CHARGE_STATE + 1)):
            dimer_charges = [1, 2, 3]
            charge_distance = hash_mz(NEUTRON / charge)

            if charge in dimer_charges:
                dimer_charge_distance = hash_mz(NEUTRON / (charge * 2))

            t_successive_combfilter = []
            teeth = charge if charge < COMB_FILTER_TEETH_CHARGE else COMB_FILTER_TEETH_CHARGE
            for z in range(-teeth, (teeth + 1)):
                t_vec = np.zeros(VECTOR_ARRAY_LENGTH)
                t_index = int(MZ_MATCH_INDEX + (z * charge_distance))

                if (t_index - ERROR_RANGE_HASHED > 0) & (t_index + ERROR_RANGE_HASHED < len(t_vec) - 1):
                    t_vec[(t_index - ERROR_RANGE_HASHED):(t_index + ERROR_RANGE_HASHED + 1)] = 1

                    t_successive_combfilter.append(t_vec)

                    if (charge in dimer_charges) & (z < teeth):
                        roll = np.roll(t_vec, dimer_charge_distance)
                        t_successive_combfilter.append(roll * -1)

            succesive_comb_filters.append(np.array(t_successive_combfilter))

        return np.array(succesive_comb_filters, dtype=object)

    def charge_preprocessor(self
                            , scan_array: np.array
                            , succesive_comb_filters: np.array) -> np.array:

        x = []
        for charge in range(MAX_CHARGE_STATE):
            t_product = (scan_array * succesive_comb_filters[charge]).sum(axis=1)
            t_product /= t_product.max()
            [x.append(xx) for xx in t_product]

        return np.array(x)



class Monoisotope_Determinator:


    def __init__(self, neural_network_models_path: str):
        self.neural_network_models_path = neural_network_models_path
        self.monoisotope_models = []
        self.successive_comb_filters = self.build_successive_comb_filter_mono()
        self.load_models()


    def load_models(self):

        for charge in range(1, MAX_CHARGE_STATE + 1):
            if not os.path.isfile(self.neural_network_models_path + ('mono_charge_%s.nn' % charge)):
                print("model not found")
            else:
                model = load_model(self.neural_network_models_path + ('mono_charge_%s.nn' % charge))
                self.monoisotope_models.append(model)


    def determine_monoisotope(self
                              , scan: np.array
                              , charge: int
                              , mz: float):

        processed_scan = self.mono_preprocessor(scan, mz, charge, self.successive_comb_filters)
        model = self.monoisotope_models[charge - 1]
        monoisotope = model.predict(np.array([processed_scan]))[0]

        try:
            mono_offset = np.where(monoisotope == monoisotope.max())[0][0], monoisotope.max()
            return mono_offset
        except:
            print(scan)




    def build_successive_comb_filter_mono(self) -> np.array:

        succesive_comb_filters = []

        for charge in range(1, (MAX_CHARGE_STATE + 1)):

            charge_distance = hash_mz(NEUTRON / charge)

            tt_successive_combfilter = []

            for rollit in range(charge + 2):
                t_successive_combfilter = []
                teeth = charge + 1
                for z in range(-teeth, 1):

                    t_vec = np.zeros(VECTOR_ARRAY_LENGTH)
                    t_index = MZ_MATCH_INDEX + (z * charge_distance) + (charge_distance * rollit)

                    if (t_index - ERROR_RANGE_HASHED > 0) & (t_index + ERROR_RANGE_HASHED < len(t_vec) - 1):
                        t_vec[(t_index - ERROR_RANGE_HASHED):(t_index + ERROR_RANGE_HASHED + 1)] = 1

                        if z == -teeth:
                            roll = np.roll(t_vec.copy(), -charge_distance)
                            t_successive_combfilter.append(roll * -1)

                        t_successive_combfilter.append(t_vec)

                tt_successive_combfilter.append(np.array(t_successive_combfilter))

            succesive_comb_filters.append(tt_successive_combfilter)

        return np.array(succesive_comb_filters, dtype=object)


    def mono_preprocessor(self
                          ,scan: np.array
                          , mz: float
                          , charge: int
                          , successive_comb_filters):

        divisor = 1 if np.max(scan[MZ_MATCH_INDEX - 2: MZ_MATCH_INDEX + 2]) == 0 else np.max(scan[MZ_MATCH_INDEX - 2: MZ_MATCH_INDEX + 2])
        scan /= divisor
        x = [(charge * mz) / 10000]

        for t_succesive_comb_filter in successive_comb_filters[charge - 1]:
            t_product = (scan * t_succesive_comb_filter)
            [x.append(xx) for xx in t_product.sum(axis=1)[1:]]

        return x



class Charge_Mono_Caller:

    def __init__(self
                 , charge_weights_path :str
                 , mono_charge_path: str):
        self.charge_determinator = Charge_Determinator(charge_weights_path)
        self.mono_determinator = Monoisotope_Determinator(mono_charge_path)


    def process(self, point2dList :list, mz :float):

        vec = self.convert_point2dlist_to_vector(point2dList)
        vec = self.extract_mz_vector(vec, mz)

        mz = float(mz)



        charge, score_charge \
            = self.charge_determinator.determine_charge(vec)

        if charge > 0:


            mono_offset, score_mono \
                = self.mono_determinator.determine_monoisotope(vec, charge, mz)



            result = {}
            result['charge'] = charge
            result['mono_offset'] = mono_offset
            result['mz_observed'] = mz
            result['mz_mono_theo'] = self.calculate_mono_mz(mz, charge, mono_offset)
            result['monoisotopic_mass'] = self.calculate_mono_mw(mz, charge, mono_offset)
            result['score_charge'] = score_charge
            result['score_mono'] = score_mono

            return result

        return {}


    def calculate_mono_mw(self, mz :float, charge :int, mono_offset :int)-> float:
        return (mz * charge) - ((charge + mono_offset) * PROTON)


    def calculate_mono_mz(self, mz :float, charge :int, mono_offset :int)-> float:
        return mz - ((mono_offset * PROTON) / charge)


    def extract_mz_vector(self, vec :np.array, mz :float) -> np.array:
        mz_hashed = hash_mz(mz)
        return vec[mz_hashed - APEX_CHARGE_CLUSTERING_HASHED:mz_hashed + APEX_CHARGE_CLUSTERING_HASHED + 1]


    def convert_point2dlist_to_vector(self, point2dList :list)->np.array:

        vec_length = int(MZ_END * VECTOR_ARRAY_LENGTH)

        vec = np.zeros(vec_length)
        for point in point2dList:
            mz_hashed = hash_mz(point[0])

            if (mz_hashed < 0) or (mz_hashed > vec_length):
                continue

            vec[mz_hashed] = point[1]

        return vec



##### Usage Example ###########
#
# charge_mono_caller = Charge_Mono_Caller(NN_SIMPLE_CHARGE_WEIGHTS, NN_MONO_WEIGHTS)
#
# scan_list = [(600,1),(601,0.5),(602,0.25),(800,1),(800.5,0.6),(801,0.3)]
# print(charge_mono_caller.process(scan_list, 600))
# print(charge_mono_caller.process(scan_list, 800.5))
