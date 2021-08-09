import pgfinder.matching as matching
import pgfinder.validation as validation

csv_filepath = "Data/RL_DB2.csv"
ftrs_filepath = "Data/20210618_RhiLeg_ndslt_TY_1.ftrs"

raw_data = matching.ftrs_reader(ftrs_filepath)
validation.validate_raw_data_df(raw_data)

theo_masses = matching.theo_masses_reader(csv_filepath)
validation.validate_theo_masses_df(theo_masses)

mod_test = ["Decay", "Multimers", "multimers_Glyco"]
validation.validate_enabled_mod_list(mod_test)

results = matching.data_analysis(raw_data, theo_masses, 0.5, mod_test, 10)

results.to_csv("Data/ms1.csv")
