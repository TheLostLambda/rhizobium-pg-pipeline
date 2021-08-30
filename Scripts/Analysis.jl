using CSV
using DataFrames
using Statistics

"Read an MS1 file into a DataFrame – dropping unresolved structures"
load_ms1(csvfile) = dropmissing(DataFrame(CSV.File(csvfile)), :inferredStructure)

"Split a String listing of structures in a Vector of Strings"
split_structures(st) = [String(m.captures[1]) for m in eachmatch(r"([A-Z-=~]+)", st)]

"Filter a DataFrame for rows containing a particular structure"
function find_structure(ms1data, structure)
  filter(:inferredStructure => st -> structure in split_structures(st), ms1data)
end

"Get all structures from a DataFrame, returning a Vector"
function all_structures(ms1data)
  Set([s for ss in ms1data[:, :inferredStructure] for s in split_structures(ss)])
end

"Write structures to a CSV file for MS2 processing"
function write_csv(file, structures)
  CSV.write(file, DataFrame(structure=collect(structures)))
end

"Summarise MS1 data"
summarise(ms1data) = select(ms1data, :ionCount, :rt, :mwMonoisotopic, :inferredStructure)

"Load MS2 scores from the filenames within a directory"
function load_ms2_scores(ms2dir)
  frag_files = filter(f -> '%' in f, readdir(ms2dir))
  [parse(Int, split(f, '%')[1]) for f in frag_files]
end

"Summarise MS2 runs using several key statistics"
function summarise_ms2(ms2dir)
  ms2_scores = load_ms2_scores(ms2dir)
  summary_file = [f for f in readdir(ms2dir, join=true) if occursin("Observed", f)][1]
  summary_match = match(r"\d+ of (\d+)", summary_file)
  summary_data = DataFrame(CSV.File(summary_file))[:, :Count]
  coverage = sum(summary_data) / (parse(Int, summary_match.captures[1]) * length(ms2_scores))
  print("""
  Matched Scans: $(length(ms2_scores))
  Average Percent Observed: $(round(mean(ms2_scores), digits=2))%
  Total Fragments Observed: $(summary_match.match)
  Coverage Percent: $(round(coverage * 100, digits=2))%
  """)
end

"Performs a simple analysis on the data"
function pipeline()
  ms = load_ms1("Data/Outputs/MS1.csv")
  selected = find_structure(ms, "GM-AEJA=GM-AEJ") |> all_structures
  write_csv("Data/Outputs/selected_structures.csv", selected)
end
