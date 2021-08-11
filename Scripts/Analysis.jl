using CSV
using DataFrames

"Read an MS1 file into a DataFrame – dropping unresolved structures"
load_ms1(csvfile) = dropmissing(DataFrame(CSV.File(csvfile)), :inferredStructure)

"Split a String listing of structures in a Vector of Strings"
split_structures(st) = [String(m.captures[1]) for m in eachmatch(r"([A-Z-]+)", st)] 

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
  CSV.write(file, DataFrame(structure = collect(structures)))
end

"Summarise MS1 data"
summarise(ms1data) = select(ms1data, :ionCount, :rt, :mwMonoisotopic, :inferredStructure)
