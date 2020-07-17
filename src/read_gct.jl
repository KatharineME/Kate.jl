using CSV
using DataFrames


function read(gct_file_path::String, column_1_name::String)

    tsv = select(
        CSV.read(gct_file_path; header = 3, delim = '\t',),
        Not(Symbol("Description")),
    )

    names!(tsv, replace(names(tsv), Symbol("Name") => Symbol(column_1_name)))

    return tsv

end
