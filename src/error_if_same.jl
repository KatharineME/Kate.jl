using Dates

include("print_and_run_cmd.jl")

function error_if_same(
    _1::String,
    _2::String,
)

    if _1 == _2

        error("$_1 and $_2 are the same.")

    end

end
