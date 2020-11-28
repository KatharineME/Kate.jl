function error_if_equal(_1::String, _2::String)

    if _1 == _2

        error("$_1 == $_2.")

    end

    return nothing

end
