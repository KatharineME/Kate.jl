using Kate

using Revise

using JSON

pase = abspath(readlink("setting.json"))

SE = JSON.parsefile(pase)

PAK = joinpath(dirname(dirname(dirname(pase))), "Kate.jl")

PAI = joinpath(PAK, "input/")

PAO = joinpath(PAK, "output/")

PASN = "/opt/snpeff/snpEff/snpEff.jar";

println("Settings are loaded.\n")

Kate.test()

println("\nEnvironment passed.")
