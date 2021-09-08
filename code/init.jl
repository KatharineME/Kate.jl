using Kate

using Revise

using JSON

pase = abspath(readlink("setting.json"))
    
se = JSON.parsefile(pase)

par = dirname(dirname(dirname(pase)))

pai = joinpath(par, "input/")

pao = joinpath(par, "output/")

pasn = "/opt/snpeff/snpEff/snpEff.jar";

Kate.test()

println("\nSettings are loaded and environment tests passed.")