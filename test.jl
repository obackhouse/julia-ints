include("src/driver.jl")

basis = parse_basis(
    "sto-3g", 
    ["O", "H", "H"],
    [[0.,0.,0.] [0.,0.,2.] [0.,2.,0.]]
)

s = driver("int1e_ovlp", basis)

show(stdout, "text/plain", s)
println("")

t = driver("int1e_kin", basis)

show(stdout, "text/plain", t)
println("")
