include("src/driver.jl")

#atoms = [("O", [0., 0., 0.]),
#         ("H", [0., 2., 0.]),
#         ("H", [0., 0., 2.])]
atoms = [("O", [ 0.,            -0.143225816552, 0.]),
         ("H", [ 1.638036840407, 1.136548822547, 0.]),
         ("H", [-1.638036840407, 1.136548822547, 0.])]

basis = parse_basis(atoms, "sto-3g")
#summarise(basis)


#s = driver("int1e_ovlp", basis)
#
#show(stdout, "text/plain", s)
#println("")
#
#t = driver("int1e_kin", basis)
#
#show(stdout, "text/plain", t)
#println("")
#
#v = driver("int1e_nuc", basis)
#
#show(stdout, "text/plain", v)
#println("")

eri = driver("int2e", basis)

show(stdout, "text/plain", sum(eri, dims=(3,4)))
println("")
