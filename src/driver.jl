# Driver for populating arrays of integrals

include("int1e_ovlp.jl")
include("int1e_kin.jl")
include("int1e_nuc.jl")


@doc """
Driver function for the module. Computes requested integral and
populates array of the suitable size.

# Arguments
- `int::String`: the type of integral to compute.
- `basis::Basis`: the basis set struct.
""" ->
function driver(int::String, basis::Basis, args...; kwargs...)
    fint, shape = _get_driver_info(int, basis)
    
    out = fint(basis, args...; kwargs...)

    out
end


function _get_driver_info(int::String, basis::Basis)
    nao = basis.size

    if int == "int1e_ovlp"
        fint = _int1e_ovlp_matrix
        shape = (nao, nao)
    elseif int == "int1e_kin"
        fint = _int1e_kin_matrix
        shape = (nao, nao)
    elseif int == "int1e_nuc"
        fint = _int1e_nuc_matrix
        shape = (nao, nao)
    else
        throw(DomainError(int, "does not support " * int))
    end

    fint, shape
end
