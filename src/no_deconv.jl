function build_dynamical_matrix_no_deconv(
    lattvecs::Matrix{Float64},
    ifc2::Array{Float64,7},
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
    qpoint_cryst::Vector{Float64},
)
    numatoms = size(basisatom2species, 1)
    dynmat = zeros(ComplexF64, (3 * numatoms, 3 * numatoms))

    # Lattvecs are given in nm we need the reclattvecs in bohr
    # Reciprocal lattice vectors have units of 1/length such that a 
    # multiplication with the Bohr-radius in nm will convert from 1/nm to 1/bohr
    reclattvecs = calc_reciprocal_lattvecs(lattvecs) * a0_nm
    qpoint_cart = reclattvecs * qpoint_cryst

    unitcell_mults = size(ifc2)[5:7]

    # Build vectors that point to unit cells
    unitpoints = Matrix{Float64}(undef, (prod(unitcell_mults), 3))
    unitpoints_ifc2 = Matrix{Int64}(undef, (prod(unitcell_mults), 3))
    counter = 1
    for i in 1:unitcell_mults[1]
        for j in 1:unitcell_mults[2]
            for k in 1:unitcell_mults[3]
                unitpoints[counter, :] .=
                    i * lattvecs[:, 1] + j * lattvecs[:, 2] + k * lattvecs[:, 3]
                unitpoints_ifc2[counter, :] .= [i, j, k]

                counter += 1
            end
        end
    end

    mass_prefactor = build_mass_prefactor(basisatom2species, species2mass)
    for iat in 1:numatoms
        for jat in 1:numatoms
            for icart in 1:3
                i = mux2to1(iat, icart)
                for jcart in 1:3
                    j = mux2to1(jat, jcart)

                    for l in 1:prod(unitcell_mults)
                        dynmat[i, j] += begin
                            ifc2[
                                icart,
                                jcart,
                                iat,
                                jat,
                                unitpoints_ifc2[l, 1],
                                unitpoints_ifc2[l, 2],
                                unitpoints_ifc2[l, 3],
                            ] *
                            exp(im * LinAlg.dot(qpoint_cart, unitpoints[l, :]))
                        end
                    end

                    dynmat[i, j] /= mass_prefactor[iat, jat]
                end
            end
        end
    end

    return dynmat
end
