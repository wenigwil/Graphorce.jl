using Random
using LinearAlgebra
using CUDA
using BenchmarkTools
using Test

function diag_dynmats(dynmats)
    #Diagonalization of the the dynamical matrix on the CPU.
    #This mutates the argument dynmats. This is bad coding practice...

    eigen(dynmats)

    return 1
end

function diag_dynmats_gpu(dynmats)
    #Diagonalization of the the dynamical matrix on the GPU.

    #Batched over the last index.
    #Works on both complex and floats according to
    #https://github.com/JuliaGPU/CUDA.jl/blob/master/lib/cusolver/dense.jl
    CUDA.@sync sols = CUDA.CUSOLVER.heevjBatched!('V', 'U', dynmats)

    #Only seems to work on floats
    #CUDA.@sync sols = CUDA.CUSOLVER.syevjBatched!('V', 'U', dynmats)
end

function calc_phases(rs, q)
    #Multithreaded CPU calculation of phases

    shape = size(rs)
    nrs = shape[2]

    phases = zeros(ComplexF32, nrs)

    Threads.@threads for ir in 1:nrs
        phases[ir] = exp(2.0f0π * 1im * dot(rs[1:3, ir], q))
    end

    return phases
end

function calc_phases_gpu!(rs, q, phases)
    #GPU calculation of phases

    CUDA.@sync phases = map(
        iθ -> exp(iθ),
        2.0f0π * 1im .* (rs[1, :] .* q[1] + rs[2, :] .* q[2] + rs[3, :] .* q[3]),
    )
end

function test_diag_dynmats()
    #Multithreaded CPU diagonalization of the dynamical matrix

    nqs = 1000000
    natoms = 1
    nbands = 3 * natoms

    dynmats = fill(rand(ComplexF32), (nbands, nbands, nqs))

    #Make Hermitian
    for iq in 1:nqs
        dynmats[:, :, iq] = hermitianpart(dynmats[:, :, iq])
    end

    Threads.@threads for iq in 1:nqs
        diag_dynmats(dynmats[:, :, iq])
    end
end

function test_calc_phases()
    #Just testing the phase calculation on the CPU for a fixed q vector.
    #Also, serves as a benchmark for CPU vs GPU computations.

    nrs = 100000000
    rs = fill(1.0f0, (3, nrs))
    q = [1.0f0 / 8.0f0, 0.0f0, 0.0f0]

    phases = calc_phases(rs, q)

    @test all(phases .≈ (1.0f0 + 1.0f0im) / sqrt(2.0f0))
end

function test_calc_phases_gpu()
    #Just testing the phase calculation on the GPU for a fixed q vector.
    #Also, serves as a benchmark for CPU vs GPU computations.

    nrs = 100000000
    rs = CUDA.fill(1.0f0, (3, nrs))
    q = [1.0f0 / 8.0f0, 0.0f0, 0.0f0]

    phases = CUDA.fill(0.0f0, nrs)
    phases = calc_phases_gpu!(rs, q, phases)

    @test all(phases .≈ (1.0f0 + 1.0f0im) / sqrt(2.0f0))
end

function test_diag_dynmats_gpu()
    nqs = 1000000 #gpu won't be able to handle more
    natoms = 1
    nbands = 3 * natoms

    #Construct the dynamical matrix. In real life, we will have
    #to construct it by Fourier transforming the 2nd order force constants tensor.
    #dynmats = fill(rand(ComplexF64), (nbands, nbands, nqs))
    dynmats = fill(rand(ComplexF32), (nbands, nbands, nqs))

    #Make Hermitian
    for iq in 1:nqs
        dynmats[:, :, iq] = hermitianpart(dynmats[:, :, iq])
    end

    #Copy to device
    dynmats_cu = CuArray(dynmats)

    solutions = diag_dynmats_gpu(dynmats_cu)
end

#@btime test_calc_phases()
#@btime test_calc_phases_gpu()

@btime test_diag_dynmats()
@btime test_diag_dynmats_gpu()
