# Notes on the `kickstart.jl` script

## Script runinfo on personal machine
>[!IMPORTANT]
> Changed `nqs=500000` in the `test_diag_dynmats`-function
> 
> Ran on a GeForce GTX 1060 with 6 GB and Intel i5-7500 with 3.8 GHz


~~~txt
2.433 s (14000010 allocations: 1.33 GiB)
178.205 ms (2000073 allocations: 175.48 MiB)
~~~


## How is GPU offloading done in julia?

- Implemented using [`CUDA.jl`](https://cuda.juliagpu.org/stable/) which wraps (?) `CUDA`-libraries to make them available in julia and easy to use.


## Diagonalization

### CPU

- Using the base implementation from `LinearAlgebra`


### GPU

- `CUDA.CUSOLVER.heevjBatched!` is **H**ermitian Matrix **J**acobi based **E**igen **V**alue **D**ecomposition which originates from the `CUSOLVER` library [here](https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-syevjbatched)

> [!CAUTION]
> What is batching over an index?

>[!CAUTION]
> What does an expression of the type `CUDA.@sync` mean in front of a variable (syntactically in general)?



## Questions

**1.** Would it be less bad coding practice if the function tells us that it mutates
    it's arguments?

~~~julia
function diag_dynmats(dynmats) 
    eigen(dynmats)
    return 1
end
~~~

