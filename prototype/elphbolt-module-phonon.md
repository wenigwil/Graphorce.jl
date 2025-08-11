# elphbolt-module-phonon
This is a place for notes on the `phonon.f90` which is located at 
`elphbolt/src/`.

## Order of calls

Here I just want to catch which subroutines are being called by the public 
subroutine `initialize()`. I imagine this to be quite helpful when I 
rebuild...

**0.** Initialize quantities: number of bands, wave-vector mesh, number of 
wave-vectors

**1.** `read_ifc2()`

**2.** `phonon_espresso_precompute()`

**3.** `calculate_phonons()`

**4.** `read_ifc3()`


## types

### phonon

- **`prefix = 'ph'`** *`string`*
    - prefix for identifying particle type
- **`scell(3)`** *`integer`*
    - number of unitcells in the supercell used for DFPT or finite 
      displacement methods
- **`ifc3(:,:,:,:)`** *`array of floats`*
    - 3rd order force constant tensor with (probably) muxed indices
- **`numrtiplets`** *`integer`*
    - number of triplets in the ifc3-file
    - triplets in the sense analogous to 2nd order force constant atom 
      pairs
- **`R_j(:,:), R_k(:,:)`** *`array of floats`*
    - positions of the 2nd and 3rd unitcell in the supercell for an ifc3 
      triplet
- **`Index_i(:), Index_j(:), Index_k(:)`** *`vector of integers`*
    - labels of primitive cell atoms in the ifc3 triplet
- **`simplex_squared_evals`** *`array of floats`*
    - simplex (???) vertices filled with squared eigenvalues
    - i don't understand any part of this yet
- **`cg_indexlist_irred(:)`** *`vector of integers`*
    - probably a list with group theory contents
- **`cg_indexlist`** *`vector of integers`*
    - probably a list with group theory contents
- **`cgset_indexlist_irred(:)`** *`vector of integers`*
    - probably a list with group theory contents
- **`cgset_indexlist`** *`vector of integers`*
    - probably a list with group theory contents
- **`ifc2(:,:,:,:,:,:,:)`** *`array of floats`*
    - 2nd order force constant tensor
- **`rws(124,0:3)`** *`array of floats`*
    - positions of supercells from stretched lattice vectors of the 
      unitcell
    - stretched by `scell(:)` 
    - supercells are part of a fixed `5x5x5` super-supercell
- **`cell_r(1:3, 0:3)`** *`array of floats`*
    - transposed `crys%lattvecs` (lattvecs are put into inputfile in 
      nanometers) and converted to units of bohr ( by *dividing* with 
      `bohr2nm`)
    - columns of `cell_r` are the cartesian axis and rows number the 
      lattice vectors
    - IMPORTANT: first columns denotes the length of the vector
- **`cell_g(1:3, 0:3)`** *`array of floats`*
    - analogous to `cell_r`
    - transposed of `crys%reclattvecs` which are computed from the 
      `crys%lattvecs`
    - converted from 1/nm to 1/bohr (by multiplying with bohr2nm)
    - first column contains the length of the vector
- **`mm(:,:)`** *`array of floats`*
    - symmetric matrix that contains all masses of all atoms in the 
      unitcell (by the counts of basis vectors *not* species) on the 
      diagonal
    - off-diagonal elements contain the mass factors (sqrt(mi * mj)) of all 
      atom combinations
- **`rr(:,:,:)`** *`array of floats`*
    - tensor that holds the difference of all atom positions in the unit 
      cell
    - can be thought of as a matrix of vectors. Diagonal is the zero-vector
    - off-diagonal elements are the difference in atom positions in the 
      unitcell (for 2 atoms in the unitcell this is of shape *2x2*x3)
- **`ws_cell(:, :)`** *`array of integers`*
    - IMPORTANT: this is *easily* (bad naming) mistaken for `wscell(3)` 
      which contains the stretched unitcell lattice vectors **This is 
      something else!**
    - allocated with size `(3, product(num_ws_cells)*crys%nu atoms**2)` 
      which is as if you save a vector of size equal to `ws_weights` in 
      each column of `ws_cell`
    - takes part in the `dynmat`-symmetry restore process (`wsweights`)
- **`ws_weight(:)`** *`array of floats`*
    - contains the computed `weight` corresponding to each `r_ws(3)` (see 
      `phonon_espresso_precompute()` for what `r_ws(3)` is)
    - first I thought this is storing a weight for every unitcell in the 
      ultracell, but it stores a weight for every atom pair combination for 
      every unitcell position in the ultracell (there are 4 atom pair 
      combinations if you got 2 atoms in the unitcell)
    - allocated with length `product(num_ws_cells)*crys%numatoms**2`. 
      `product(num_ws_cells)` is the number of all unitcell positions in 
      the ultracell (check the definition of `num_ws_cells` at 
      `phonon_espresso_precompute()`).

## `initialize`

This subroutine is said to initialize the phonon data type, calculate 
*ground state* phonon properties and read ifc3 data.
1. The subroutine begins with type declarations (?) of class-variables for 
   the modules crystal, symmetry, numerics and wannier! So far so good. 
  
~~~fortran
class(phonon), intent(out) :: self
type(crystal), intent(in) :: crys
type(symmetry), intent(in) :: sym
type(numerics), intent(in) :: num
type(wannier), intent(in), optional :: wann
~~~

2. Interestingly some variables are now defined that look like they are 
   type inferred

~~~fortran
!Set phonon branches
self%numbands = crys%numatoms*3
!Set wave vector mesh
self%wvmesh = num%qmesh
!Set number of phonon wave vectors
self%nwv = product(self%wvmesh(:))
~~~

Pretty self-explanatory what those variables mean, but why can you just 
define them like that. Didn't every variable in fortran need to have an 
explicitly defined type associated with itself?

3. Checking whether to use ifc2s that stem from a Wannier calculation. If 
   Wannier-calculated ifc2s should be used this is skipped.

~~~fortran
if(.not. num%use_Wannier_ifc2s) then
   !Read ifc2 and related quantities
   call read_ifc2(self, crys)

   !Precompute dynamical matrix related quantities
   call phonon_espresso_precompute(self, crys)
end if
~~~

4. If wannier object is present, call `calculate_phonons()` with this 
   `wann`-object. If not do not.

~~~fortran
if(present(wann)) then
   call calculate_phonons(self, crys, sym, num, wann)
else
   call calculate_phonons(self, crys, sym, num)
end if
~~~

5. Only do `read_ifc3()` if `num%onlyebte` is false and `num%runlevel` is 
   3. 

~~~fortran
if(.not. num%onlyebte .and. num%runlevel /= 3) then
   !Read ifc3s and related quantities
   call read_ifc3(self, crys)
end if
~~~


## `phonon_espresso_precompute`

This subroutine will compute quantities for constructing the dynamical 
matrix. These quantities do not need $\bar{q}$-vectors as input

### Variables

- **`i, j`** *`integer`*
    - usually cartesian counting indeces
- **`iat, jat`** *`integer`*
    - dummies for looping through the number of atoms
- **`m1 ,m2, m3`** *`integer`*
    - dummies for looping through the number of unitcells in the 
      super-supercell
- **`num_ws_cell(3)`** *`integer`*
    - `self%scell` has the *number of unitcells* that makes up a supercell 
      in each lattice vector direction. `num_ws_cell(:) = 4*self%scell(:) + 
      1` converts this to the number of positions of unitcells contained in 
      the ultracell in each lattice vector direction
    - number of vectors that point to unitcell origins
    - example for a 6x6x6 supercell: for all `i` we have 
      `num_ws_cell(i)=4*6+1`
    - good name for the super-supercell would be **ultracell**
    - **IMPORTANT** this is not to be confused with number of unitcells in 
      the ultracell which would be just `4*self%scell(i)` for a given 
      lattice vector direction `i`. This is the famous duality between the 
      number of connections between points and the number of points.
- **`counter`** *`integer`*
    - Counter for the muxed unitcell positions inside the **ultracell**
- **`ir`** *`integer`*
    - Dummy for looping through the muxed supercell positions inside the 
      ultracell analogous to `counter`
- **`deg`** *`integer`*
    - Counter for the "degeneracy" of how many voronoi planes a unitcell 
      position  is an element of
- **`distance`** *`float`*
    - explicit value of the Hesse Normalform of the voronoi plane that is 
      checked against, whether a point is on i) the side of the origin, ii) 
      opposite side of the origin or iii) on the voronoi plane
- **`weight`** *`float`*
    - counter defined as the inverse unit-incremented `deg` (`1/(deg+1)`)
- **`r_ws`** *`vector of floats`*
    - vector that points to a unitcell inside the ultracell (see `t(3)`) 
      but add the difference of two atoms `i` and `j` to it (see `rr(:,:)`)
    - depending on the order of `i` and `j` this points somewhere between 
      the atoms
    - `r_ws` only exists temporary in each computation for one atom pair 
      combination for one unitcell at a time and is overwritten every 
      iteration of the loop it takes part in
- **`t(3)`** *`vector of floats`*
    - vector that points to one unitcell position inside the ultracell

### Operations

1. Initialize `ws_weight` and `ws_cell` (which has nothing to do with 
   `wscell`)

~~~fortran
allocate(self%ws_weight(product(num_ws_cells)*crys%numatoms**2))
allocate(self%ws_cell(3, product(num_ws_cells)*crys%numatoms**2))
~~~

2. Set all values (?) of `ws_weight` and `ws_cell` to zero as well as 
   `counter`

~~~fortran
self%ws_weight = 0
self%ws_cell = 0.0_r64
counter = 0
~~~

3. Set up loops for dummies `iat, jat, m1, m2, m3` and *wire in* `counter`

~~~fortran
do iat = 1, crys%numatoms
   do jat = 1, crys%numatoms
      do m1 = -2*self%scell(1), 2*self%scell(1)
         do m2 = -2*self%scell(2), 2*self%scell(2)
            do m3 = -2*self%scell(3), 2*self%scell(3)
               counter = counter + 1
~~~

4. Compute the values for `r_ws(3)`. `r_ws` is a temporary variable which 
   is only consistent with itself inside one loop iteration.

~~~fortran
do i = 1, 3
  t(i) = m1*self%cell_r(1, i) + m2*self%cell_r(2, i) + 
  m3*self%cell_r(3, i)
  r_ws(i) = t(i) + self%rr(iat, jat, i)
end do
~~~

5. This part is mostly identical to a function found at Quantum Espresso is 
   `PW/src/wsweights.f90`. The current `r_ws(3)` is checked against postion 
   relative to a voronoi (wigner-seitz) plane between the origin and `rws` 
   (which is the position of a supercell inside the ultracell). `distance` 
   equates to the plane condition with three cases checked corresponding to 
   the 3 regions the plane is dividing the space into (Nullset included). 
   
    - Case i) `distance > 0` which means `r_ws` is on the side of the plane 
      that `rws` points towards (opposite of the origin).
    - Case ii) `abs(distance) < 1.0e-6_r64` which mean that `distance` is 
      close to 0. `r_ws` most certainly very close or on the plane.
    - Case iii) `distance < 0` which means that the point `r_ws` is not on 
      the plane but on the side of the origin.

    In Case i) `weight` will stay 0 but `j=1` (`j` is just for sorting out 
    this case). For Case ii) the number `deg = deg + 1` and `weight = 
    1/deg` (this weight will be added to `ws_weight` which is the total 
    weight corresponding to a particular atom pair combination and unitcell 
    in the ultracell). In the *hidden* Case iii) `j` will stay 0 and 
    `deg=1` such that `weight=1`.

~~~fortran
weight = 0.0_r64
deg = 1
j = 0
do ir = 1, 124
  distance = dot_product(r_ws, self%rws(ir, 1:3)) - self%rws(ir, 0)
  if(distance > 1.0e-6_r64) then
     j = 1
     cycle
  end if
  if(abs(distance) < 1.0e-6_r64) then
     deg = deg + 1
  end if
end do

if(j == 0) then
  weight = 1.0_r64/dble(deg)
end if
self%ws_weight(counter) = self%ws_weight(counter) + weight
~~~

6. The specifics about the `mod(A,P)` function are 
   [here](https://gcc.gnu.org/onlinedocs/gfortran/MOD.html). The 
   `mod`-function in Fortran is basically returning `A-int(A/P)*P` where 
   `int(A/P)` is [**rounding towards 
   zero**](https://en.wikipedia.org/wiki/Rounding#Rounding_toward_zero) 
   (mathematically
   $\mathrm{trunc}(x) = \mathrm{sgn}(x) \lfloor | x | \rfloor$). 
   `mod`-functions in every programming language are subject to 
   implementation which is due to the freedom of convention. In Fortran the 
   `mod`-function is not implemented *clock-like*. A clock-like 
   `mod`-function will always return values between `0` and `P`, whereas 
   Fortrans `mod`-functions' return values have `-P < 0` for `A<0` and `0 > 
   P` for `A>0`. The dummy variables `m1, m2, m3` are non-cartesian 
   postions of unitcells in the ultracell. `self%scell` is a copy from the 
   ifc2-file of the number of unit cells in a supercell in each lattice 
   vector direction.

   This code is inside the loop over all atom pair combinations for all 
   unitcell positions and only takes place for `r_ws` on or inside some 
   voronoi plane. We muxed the tuple `(iat, jat, m1, m2, m3)` into a number 
   which is `counter`. Know after we have gone through a `weight` 
   assignement for an `r_ws` we will save a new unitcell position in 
   crystal coordinates. This new unitcell position is saved at 
   `self%ws_cell(:,counter)` and is calculated modulo the *supercell*. The 
   Fortran `mod`-function works not clock-like but is converted into 
   clock-like modulo with an if-statement below that adds (`mod(A,P)`) `P` 
   in the case for `A<=0`.

   For some intellectual reason we shift Fortrans `mod`-function by `1` 
   (e.g. `mod(m1 + 1, self%scell(1))`). This causes the new vector for each 
   `r_ws` to never fold back into the origin unitcell. This will probably 
   clear up in the `phonon_espresso()`-subroutine.


~~~fortran
if(weight > 0.0_r64) then
  self%ws_cell(1, counter) = mod(m1 + 1, self%scell(1))
  if(self%ws_cell(1, counter) <= 0) then
     self%ws_cell(1, counter) = self%ws_cell(1, counter) + self%scell(1)
  end if
  
  self%ws_cell(2, counter) = mod(m2 + 1,self%scell(2))
  if(self%ws_cell(2, counter) <= 0) then
     self%ws_cell(2, counter) = self%ws_cell(2, counter) + self%scell(2)
  end if

  self%ws_cell(3, counter) = mod(m3 + 1,self%scell(3))
  if(self%ws_cell(3, counter) <= 0) then
     self%ws_cell(3, counter) = self%ws_cell(3, counter) + self%scell(3)
  end if
end if
~~~
   
## `phonon_espresso`

This subroutine is called by the `calculate_phonons()`-subroutine which 
calls `phonon_espresso()` *after* `phonon_espresso_precompute`.

### Variables

- **`nq`** *`integer`*
    - probably the number of qpoints
    - this is passed as an argument of `phonon_espresso`
- **`qpoints(nq,3)`** *`vector of floats`*
    - probably the for which the dynamical matrix is constructed
- **`omegas(nq, self%numbands)`** *`matrix of floats`*
    - frequencies calculated from the dynamical matrix
- **`velocities(nq, self%numbands, 3)`** *`array of floats`*
    - phonon group velocities calculated from the eigenvectors (?)
- **`eigenvect(nq, self%numbands, self%numbands)`** *`array of complex`*
    - eigenvectors obtained from diagonalization of dynamical matrix
    - this is an argument of `phonon_espresso` because the container has to 
      be passed into which the solutions are wrote
- **`toTHz=20670.687_r64`** *`float`*
    - NOT the conversion from energy in `Ry=13.605 eV`
    - the ifc2-file of quantum espresso was given in a Rydberg energy scale 
      and so are the eigenvectors computed from the formed dynamical matrix 
      (probably...). When the velocities are formed `toTHz` is used to 
      convert them into `km/s` but I dont know how yet.
- **`counter`** *`integer`*
    - muxes the specific `r_ws` into a number
- **`ntype=crys%numelements`** *`integer`*
    - number of species
- **`nat`** *`integer`*
    - number of atoms in a unitcell
- **`nbranches=3*nat`** *`integer`*
    - number of solutions to the equation of motion for one q point
- **`i, j`** *`integer`*
    - DESCRIPTION
- **`ipol, jpol`** *`integer`*
    - dummies that go through the atom displacement cartesian index of the 
      `self%ifc2`.
- **`iat, jat`** *`integer`*
    - just as usual the dummies that go from 1 to every atom number until 
      `nat`
- **`idim, jdim`** *`integer`*
    - encodes the type of atom **into** the structure of the dynamical 
      matrix. Usually in theory the dynamical matrix is an 
      `nat`x`nat`-matrix for every q-point. A field of matrices. But here 
      it is a big `3*nat`x`3*nat`-matrix, which is also apparent from the 
      definition of `ndim` and the allocation of `dyn_s`
- **`t1, t2, t3`** *`integer`*
    - components of the newly mapped `r_ws` vector at a specific `counter`.
- **`m1, m2, m3`** *`integer`* usual dummies for running through all unitcell positions
- **`iq`** *`integer`*
    - dummy for running from 1 to the number of qpoints `nq`
- **`ndim = 3*nat`** *`integer`*
    - the dynamical matrix will be of shape `ndim`x`ndim`
- **`nwork`** *`integer`*
    - something technical related to the `zheev()` diagonalization...
- **`ncell_g(3)`** *`vector of integers`*
    - related to the nonanalytic correction of the dynamical matrix
- **`total_weight`** *`float`*
    - it is not used somehow. Just set to `0`
- **`weight`** *`float`*
    - copied from `self%ws_weight(counter)` which is the weight assigned to 
      a particular `r_ws`
- **`alpha`** *`float`*
    - related to the nonanalytic correction of the dynamical matrix
- **`geg`** *`float`*
    - related to the nonanalytic correction of the dynamical matrix
- **`gmax`** *`float`*
    - related to the nonanalytic correction of the dynamical matrix
- **`qt`** *`real`*
    - temp for the dot-product of a qpoint with the vector saved by 
      `t1,t2,t3`
- **`t(0:3)`** *`vector of floats`*
    - cartesian unitcell positions
- **`omega2(:)`** *`vector of floats`*
    - phonon frequency squared. Is allocated with length of `nbranches`
- **`rwork(:)`** *`vector of floats`*
    - something technical related to the `zheev`-diagonalization
- **`q`** *`matrix of floats`*
    - matrix of qpoints in the cartesian basis of the reciprocal lattice 
      vectors
    - units are converted from 1/nm to 1/bohr by *multiplying* with 
      `bohr2nm`
- **`vels`** *`array of floats`*
    - don't know yet. Apparently passed also appears in 
      `calculate_phonons()`
- **`eigenvectors`** *`matrix of complex`*
    - probably eigenvectors from `zheev` (?)
- **`work`** *`array of complex`*
    - technical, don't know yet
- **`dyn(:,:)`** *`array of complex`*
    - muxed dynamical matrix with long and short range parts
- **`dyn_s(:,:,:)`** *`array of complex`*
    - short range part of `dyn(:,:)`
- **`dyn_l(:,:,:)`** *`array of complex`*
    - long range part of `dyn(:,:)`
- **`ddyn(:,:,:)`** *`array of complex`*
    - don't know yet
- **`ddyn_s(:,:,:,:)`** *`array of complex`*
    - don't know yet
- **`ddyn_l(:,:,:,:)`** *`array of complex`*
    - don't know yet

### Operations

I'm not going to list the allocations of the variables.

1. copy some necessary quantities

~~~fortran
nwork = 1
ntype = crys%numelements
nat = crys%numatoms
ndim = 3*nat
~~~

2. This is already hinted at in the definition for the matrix `q(:,:)` 
   above. The `qpoints` are transformed into cartesian coordinates using 
   the reciprocal lattice vectors `crys%reclattvecs` and converted to units 
   of 1/bohr from 1/nm.
>[!WARNING]
> `volume_r` is defined and set but never used lol

~~~fortran
do iq = 1, nq
   q(iq, :) = matmul(crys%reclattvecs, qpoints(iq, :))
end do
q = q*bohr2nm
volume_r = crys%volume/bohr2nm**3
~~~

3. Some technical things regarding the non-analytical correction are 
   happening

~~~fortran
gmax = 14.0_r64
alpha = (twopi*bohr2nm/dnrm2(3, crys%lattvecs(:,1), 1))**2
geg = gmax*4.0_r64*alpha
ncell_g = int(sqrt(geg)/self%cell_g(:, 0)) + 1
~~~

4. The usual loop over all atoms in all supercells in the ultra is done

~~~fortran
do iat = 1, nat
   do jat = 1, nat
      do m1 = -2*self%scell(1), 2*self%scell(1)
         do m2 = -2*self%scell(2), 2*self%scell(2)
            do m3 = -2*self%scell(3), 2*self%scell(3)
~~~





## read_ifc2

This subroutine reads in the 2nd order force constants from the Quantum 
Espresso format. The documentation of QEs PHonon module is 
[here](https://www.quantum-espresso.org/Doc/user_guide_PDF/ph_user_guide.pdf). 
Strangely this is labeled as a user guide but assumes how the phonon module 
of them works.


### Parameters and Variables

The subroutine itself only has two parameters which are the `inout` phonon 
derived type (or class) and an `in` variable of the crystal type.

Quite a number of local variable are declared which I will try to fill up 
with meaning while going through the script

#### 64 bit integers

- `qscell(3)`: qgrid used in the DFPT calculation
- `tipo(crys%numatoms)`: vector which holds the species type for each atom
- `t1`: temp `ifc2` tensor location index 5
- `t2`: temp `ifc2` tensor location index 6
- `t3`: temp `ifc2` tensor location index 7
- `i`: dummy index for counting
- `j`: dummy index for counting
- `iat`: temp `ifc2` tensor location index 3 
- `jat`: temp `ifc2` tensor location index 4
- `ibrav`: ???
- `ipol`: temp `ifc2` tensor location index 1
- `jpol`: temp `ifc2` tensor location index 2
- `m1`: iterator for the number of supercells positions in the ultracell
- `m2`: iterator for the number of supercells positions in the ultracell
- `m3`: iterator for the number of supercells positions in the ultracell
- `ntype`: number of atomic species
- `nat`: number of atoms in the basis
- `nfc2`: number of force constants `3*3*nat*nat`

#### 64 bit floating

- `r(crys%numatoms, 3)`: relative coordinates for each atom in basis
- `wscell(3,0:3)`: stretched unitcell lattice vectors. Supercell lattice 
  vectors
- `celldm(6)`: vector of cell dimensions (why 6d?)
- `at(3,3)`: ???
- `mass(crys%numelements)`: mass for each atomic species
- `zeff(crys%numatoms, 3, 3)`: born effective charges
- `eps(3, 3)`: dielectric tensor
- `dnrm2`: LAPACK for euclidean norm of a vector
- `massfactor`: species mass conversion factor

#### Characters

- `polar_key`: polar material or not? probably eiter 'T' or 'F'
- `label(crys%numelements)`: label of atomic species

### Routine analysis

Lets go over what the routine does and how it uses the parameters and 
variables.

1. [Allocate 
   arrays](https://fortran-lang.org/learn/best_practices/allocatable_arrays/) 
   `mm` and `rr`. These arrays are declared in the derived type of the 
   `phonon` class. By best practice these *should* be scratch arrays. But I 
   am [reading 
   now](https://fortran-lang.org/learn/quickstart/arrays_strings/) that 
   they are just dynamic arrays.

~~~fortran
allocate(self%mm(crys%numatoms, crys%numatoms))
allocate(self%rr(crys%numatoms, crys%numatoms, 3))
~~~

2. [List-directed 
   I/O](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc5/index.html) 
   is a free-form I/O and it follows the rules given in source on the link. 
   It reads the first row of the ifc2-file with the number of species 
   `ntype`, number of atoms in the basis `nat`, `ibrav` ?? and 
   `celldm`-array

~~~fortran
read(1,*) ntype, nat, ibrav, celldm(1:6)
~~~

3. Special case with short-hand [implied 
   do](https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/io.html). I 
   can only guess that `ibrav` controls the specific lattice or something 
   similar. The implied do is used to read all 9 components of the 
   `at`-array.

~~~fortran
if (ibrav==0) then
   read(1,*) ((at(i,j),i=1,3),j=1,3)
end if
~~~

>[!CAUTION]
> What is `ibrav` and what denotes the special case for `ibrav=0`?

4. Read all species type lines which contain the species `label` and 
   species `mass`

~~~fortran
    do i = 1, ntype
       read(1, *) j, label(i), mass(i)
    end do
    mass = crys%masses/massfactor
~~~

5. For the number of atoms in the basis read `tipo` ?? and the relative 
   coordinates to each other

~~~fortran
do i = 1, nat
   read(1, *) j, tipo(i), r(i, 1:3)
end do
r = transpose(matmul(crys%lattvecs, crys%basis))/bohr2nm
~~~

>[!CAUTION]
> Why does this step override the `r`-array anyways? Why read it in the 
> first place?

6. If it is a polar material read the $\varepsilon$-tensor and for each 
   atom in the basis read the Born-effective charges

~~~fortran
read(1, *) polar_key
if(polar_key == "T") then
   do i = 1, 3
      read(1, *) eps(i, 1:3)
   end do
   do i = 1, nat
      read(1, *)
      do j = 1, 3
         read(1, *) zeff(i, j, 1:3)
      end do
   end do
end if
~~~

7. Which q-grid was used and copy it to `phonon%scell`

~~~fortran
read(1,*) qscell(1:3)
self%scell = qscell
~~~

>[!CAUTION]
> Why copy `qscell` to `scell`? **ANSWER: `scell` is the number of unit 
> cells, the supercell is made up of! `scell` is used later at step 12**

8. Allocate the `ifc2` tensor of rank 7 and store the number of force 
   constants

~~~fortran
allocate(self%ifc2(3, 3, nat, nat, self%scell(1), self%scell(2), 
self%scell(3)))
nfc2 = 3*3*nat*nat
~~~

9. Read the actual force constants by adresses. Every lines' last value 
   represents one element in the tensor. The rest is just adresses. After 
   this the actual reading of data is completed.
   The tensor has rank 7. 

$$U = U_0 + \frac{1}{2} \sum_{n,n',\alpha,\alpha'} \frac{\partial^2 
U}{\partial \xi^\alpha_n \,\partial\xi^{\alpha'}_{n'}}
\xi^{\alpha}_{n}\xi^{\alpha'}_{n'} + \dots$$


    

~~~fortran
do i = 1, nfc2
   read(1, *) ipol, jpol, iat, jat
   do j = 1, self%scell(1)*self%scell(2)*self%scell(3)
      read(1, *) t1, t2, t3, &
           self%ifc2(ipol, jpol, iat, jat, t1, t2, t3)
   end do
end do
~~~

10. Forcing the conservation of momentum into the `ifc2`-tensor.

~~~fortran
do i = 1, 3
   do j = 1, 3
      do iat = 1, nat
         self%ifc2(i, j, iat, iat, 1, 1, 1) = self%ifc2(i, j, iat, iat, 1, 1, 1) - &
              sum(self%ifc2(i, j, iat, :, :, :, :))
      end do
   end do
end do
~~~


>[!CAUTION]
> The ordering of this loop goes against the column-major nature of Fortran 
> as far as I can tell. It is probably not that import in this part because 
> the subroutine should just read in the file. But still weird.


11. Grab reciprocal and direct lattice vectors, convert them to nanometers 
    and save them into arrays. The lattice vectors are stored in columns. 
    Rows goes through the 3 dimensions. The 0-column denotes the length of 
    the lattice vectors

~~~fortran
self%cell_r(:, 1:3) = transpose(crys%lattvecs)/bohr2nm
do i = 1, 3
   self%cell_r(i, 0) = dnrm2(3, self%cell_r(i, 1:3), 1)
end do
self%cell_g(:, 1:3) = transpose(crys%reclattvecs)*bohr2nm
do i = 1, 3
   self%cell_g(i, 0) = dnrm2(3, self%cell_g(i, 1:3), 1)
end do
~~~

12. stretch the unitcell lattice vector by the number of unitcells in each 
    direction which make up the supercell (`scell` contains the number of 
    unit cells in each of the lattice vector directions)

~~~fortran
wscell(1,1:3) = self%cell_r(1,1:3)*self%scell(1)
wscell(2,1:3) = self%cell_r(2,1:3)*self%scell(2)
wscell(3,1:3) = self%cell_r(3,1:3)*self%scell(3)
~~~

13. The supercell is made up of unitcells. Now we build data corresponding 
    to a super-supercell that is made up of 5 supercells in each lattice 
    vector (stretched and contained in `wscell`-array) direction. We are 
    storing all vectors that point to a specific supercell in the 
    super-supercell (columns `1:3` of `self%rws`) AND *half of their 
    squared lengths* (column `0`).

~~~fortran
j = 1
do m1 = -2, 2
   do m2 = -2, 2
      do m3 = -2, 2
         if(all((/m1, m2, m3/).eq.0)) then
            cycle
         end if
         do i = 1, 3
            self%rws(j, i) = wscell(1, i)*m1 + wscell(2, i)*m2 + wscell(3, i)*m3
         end do
         self%rws(j, 0) = 0.5*dot_product(self%rws(j, 1:3), self%rws(j, 1:3))
         j = j + 1
      end do
   end do
end do
~~~

14. We are building symmetric matrices. They hold **a)** Masses and Mass 
    factors (`self%mm`). Diagonal is the mass of atom `i`. All other 
    elements are the mass factors of atoms `i` and `j`. **b)** Distances 
    between atom `i` and `j` (`self%rr`). Diagonal is `0`.

~~~fortran
do i = 1, nat
   self%mm(i, i) = mass(tipo(i))
   self%rr(i, i, :) = 0
   do j = i + 1, nat
      self%mm(i, j) = sqrt(mass(tipo(i))*mass(tipo(j)))
      self%rr(i, j, 1:3) = r(i, 1:3) - r(j, 1:3)
      self%mm(j, i) = self%mm(i, j)
      self%rr(j, i, 1:3) = -self%rr(i, j, 1:3)
   end do
end do
~~~

## read_ifc3

## deallocate_phonon_quantities

## save_xmassvar

## allocate_xmassvar

## clean_xmassvar

