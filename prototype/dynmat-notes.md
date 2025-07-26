# dynmat-notes

These notes capture what's happening inside the `dynmat.py`-script in this 
directory.

## readqeifc2-function

Pretty sure this is just reading a qe-dfpt-output file. Yet notable is that 
on top of reading the file, does unit conversion, reshaping and enforcement 
of the acustic sum rule.

- conversion of the force constants from ryd/bohr^2 to ev/nm^2

~~~python
ryd2ev = 13.605698066
#bohr2nm = 0.052917725
bohr2nm = 0.052917721092 #shengbte
fc *= ryd2ev/bohr2nm/bohr2nm
~~~

- reshaping converts from a muxed column-major design of the ifc2-tensor to 
  a row-major de-muxed one

~~~python
#Reshape ifc2 file from (m,lx,ly,lz) to (lx,ly,lz,j,a,jp,b) (standard) form
# lx,ly,lz are supcercell indices
fc_std = np.zeros((scell[0],scell[1],scell[2],natoms,3,natoms,3))
for m in xrange(numfc):
    a = ind[m][0] - 1 #pol
    b = ind[m][1] - 1 #pol
    j = ind[m][2] - 1 #atom
    jp = ind[m][3] - 1 #atom

    for lx in xrange(scell[0]):
        for ly in xrange(scell[1]):
            for lz in xrange(scell[2]):
                fc_std[lx][ly][lz][j][a][jp][b] = fc[m][lx][ly][lz]
~~~

- enforcing the sum rule is done similiarly as in `elphbolt`

~~~python
for a in xrange(3):
    for b in xrange(3):
        for j in xrange(natoms):
            fc_std[0][0][0][j][a][j][b] -= np.sum(fc_std[:,:,:,j,a,:,b])
~~~


## formdynmat-function

### Variables

- `cell`: contains all sorts of information about the crystal
    - `len(cell[1])` is the number of atoms
    - `a_i = cell[0][i]` probably lattice constants
    - `cell[0][iat]` position of atom `iat` 
- `numq`: number of q-points
- `natoms`: number of atoms
- `scell`: probably super-cell information container
    - `np.prod(scell)` is number of super-cells
    - `sc_i = scell[i]` either numbers or vectors

- *non-analytical* correction variables come here

- `dynmat`: complex matrix of size `(numq, natoms*3, natoms*3)`
- `a1`, `a2`, `a3`: lattice constants (?)
- `sc1`, `sc2`, `sc3`: unit-cell repetetions that make up the supercell
- `sclvec`: super-cell lattice vectors from 
  `np.multiply(system.scell,system.lvec)`
    - columns are the cartesian components and rows label the different 
      lattice vectors
- `rssc`: `(5*5*5, 4)`-array for the 5x5x5 super-supercell position vectors
- `cnt`: counting variable (muxing)

- `r_ssc`: (great naming there) 3-vector
- `t`: 3-vector
- `eps`: small number
- `totwgth`: weight-factor
- `massfac`: mass factor

### Build Vector To Supercells in the Super-Supercell

~~~python
cnt = 0
for i,j,k in it.product(xrange(-2,3),xrange(-2,3),xrange(-2,3)):
    if i == 0 and j == 0 and k == 0:
        continue
    for a in xrange(3):
        rssc[cnt][a] = i*sclvec[0][a] + j*sclvec[1][a] + k*sclvec[2][a]
    rssc[cnt][3] = 0.5*np.dot(rssc[cnt][0:3],rssc[cnt][0:3])
    cnt += 1
~~~

1. `it.product(xrange(-2,3),xrange(-2,3),xrange(-2,3))` builds an iterable 
   object of 3-tuples that contain all combinations of numbers in 
   `(-2,-1,0,1,2)`

2. loop skips the tuple `(0,0,0)`

3. fill the `rssc`-arrays' first 3 columns with cartesian components. Last 
   column is half of the squared length of the aforementioned cartesian 
   components

4. `rssc` is filled by the muxed tuples from step 1 above (by numbers 0 to 
   124)

### Populate The Dynamical Matrix

1. iterate `(iat, jat)` through the all 2-tuples for combinations of 
   numbers from 0 to `natoms`

~~~python
for iat, jat in it.product(xrange(natoms),xrange(natoms)):
~~~

2. calculate mass factor

~~~python
massfac = 1.0/np.sqrt(
    system.masses[system.types[iat]-1]*system.masses[system.types[jat]-1])
~~~

3. iterate `(m1, m2, m3)` through 3-tuples. The 3-tuples are all unitcell 
   adresses in the super-supercell, because the `sc1`, `sc2`, `sc3` are the 
   numbers of unitcells in each lattice vector (?) direction that make up 
   the supercell

~~~python
for m1,m2,m3 in it.product(
    xrange(-2*sc1,2*sc1+1),
    xrange(-2*sc2,2*sc2+1),
    xrange(-2*sc3,2*sc3+1)
    ):
~~~

4. convert the iteration integers `m1`, `m2`, `m3` into crystal coordinates

~~~python
t = np.dot([m1,m2,m3],system.lvec)
~~~

5. translate the vector that points from atom j to atom i in the unit cell 
   into 









