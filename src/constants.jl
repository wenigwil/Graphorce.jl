# Bohr Radius (nm)
const a0_nm = 5.29177210544e-2

const ev_joule = 1.602176634e-19

const hartree_ev = 27.211386245981
const rydberg_ev = 0.5 * hartree_ev
const rydberg_joule = rydberg_ev * ev_joule

# hbar in J/THz
const hbar_Thz = 1.05457172647e-22

# angular frequency in rydberg to ang. THz
const RydtoTHz = (rydberg_joule / hbar_Thz) / (2 * pi)

# Mass of an electron in kg
const m_e = 9.10938291e-31

# Dalton (Unified atomic mass unit) in kg
const m_u = 1.660538921e-27
