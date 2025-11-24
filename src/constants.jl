# Bohr Radius (nm)
const a0_nm = 5.29177210544e-2

const ev_joule = 1.602176634e-19

const hartree_ev = 27.211386245981
const rydberg_ev = 0.5 * hartree_ev
const rydberg_joule = rydberg_ev * ev_joule

# Planksches Wirkungsquantum in various units and versions
const h_Js = 6.62607015e-34
const hbar_J_over_THz = 1e12 * h_Js / (2 * pi)
const hbar_eV_over_THz = hbar_J_over_THz / ev_joule

# angular frequency in rydberg to ang. THz
const RydtoTHz = (rydberg_joule / hbar_J_over_THz) / (2 * pi)

# Mass of an electron in kg
const m_e = 9.10938291e-31

# Dalton (Unified atomic mass unit) in kg
const m_u = 1.660538921e-27
