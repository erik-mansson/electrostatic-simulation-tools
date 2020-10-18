% Defines a set of physicsl constants in SI units (and some other units)
% in the current workspace, and most of them also as globals.
%
% SEE ALSO
% dimension_analysis for conversion to/from atomic units

global c_0 permittivity_0 mu_0 impedance_0 q hbar h_Planck electron_mass atomic_mass_unit k_Boltzmann fs eV N_Avogrado kcal atomic_momentum_unit atomic_time_unit Hartree elementary_charge Rydberg_inf Rydberg_inf_hc_eV
c_0 = 299792458; %[m/s] speed of light
permittivity_0 = 8.854187817e-12; %[A^2 s^4 kg^-1 m^-3]=[C^2 N^-1 m^-2] electric constant = permittivity of vacuum (epsilon_0)
mu_0 = 4*pi*1e-7;
impedance_0 = mu_0 * c_0; %[Ohm?] (eta_0)
% elementary_charge = 1.602176487e-19; %[C] elementary charge, previous value, CODATA2006
elementary_charge = 1.602176565e-19; %[C] elementary charge, NIST CODATA2010 retrieved 2013-08-07
q = elementary_charge; %synonym for elementary charge
%hbar = 1.054571628e-34; %[J s] or [J s / rad] Planck's constant over 2 pi, previous value, CODATA2006
hbar = 1.054571726e-34; %[J s] or [J s / rad] Planck's constant over 2 pi, NIST CODATA2010 retrieved 2013-08-07
%h_Planck = 6.62606896e-34; %[J s] Planck's constant, previous value, CODATA2006
h_Planck = 6.62606957e-34; %[J s] Planck's constant, NIST CODATA2010 retrieved 2013-08-07
%electron_mass = 9.10938215e-31; %[kg] previous value, CODATA2006
electron_mass = 9.10938291e-31; %[kg] NIST CODATA2010 retrieved 2013-08-07
%atomic_mass_unit = 1.660538782e-27;%[kg] previous value, CODATA2006
atomic_mass_unit = 1.660538921e-27;%[kg] CODATA2010 retrieved 2013-08-07
%k_Boltzmann = 1.3806504e-23; %[J K^-1] Boltzmann constant, previous value, CODATA2006
k_Boltzmann = 1.3806488e-23; %[J K^-1] Boltzmann constant, CODATA2010 retrieved 2013-08-07
%atomic_momentum_unit = 1.992851565E-24; % Atomic unit of momentum [kg m/s] = hbar/a_0 = m_e*c_0*alpha (CODATA 2006)
atomic_momentum_unit = 1.992851740E-24; % Atomic unit of momentum [kg m/s] = hbar/a_0 = m_e*c_0*alpha (CODATA 2010, retrieved 2013-08-07)
Rydberg_inf = 10973731.568539; %[m^-1] Rydberg constant CODATA2010
Rydberg_inf_hc_eV = 13.60569253; %[eV] Rydberg_inf*h_Planck*c_0/eV, from CODATA to avoid combining rounding errors
fs = 1e-15; %[s] a femtosecond
%eV = 1.602176487e-19; %[J] an electron volt (1eV / 1J = 1 q / 1 C), previous value, CODATA2006
eV = 1.602176565e-19; %[J] an electron volt (1eV / 1J = 1 q / 1 C), NIST CODATA2010 retrieved 2013-08-07
%N_Avogrado = 6.02214179e23; % [1/mol] Avogadro constant, previous value, CODATA2006
N_Avogrado = 6.02214129e23; % [1/mol] Avogadro constant, NIST CODATA2010 retrieved 2013-08-07
kcal = 4184; %[J] 1 kilocalorie = 1 "large Calorie"
T_0 = 273.15; %[degC] the absolute zero temperature in degrees Celcius, or the freezing point of water
alpha_finestructure = 7.2973525698E-3; %[1] fine structure constant = q^2/(4*pi*permittivity_0*hbar*c_0), CODATA2010
a_0 =  5.2917721092E-11; %[m] Bohr radius = 4*pi*permittivity_0*(hbar/q)^2/electron_mass, CODATA2010
Hartree = 4.35974434e-18; %[J] corresponding to 1 Hartree = 27.21138505(60) eV [NIST CODATA2010], a.k.a. atomic_energy_unit
atomic_time_unit = 2.418884326502E-17; %[s] = hbar/Hartree
mbar = 100; % [Pa] 1 mbar = 100 Pa

% IMPROVEMENT: allow some argument or global variable to switch to older
% CODATA2006 version, to maintain/check old behaviour even after upgrade of
% the constants.

% Conversion factors
global inverse_cm cm_minus_1 kcal_per_mol
inverse_cm = 2*pi*hbar*c_0*100; % [J] corresponding to 1 cm^-1 wavenumber
cm_minus_1 = inverse_cm;
kcal_per_mol = kcal/N_Avogrado; %[J] corresponding to 1 kcal/mol = 4184 J/mol

hf_800nm_in_eV = h_Planck*c_0/800e-9/eV; % [eV] energy of 800 nm photon (circa 1.55 eV)
Debye_from_SI = 1E-21 / c_0; %= 3.33564E-30; %[Coulomb m]
Debye_from_CGS = 1E-21 / c_0 / sqrt(4*pi*permittivity_0); %= 3.33564E-30; %[Coulomb m]
