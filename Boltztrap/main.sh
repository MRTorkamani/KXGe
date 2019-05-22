for name in K*
do
cd $name
for spin in dn up
do
rm -r $name$spin

mkdir $name$spin

sed '1,7d' $name.energy$spin*	>>	$name$spin/$name.energy
cat $name.scf			>>	$name$spin/$name.scf
cat $name.struct		>>	$name$spin/$name.struct

cd $name$spin
cat > $name.py << EOF
import sys
import copy
import os.path
import contextlib
import os
import pickle
import logging

import numpy as np
import matplotlib.pylab as pl

#from environment import data_dir

from BoltzTraP2 import dft as BTP
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib
from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import serialization
from BoltzTraP2 import units

logging.basicConfig(
    level=logging.DEBUG, format="{levelname:8s}│ {message:s}", style="{")
data_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), ".."))
dirname = os.path.join(data_dir, "$name$spin")
bt2file = "$name$spin.bt2"

# If a ready-made file with the interpolation results is available, use it
# Otherwise, create the file.
if not os.path.exists(bt2file):
    # Load the input
    data = BTP.DFTData(os.path.join(data_dir, dirname))
    # Select the interesting bands
    data.bandana(emin=data.fermi - .2, emax=data.fermi + .2)
    # Set up a k point grid with roughly five times the density of the input
    equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                           len(data.kpoints) * 5)
    # Perform the interpolation
    coeffs = fite.fitde3D(data, equivalences)
    # Save the result
    serialization.save_calculation(bt2file, data, equivalences, coeffs,
                                   serialization.gen_bt2_metadata(
                                       data, data.mommat is not None))

# Load the interpolation results
data, equivalences, coeffs, metadata = serialization.load_calculation(bt2file)

# Reconstruct the bands
lattvec = data.get_lattvec()
eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)

# Obtain the Fermi integrals for different chemical potentials at
# room temperature.
TEMP = np.array([300.])
epsilon, dos, vvdos, cdos = BL.BTPDOS(eband, vvband, npts=4000)
margin = 9. * units.BOLTZMANN * TEMP.max()
mur_indices = np.logical_and(epsilon > epsilon.min() + margin,
                             epsilon < epsilon.max() - margin)
mur = epsilon[mur_indices]
N, L0, L1, L2, Lm11 = BL.fermiintegrals(
    epsilon, dos, vvdos, mur=mur, Tr=TEMP, dosweight=data.dosweight)

# Compute the Onsager coefficients from those Fermi integrals
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, TEMP, UCvol)

fermi = BL.solve_for_mu(epsilon, dos, data.nelect, 300, data.dosweight)

# Plot the results
fig1, ax1 = pl.subplots(1, figsize=(6, 3))
ax1.set_xlim([-1, 1])
ax1.set_ylim([-300, 300])
ax1.plot((mur - fermi) / BL.eV, seebeck[0, :, 0, 0] * 1e6, "k-", label="xx")
ax1.plot((mur - fermi) / BL.eV, seebeck[0, :, 2, 2] * 1e6, "k--", label="zz")
ax1.set_xlabel("\$\mu\$ [eV]")
ax1.set_ylabel("S [\$\mu\$ V/K]")
ax1.legend()
fig1.tight_layout(pad=1.)
fig1.savefig("seebeck$name$spin.pdf")

fig2, ax2 = pl.subplots(1, figsize=(6, 3))
ax2.set_xlim([-1, 1])
ax2.set_ylim([0, 70])
ax2.plot(
    (mur - fermi) / BL.eV,
    seebeck[0, :, 0, 0]**2 * sigma[0, :, 0, 0] * 1e-10,
    "k-",
    label="xx")
ax2.plot(
    (mur - fermi) / BL.eV,
    seebeck[0, :, 2, 2]**2 * sigma[0, :, 2, 2] * 1e-10,
    "k--",
    label="zz")
ax2.set_xlabel("\$\mu\$ [eV]")
ax2.set_ylabel(r"S\$^2\$ \$\sigma/\$ \$\tau\$ [10\$^{14}\$ \$\mu\$W/(cm K\$^2\$) s]")
ax2.legend()
fig2.tight_layout(pad=1.)
fig2.savefig("PF$name$spin.pdf")

fig3, ax3 = pl.subplots(1, figsize=(6, 3))
ax3.set_xlim([-1, 1])
ax3.set_ylim([0, 700])
ax3.plot(
    (mur - fermi) / BL.eV,
    sigma[0, :, 0, 0] * 1e-18,
    "k-",
    label="xx")
ax3.plot(
    (mur - fermi) / BL.eV,
    sigma[0, :, 2, 2] * 1e-18,
    "k--",
    label="zz")
ax3.set_xlabel('\$\mu\$ [eV]')
ax3.set_ylabel(r'\$\sigma/ \tau\$ [\$(\Omega cm)^{-1}\$  s]')
ax3.legend()
fig3.tight_layout(pad=1.)
fig3.savefig("sigma")
######################################################


# Use the Fermi integrals to obtain the Onsager coefficients
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = bandlib.calc_Onsager_coefficients(
    L0, L1, L2, mur, TEMP, UCvol)

# Plot the results
fig1, ax1 = pl.subplots(1, figsize=(6, 4))
kappa = (kappa[0, :, 0, 0] + kappa[0, :, 1, 1] + kappa[0, :, 2, 2]) / 3.
kappaWF = (sigma[0, :, 0, 0] + sigma[0, :, 1, 1] + sigma[0, :, 2, 2]
           ) / 3. * 2.44E-8 * TEMP[0]
kappaD = (seebeck[0, :, 0, 0]**2 * sigma[0, :, 0, 0] +
          seebeck[0, :, 1, 1]**2 * sigma[0, :, 1, 1] +
          seebeck[0, :, 2, 2]**2 * sigma[0, :, 2, 2]) / 3. * TEMP[0]

ax1.plot((mur - fermi) / bandlib.eV, kappa * 1E-14, "k-", label='\$\kappa_e\$')
ax1.plot(
    (mur - fermi) / bandlib.eV, kappaWF * 1E-14, "k:", label="\$L\sigma T\$")
ax1.set_ylim([0, 60])
ax1.set_xlim([-1, 1])
ax1.set_xlabel('\$\mu\$ [eV]')
ax1.set_ylabel(
    r'\$\kappa/\tau \;\left(10^{14}\; \mathrm{W m}^{-1} \mathrm{K}^{-1}\right)\$'
)
ax1.legend()
fig1.savefig("kappa$name$spin.pdf")
#########################
#pl.show()

EOF

python3.7 $name.py
cd ..
done
cd ..
done
