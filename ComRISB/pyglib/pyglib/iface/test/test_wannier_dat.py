from pyglib.iface.ifwannier import get_wannier_dat


reals_latt, recip_latt, kpts, include_bands, wfwannier_list, bnd_es = \
        get_wannier_dat()

print(reals_latt)
print(recip_latt)

print(reals_latt.dot(recip_latt.T))
