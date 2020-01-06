import h5py, glob, numpy
from pyglib.estructure.fermi import get_fermi_level
from pyglib.iface.ifwannier import get_risb_bndes


with h5py.File("GPARAMBANDS.h5", "r") as f:
    num_elec = f['/nelectron'][0]
    ismear = f["/ismear"][0]
    delta = f["/delta"][0]
print("num of electrons = {}".format(num_elec))

bnd_es = get_risb_bndes()
nk = bnd_es.shape[1]
wklist = [1./nk for k in range(nk)]

print(bnd_es[0][0])
with h5py.File("gband_es.h5", "w") as f:
    f["/gbnd_es"] = bnd_es

efermi = get_fermi_level(bnd_es, wklist, num_elec,
        delta=delta, ismear=ismear)
print("lda fermi level = {}".format(efermi))
