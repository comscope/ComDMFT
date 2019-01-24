import h5py

global f
f = None

def open(fname):
    global f
    f = h5py.File(fname, 'w')


def close():
    f.close()
