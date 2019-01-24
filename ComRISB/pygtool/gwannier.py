#!/usr/bin/env python
from pyglib.iface.ifwannier import if_gwannier
import pickle


# get if_gwannier parameters.
with open("gwannier_params.pkl", "rb") as f:
    params = pickle.load(f)

assert "corbs_list" in params, "missing corbs_list in gwannier_params.pkl!"
params["delta_charge"] = params.get("delta_charge", 0.)
params["wpath"] = params.get("wpath", "../wannier")
params["lpath"] = params.get("lpath", "../lattice")
params["wprefix"] = params.get("wprefix", "wannier")
params["lprefix"] = params.get("lprefix", "mdl")
params["lrot_list"] = params.get("lrot_list", None)
params["iso"] = params.get("iso", 1)
params["ispin"] = params.get("ispin", 1)
params["ismear"] = params.get("ismear", 0)
params["delta"] = params.get("delta", 0.005)
params["icycle"] = params.get("icycle", 0)

if_gwannier(corbs_list=params["corbs_list"],
        delta_charge=params["delta_charge"], wpath=params["wpath"],
        lpath=params["lpath"],
        wprefix=params["wprefix"], lprefix=params["lprefix"],
        lrot_list=params["lrot_list"], iso=params["iso"],
        ispin=params["ispin"], ismear=params["ismear"],
        delta=params["delta"], icycle=params["icycle"])
