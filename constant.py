from pyteomics import mass

ion_cap_formula = {
    "a": "H-2O-1" + "C-1O-1",
    "b": "H-2O-1",
    "x": "H-2O-1" + "CO2",
    "y": "",
    "c-1": "H-2O-1" + "NH3" + "H-1",
    "c": "H-2O-1" + "NH3",
    "cdot": "H-2O-1" + "NH3" + "H1",
    "c+1": "H-2O-1" + "NH5",
    "zdot": "H-2O-1" + "N-1" + "O",
    "z+1": "H-2O-1" + "N-1" + "OH",
    "z+2": "H-2O-1" + "N-1" + "OH2",
    "z+3": "H-2O-1" + "N-1" + "OH3",
    "t": "",
}

ion_cap_delta_mass = {
    name: mass.calculate_mass(formula, absolute=False) for (name, formula) in ion_cap_formula.items()
}


ion_direction = {
    "a": "n-term",
    "b": "n-term",
    "x": "c-term",
    "y": "c-term",
    "cdot": "n-term",
    "c": "n-term",
    "c-1": "n-term",
    "c+1": "n-term",
    "zdot": "c-term",
    "z+1": "c-term",
    "z+2": "c-term",
    "z+3": "c-term",
}


# ---------------------------------------------------------------------------- #
#                               For visualization                              #
# ---------------------------------------------------------------------------- #


colors = [
    "#9b2226",
    "#005f73",
    "#ee9b00",
    "#0a9396",
    "#94d2bd",
    "#ca6702",
    "#e9d8a6",
    "#bb3e03",
    "#001219",
    "#006BA6",
    "#35A7FF",
    "#EFA8B8",
    "#BFACC8",
    "#476A6F",
    "#7067CF",
    "#364156",
    "#98FB98",
    "#8A2BE2",
    "#35682D",
    "#252850",
    "#7E7B52",
]

colors = colors + (["#808080"] * 1000)
