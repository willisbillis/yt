import yt
from yt import YTQuantity
from yt.units import dimensions

def partcntH(field,data):
    molmass = YTQuantity(1.00794, "g")
    partcnt = data['flash','h   ']*data['gas','cell_mass']/molmass*6.0221409e23
    return partcnt

yt.add_field(('gas','partcntH'), function = partcntH, units = "")

def partcntO(field,data):
    molmass = YTQuantity(15.9994, "g")
    partcnt = (data['flash','o   '] +
                data['flash','o1  '] +
                data['flash','o2  '] +
                data['flash','o3  '] +
                data['flash','o4  '] +
                data['flash','o5  '] +
                data['flash','o6  '] +
                data['flash','o7  '] +
                data['flash','o8  '] )*data['gas','cell_mass']/molmass*6.0221409e23
    return partcnt

yt.add_field(('gas','partcntO'), function = partcntO, units = "")

def met_O(field,data):
    met = data['gas','partcntO']/data['gas','partcntH']
    solar_met = 6.61e-4
    met = met/solar_met
    return met

yt.add_field(('gas','met_O'), function = met_O, units = "")

def thermal_energy(field,data):
    k_boltzmann = YTQuantity(1.38064852e-16, "erg/K")
    return 3/2*(data['temperature'])*(k_boltzmann)

yt.add_field(('gas','thermal_energy'), function = thermal_energy, units="cm**2/s**2", force_override=True)
