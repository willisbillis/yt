from __future__ import division
import yt
from yt import YTQuantity, ValidateDataField
from yt.units import dimensions
import numpy as np

################################################################################
def timer(t):
    # outputs time in more readable units
    if t <= 60:
        return "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        return "Elapsed time: %.1f min" % t
################################################################################
    # Create new derived fields

def kmsvel_z(field,data):
    kmsvel_z = data["velocity_z"].in_units("km/s")
    return kmsvel_z
yt.add_field(("gas","velocity_Z"),
              function = kmsvel_z,
              units = "km/s",
              validators=[ValidateDataField(["velocity_z"])])

def kmsvel_x(field,data):
    kmsvel_x = data["velocity_x"].in_units("km/s")
    return kmsvel_x
yt.add_field(("gas","velocity_X"),
              function = kmsvel_x,
              units = "km/s",
              validators=[ValidateDataField(["velocity_x"])])

def partcntH(field,data):
    molmass = YTQuantity(1.00794, "g")
    partcnt = data['flash','h   ']*data['gas','cell_mass']/molmass*6.0221409e23
    return partcnt
yt.add_field(('gas','partcntH'),
              function = partcntH,
              units = "",
              validators=[ValidateDataField([['flash','h   '],['gas','cell_mass']])])

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
yt.add_field(('gas','partcntO'),
              function = partcntO,
              units = "",
              validators=[ValidateDataField([['flash','o   '],
                                            ['flash','o1  '],
                                            ['flash','o2  '],
                                            ['flash','o3  '],
                                            ['flash','o4  '],
                                            ['flash','o5  '],
                                            ['flash','o6  '],
                                            ['flash','o7  '],
                                            ['flash','o8  '],
                                            ['gas','cell_mass']])])

def met_O(field,data):
    met = data['gas','partcntO']/data['gas','partcntH']
    solar_met = 6.61e-4
    met = met/solar_met
    return met
yt.add_field(('gas','met_O'),
              function = met_O,
              units = "",
              validators=[ValidateDataField([['gas','partcntO'],['gas','partcntH']])])

def thermal_energy(field,data):
    k_boltzmann = YTQuantity(1.38064852e-16, "erg/K")
    return 3/2*(data['temperature'])*(k_boltzmann)
yt.add_field(('gas','thermal_energy'),
              function = thermal_energy,
              units="cm**2/s**2",
              validators=[ValidateDataField(['temperature'])])

def cloud_velz(field,data):
    bulk_vel=YTQuantity(140e5, 'cm/s')
    return data['flash', 'velz'].in_units('cm/s')-bulk_vel
yt.add_field(("gas", "cloud_velz"),
              units="cm/s",
              function=cloud_velz,
              validators=[ValidateDataField(['flash', 'velz'])])

def dyn_pressure(field,data):
    dyn_press=data['gas', 'kinetic_energy']/data['gas', 'cell_volume']
    return dyn_press
yt.add_field(("gas", "dyn_pressure"),
              units="auto",
              dimensions=dimensions.pressure,
              function=dyn_pressure,
              validators=[ValidateDataField([['gas', 'kinetic_energy'],['gas', 'cell_volume']])])

def pressure(field,data):
    R = YTQuantity(8.3144598e6, 'cm**3*Pa/K')
    pressure = (data['gas','partcntO']+data['gas','partcntH'])/6.0221409e23*R*data['gas','temperature'].in_units("K")/data['gas','cell_volume'].in_units("cm**3")
    return pressure.convert_to_units("Pa")
yt.add_field(("gas","pressure"),
              units = "Pa",
              function=pressure,
              validators=[ValidateDataField([['gas','partcntO'],['gas','partcntH'],['gas','temperature'],['gas','cell_volume']])])

def ram_pressure(field,data):
    ram_press=0.5*data['gas','density']*((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))
    return ram_press.convert_to_units("Pa")
yt.add_field(("gas", "ram_pressure"),
              units="Pa",
              function=ram_pressure,
              validators=[ValidateDataField([['gas','density'],['gas', 'velocity_x'],['gas', 'velocity_y'],['gas', 'velocity_z']])])

def mach_speed(field,data):
    gamma = 5.0/3.0
    sound_speed = (gamma*data['gas','pressure']/data['gas','density'].in_units('kg/m**3'))**0.5
    mach_speed = (((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))**(0.5)).in_units("m/s")/sound_speed
    return mach_speed
yt.add_field(('gas','mach_speed'),
              units="",
              function=mach_speed,
              validators=[ValidateDataField([['gas','pressure'],
                                            ['gas','density'],
                                            ['gas', 'velocity_x'],
                                            ['gas', 'velocity_y'],
                                            ['gas', 'velocity_z']])])

def specific_KE(field,data):
    ske = 0.5*data["velocity_z"]**2
    return ske.in_units("J/g")
yt.add_field(('gas','specific_KE'),
              units="J/g",
              function=specific_KE,
              validators=[ValidateDataField(["velocity_z"])])
