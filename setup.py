import yt
from yt import YTQuantity
from yt.units import dimensions

################################################################################
    # insert into script to print full callstack (for debugging)
def trace():
    for line in traceback.format_stack():
        print line.strip()

def time_fix(t):
    # outputs time in more readable units
    if t <= 60:
        print "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        print "Elapsed time: %.1f min" % t
    return
################################################################################
    # Create new derived fields

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

def cloud_velz(field,data):
    bulk_vel=YTQuantity(140e5, 'cm/s')
    return data['flash', 'velz'].in_units('cm/s')-bulk_vel
yt.add_field(("gas", "cloud_velz"), units="cm/s", function=cloud_velz)

def dyn_pressure(field,data):
    dyn_press=data['gas', 'kinetic_energy']/data['gas', 'cell_volume']
    return dyn_press
yt.add_field(("gas", "dyn_pressure"), units="auto", dimensions=dimensions.pressure, function=dyn_pressure)

def pressure(field,data):
    R = YTQuantity(8.3144598e6, 'cm**3*Pa/K')
    pressure = (data['gas','partcntO']+data['gas','partcntH'])/6.0221409e23*R*data['gas','temperature'].in_units("K")/data['gas','cell_volume'].in_units("cm**3")
    return pressure.convert_to_units("Pa")
yt.add_field(("gas","pressure"), units = "Pa", function=pressure, force_override=True)

def ram_pressure(field,data):
    ram_press=0.5*data['gas','density']*((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))
    return ram_press.convert_to_units("Pa")
yt.add_field(("gas", "ram_pressure"), units="Pa", function=ram_pressure, force_override=True)

def mach_speed(field,data):
    sound_speed = (data['gas','pressure']/data['gas','density'].in_units('kg/m**3'))**0.5
    mach_speed = (((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))**(0.5)).in_units("m/s")/sound_speed
    return mach_speed
yt.add_field(('gas','mach_speed'), units="", function=mach_speed)

