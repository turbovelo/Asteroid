"""
Author: Konrad Barboutie
Current date: 2025-11-10
"""
import requests
import re
import plotly.io as pio
pio.renderers.default = "browser"
import webbrowser
import os
#import matplotlib.pyplot as plt
#from pprint import pprint
#import json
import sys
import numpy as np
#from jplephem.spk import SPK
from jplephem.calendar import compute_julian_date
from poliastro.bodies import Earth, Sun
#from poliastro.twobody import Orbit
from poliastro.twobody.orbit import Orbit
#from poliastro.plotting import OrbitPlotter3D
from astropy import units as u
#from poliastro.core.elements import coe2rv
from astropy.time import Time
from poliastro.plotting.interactive import OrbitPlotter3D
from astropy.coordinates import get_body_barycentric_posvel


print("hello world")

# Initial date setting
current_date = compute_julian_date(2015, 2, 8)
launch_date = compute_julian_date(2030, 6, 1)
ephemeris_start = compute_julian_date(2030,6,1)
ephemeris_stop = compute_julian_date(2030,6,2)

# Request the ephemeris from Horizon API

url_2008_ev5 = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='DES=2008 EV5'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_earth = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='399'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_venus = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='299'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_mars = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='499'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_mercury = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='199'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
def get_elements(url):
    # Acquiring and Storing the data as a json file
    response = requests.get(url)
    data = response.json()

    # Verify connection was smooth
    if response.status_code != 200:
        print("ERROR GETTING THE DATA")
        sys.exit()
    else:
        print("Status Code: ", response.status_code) # 200 means everything ok

    # Extracting orbital elements code by Gemini AI, improved by me
    start_marker = '$$SOE'
    end_marker = '$$EOE'
    start_index = data['result'].find(start_marker)
    end_index = data['result'].find(end_marker)
    if start_index == -1 or end_index == -1:
        print("Error: Could not find the $$SOE / $$EOE markers in the response.")
        exit()
    # Extract the block (start_index + length of '$$SOE' marker)
    ephemeris_block = data['result'][start_index + len(start_marker): end_index].strip()

    # Focus on the data for the first date (June 1, 2030)
    # This is the block starting after $$SOE and ending before the second date
    first_date_block = ephemeris_block.split('2462654.500000000')[0].strip()

    # Extract the key-value pairs using regular expressions
    # This regex looks for: two capital letters (the element code, e.g., EC, QR),
    # followed by an '=' sign and then the floating-point number.
    elements_dictionnary = {}
    matches = re.findall(r'(\s[A-Z]{1,2})\s*=\s*([0-9.\-E+]+)', first_date_block)

    for key, value in matches:
        # Clean the key (remove leading/trailing spaces) and store as a float
        elements_dictionnary[key.strip()] = float(value)
    #   EC      Eccentricity
    #   QR      Periapsis distance
    #   IN      Inclination w.r.t. xy-plane (degrees)
    #   OM      Longitude of Ascending Node (degrees)
    #   W       Argument of Perifocus (degrees)
    #   Tp      Periapsis time (user specifies absolute or relative date), not available for 2008 EV5 so not included in the key list
    #   N       Mean motion (degrees/DU)
    #   MA      Mean anomaly (degrees)
    #   TA      True anomaly (degrees)
    #   A       Semi-major axis
    #   AD      Apoapsis distance
    #   PER     Orbital Period or PR ?

    KEY_ORDER = ['A', 'EC', 'IN', 'OM', 'W', 'TA', 'QR', 'N', 'MA', 'AD', 'PR'] # a, e, i, RAAN, argp, nu
    elements_sorted_values = [elements_dictionnary[key] for key in KEY_ORDER]
    # Storing the values in a numpy array for future use
    elements = np.array(elements_sorted_values)
    return elements

########################################################################################################################
# Poliastro Orbit calculations
########################################################################################################################
mu_sun = 1.32712440042e11
def radius_from_elements_pqw(a, e, nu):
    r = a*(1-pow(e,2))/(1+e*np.cos(nu))
    p = r*np.cos(nu)
    q = r*np.sin(nu)
    w = 0
    return [p, q, w]
def velocity_from_elements_pqw(a, e, nu, mu):
    semi_latus_rectum = a*(1-pow(e,2))
    r = np.sqrt(mu/semi_latus_rectum)
    p = -r*np.sin(nu)
    q = r*(e+np.cos(nu))
    w = 0
    return[p, q, w]
def rv_calculation(elements):
    r = radius_from_elements_pqw(elements[0], elements[1], elements[6]) * u.km
    v = velocity_from_elements_pqw(elements[0], elements[1], elements[6], mu_sun) * u.km/u.s
    return r, v

elements_mercury = get_elements(url_mercury)
elements_venus = get_elements(url_venus)
elements_earth = get_elements(url_earth)
elements_mars = get_elements(url_mars)
elements_asteroid = get_elements(url_2008_ev5)


epoch = Time(launch_date, format='jd', scale ='tdb')
#r_earth, v_earth = get_body_barycentric_posvel(Earth.name, epoch)
#earth_orbit = Orbit.from_vectors(Sun, r_earth.xyz, v_earth.xyz, epoch=epoch)

r_mercury, v_mercury = rv_calculation(elements_mercury)
r_venus, v_venus = rv_calculation(elements_venus)
r_earth, v_earth = rv_calculation(elements_earth)
r_mars, v_mars = rv_calculation(elements_mars)
r_asteroid, v_asteroid = rv_calculation(elements_asteroid)


#orb_mercury = Orbit.from_vectors(Sun, r_mercury, v_mercury)
#orb_venus = Orbit.from_vectors(Sun, r_venus, v_venus)
#orb_earth = Orbit.from_vectors(Sun, r_earth, v_earth)
#orb_mars = Orbit.from_vectors(Sun, r_mars, v_mars)
#orb_asteroid = Orbit.from_vectors(Sun, r_asteroid, v_asteroid)


def define_orbit(origin, classical_elements):
    orb = Orbit.from_classical(origin, classical_elements[0] * u.km, classical_elements[1] * u.one, classical_elements[2] * u.deg, classical_elements[3] * u.deg, classical_elements[4] * u.deg, classical_elements[5] *u.deg)
    return orb

orb_mercury = define_orbit(Sun, elements_mercury)
orb_venus = define_orbit(Sun, elements_venus)
orb_earth = define_orbit(Sun, elements_earth)
orb_mars = define_orbit(Sun, elements_mars)
orb_asteroid = define_orbit(Sun, elements_asteroid)


#orb = Orbit.from_vectors(Sun, r * u.km, v * u.km / u.s)
plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)
plotter.plot(orb_mercury, label = 'Mercury', color = '#B1B1B1')
plotter.plot(orb_venus, label = 'Venus', color = '#EEDC82')
plotter.plot(orb_earth, label = 'Earth', color = '#2E8B57')
plotter.plot(orb_mars, label = 'Mars', color = '#C1440E')
plotter.plot(orb_asteroid, label = '2008 EV5', color = 'purple')
fig = plotter._figure
fig.write_html("orbit.html")
webbrowser.open('file://' + os.path.realpath("orbit.html"))