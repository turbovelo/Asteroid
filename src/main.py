"""
Author: Konrad Barboutie
Current date: 2025-11-10
"""

from math import floor
import webbrowser, os, sys, re, requests, numpy as np, plotly.io as pio, json

pio.renderers.default = "browser"
from jplephem.calendar import compute_julian_date
import tkinter as tk
from poliastro.maneuver import Maneuver
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range
from poliastro.plotting import OrbitPlotter3D
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from poliastro.plotting.interactive import OrbitPlotter3D
from poliastro import iod
from poliastro.bodies import Earth, Mars, Sun
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.util import time_range
from poliastro.core.elements import coe_rotation_matrix, coe2rv
from poliastro.plotting._base import Trajectory
import plotly.io as pio

########################################################################################################################
# Function definitions
########################################################################################################################


def get_data(url, year):
    element = [0, 0, 0, 0, 0, 0]
    response = requests.get(url)
    raw_data = response.json()["result"]
    with open("data.txt", "w") as file:
        file.write(raw_data)
    file.close()
    with open("data.txt", "r") as f:
        lines = f.readlines()
        start_index = 0
        for i in range(len(lines)):
            if "$$SOE" in lines[i]:
                start_index = i
        if start_index == 0:
            print("Fatal Error $$SOE not found")
            return 0

        count = 0
        for j in range(start_index + 1, len(lines)):
            if count < 2:
                lines[j] = lines[j].replace("=", " ")
                if "EC" in lines[j]:
                    match = re.search(r"EC\s*([0-9.E+-]+)", lines[j])
                    element[1] = float(match.group(1))
                if "IN" in lines[j]:
                    match = re.search(r"IN\s*([0-9.E+-]+)", lines[j])
                    element[2] = float(match.group(1))

                if "OM" in lines[j]:
                    match = re.search(r"OM\s*([0-9.E+-]+)", lines[j])
                    element[4] = float(match.group(1))

                if "W" in lines[j]:
                    match = re.search(r"W\s*([0-9.E+-]+)", lines[j])
                    element[3] = float(match.group(1))

                if "TA" in lines[j]:
                    match = re.search(r"TA\s*([0-9.E+-]+)", lines[j])
                    element[5] = float(match.group(1))

                if " A " in lines[j]:
                    match = re.search(r"A\s*([0-9.E+-]+)", lines[j])
                    element[0] = float(match.group(1))

                if str(year) in lines[j]:
                    count += 1

                # print(f'j: {j}')
                # print(lines[j])
                # print(element)
            else:
                return element


def get_elements(url):
    # Acquiring and Storing the data as a json file
    response = requests.get(url)
    data = response.json()
    # Verify connection was smooth
    if response.status_code != 200:
        print("ERROR GETTING THE DATA")
        sys.exit()
    else:
        print(
            f"Status Code: {response.status_code}, successful"
        )  # 200 means everything ok
    start_marker = "$$SOE"
    end_marker = "$$EOE"
    pattern_block = r"\$\$SOE\n.*?\n(.*?)(?=\n\d+\.\d+\s=)"
    match_block = re.search(pattern_block, data, re.DOTALL)

    if match_block:
        raw_text = match_block.group(1)

        # 2. Parse the key-value pairs from that specific block
        # Matches things like "EC= 1.23E-02"
        pattern_values = r"([A-Za-z]+)\s*=\s*([-+]?\d*\.\d+(?:[Ee][-+]?\d+)?)"
        extracted_data = dict(re.findall(pattern_values, raw_text))

        print("Extracted Data:")
        for key, value in extracted_data.items():
            print(f"{key}: {value}")
    else:
        print("Pattern not found.")


def normalize(v):
    norm = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    v_norm = v / norm
    return v_norm


def length(v):
    v = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v


def node_line(x, y):
    return np.cross(x, y)


def eccentricity_vector(mu, r, v):
    rl = length(r)
    vl = length(v)
    e = 1 / mu * ((vl**2 - mu / rl) * r - (np.dot(r, v)) * v)
    return e


def radius_from_elements_pqw(a, e, nu):
    r = a * (1 - pow(e, 2)) / (1 + e * np.cos(nu))
    p = r * np.cos(nu)
    q = r * np.sin(nu)
    w = 0
    return [p, q, w]


def velocity_from_elements_pqw(a, e, nu, mu):
    semi_latus_rectum = a * (1 - pow(e, 2))
    r = np.sqrt(mu / semi_latus_rectum)
    p = -r * np.sin(nu)
    q = r * (e + np.cos(nu))
    w = 0
    return [p, q, w]


def rv_calculation(elements):
    r = radius_from_elements_pqw(elements[0], elements[1], elements[6]) * u.km
    v = (
        velocity_from_elements_pqw(elements[0], elements[1], elements[6], mu_sun)
        * u.km
        / u.s
    )
    return r, v


def define_orbit(origin, classical_elements, epoch):
    orb = Orbit.from_classical(
        origin,
        classical_elements[0] * u.km,
        classical_elements[1] * u.one,
        classical_elements[2] * u.deg,
        classical_elements[3] * u.deg,
        classical_elements[4] * u.deg,
        classical_elements[5] * u.deg,
        epoch,
    )
    return orb


def define_ephem(origin, classical_elements):
    orb = Ephem.from_orbit(
        origin,
        classical_elements[0] * u.km,
        classical_elements[1] * u.one,
        classical_elements[2] * u.deg,
        classical_elements[3] * u.deg,
        classical_elements[4] * u.deg,
        classical_elements[5] * u.deg,
    )
    return orb


def sphere_of_influence(m1, r):
    mu = 1.32712440042e11  # mu of sun
    return r * pow(m1 / mu, 2 / 5)


def angular_momentum_normalized(i, raan):
    h = [np.sin(i) * np.sin(raan), -np.sin(i) * np.cos(raan), np.cos(i)]
    return h


def nu_to_eccentric_anomaly(e, nu):
    E = 2 * np.arctan(np.tan(nu / 2) * np.sqrt((1 - e) / (1 + e)))
    return E


def time_to_node_line(a, e, w, nu, mu):
    # nu corresponds to current true anomaly
    E0 = nu_to_eccentric_anomaly(e, nu)  # Initial E
    E_asc = nu_to_eccentric_anomaly(e, 2 * np.pi - w)  # Final E
    E_des = nu_to_eccentric_anomaly(e, np.pi - w)
    dt_asc = np.sqrt(a**3 / mu) * ((E_asc - e * np.sin(E_asc)) - (E0 - e * np.sin(E0)))
    dt_des = np.sqrt(a**3 / mu) * ((E_des - e * np.sin(E_des)) - (E0 - e * np.sin(E0)))
    if dt_asc < 0:
        dt_asc += np.sqrt(a**3 / mu) * (
            2 * np.pi + (E_asc - e * np.sin(E_asc)) - (E0 - e * np.sin(E0))
        )
    if dt_des < 0:
        dt_des = np.sqrt(a**3 / mu) * (
            2 * np.pi + (E_des - e * np.sin(E_des)) - (E0 - e * np.sin(E0))
        )
    return dt_asc, dt_des


def julian_to_current(jd):
    j = jd + 0.5
    z = floor(j)
    alpha = floor((z - 1867216.25) / 36524.25)
    A = z + 1 + alpha - floor(alpha / 4)
    B = A + 1524
    C = floor((B - 122.1) / 365.25)
    D = floor(365.25 * C)
    E = floor((B - D) / 30.6001)
    day = B - D - floor(30.6001 * E)
    if E < 13:
        month = E - 1
    else:
        month = E - 13
    if month == 1 or month == 2:
        year = C - 4716
    else:
        year = C - 4715
    day_int = round(day)
    date_str = f"{year}-{month:02d}-{day_int:02d}"
    return date_str


def argument_perigee(n, e):
    nl = length(n)
    el = length(e)
    w = np.arccos(np.dot(n, e) / (el * nl))
    if e[2] < 0:
        w = w + np.pi
        return w
    else:
        return w


def angular_momentum(r, v):
    return np.cross(r, v)


########################################################################################################################
# Data Collection & Initialization
########################################################################################################################
# Write Horizon ephemeris data to file for faster execution

figure_name = "orbit"
# Initial date setting
start_year = 2037
start_month = 6
start_day = 24
end_year = start_year + 1
end_month = 6
end_day = 25

# 2038-02-04 arrival time
observation_time = 1
fast_forward = False  # See fast forward in time by latter specified amount of years

ti = Time(f"{start_year}-{start_month}-{start_day}", scale="utc")  # initial time
tf = Time(
    f"{start_year+observation_time}-{start_month}-{start_day}", scale="utc"
)  # final time
epochs = time_range(start=ti, end=tf)
# Url's to connect to jpl Horizon API and retrieve data
url_2008_ev5 = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='DES=2008 EV5'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start_year}-{start_month}-{start_day}'&STOP_TIME='{end_year}-{end_month}-{end_day}"
url_earth = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='399'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start_year}-{start_month}-{start_day}'&STOP_TIME='{end_year}-{end_month}-{end_day}"


########################################################################################################################
# Poliastro Orbit calculations
########################################################################################################################
mu_sun = 1.32712440042e11
mu_earth = 3.986e5

# Orbital Elements lists

elements_earth = get_data(url_earth, start_year)
elements_asteroid = get_data(url_2008_ev5, start_year)
print(elements_earth, "\n", elements_asteroid)
elements_earth_radians = elements_earth.copy()
elements_asteroid_radians = elements_asteroid.copy()
elements_earth_radians[2:6] = np.radians(elements_earth[2:6])
elements_asteroid_radians[2:6] = np.radians(elements_asteroid[2:6])


# Define orbits
earth_initial = define_orbit(Sun, elements_earth, ti)
asteroid_initial = define_orbit(Sun, elements_asteroid, ti)

# Define The orbit that will be seen in the plot after the elapsed ephem time chosen above
earth_ephem = Ephem.from_orbit(orbit=earth_initial, epochs=epochs)
asteroid_ephem = Ephem.from_orbit(orbit=asteroid_initial, epochs=epochs)


########################################################################################################################
# Maneuvering
########################################################################################################################
# Determine line of nodes between the two inclined orbits for Bi-Elliptic Transfer calculations
r_earth, v_earth = earth_initial.rv()
r_earth, v_earth = r_earth.value, v_earth.value
r_asteroid, v_asteroid = asteroid_initial.rv()
r_asteroid, v_asteroid = r_asteroid.value, v_asteroid.value
h_earth = angular_momentum(r_earth, v_earth)
h_asteroid = angular_momentum(r_asteroid, v_asteroid)
n = node_line(h_earth, h_asteroid)
print(f"node line {normalize(n)}")
e_earth = eccentricity_vector(mu_sun, r_earth, v_earth)
w_earth = argument_perigee(n, e_earth)
e_asteroid = eccentricity_vector(mu_sun, r_asteroid, v_asteroid)
w_asteroid = argument_perigee(n, e_asteroid)
dt_ascending_earth, dt_descending_earth = time_to_node_line(
    elements_earth[0], elements_earth[1], w_earth, elements_earth_radians[5], mu_sun
)  # elements_earth_radians[5]
dt_ascending_asteroid, dt_descending_asteroid = time_to_node_line(
    elements_asteroid[0],
    elements_asteroid[1],
    w_asteroid,
    elements_asteroid_radians[5],
    mu_sun,
)
print(f"Time to ascending node Earth: {dt_ascending_earth/86400} days")
print(f"Time to descending node Earth: {dt_descending_earth/86400} days")
print(f"Time to ascending node Asteroid: {dt_ascending_asteroid/86400} days")
print(f"Time to descending node Asteroid: {dt_descending_asteroid/86400} days")
dt = dt_ascending_earth / 86400
transfer_date = ti + dt
dt /= dt * u.day
ttf = Time(transfer_date, format="jd")
print(f"Julian transfer date: {ttf.jd}")
print(f"Date of transfer @ RAAN: {julian_to_current(transfer_date.jd)}")

########################################################################################################################
# Plotting
########################################################################################################################


plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)

plotter.plot_ephem(
    earth_ephem,
    ti,
    label=f"Earth @ Launch Position {start_year, start_month, start_day}",
    color="blue",
)
plotter.plot_ephem(
    asteroid_ephem, ti, label="Asteroid @ Launch Position", color="orange"
)


fig = plotter._figure
fig.write_html(f"figures/{figure_name}.html")
webbrowser.open("file://" + os.path.realpath(f"figures/{figure_name}.html"))
