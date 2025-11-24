"""
Author: Konrad Barboutie
Current date: 2025-11-10
"""

from copy import deepcopy
from datetime import datetime, timedelta
import dateutil
from erfa import ErfaWarning
import logging
from pathlib import Path
import warnings
import webbrowser

from astropy import units as u
from astropy.time import Time, TimeDelta
import plotly
from poliastro.bodies import Sun
from poliastro.ephem import Ephem
from poliastro.plotting.interactive import OrbitPlotter3D
from poliastro.util import time_range
from rich.logging import RichHandler

import cli
import data
import orbitmath as om

plotly.io.renderers.default = "browser"

PLANET_IDS = {
    "199": "Mercury",
    "299": "Venus",
    "399": "Earth",
    "499": "Mars",
    "599": "Jupiter",
    "699": "Saturn",
    "799": "Uranus",
    "899": "Neptune",
    "999": "Pluto",
}

args = cli.get_inputs()

if args.out:
    OUT_DIR = Path(args.out).resolve()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
else:
    BASE_DIR = Path(__file__).resolve().parent
    ROOT_DIR = BASE_DIR.parent
    OUT_DIR = ROOT_DIR / "out"

if args.debug is True:
    level = logging.DEBUG
else:
    level = logging.INFO

# Set up the logger.
logging.basicConfig(
    level=level,
    format="{message}",
    style="{",
    datefmt="[%X]",
    handlers=[
        RichHandler(),
        logging.FileHandler(OUT_DIR / "log.txt", mode="w", encoding="utf-8"),
    ],
)

# Capture all warnings from other modules.
logging.captureWarnings(True)

# Ignore the ErfaWarning for now.
warnings.filterwarnings("ignore", category=ErfaWarning)

logging.info(f"{datetime.now().strftime('%c')}")
logging.info(f"Selected output directory is located at {OUT_DIR}")

########################################################################################################################
# Data Collection & Initialization
########################################################################################################################
# Write Horizon ephemeris data to file for faster execution

figure_name = "orbit"
logging.info(f"{dateutil.parser.parse(args.start)}")

ID_ORIGIN = args.origin
ID_TARGET = args.target

START_TIME = Time(dateutil.parser.parse(args.start), format="datetime")
STOP_TIME = Time(dateutil.parser.parse(args.stop), format="datetime")


query_earth = data.get(data.Query(ID_ORIGIN, START_TIME, STOP_TIME))
query_2008_ev5 = data.get(data.Query(ID_TARGET, START_TIME, STOP_TIME))

data_earth = data.parse(query_earth)
data_asteroid = data.parse(query_2008_ev5)

data.save(data_earth, OUT_DIR / "data_earth.json")
data.save(data_asteroid, OUT_DIR / "data_asteroid.json")

logging.info("Combobulating the discombobulator...")

elements_earth = data.extract(START_TIME, data_earth)
elements_asteroid = data.extract(START_TIME, data_asteroid)

logging.info(f"{elements_earth=}")
logging.info(f"{elements_asteroid=}")

# 2038-02-04 arrival time
observation_time = 1
fast_forward = False  # See fast forward in time by latter specified amount of years


ti = START_TIME
tf = ti + TimeDelta(timedelta(days=365), format="datetime")


epochs = time_range(start=ti, end=tf)

########################################################################################################################
# Poliastro Orbit calculations
########################################################################################################################

logging.info("Determining orbits...")

mu_sun = 1.32712440042e11
mu_earth = 3.986e5

# Orbital Elements lists
elements_earth_radians = deepcopy(elements_earth)
elements_asteroid_radians = deepcopy(elements_asteroid)
elements_earth_radians.to_radians()
elements_asteroid_radians.to_radians()

# Define orbits
earth_initial = om.define_orbit(Sun, elements_earth, ti)
asteroid_initial = om.define_orbit(Sun, elements_asteroid, ti)

# Define The orbit that will be seen in the plot after the elapsed ephem time chosen above
earth_ephem = Ephem.from_orbit(orbit=earth_initial, epochs=epochs)
asteroid_ephem = Ephem.from_orbit(orbit=asteroid_initial, epochs=epochs)

########################################################################################################################
# Maneuvering
########################################################################################################################

logging.info("Calculating manoeuver...")

# Determine line of nodes between the two inclined orbits for Bi-Elliptic Transfer calculations
r_earth, v_earth = earth_initial.rv()
r_earth, v_earth = r_earth.value, v_earth.value
r_asteroid, v_asteroid = asteroid_initial.rv()
r_asteroid, v_asteroid = r_asteroid.value, v_asteroid.value
h_earth = om.angular_momentum(r_earth, v_earth)
h_asteroid = om.angular_momentum(r_asteroid, v_asteroid)
n = om.node_line(h_earth, h_asteroid)
logging.info(f"node line {om.normalize(n)}")
e_earth = om.eccentricity_vector(mu_sun, r_earth, v_earth)
w_earth = om.argument_perigee(n, e_earth)
e_asteroid = om.eccentricity_vector(mu_sun, r_asteroid, v_asteroid)
w_asteroid = om.argument_perigee(n, e_asteroid)
dt_ascending_earth, dt_descending_earth = om.time_to_node_line(
    elements_earth.A,
    elements_earth.EC,
    w_earth,
    elements_earth_radians.TA,
    mu_sun,
)  # elements_earth_radians[5]
dt_ascending_asteroid, dt_descending_asteroid = om.time_to_node_line(
    elements_asteroid.A,
    elements_asteroid.EC,
    w_asteroid,
    elements_asteroid_radians.TA,
    mu_sun,
)
logging.info(f"Time to ascending node Earth: {dt_ascending_earth/86400} days")
logging.info(f"Time to descending node Earth: {dt_descending_earth/86400} days")
logging.info(f"Time to ascending node Asteroid: {dt_ascending_asteroid/86400} days")
logging.info(f"Time to descending node Asteroid: {dt_descending_asteroid/86400} days")
dt = dt_ascending_earth / 86400
transfer_date = ti + dt
dt /= dt * u.day
ttf = Time(transfer_date, format="jd")
logging.info(f"Julian transfer date: {ttf.jd}")
logging.info(f"Date of transfer @ RAAN: {om.julian_to_current(transfer_date.jd)}")

########################################################################################################################
# Plotting
########################################################################################################################

logging.info("Plotting the orbits...")
plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)

if ID_ORIGIN in PLANET_IDS.keys():
    origin_name = PLANET_IDS[ID_ORIGIN]
else:
    origin_name = ID_ORIGIN

if ID_TARGET in PLANET_IDS.keys():
    target_name = PLANET_IDS[ID_ORIGIN]
else:
    target_name = ID_ORIGIN

plotter.plot_ephem(
    earth_ephem,
    ti,
    label=f"{origin_name} @ Launch Position {ti.to_datetime().strftime('%Y %M %d')}",
    color="blue",
)
plotter.plot_ephem(
    asteroid_ephem, ti, label=f"{target_name} @ Launch Position", color="orange"
)

fig = plotter._figure
path = OUT_DIR / f"{figure_name}.html"
fig.write_html(path)
logging.info(f"Saved plot to {path}")

webbrowser.open("file://" / path)
