"""
Author: Konrad Barboutie
Current date: 2025-11-10
"""

import webbrowser
import os
import numpy as np
import logging
from erfa import ErfaWarning
import warnings

from astropy import units as u
from astropy.time import Time
import plotly
from poliastro.bodies import Sun
from poliastro.ephem import Ephem
from poliastro.plotting.interactive import OrbitPlotter3D
from poliastro.util import time_range

import orbitmath as om
from data import get_data

plotly.io.renderers.default = "browser"

########################################################################################################################
# Function definitions
########################################################################################################################


if __name__ == "__main__":

    # Set up the logger.
    logging.basicConfig(
        level=logging.INFO,
        format="{asctime} {levelname:<8}- {message}",
        datefmt="[%H:%M:%S]",
        style="{",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("log.txt", mode="w", encoding="utf-8"),
        ],
    )

    # Capture all warnings from other modules.
    logging.captureWarnings(True)

    # Ignore the ErfaWarning for now.
    warnings.filterwarnings("ignore", category=ErfaWarning)

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
    logging.info(f"{elements_earth=}")
    logging.info(f"{elements_asteroid=}")
    elements_earth_radians = elements_earth.copy()
    elements_asteroid_radians = elements_asteroid.copy()
    elements_earth_radians[2:6] = np.radians(elements_earth[2:6])
    elements_asteroid_radians[2:6] = np.radians(elements_asteroid[2:6])

    # Define orbits
    earth_initial = om.define_orbit(Sun, elements_earth, ti)
    asteroid_initial = om.define_orbit(Sun, elements_asteroid, ti)

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
    h_earth = om.angular_momentum(r_earth, v_earth)
    h_asteroid = om.angular_momentum(r_asteroid, v_asteroid)
    n = om.node_line(h_earth, h_asteroid)
    logging.info(f"node line {om.normalize(n)}")
    e_earth = om.eccentricity_vector(mu_sun, r_earth, v_earth)
    w_earth = om.argument_perigee(n, e_earth)
    e_asteroid = om.eccentricity_vector(mu_sun, r_asteroid, v_asteroid)
    w_asteroid = om.argument_perigee(n, e_asteroid)
    dt_ascending_earth, dt_descending_earth = om.time_to_node_line(
        elements_earth[0], elements_earth[1], w_earth, elements_earth_radians[5], mu_sun
    )  # elements_earth_radians[5]
    dt_ascending_asteroid, dt_descending_asteroid = om.time_to_node_line(
        elements_asteroid[0],
        elements_asteroid[1],
        w_asteroid,
        elements_asteroid_radians[5],
        mu_sun,
    )
    logging.info(f"Time to ascending node Earth: {dt_ascending_earth/86400} days")
    logging.info(f"Time to descending node Earth: {dt_descending_earth/86400} days")
    logging.info(f"Time to ascending node Asteroid: {dt_ascending_asteroid/86400} days")
    logging.info(
        f"Time to descending node Asteroid: {dt_descending_asteroid/86400} days"
    )
    dt = dt_ascending_earth / 86400
    transfer_date = ti + dt
    dt /= dt * u.day
    ttf = Time(transfer_date, format="jd")
    logging.info(f"Julian transfer date: {ttf.jd}")
    logging.info(f"Date of transfer @ RAAN: {om.julian_to_current(transfer_date.jd)}")

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
