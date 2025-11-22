from dataclasses import dataclass, field
from math import floor
from typing import List, Tuple

from astropy import units
from astropy.time import Time
import numpy as np
from poliastro.bodies import Body
from poliastro.ephem import Ephem
from poliastro.twobody import Orbit


@dataclass
class OrbitalElements:
    A: float = field(default=0)
    EC: float = field(default=0)
    IN: float = field(default=0)
    W: float = field(default=0)
    OM: float = field(default=0)
    TA: float = field(default=0)
    NU: float = field(default=0)

    def __post_init__(self):
        self.degrees = True

    def to_radians(self) -> None:
        if self.degrees is False:
            return

        self.A = np.deg2rad(self.A)
        self.EC = np.deg2rad(self.EC)
        self.IN = np.deg2rad(self.IN)
        self.W = np.deg2rad(self.W)
        self.OM = np.deg2rad(self.OM)
        self.TA = np.deg2rad(self.TA)
        self.NU = np.deg2rad(self.NU)

        self.degrees = False

    def to_degrees(self) -> None:
        if self.degrees is True:
            return

        self.A = np.rad2deg(self.A)
        self.EC = np.rad2deg(self.EC)
        self.IN = np.rad2deg(self.IN)
        self.W = np.rad2deg(self.W)
        self.OM = np.rad2deg(self.OM)
        self.TA = np.rad2deg(self.TA)
        self.NU = np.rad2deg(self.NU)

        self.degrees = True


def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    v_norm = v / norm
    return v_norm


def length(v: np.ndarray) -> float:
    v = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v


def node_line(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.cross(x, y)


def eccentricity_vector(mu: float, r: np.ndarray, v: np.ndarray) -> np.ndarray:
    rl = length(r)
    vl = length(v)
    e = 1 / mu * ((vl**2 - mu / rl) * r - (np.dot(r, v)) * v)
    return e


def radius_from_elements_pqw(a: float, e: float, nu: float) -> List[float]:
    r = a * (1 - pow(e, 2)) / (1 + e * np.cos(nu))
    p = r * np.cos(nu)
    q = r * np.sin(nu)
    w = 0
    return [p, q, w]


def velocity_from_elements_pqw(a: float, e: float, nu: float, mu: float) -> List[float]:
    semi_latus_rectum = a * (1 - pow(e, 2))
    r = np.sqrt(mu / semi_latus_rectum)
    p = -r * np.sin(nu)
    q = r * (e + np.cos(nu))
    w = 0
    return [p, q, w]


def rv_calculation(elements: OrbitalElements) -> Tuple[np.ndarray, np.ndarray]:
    r = radius_from_elements_pqw(elements.A, elements.EC, elements.TA) * units.km
    v = (
        velocity_from_elements_pqw(elements.A, elements.EC, elements.TA, mu_sun)
        * units.km
        / units.s
    )
    return r, v


def define_orbit(origin: Body, classical_elements: OrbitalElements, epoch: Time):
    orb = Orbit.from_classical(
        origin,
        classical_elements.A * units.km,
        classical_elements.EC * units.one,
        classical_elements.IN * units.deg,
        classical_elements.W * units.deg,
        classical_elements.OM * units.deg,
        classical_elements.TA * units.deg,
        epoch,
    )
    return orb


def define_ephem(origin: Body, classical_elements: OrbitalElements):
    orb = Ephem.from_orbit(
        origin,
        classical_elements.A * units.km,
        classical_elements.EC * units.one,
        classical_elements.IN * units.deg,
        classical_elements.W * units.deg,
        classical_elements.OM * units.deg,
        classical_elements.TA * units.deg,
    )
    return orb


def sphere_of_influence(m1: float, r: float) -> float:
    mu = 1.32712440042e11  # mu of sun
    return r * (m1 / mu) ** (2 / 5)


def angular_momentum_normalized(i: float, raan: float) -> np.ndarray:
    h = [np.sin(i) * np.sin(raan), -np.sin(i) * np.cos(raan), np.cos(i)]
    return np.array(h)


def nu_to_eccentric_anomaly(e: float, nu: float) -> float:
    E = 2 * np.arctan(np.tan(nu / 2) * np.sqrt((1 - e) / (1 + e)))
    return E


def time_to_node_line(a: float, e: float, w: float, nu: float, mu: float):
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


def julian_to_current(jd: float):
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


def argument_perigee(n: np.ndarray, e: np.ndarray) -> float:
    nl = length(n)
    el = length(e)
    w = np.arccos(np.dot(n, e) / (el * nl))
    if e[2] < 0:
        w = w + np.pi
        return w
    else:
        return w


def angular_momentum(r: np.ndarray, v: np.ndarray):
    return np.cross(r, v)
