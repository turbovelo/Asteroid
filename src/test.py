from poliastro.twobody import Orbit
from jplephem.calendar import compute_julian_date
from poliastro.bodies import Earth, Mars, Sun
from poliastro.ephem import Ephem
from poliastro.plotting.interactive import OrbitPlotter3D
from astropy import units as u
from astropy.time import Time, TimeDelta
import webbrowser, os
figure_name = 'test'
epoch = Time("2037-06-24 12:00:00", scale="tdb")
A = 1.516661905270592E+08 * u.km
EC= 1.084686128924439E-02 * u.one
IN= 3.183372510637163E-03 * u.deg
W = 3.043798594895089E+02 * u.deg
OM= 1.710197739423058E+02 * u.deg
TA= 1.570207833379229E+02 * u.deg
#X =-1.378249046368877E+02 Y = 1.675263270910506E+03 Z =-6.139195868505277E+03\n VX=-2.901410011378197E-01 VY=-7.604139085540065E-03 VZ= 4.438646604511316E-03
re = [-1.378249046368877E+02, 1.675263270910506E+03,-6.139195868505277E+03] * u.km
ve= [-2.901410011378197E-01, -7.604139085540065E-03, 4.438646604511316E-03] * u.km/u.s
ra = [-1.449391930783899e08, 9.746148285124145e07, 1.848501659581006e07] * u.km
va= [-1.728374960605640E+01, -2.738421095798672E+01, -1.409369092272017E+00] * u.km/u.s

orb_e = Orbit.from_classical(Sun,A, EC, IN, W, OM, TA)
orb_a = Orbit.from_vectors(Sun,ra,va)
plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)
plotter.plot(
    orb_e
)
plotter.plot(
    orb_a
)

fig = plotter._figure
fig.write_html(f"figures/{figure_name}.html")
webbrowser.open('file://' + os.path.realpath(f"figures/{figure_name}.html"))