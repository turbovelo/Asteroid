from matplotlib.figure import Figure
import time



# Load the data from the ephemeris
kernel = SPK.open('data/de421.bsp')
# Make vector output attractive
np.set_printoptions(precision=3)


# Position in km, velocity in km/day
mars_position, mars_velocity = kernel[0,4].compute_and_differentiate(date)
mars_velocity_seconds = mars_velocity/86400.0

#print(mars_position, mars_velocity, mars_velocity_seconds)

kernel.close()
try:
    element_data = data['result']['table'][0]
    eccentricity = element_data['EC']
    semimajor_axis = element_data['AC']
    inclination = element_data['IN']
    print(f"Success")
    print(f"Eccentricity: {eccentricity}")
    print(f"Semimajor-axis: {semimajor_axis}")
    print(f"Inclination: {inclination}")
except KeyError as e:
    print(f"error extracting data")
    if 'error' in data:
        print(f"API error message: {data['error']}")
try:
    # Get the controller for 'firefox'
    browser = webbrowser.get('firefox')

except webbrowser.Error:
    print(f"Could not find Firefox. Plot saved to {file_path}")
    # Fallback to default browser if Firefox isn't found
    # webbrowser.open('file://' + os.path.realpath(file_path))
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

elements_mercury = get_elements(url_mercury)
elements_venus = get_elements(url_venus)
elements_mars = get_elements(url_mars)
orb_mercury = define_orbit(Sun, elements_mercury)
orb_venus = define_orbit(Sun, elements_venus)

orb_mars = define_orbit(Sun, elements_mars)


#plotter.plot(orb_mercury, label = 'Mercury', color = '#B1B1B1')
#plotter.plot(orb_venus, label = 'Venus', color = '#EEDC82')
#plotter.plot(orb_mars, label = 'Mars', color = '#C1440E')

plotter.plot_ephem(
    orb_earth,
    label = 'Earth',
    color = '#2E8B57'
)

h_earth = angular_momentum_normalized(elements_earth[2],elements_earth[3])
h_asteroid = angular_momentum_normalized(elements_asteroid[2],elements_asteroid[3])
relative_raan_vector = np.cross(h_earth, h_asteroid)

relative_raan_vector_normalized = relative_raan_vector/np.sqrt(relative_raan_vector[0]**2+relative_raan_vector[1]**2+relative_raan_vector[2]**2)
relative_raan_angle = np.arctan2(relative_raan_vector[0], relative_raan_vector[1])*180/np.pi
print(relative_raan_vector_normalized)
print(relative_raan_angle)
print("Line of nodes angle between the two orbits:", relative_raan_angle)
# Determine transfer date for impulsive maneuver
dt = time_to_raan(elements_asteroid_radians[0], elements_asteroid_radians[1], elements_asteroid_radians[4], elements_asteroid_radians[5], mu_sun)



k = np.linspace(-149e6, 149e6, 100)
r_coordinates = np.outer(k, relative_raan_vector_normalized)
r_coordinates_quantity = r_coordinates * u.km
x,y,z = r_coordinates_quantity[0], r_coordinates_quantity[1], r_coordinates_quantity[2]
relative_raan_coordinates = CartesianRepresentation(x,y,z)

else:
    with open(file, 'r') as file:
        data = json.load(file)
        elements_earth = np.array(data['Earth'])
        elements_asteroid = np.array(data['Asteroid'])
        elements_asteroid_radians = elements_asteroid.copy()
        elements_earth_radians = elements_earth.copy()
        elements_earth_radians[2:6] = np.radians(elements_earth[2:6])
        elements_asteroid_radians[2:6] = np.radians(elements_asteroid[2:6])

if flag != False:

    if fast_forward != False:
        plotter.plot(
            asteroid_after,
            label=f"asteroid after {years} years",
            color="orange"
        )
        plotter.plot(
            earth_after,
            label=f"earth after {years} years",
            color="blue"
        )

        plotter.plot_maneuver(
            initial_orbit=earth_maneuver_1,
            maneuver=impulse1,
            color="purple",
            label="Transfer Orbit 1"
        )
        plotter.plot_maneuver(
            initial_orbit=transfer_1_end,
            maneuver=impulse2,
            color="green",
            label="Transfer Orbit 2"
        )

        plotter.plot_trajectory(
            transfer_1.sample(max_anomaly=180 * u.deg),
            color='purple',
            label='Transfer Orbit 1'
        )
        plotter.plot_trajectory(
            transfer_2.sample(max_anomaly=360 * u.deg),
            color='green',
            label='Transfer Orbit 2'
        )

# Fast-forward to maneu}ver date
#earth_maneuver_1 = earth_initial.propagate(ttf)
#asteroid_maneuver_1 = asteroid_initial.propagate(ttf)
#r_earth, v_earth = earth_maneuver_1.rv()
#r_asteroid, v_asteroid = asteroid_maneuver_1.rv()
#v = normalize(v_earth.value)
print("v=", v)
dv_tnz = [5, 0.0, 0.0] * u.km/u.s # to calculate! tangential normal -nadir thrust component then rotate it to the pqw
#dv1 = 5*v * u.km/u.s
#dv2 = (0.75*v) * u.km/u.s + [0, 0, 1.75] * u.km/u.s
#impulse1 = Maneuver.impulse(dv1)
#transfer_1 = earth_maneuver_1.apply_maneuver(impulse1)
#a = transfer_1.a
#print(a.value)
#dt2 = np.sqrt(a.value**3/mu_sun)*np.pi/2
transfer_date += dt2
ttf = Time(transfer_date, format = 'jd')
transfer_1_end = transfer_1.propagate_to_anomaly(180 * u.deg)

impulse2 = Maneuver.impulse(dv2)
transfer_2 = transfer_1_end.apply_maneuver(impulse2)

years = 8
asteroid_after = asteroid_initial.propagate(years * u.year)
earth_after = earth_initial.propagate(years * u.year)

phase_angle = np.pi*(1-np.sqrt(((1+r_earth/r_asteroid)/2)**3))

def get_elements(url):
    # Acquiring and Storing the data as a json file
    response = requests.get(url)
    data = response.json()

    # Verify connection was smooth
    if response.status_code != 200:
        print("ERROR GETTING THE DATA")
        sys.exit()
    else:
        print(f"Status Code: {response.status_code}, successful") # 200 means everything ok

    # Extracting orbital elements code by Gemini AI, improved by me
    start_marker = '$$SOE'
    end_marker = '$$EOE'
    start_index = data['result'].find(start_marker)
    end_index = data['result'].find(end_marker)
    if start_index == -1 or end_index == -1:
        print("Error: Could not find the $$SOE / $$EOE markers in the response.")
        exit()
    ephemeris_block = data['result'][start_index + len(start_marker): end_index].strip()
    # Extract the block (start_index + length of '$$SOE' marker)

    # Focus on the data for the first date (June 1, 2030)
    # This is the block starting after $$SOE and ending before the second date
    first_date_block = ephemeris_block.split('\n\n')[0].strip()

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
match lines[j]:
                    case line if 'A=' in lines[j]:
                        print(line)
                    case line if 'EC' in lines[j]:
                        print(line)
                    case line if 'IN' in lines[j]:
                        print(line)
                    case line if 'W' in lines[j]:
                        print(line)
                    case line if 'OM' in lines[j]:
                        print(line)
                    case line if 'TA' in lines[j]:
                        print(line)
                    case _:
                        print('No match found')