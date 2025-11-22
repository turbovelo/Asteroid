########################################################################################################################
# Function definitions
########################################################################################################################


def get_elements(url):
    # Acquiring and Storing the data as a json file
    response = requests.get(url)
    data = response.json()

    # Verify connection was smooth
    if response.status_code != 200:
        print("ERROR GETTING THE DATA")
        sys.exit()
    else:
        print("Status Code: ", response.status_code)  # 200 means everything ok

    # Extracting orbital elements code by Gemini AI, improved by me
    start_marker = "$$SOE"
    end_marker = "$$EOE"
    start_index = data["result"].find(start_marker)
    end_index = data["result"].find(end_marker)
    if start_index == -1 or end_index == -1:
        print("Error: Could not find the $$SOE / $$EOE markers in the response.")
        exit()
    # Extract the block (start_index + length of '$$SOE' marker)
    ephemeris_block = data["result"][
        start_index + len(start_marker) : end_index
    ].strip()

    # Focus on the data for the first date (June 1, 2030)
    # This is the block starting after $$SOE and ending before the second date
    first_date_block = ephemeris_block.split("2462654.500000000")[0].strip()

    # Extract the key-value pairs using regular expressions
    # This regex looks for: two capital letters (the element code, e.g., EC, QR),
    # followed by an '=' sign and then the floating-point number.
    elements_dictionnary = {}
    matches = re.findall(r"(\s[A-Z]{1,2})\s*=\s*([0-9.\-E+]+)", first_date_block)

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

    KEY_ORDER = [
        "A",
        "EC",
        "IN",
        "OM",
        "W",
        "TA",
        "QR",
        "N",
        "MA",
        "AD",
        "PR",
    ]  # a, e, i, RAAN, argp, nu
    elements_sorted_values = [elements_dictionnary[key] for key in KEY_ORDER]
    # Storing the values in a numpy array for future use
    elements = np.array(elements_sorted_values)
    return elements


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
