import re, requests, json

start_year = 2037
start_month = 6
start_day = 24
end_year = 2037
end_month = 6
end_day = 25

url_ast = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='DES=2008 EV5'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start_year}-{start_month}-{start_day}'&STOP_TIME='{end_year}-{end_month}-{end_day}"
url_earth = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='399'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start_year}-{start_month}-{start_day}'&STOP_TIME='{end_year}-{end_month}-{end_day}"


def get_data(url, year):
    element = [0, 0, 0, 0, 0, 0]
    response = requests.get(url)
    raw_data = response.json()["result"]
    with open("data.txt", "w") as file:
        file.write(raw_data)
    file.close()
    with open("data.txt", "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if "$$SOE" in lines[i]:
                start_index = i

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

                print(f"j: {j}")
                print(lines[j])
                print(element)
            else:
                return element


get_data(url_earth, start_year)
