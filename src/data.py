from dataclasses import dataclass, field, fields, InitVar, asdict
import json
import logging
import re
import requests
from typing import List, Dict

from astropy.time.core import Time

from orbitmath import OrbitalElements


@dataclass
class Query:
    body: InitVar[str]
    start_year: InitVar[int]
    start_month: InitVar[int]
    start_day: InitVar[int]
    end_year: InitVar[int]
    end_month: InitVar[int]
    end_day: InitVar[int]
    url: str = field(init=False)

    def __post_init__(
        self,
        body: int,
        start_year: int,
        start_month: int,
        start_day: int,
        end_year: int,
        end_month: int,
        end_day: int,
    ):
        """Constructs the query url from the initial inputs."""

        # Url's to connect to jpl Horizon API and retrieve data
        self.url = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='{body}'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start_year}-{start_month}-{start_day}'&STOP_TIME='{end_year}-{end_month}-{end_day}"


class MarkerNotFoundError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def get(query: Query) -> dict:
    logging.info(f'Requesting data from "{query.url}"...')
    response = requests.get(query.url)
    return response.json()


def parse(response: dict, year: int) -> Dict[str, OrbitalElements]:
    if "result" not in response:
        message = "Key 'result' not found in the data!"
        logging.error(message)
        raise KeyError(message)

    contents: str = response["result"]
    MARKER_START = "$$SOE"
    MARKER_STOP = "$$EOE"

    start = contents.find(MARKER_START)
    if start == -1:
        message = f"Marker '{MARKER_START}' not found in the file!"
        logging.error(message)
        raise MarkerNotFoundError(message)

    # Ignore the marker itself.
    start += len(MARKER_START)

    stop = contents.find(MARKER_STOP)
    if stop == -1:
        message = f"Marker '{MARKER_STOP}' not found in the file!"
        logging.error(message)
        raise MarkerNotFoundError(message)

    # Keep only the data section
    lines = contents[start:stop].strip().split("\n")

    STEP = 5  # Size of a full block, containing the date and the orbital elements
    extracted_data = {}
    for block_index in range(0, len(lines), STEP):
        # The date is on the first line, and the rest of the block contains the data
        date = lines[block_index].split("=")[0]
        extracted_data[date] = asdict(OrbitalElements())
        data_block = lines[block_index + 1 : block_index + STEP]

        for line in data_block:
            # Spacing is inconsistent in the file.
            line = line.replace("= ", "=").replace("= ", "=").replace(" =", "=")
            values = line.strip().split(" ")
            values = [val for value in values for val in value.split("=")]

            for var, val in zip(values[::2], values[1::2]):

                # Skip keys that are not required.
                if var not in fields(OrbitalElements):
                    continue

                # Sanity check in case parsing screws up
                if not var.isalpha():
                    message = f"Error while parsing the key: '{var}' is not an alphabetic string"
                    logging.error(message)
                    raise KeyError(message)

                extracted_data[date][var] = float(val)

    # Convert the extracted data to dataclasses
    for key in extracted_data.keys():
        extracted_data[key] = OrbitalElements(**extracted_data[key])

    return extracted_data


def save(elements: dict[str, OrbitalElements], file_name: str) -> None:
    logging.info(f"{elements=}")
    # Convert the dataclasses to dicts
    for key in elements.keys():
        elements[key] = asdict(elements[key])

    with open(file_name, "w") as file:
        json.dump(elements, file, indent=4)


def load(file_name: str) -> dict[str, OrbitalElements]:
    with open(file_name, "r") as file:
        data = json.load(file)

    # Convert to OrbitalElements
    for key in data.keys():
        data[key] = OrbitalElements(**data[key])

    return data


def get_data(url: str, year: int) -> List[float] | None:
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
            logging.info("Fatal Error $$SOE not found")
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
            else:
                return element
