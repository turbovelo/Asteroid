from copy import deepcopy
from dataclasses import dataclass, field, InitVar, asdict
import json
import logging
import requests
from typing import Dict

from astropy.time.core import Time

from orbitmath import OrbitalElements


@dataclass
class Query:
    body: InitVar[str]
    start_time: InitVar[Time]
    stop_time: InitVar[Time]
    url: str = field(init=False)

    def __post_init__(self, body: int, start_time: Time, stop_time: Time):
        """Constructs the query url from the initial inputs."""

        # For easy access to the year, month and day
        start = start_time.to_datetime()
        stop = stop_time.to_datetime()

        # Url's to connect to jpl Horizon API and retrieve data
        self.url = f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='{body}'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&CENTER='500@10'&EPHEM_TYPE='ELEMENTS'&START_TIME='{start.year}-{start.month}-{start.day}'&STOP_TIME='{stop.year}-{stop.month}-{stop.day}"


class MarkerNotFoundError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def get(query: Query) -> dict:
    logging.info(f'Requesting data from "{query.url}"...')
    response = requests.get(query.url, timeout=10)

    if response.status_code != requests.codes.ok:
        logging.error("Error while getting the data")
        response.raise_for_status()

    return response.json()


def parse(response: dict) -> Dict[str, OrbitalElements]:
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
                if var not in extracted_data[date].keys():
                    continue

                # Sanity check in case parsing screws up
                if not var.isalpha():
                    message = f"Error while parsing the key: '{var}' is not an alphabetic string"
                    logging.error(message)
                    raise KeyError(message)

                extracted_data[date][var] = float(val)

    # Convert the extracted data to dataclasses
    for key, value in extracted_data.items():
        extracted_data[key] = OrbitalElements(**value)

    return extracted_data


def save(elements: dict[str, OrbitalElements], file_name: str) -> None:
    # Convert the dataclasses to dicts
    data = deepcopy(elements)
    for key in elements.keys():
        data[key] = asdict(data[key])

    with open(file_name, "w") as file:
        json.dump(data, file, indent=4)

    logging.info(f"Saved data to {file_name}")


def load(file_name: str) -> dict[str, OrbitalElements]:
    with open(file_name, "r") as file:
        data = json.load(file)

    # Convert to OrbitalElements
    for key in data.keys():
        data[key] = OrbitalElements(**data[key])

    logging.info(f"Loaded data from {file_name}")

    return data


def extract(time: Time, data: Dict[str, OrbitalElements]) -> OrbitalElements:
    FORMAT = "ymdhms"
    time.format = FORMAT  # ensure all formats are the same

    for date, elements in data.items():
        # Prepare the date format for comparison
        date = Time(date, format="jd")
        date.format = FORMAT

        if date.value != time:
            continue

        logging.info(f"Found match for the date: {time.value}")
        return elements

    logging.warning(
        f"No match for the date {time.value} was found in the data, returning empty OrbitalElements"
    )

    return OrbitalElements()
