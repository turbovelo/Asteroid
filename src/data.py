from dataclasses import dataclass, field, InitVar
import logging
import re
import requests
from typing import List


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
