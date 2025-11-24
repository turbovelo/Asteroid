import argparse


def get_inputs() -> argparse.Namespace:
    parser = argparse.ArgumentParser("asteroid", add_help=True)
    parser.add_argument(
        "-origin",
        "-o",
        default="399",
        help="The ID of the body to transfer from",
    )
    parser.add_argument(
        "-target",
        "-t",
        default="DES=2008 EV5",
        help="The ID of the body to transfer to",
    )
    parser.add_argument(
        "-start",
        default="2037-6-24",
        help="The start date for the transfer",
    )
    parser.add_argument(
        "-stop",
        default="2038-6-25",
        help="The stop date for the transfer",
    )
    parser.add_argument(
        "-out",
        required=False,
        help="Output directory for the generated files",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        help="Enable debug level for logging",
    )

    return parser.parse_args()
