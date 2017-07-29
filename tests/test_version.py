import cubric_mrs

import requests

pypi_url = 'https://pypi.python.org/pypi/cubric_mrs/json'


def test_version():
    r = requests.get(pypi_url)
    assert cubric_mrs.__version__ not in r.json()["releases"].keys()
