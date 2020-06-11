import requests
import time
from requests.adapters import HTTPAdapter
from urllib3 import Retry


def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session


def make_request(url):
    t0 = time.time()
    request = None;
    try:
        time.sleep(1)
        request = requests.get(url=url)
    except Exception as x:
        time.sleep(1)
        print('request failed :(', x.__class__.__name__)
    else:
        print('request is successful', request.status_code)
    finally:
        t1 = time.time()
        print('Took', t1 - t0, 'seconds')
    return request

