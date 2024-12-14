import requests
import keyring
import urllib
import json
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query ADS by title and return citation in custom format."
    )
    parser.add_argument("title", nargs=1, type=str, help="Paper title.")
    parser.add_argument(
        "--fmt",
        nargs=1,
        type=str,
        help="Citation format string.",
    )

    args = parser.parse_args()
    title = args.title[0]

    if args.fmt is not None:
        fmt = args.fmt[0]
    else:
        fmt = '%l (%Y) "%T", %J, %V(%p), http://dx.doi.org/%d'

    token = keyring.get_password("ADS_API_TOKEN", "ADS_API_TOKEN")

    titleenc = urllib.parse.quote(title)

    # the query parameters can be included as part of the URL
    r = requests.get(
        rf'https://api.adsabs.harvard.edu/v1/search/query?q=title%3A"{titleenc}"&fl=bibcode&rows=1',
        headers={"Authorization": "Bearer " + token},
    )

    assert r.status_code == 200, "Request failed."
    bibcode = r.json()["response"]["docs"][0]["bibcode"]

    r2 = requests.post(
        "https://api.adsabs.harvard.edu/v1/export/custom",
        headers={
            "Authorization": "Bearer " + token,
            "Content-type": "application/json",
        },
        data=json.dumps({"bibcode": [bibcode], "format": fmt}),
    )

    assert r2.status_code == 200, "Second request failed."

    print(r2.json()["export"])
