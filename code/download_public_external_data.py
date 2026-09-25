"""Download only public, URL-confirmed external inputs; never bypass licenses."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

PUBLIC = {
    "hawrylycz_table6": (
        "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4171/MediaObjects/41593_2015_BFnn4171_MOESM97_ESM.csv",
        "Hawrylycz_2015_Supplementary_Table_6.csv",
    )
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("resource", choices=sorted(PUBLIC))
    parser.add_argument("--destination", required=True)
    args = parser.parse_args()
    url, filename = PUBLIC[args.resource]
    out = Path(args.destination).resolve() / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise SystemExit(f"refusing to overwrite existing file: {out}")
    urllib.request.urlretrieve(url, out)
    print(out)


if __name__ == "__main__":
    main()
