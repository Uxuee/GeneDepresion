"""Verify recovered Tier4 CSVs directly from immutable Git blobs."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance/tier4_ranked_signature_recovery_manifest.json"
REPORT = ROOT / "audit/tier4_git_artifact_verification.json"


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks = []
    ok = True
    for item in manifest["files"]:
        payload = subprocess.check_output(["git", "cat-file", "blob", item["blob"]], cwd=ROOT)
        digest = hashlib.sha256(payload).hexdigest()
        rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig"))))
        ids = [row["gene_id"] for row in rows]
        record = {
            "region": item["region"],
            "blob": item["blob"],
            "sha256": digest,
            "rows": len(rows),
            "unique_gene_ids": len(set(ids)),
            "hash_matches": digest == item["sha256"],
            "row_count_matches": len(rows) == item["rows"],
            "gene_ids_unique": len(ids) == len(set(ids)),
        }
        record["ok"] = record["hash_matches"] and record["row_count_matches"] and record["gene_ids_unique"]
        ok &= record["ok"]
        checks.append(record)
    REPORT.write_text(json.dumps({"ok": ok, "checks": checks}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"ok": ok, "artifacts_verified": len(checks)}, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
