#!/usr/bin/env python3
"""
get_rsid_positions.py

For each rsID in alleles.csv, fetch its position in both GRCh37 and GRCh38 using Ensembl REST API.
Outputs: rsid_positions.csv with columns: rsID,chr37,pos37,chr38,pos38
"""
import csv
import requests
import time

ENSEMBL_REST_37 = "https://grch37.rest.ensembl.org/variation/human/{}?content-type=application/json"
ENSEMBL_REST_38 = "https://rest.ensembl.org/variation/human/{}?content-type=application/json"

input_file = "alleles.csv"

output_file = "rsid_positions.csv"

rsids = []
with open(input_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        rsid = row["SNP"].strip()
        if rsid:
            rsids.append(rsid)

results = []
total = len(rsids)
for idx, rsid in enumerate(rsids, 1):
    print(f"[{idx}/{total}] Processing {rsid}...", flush=True)
    row = {"rsID": rsid, "chr37": "", "pos37": "",
           "chr38": "", "pos38": ""}
    # GRCh37
    try:
        r37 = requests.get(ENSEMBL_REST_37.format(rsid), timeout=10)
        if r37.ok:
            data = r37.json()
            print(f"  [GRCh37 API response] {data}", flush=True)
            mappings = data.get("mappings", [])
            for m in mappings:
                if m.get("assembly_name") == "GRCh37":
                    row["chr37"] = m.get("seq_region_name", "")
                    row["pos37"] = m.get("start", "")
                    break
        time.sleep(0.1)
    except Exception as e:
        print(f"  [Warning] GRCh37 lookup failed for {rsid}: {e}", flush=True)
    # GRCh38
    try:
        r38 = requests.get(ENSEMBL_REST_38.format(rsid), timeout=10)
        if r38.ok:
            data = r38.json()
            print(f"  [GRCh38 API response] {data}", flush=True)
            mappings = data.get("mappings", [])
            for m in mappings:
                if m.get("assembly_name") == "GRCh38":
                    row["chr38"] = m.get("seq_region_name", "")
                    row["pos38"] = m.get("start", "")
                    break
        time.sleep(0.1)
    except Exception as e:
        print(f"  [Warning] GRCh38 lookup failed for {rsid}: {e}", flush=True)
    results.append(row)


with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=["rsID", "chr37", "pos37", "chr38", "pos38"])
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"Wrote {len(results)} rsIDs to {output_file}")
