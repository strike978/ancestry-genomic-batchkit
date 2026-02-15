#!/usr/bin/env python3
from collections import defaultdict
import csv
import os
import pysam
import urllib.request
import urllib.parse
import json
import time


# HGDP VCF URLs from Sanger Institute
HGDP_VCF_URLS = {
    str(i): f"https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr{i}.vcf.gz"
    for i in list(range(1, 23)) + ["X", "Y"]
}
ALLELES_CSV = "andrei_snp.csv"


def download_hgdp_metadata():
    """Download HGDP sample metadata with population information"""
    metadata = {}
    try:
        url = "https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/hgdp_wgs.20190516.metadata.txt"
        print("[INFO] Downloading HGDP metadata...")

        with urllib.request.urlopen(url) as response:
            content = response.read().decode('utf-8')

        lines = content.strip().split('\n')
        header = lines[0].split('\t')

        # Find column indices
        sample_idx = header.index('sample')
        population_idx = header.index('population')
        region_idx = header.index('region')
        latitude_idx = header.index('latitude')
        longitude_idx = header.index('longitude')
        sex_idx = header.index('sex')

        for line in lines[1:]:
            fields = line.split('\t')
            if len(fields) > max(sample_idx, population_idx, region_idx, latitude_idx, longitude_idx, sex_idx):
                sample_id = fields[sample_idx]
                metadata[sample_id] = {
                    'population': fields[population_idx],
                    'region': fields[region_idx],
                    'latitude': fields[latitude_idx],
                    'longitude': fields[longitude_idx],
                    'sex': fields[sex_idx]
                }

        print(f"[INFO] Downloaded metadata for {len(metadata)} HGDP samples")
        return metadata

    except Exception as e:
        print(f"[ERROR] Failed to download HGDP metadata: {e}")
        return {}


def collect_hgdp_snp_data(snps_by_chrom, snp_info):
    """Collect all SNP data across chromosomes and organize by sample"""
    print("[INFO] Collecting all HGDP SNP data...")

    # sample_id -> [(rsid, chrom, pos, genotype), ...]
    all_sample_data = defaultdict(list)

    for idx, (chrom, positions) in enumerate(snps_by_chrom.items(), 1):
        print(
            f"\n[INFO] Processing chromosome {chrom} ({idx}/{len(snps_by_chrom)}) with {len(positions)} SNPs...")

        vcf_url = HGDP_VCF_URLS.get(chrom)
        if not vcf_url:
            print(f"[WARNING] No VCF URL for chromosome {chrom}. Skipping.")
            continue
        print(
            f"    [INFO] Fetching data for chr{chrom} from remote HGDP VCF...")
        try:
            vcf = pysam.VariantFile(vcf_url)
        except Exception as e:
            print(f"[ERROR] Could not open VCF for chr{chrom}: {e}")
            continue

        for pos_idx, pos in enumerate(positions, 1):
            print(
                f"      [PROGRESS] Fetching variant {pos_idx}/{len(positions)} at {chrom}:{pos}")
            try:
                # Try both with and without 'chr' prefix
                found_records = False
                for chr_format in [f"chr{chrom}", chrom]:
                    try:
                        for record in vcf.fetch(chr_format, int(pos)-1, int(pos)):
                            found_records = True
                            chrom_col = record.chrom.replace("chr", "")
                            pos_col = str(record.pos)
                            ref = record.ref
                            alt = record.alts[0] if record.alts else ''

                            key = (chrom_col, pos_col)
                            if key in snp_info:
                                current_variant = {
                                    'rsid': snp_info[key]['rsID'],
                                    'chrom': snp_info[key]['chrom'],
                                    'pos': snp_info[key]['pos'],
                                    'ref': ref,
                                    'alt': alt
                                }
                                print(
                                    f"        [PROGRESS] Variant {current_variant['rsid']} at {chrom_col}:{pos_col}")

                                for sample_idx, sample in enumerate(record.samples, 1):
                                    gt = record.samples[sample].get('GT')
                                    if gt is None or None in gt:
                                        genotype = '--'
                                    else:
                                        gt0, gt1 = gt
                                        if gt0 == 0 and gt1 == 0:
                                            genotype = ref + ref
                                        elif gt0 == 1 and gt1 == 1:
                                            genotype = alt + alt
                                        elif (gt0 == 0 and gt1 == 1) or (gt0 == 1 and gt1 == 0):
                                            genotype = ref + alt
                                        else:
                                            genotype = '--'

                                    all_sample_data[sample].append((
                                        current_variant['rsid'],
                                        current_variant['chrom'],
                                        current_variant['pos'],
                                        genotype
                                    ))
                                    if sample_idx % 500 == 0:
                                        print(
                                            f"          [PROGRESS] {sample_idx} samples processed for variant {current_variant['rsid']}")
                        if found_records:
                            break
                    except:
                        continue
            except Exception as e:
                print(f"      [WARN] Error processing position {pos}: {e}")
                continue

        vcf.close()
        print(f"    [INFO] Chromosome {chrom} processed.")

    return all_sample_data


def write_hgdp_sample_files(sample_data, metadata):
    """Write HGDP individual sample files (one per sample) with population information"""
    print("[INFO] Writing HGDP individual sample files...")

    output_dir = "hgdp"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for sample_id in sample_data.keys():
        sample_info = metadata.get(sample_id, {})
        population = sample_info.get('population', 'Unknown')
        region = sample_info.get('region', 'Unknown')
        sex = sample_info.get('sex', 'Unknown')
        latitude = sample_info.get('latitude', 'Unknown')
        longitude = sample_info.get('longitude', 'Unknown')

        # Sanitize population name for folder
        pop_folder = population.replace(' ', '_').replace(
            '/', '_') if population else 'Unknown'
        pop_dir = os.path.join(output_dir, pop_folder)
        if not os.path.exists(pop_dir):
            os.makedirs(pop_dir)

        sample_file = os.path.join(pop_dir, f"{sample_id}.txt")
        with open(sample_file, 'w') as f:
            f.write(f"# Sample: {sample_id}\n")
            f.write(f"# Population: {population}\n")
            f.write(f"# Region: {region}\n")
            f.write(f"# Sex: {sex}\n")
            f.write(f"# Coordinates: {latitude}, {longitude}\n")
            f.write("# Data source: Human Genome Diversity Project (HGDP)\n")
            f.write("# rsid\tchromosome\tposition\tgenotype\n")

            snp_data = sorted(sample_data[sample_id], key=lambda x: (
                int(x[1]) if x[1].isdigit() else 999, int(x[2])))
            for rsid, chrom, pos, genotype in snp_data:
                f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")

    print(
        f"[INFO] HGDP sample files with population information written to {output_dir}/<Population>/")
    return output_dir


def main():
    print("[INFO] Script started. Creating HGDP individual sample files from Human Origins data.")

    # Download HGDP metadata first
    hgdp_metadata = download_hgdp_metadata()
    if not hgdp_metadata:
        print("[WARNING] Could not download HGDP metadata. Files will not include population information.")

    print(
        f"[INFO] Reading SNP list from {ALLELES_CSV} and grouping by chromosome (using GRCh38 positions)...")
    snps_by_chrom = defaultdict(list)
    snp_info = {}

    with open(ALLELES_CSV, 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if 'chr38' in header and 'pos38' in header:
            chrom_col = 'chr38'
            pos_col = 'pos38'
            print("[INFO] Detected chr38/pos38 columns. Using GRCh38 coordinates.")
        else:
            raise ValueError(
                "Could not detect chr38 columns in the CSV header.")

        for row in reader:
            rsid = row['rsID']
            chrom = str(row[chrom_col])
            pos = str(row[pos_col])
            # Accept chromosomes 1-22, X, Y only
            if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                continue
            snps_by_chrom[chrom].append(pos)
            snp_info[(chrom, pos)] = {
                'rsID': rsid,
                'chrom': chrom,
                'pos': pos
            }

    print(
        f"[INFO] Loaded {len(snp_info)} SNPs across {len(snps_by_chrom)} chromosomes")

    # Collect HGDP SNP data
    all_sample_data = collect_hgdp_snp_data(snps_by_chrom, snp_info)
    print(f"[INFO] Collected data for {len(all_sample_data)} samples")

    # Write HGDP individual sample files with population information
    output_dir = write_hgdp_sample_files(all_sample_data, hgdp_metadata)

    print(
        f"\n[INFO] Done! HGDP individual sample files with population groups written to {output_dir}/")
    print(f"[INFO] Each file contains SNPs for one HGDP sample in 23andMe-like format with population information.")


if __name__ == "__main__":
    main()
