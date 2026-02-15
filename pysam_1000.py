#!/usr/bin/env python3
from collections import defaultdict
import csv
import os
import pysam
import urllib.request

PANEL_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"


# Map chromosome number to VCF URL
VCF_URLS = {
    str(i): f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    for i in list(range(1, 23)) + ["X", "Y"]
}

ALLELES_CSV = "andrei_snp.csv"

# Download panel file if not exists


def download_file(url, filename):
    if not os.path.isfile(filename):
        print(f"[INFO] Downloading {filename} from {url} ...")
        urllib.request.urlretrieve(url, filename)
        print(f"[INFO] Downloaded {filename}.")
    else:
        print(f"[INFO] {filename} already exists. Skipping download.")

# Download VCF and index if not exists


def collect_all_snp_data(snps_by_chrom, snp_info):
    """Collect all SNP data across chromosomes and organize by sample"""
    print("[INFO] Collecting all SNP data...")

    # sample_id -> [(rsid, chrom, pos, genotype), ...]
    all_sample_data = defaultdict(list)

    for idx, (chrom, positions) in enumerate(snps_by_chrom.items(), 1):
        print(
            f"\n[INFO] Processing chromosome {chrom} ({idx}/{len(snps_by_chrom)}) with {len(positions)} SNPs...")

        vcf_url = VCF_URLS[chrom]

        print(f"    [INFO] Fetching data for chr{chrom} from remote VCF...")
        vcf = pysam.VariantFile(vcf_url)

        for pos_idx, pos in enumerate(positions, 1):
            print(
                f"      [PROGRESS] Fetching variant {pos_idx}/{len(positions)} at {chrom}:{pos}")
            for record in vcf.fetch(chrom, int(pos)-1, int(pos)):
                chrom_col = record.chrom
                pos_col = str(record.pos)
                id_col = record.id
                ref = record.ref
                alt = record.alts[0] if record.alts else ''
                af = record.info.get('AF', [None])[0]

                key = (chrom_col, pos_col)
                if key in snp_info:
                    current_variant = {
                        'rsid': snp_info[key]['SNP'],
                        'chrom': chrom_col,
                        'pos': pos_col,
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

        vcf.close()
        print(f"    [INFO] Chromosome {chrom} processed.")

    return all_sample_data, []


def write_sample_files(sample_data, panel_file):
    """Write individual sample files organized by population"""
    print("[INFO] Writing individual sample files...")

    # Load population data
    pop_info = {}
    with open(panel_file) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                # Map columns by header
                entry = dict(zip(header, parts))
                sample_id = entry.get('sample', '')
                pop = entry.get('pop', '')
                superpop = entry.get('super_pop', '')
                pop_info[sample_id] = {'pop': pop, 'superpop': superpop}

    # Create output directory structure
    output_dir = "raw_genotype_data"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Group samples by population
    samples_by_pop = defaultdict(list)
    for sample_id in sample_data.keys():
        if sample_id in pop_info:
            pop = pop_info[sample_id]['pop']
            samples_by_pop[pop].append(sample_id)

    # Write files
    for pop, sample_ids in samples_by_pop.items():
        pop_dir = os.path.join(output_dir, pop)
        if not os.path.exists(pop_dir):
            os.makedirs(pop_dir)

        print(
            f"    [INFO] Writing {len(sample_ids)} samples for population {pop}...")

        for idx, sample_id in enumerate(sample_ids, 1):
            sample_file = os.path.join(pop_dir, f"{sample_id}.txt")
            with open(sample_file, 'w') as f:
                # Write sample details at the top
                pop = pop_info.get(sample_id, {}).get('pop', '')
                superpop = pop_info.get(sample_id, {}).get('superpop', '')
                f.write(f"# Sample : {sample_id}\n")
                f.write(f"# Population : {pop}\n")
                f.write(f"# Region : {superpop}\n")
                # Write header
                f.write("# rsid\tchromosome\tposition\tgenotype\n")

                # Sort SNPs by chromosome and position
                snp_data = sorted(sample_data[sample_id], key=lambda x: (
                    int(x[1]) if x[1].isdigit() else 999, int(x[2])))

                # Write SNP data
                for rsid, chrom, pos, genotype in snp_data:
                    f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            if idx % 100 == 0 or idx == len(sample_ids):
                print(
                    f"      [PROGRESS] {idx}/{len(sample_ids)} samples written for population {pop}")

    print(
        f"[INFO] Individual sample files written to {output_dir}/ organized by population")
    return output_dir


def main():
    print("[INFO] Script started. Creating individual sample files like 23andMe format.")
    print("[INFO] Downloading panel file if needed...")
    download_file(PANEL_URL, PANEL_FILE)

    print(
        f"[INFO] Reading SNP list from {ALLELES_CSV} and grouping by chromosome...")
    snps_by_chrom = defaultdict(list)
    snp_info = {}

    # Detect which columns to use for chromosome and position (chr37 or chr38)
    with open(ALLELES_CSV, 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if 'chr37' in header and 'pos37' in header:
            chrom_col = 'chr37'
            pos_col = 'pos37'
            print("[INFO] Detected chr37/pos37 columns. Using GRCh37 coordinates.")
        elif 'chr38' in header and 'pos38' in header:
            chrom_col = 'chr38'
            pos_col = 'pos38'
            print("[INFO] Detected chr38/pos38 columns. Using GRCh38 coordinates.")
        else:
            raise ValueError(
                "Could not detect chr37/chr38 columns in the CSV header.")

        for row in reader:
            snp_id = row.get('rsID', row.get('SNP', ''))
            chrom = str(row[chrom_col])
            pos = str(row[pos_col])
            snps_by_chrom[chrom].append(pos)
            snp_info[(chrom, pos)] = {
                'SNP': snp_id,
                'Chromosome': chrom,
                'Position': pos,
                'Gene': row.get('Gene', ''),
                'Ancestral_Allele': row.get('Ancestral_Allele', ''),
                'Derived_Allele': row.get('Derived_Allele', '')
            }

    print(
        f"[INFO] Loaded {len(snp_info)} SNPs across {len(snps_by_chrom)} chromosomes")

    # Collect all SNP data organized by sample
    all_sample_data, tempfiles = collect_all_snp_data(snps_by_chrom, snp_info)
    print(f"[INFO] Collected data for {len(all_sample_data)} samples")

    # Write individual sample files organized by population
    output_dir = write_sample_files(all_sample_data, PANEL_FILE)

    print(f"\n[INFO] Done! Individual sample files written to {output_dir}/")
    print(f"[INFO] Each population folder contains individual sample files in 23andMe format.")


if __name__ == "__main__":
    main()
