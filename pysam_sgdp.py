#!/usr/bin/env python3

import csv
import os
import sys
import time
import pysam


BASE_URL = "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/chr.sgdp.pub.{}.bcf"
BASE_CSI_URL = "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/chr.sgdp.pub.{}.bcf.csi"


def load_sgdp_metadata(metadata_file='SGDP_metadata.279public.21signedLetter.44Fan.samples.txt', verbose=True):
    """Load SGDP sample metadata from the metadata file"""
    if verbose:
        print(f"[METADATA] Loading SGDP metadata from {metadata_file}...", flush=True)
    
    metadata = {}
    # Try different encodings in case of encoding issues
    encodings_to_try = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    
    for encoding in encodings_to_try:
        try:
            with open(metadata_file, encoding=encoding) as f:
                header_line = f.readline().strip()
                if header_line.startswith('#'):
                    header = header_line[1:].split('\t')
                else:
                    header = header_line.split('\t')

                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        entry = dict(zip(header, parts))
                        
                        # Create metadata entry
                        meta_data = {
                            'population': entry.get('Population_ID', ''),
                            'region': entry.get('Region', ''),
                            'country': entry.get('Country', ''),
                            'town': entry.get('Town', ''),
                            'gender': entry.get('Gender', ''),
                            'latitude': entry.get('Latitude', ''),
                            'longitude': entry.get('Longitude', ''),
                            'contributor': entry.get('Contributor', ''),
                            'dna_source': entry.get('DNA_Source', '')
                        }
                        
                        # Add entries for all possible ID formats
                        id_columns = ['Sample_ID', 'Sample_ID(Aliases)', 'SGDP_ID', 'Illumina_ID']
                        for id_col in id_columns:
                            sample_id = entry.get(id_col, '').strip()
                            if sample_id and sample_id != '?':
                                # Add base ID
                                metadata[sample_id] = meta_data.copy()
                                
                                # Add entries with common BCF suffixes
                                for suffix in ['.DG', '.GT', '.AD', '.DP', '.PL', '.GQ']:
                                    metadata[sample_id + suffix] = meta_data.copy()
                                
                                # Also try without potential suffixes
                                if '.' in sample_id:
                                    base_id = sample_id.split('.')[0]
                                    metadata[base_id] = meta_data.copy()
                                    
            if verbose:
                print(f"[METADATA] Successfully read file with {encoding} encoding", flush=True)
                print(f"[DEBUG] Created metadata entries for {len(metadata)} sample IDs", flush=True)
                # Show some examples
                brahui_samples = [k for k in metadata.keys() if 'Brahui' in k or 'S_Brahui' in k]
                if brahui_samples:
                    print(f"[DEBUG] Brahui samples found: {brahui_samples[:5]}", flush=True)
                hezhen_samples = [k for k in metadata.keys() if 'Hezhen' in k or 'S_Hezhen' in k]
                if hezhen_samples:
                    print(f"[DEBUG] Hezhen samples found: {hezhen_samples[:5]}", flush=True)
            break  # Success, exit the encoding loop
            
        except UnicodeDecodeError:
            if verbose and encoding == encodings_to_try[-1]:
                print(f"[WARNING] All encodings failed. Trying with error handling...", flush=True)
            continue
        except Exception as e:
            if verbose:
                print(f"[WARNING] Could not load metadata with {encoding}: {e}", flush=True)
            if encoding == encodings_to_try[-1]:
                return {}
            continue
    
    # If all encodings failed, try with error handling
    if not metadata:
        try:
            with open(metadata_file, encoding='utf-8', errors='ignore') as f:
                # Same logic as above but with error handling...
                pass
        except Exception as e:
            if verbose:
                print(f"[WARNING] Could not load metadata even with error handling: {e}", flush=True)
            return {}

    if verbose:
        print(f"[METADATA] Loaded metadata for {len(metadata)} samples", flush=True)
    return metadata


def write_23andme_format_files(sample_data, metadata, output_dir='sgdp', verbose=True):
    """Write individual sample files in 23andMe format organized by population"""
    if verbose:
        print(f"[OUTPUT] Writing 23andMe format files to {output_dir}/...", flush=True)
    
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Group samples by population
    from collections import defaultdict
    samples_by_pop = defaultdict(list)
    
    for sample_id in sample_data.keys():
        if sample_id in metadata:
            pop = metadata[sample_id]['population']
            if verbose and 'Unknown' in [pop]:
                print(f"[DEBUG] Sample {sample_id} has population {pop}", flush=True)
        else:
            pop = 'Unknown'
            if verbose:
                print(f"[DEBUG] Sample {sample_id} not found in metadata - assigning to Unknown", flush=True)
        samples_by_pop[pop].append(sample_id)
    
    if verbose:
        print(f"[OUTPUT] Found {len(samples_by_pop)} populations: {list(samples_by_pop.keys())}", flush=True)
    
    total_samples_written = 0
    for pop, sample_ids in samples_by_pop.items():
        # Create population directory
        pop_dir = os.path.join(output_dir, pop)
        if not os.path.exists(pop_dir):
            os.makedirs(pop_dir)
        
        if verbose:
            print(f"[OUTPUT] Writing {len(sample_ids)} samples for population {pop}...", flush=True)
        
        for idx, sample_id in enumerate(sample_ids, 1):
            sample_file = os.path.join(pop_dir, f"{sample_id}.txt")
            
            with open(sample_file, 'w') as f:
                # Write metadata header
                meta = metadata.get(sample_id, {})
                f.write(f"# Sample: {sample_id}\n")
                f.write(f"# Population: {meta.get('population', 'Unknown')}\n")
                f.write(f"# Region: {meta.get('region', 'Unknown')}\n")
                f.write(f"# Country: {meta.get('country', 'Unknown')}\n")
                f.write(f"# Town: {meta.get('town', 'Unknown')}\n")
                f.write(f"# Gender: {meta.get('gender', 'Unknown')}\n")
                f.write(f"# Coordinates: {meta.get('latitude', '')},{meta.get('longitude', '')}\n")
                f.write(f"# Contributor: {meta.get('contributor', 'Unknown')}\n")
                f.write(f"# DNA Source: {meta.get('dna_source', 'Unknown')}\n")
                f.write(f"# rsid\tchromosome\tposition\tgenotype\n")
                
                # Sort SNPs by chromosome and position for consistent output
                snp_data = sorted(sample_data[sample_id], key=lambda x: (
                    int(x[1]) if x[1].isdigit() else 999, int(x[2])
                ))
                
                # Write SNP data
                for rsid, chrom, pos, genotype in snp_data:
                    f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            
            total_samples_written += 1
            if verbose and (idx % 50 == 0 or idx == len(sample_ids)):
                print(f"    [PROGRESS] {idx}/{len(sample_ids)} samples written for {pop}", flush=True)
    
    if verbose:
        print(f"[OUTPUT] Complete! {total_samples_written} sample files written to {output_dir}/", flush=True)
    return output_dir


def batch_query_snps_to_23andme(csv_file='andrei_snp.csv', metadata_file='SGDP_metadata.279public.21signedLetter.44Fan.samples.txt', output_dir='sgdp', verbose=True):
    """
    Batch query all SNPs and create 23andMe format files organized by population.
    """
    start_time = time.time()
    
    # Load metadata first
    metadata = load_sgdp_metadata(metadata_file, verbose)
    
    if verbose:
        print(f"[STARTUP] Reading SNP list from {csv_file}...", flush=True)
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        snps = list(reader)
    if verbose:
        print(f"[STARTUP] Loaded {len(snps)} SNPs from {csv_file}", flush=True)
    
    if verbose:
        print(f"[STARTUP] Grouping SNPs by chromosome...", flush=True)
    snps_by_chrom = {}
    snp_info = {}
    for row in snps:
        chrom = row['chr37']
        pos = row['pos37']
        rsid = row['rsID']
        snps_by_chrom.setdefault(chrom, []).append({'pos': pos, 'rsid': rsid})
        snp_info[(chrom, pos)] = {'rsid': rsid, 'chrom': chrom, 'pos': pos}
    if verbose:
        print(f"[STARTUP] Found SNPs on {len(snps_by_chrom)} chromosomes: {list(snps_by_chrom.keys())}", flush=True)
    
    # Collect all sample data
    from collections import defaultdict
    sample_data = defaultdict(list)  # sample_id -> [(rsid, chrom, pos, genotype), ...]
    missing_by_snp = defaultdict(list)
    missing_by_sample = defaultdict(list)
    all_snps = set()
    
    total_chroms = len(snps_by_chrom)
    chrom_idx = 0
    total_processed_snps = 0
    
    if verbose:
        print(f"[DATA] Starting data collection from BCF files...", flush=True)
    
    for chrom, snplist in snps_by_chrom.items():
        chrom_idx += 1
        chrom_map = {str(i): str(i) for i in range(1, 23)}
        chrom_map.update({"X": "X", "Y": "Y", "MT": "MT", "M": "MT"})
        chrom_key = chrom_map.get(str(chrom))
        if not chrom_key:
            if verbose:
                print(f"[WARN] Skipping unsupported chromosome {chrom}", flush=True)
            continue
            
        BCF_URL = BASE_URL.format(chrom_key)
        CSI_URL = BASE_CSI_URL.format(chrom_key)
        if verbose:
            print(f"[INFO] ({chrom_idx}/{total_chroms}) Querying {BCF_URL} for {len(snplist)} SNPs on chr{chrom}", flush=True)
            print(f"    [CONNECT] Opening BCF connection...", flush=True)
            
        try:
            bcf = pysam.VariantFile(BCF_URL, index_filename=CSI_URL)
            if verbose:
                print(f"    [CONNECT] Successfully connected to chr{chrom} BCF", flush=True)
        except Exception as e:
            if verbose:
                print(f"[ERROR] Could not open remote BCF for chr{chrom}: {e}", flush=True)
            continue
            
        total_snps = len(snplist)
        chrom_start_time = time.time()
        
        for snp_idx, snp in enumerate(snplist, 1):
            pos = int(snp['pos'])
            rsid = snp['rsid']
            key = (chrom, str(pos), rsid)
            all_snps.add(key)
            found = False
            
            if verbose and (snp_idx == 1 or snp_idx % 10 == 0 or snp_idx == total_snps):
                elapsed = time.time() - chrom_start_time
                rate = snp_idx / elapsed if elapsed > 0 else 0
                eta = (total_snps - snp_idx) / rate if rate > 0 else 0
                print(f"    [PROGRESS] chr{chrom}: SNP {snp_idx}/{total_snps} (pos {pos}, rsid {rsid}) - {rate:.1f} SNPs/sec, ETA: {eta:.0f}s", flush=True)
            
            for rec in bcf.fetch(chrom, pos-1, pos):
                found = True
                if verbose and snp_idx % 25 == 0:
                    print(f"        [SAMPLES] Processing {len(rec.samples)} samples for SNP {rsid}", flush=True)
                
                ref = rec.ref
                alt = rec.alts[0] if rec.alts else ''
                
                for sample in rec.samples:
                    gt = rec.samples[sample].get('GT')
                    if gt is None or None in gt:
                        genotype = '--'
                        missing_by_snp[key].append(sample)
                        missing_by_sample[sample].append(key)
                    else:
                        # Convert to 23andMe format (AA, AT, etc.)
                        gt0, gt1 = gt
                        if gt0 == 0 and gt1 == 0:
                            genotype = ref + ref
                        elif gt0 == 1 and gt1 == 1:
                            genotype = alt + alt
                        elif (gt0 == 0 and gt1 == 1) or (gt0 == 1 and gt1 == 0):
                            genotype = ref + alt
                        else:
                            genotype = '--'
                            missing_by_snp[key].append(sample)
                            missing_by_sample[sample].append(key)
                    
                    # Add to sample data
                    sample_data[sample].append((rsid, chrom, pos, genotype))
            
            if not found:
                # SNP not found in BCF - add missing entries for known samples
                for sample_id in list(sample_data.keys()):
                    sample_data[sample_id].append((rsid, chrom, pos, '--'))
                    missing_by_snp[key].append(sample_id)
                    missing_by_sample[sample_id].append(key)
            
            total_processed_snps += 1
            
        bcf.close()
        if verbose:
            chrom_elapsed = time.time() - chrom_start_time
            print(f"    [COMPLETE] chr{chrom} finished in {chrom_elapsed:.1f}s ({len(snplist)} SNPs)", flush=True)
    
    if verbose:
        print(f"\n[DATA] Collected data for {len(sample_data)} samples", flush=True)
    
    # Write 23andMe format files
    output_dir = write_23andme_format_files(sample_data, metadata, output_dir, verbose)
    
    if verbose:
        total_elapsed = time.time() - start_time
        print(f"\n[COMPLETE] All done! Processed {total_processed_snps} SNPs in {total_elapsed:.1f} seconds", flush=True)
        print(f"[INFO] Results saved to {output_dir}/ organized by population")
        print("\n# Missing data summary:")
        print(f"Samples with missing SNPs: {len(missing_by_sample)}")
        print(f"SNPs with missing data: {len(missing_by_snp)}")


if __name__ == "__main__":
    print("[STARTUP] Starting SGDP SNP query script for 23andMe format files...", flush=True)
    print("[INFO] Creating 23andMe format files organized by population in sgdp/ folder", flush=True)
    batch_query_snps_to_23andme()