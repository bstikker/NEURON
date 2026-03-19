#!/usr/bin/env python3

import os
import csv
import argparse

parser = argparse.ArgumentParser(
    description="Generate cumulative read lists per time bin from sequencing summary"
)
parser.add_argument("--summary", required=True, help="Path to sequencing summary TXT/TSV file")
parser.add_argument("--sample", required=True, help="Sample name, e.g. TVU19-13652")
parser.add_argument("--barcode-label", required=True,
                    help='Exact barcode field match, e.g. "1-TVU19-13652-barcode01"')
parser.add_argument("--output", required=True, help="Output directory for read lists")
parser.add_argument("--bins", nargs="+", type=int,
                    default=[15, 30, 45, 60, 90, 120, 180, 240, 360, 480],
                    help="Cumulative time bins in minutes")
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

def detect_delimiter(header_line: str):
    return "\t" if "\t" in header_line else None

records = []

with open(args.summary, "r", newline="") as f:
    first_line = f.readline()
    if not first_line:
        raise ValueError("Sequencing summary file is empty.")

    delimiter = detect_delimiter(first_line)
    f.seek(0)

    if delimiter is not None:
        reader = csv.DictReader(f, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        required = ["read_id", "start_time", "barcode"]
        missing = [c for c in required if c not in fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        for row in reader:
            if row["barcode"] != args.barcode_label:
                continue
            try:
                start_time = float(row["start_time"])
            except (TypeError, ValueError):
                continue
            records.append((row["read_id"], start_time))
    else:
        header = first_line.strip().split()
        col_index = {name: i for i, name in enumerate(header)}

        required = ["read_id", "start_time", "barcode"]
        missing = [c for c in required if c not in col_index]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        read_id_idx = col_index["read_id"]
        start_time_idx = col_index["start_time"]
        barcode_idx = col_index["barcode"]

        for line in f:
            parts = line.strip().split()
            if len(parts) <= max(read_id_idx, start_time_idx, barcode_idx):
                continue
            if parts[barcode_idx] != args.barcode_label:
                continue
            try:
                start_time = float(parts[start_time_idx])
            except ValueError:
                continue
            records.append((parts[read_id_idx], start_time))

if not records:
    raise ValueError(f'No reads found for barcode label "{args.barcode_label}".')

records.sort(key=lambda x: x[1])
t0 = records[0][1]
records = [(rid, (st - t0) / 60.0) for rid, st in records]

counts_rows = []

for t in args.bins:
    out_file = os.path.join(args.output, f"{args.sample}_read_ids_{t}min.txt")
    n = 0
    with open(out_file, "w") as out:
        for rid, time_min in records:
            if time_min <= t:
                out.write(rid + "\n")
                n += 1
            else:
                break
    counts_rows.append((args.sample, t, n, out_file))
    print(f"Wrote {n} reads for {t} min -> {out_file}")

counts_tsv = os.path.join(args.output, f"{args.sample}_read_counts.tsv")
with open(counts_tsv, "w", newline="") as out:
    out.write("sample\ttime_min\testimated_reads\tread_list\n")
    for sample, time_min, nreads, read_list in counts_rows:
        out.write(f"{sample}\t{time_min}\t{nreads}\t{read_list}\n")

print(f"Wrote read counts summary -> {counts_tsv}")
