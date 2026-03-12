#!/usr/bin/env python3
import argparse, csv, math

def cum_sums(xs):
    s=0; out=[]
    for x in xs:
        s+=x; out.append(s)
    return out

def build_predictor():
    cpg_intervals = [51924,104073,124078,149111,173504,194399,207456,217193,232101,241278,247600,258197]
    time_intervals = [5,10,15,20,25,30,35,40,45,50,55,60]
    cum_cpg = cum_sums(cpg_intervals)
    times   = time_intervals

    def predict(cpg):
        if cpg <= 0:
            return 0.0
        if cpg <= cum_cpg[0]:
            slope = times[0] / cum_cpg[0]
            return slope * cpg
        for i in range(1, len(cum_cpg)):
            if cpg <= cum_cpg[i]:
                x0, x1 = cum_cpg[i-1], cum_cpg[i]
                y0, y1 = times[i-1],   times[i]
                return y0 + (y1 - y0) * ( (cpg - x0) / (x1 - x0) )
        x0, x1 = cum_cpg[-2], cum_cpg[-1]
        y0, y1 = times[-2],   times[-1]
        slope = (y1 - y0) / (x1 - x0)
        return y1 + slope * (cpg - x1)

    return predict

def read_rows(path):
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            yield row

def write_out(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--threshold", type=float, default=0.8)
    ap.add_argument("--stable_n", type=int, default=4, help="consecutive points")
    args = ap.parse_args()

    predict_time = build_predictor()

    # Collect confident rows in order
    conf_rows = []
    for r in read_rows(args.in_csv):
        try:
            if float(r["confidence"]) >= args.threshold:
                conf_rows.append(r)
        except (KeyError, ValueError):
            continue

    out = []
    if len(conf_rows) >= args.stable_n:
        # First window of N confident points
        idx = args.stable_n - 1
        covered = float(conf_rows[idx]["covered_cpg"])
        probes  = float(conf_rows[idx].get("number_probes", "nan"))
        tmin    = predict_time(covered)
        out.append({
            "label": args.label,
            "covered_cpg": covered,
            "number_probes": probes,
            "time_min": round(tmin, 3),
        })

    write_out(args.out_csv, out, ["label","covered_cpg","number_probes","time_min"])

if __name__ == "__main__":
    main()
