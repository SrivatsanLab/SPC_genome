"""Count polyA/T reads per cell and (optionally) inject as obs column into an h5ad.

polyA/T rule: SEQ contains a run of >= --min-run A or T (default 10). Matches
CapGTA_dev/scripts/filter_polyat_reads.py, which is the RNA-read rule the
splicedpolyat GEX pipeline uses. SEQ includes soft-clipped bases, so
soft-clipped polyA tails are captured.

Two-step so the 6 GB h5ad is never fully rewritten:
  1. counting  — per-cell samtools calls in a process pool, results streamed to CSV
  2. injection — with --inject, write the counts as a new obs dataset via h5py
                 (only ~20 KB touched in the h5ad)

The CSV is resumable: on rerun, cells already in the CSV are skipped.

Usage:
    python count_polyat_per_cell.py <adata_h5ad> <out_csv> \\
        [--min-run 10] [--threads 16] [--samtools-threads 2] [--inject]
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


def count_one(args: tuple[str, str, int, int]) -> tuple[str, int]:
    cell_id, bam_path, min_run, sam_threads = args
    expr = f'seq =~ "A{{{min_run},}}" || seq =~ "T{{{min_run},}}"'
    cmd = ["samtools", "view", "-c", "-@", str(sam_threads), "-e", expr, bam_path]
    out = subprocess.run(cmd, capture_output=True, check=True, text=True)
    return cell_id, int(out.stdout.strip())


def load_already_done(csv_path: Path) -> dict[str, int]:
    if not csv_path.exists():
        return {}
    done: dict[str, int] = {}
    with csv_path.open() as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            done[row["cell_id"]] = int(row["polyat_reads"])
    return done


def inject_into_h5ad(h5ad: Path, counts_by_cell: dict[str, int]) -> None:
    """Append polyat_reads dataset + update column-order attr, in place."""
    with h5py.File(str(h5ad), "r+") as f:
        obs = f["obs"]
        # obs uses '_index' as its own index dataset; entries are bytes.
        idx_name = obs.attrs.get("_index", "_index")
        index = [x.decode() if isinstance(x, bytes) else x for x in obs[idx_name][()]]
        arr = np.array([counts_by_cell[c] for c in index], dtype=np.int64)

        if "polyat_reads" in obs:
            del obs["polyat_reads"]
        d = obs.create_dataset("polyat_reads", data=arr)
        d.attrs["encoding-type"] = "array"
        d.attrs["encoding-version"] = "0.2.0"

        col_order = list(obs.attrs.get("column-order", []))
        col_order = [c.decode() if isinstance(c, bytes) else c for c in col_order]
        if "polyat_reads" not in col_order:
            col_order.append("polyat_reads")
        obs.attrs["column-order"] = np.array(col_order, dtype=object)
    print(f"  injected polyat_reads into {h5ad}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("adata_h5ad")
    p.add_argument("out_csv")
    p.add_argument("--min-run", type=int, default=10)
    p.add_argument("--threads", type=int, default=16,
                   help="Number of BAMs processed in parallel")
    p.add_argument("--samtools-threads", type=int, default=2,
                   help="Threads per samtools invocation")
    p.add_argument("--inject", action="store_true",
                   help="After counting, add polyat_reads to h5ad obs in-place")
    args = p.parse_args()

    h5ad = Path(args.adata_h5ad)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading obs (backed) from {h5ad}…")
    a = ad.read_h5ad(str(h5ad), backed="r")
    if "bam_path" not in a.obs.columns:
        print("Error: adata.obs missing 'bam_path'", file=sys.stderr)
        return 1
    cells = list(a.obs.index)
    bams = list(a.obs["bam_path"])
    print(f"  {len(cells)} cells")

    done = load_already_done(out_csv)
    if done:
        print(f"  {len(done)} cells already in {out_csv} — skipping")

    todo = [(c, b, args.min_run, args.samtools_threads)
            for c, b in zip(cells, bams) if c not in done]
    print(f"  {len(todo)} cells to count")

    # Append mode; write header only if starting fresh.
    write_header = not out_csv.exists()
    with out_csv.open("a", newline="") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(["cell_id", "polyat_reads"])
        if todo:
            with ProcessPoolExecutor(max_workers=args.threads) as ex:
                futs = [ex.submit(count_one, t) for t in todo]
                for i, fut in enumerate(as_completed(futs), 1):
                    cell_id, n = fut.result()
                    done[cell_id] = n
                    w.writerow([cell_id, n])
                    if i % 100 == 0 or i == len(todo):
                        fh.flush()
                        print(f"  {i}/{len(todo)} cells counted", flush=True)

    print(f"\nWrote {out_csv} ({len(done)} rows)")
    total = sum(done.values())
    print(f"  total polyA/T reads across cells: {total:,}")
    print(f"  mean per cell: {total/len(done):,.0f}")

    if args.inject:
        missing = [c for c in cells if c not in done]
        if missing:
            print(f"Error: {len(missing)} cells missing counts — refusing to inject",
                  file=sys.stderr)
            return 1
        print(f"\nInjecting polyat_reads into {h5ad}…")
        inject_into_h5ad(h5ad, done)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
