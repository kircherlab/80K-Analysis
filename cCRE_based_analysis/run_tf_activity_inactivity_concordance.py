#!/usr/bin/env python
"""
Concordance pipeline for TFBS analyses.

Ridge logistic regression per TF on all ~24k tested cCREs:
    TF_bound (binary) ~ normalized_activity_NGN2 + tf_density_excluding_TF + gc_content

Outputs:
  - tf_activity_ridge_full.tsv.gz         (every tested TF, ridge bootstrap results)
  - tf_activity_fisher_full.tsv.gz        (every TF, Fisher active vs inactive)
  - supplementary_table_concordance.tsv   (the focused subset for the figure)
  - figure_concordance.svg / .png         (Extended Data Fig. 4 replacement)

Interpretation of ridge output:
  beta_activity > 0, p_adj < threshold → TF more likely bound as activity rises
  beta_activity < 0, p_adj < threshold → TF less likely bound as activity rises
"""

import argparse
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from joblib import Parallel, delayed
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests


# ----------------------------- CLI ----------------------------- #

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--active_counts", required=True,
                   help="TFBS counts for active cCREs. Long format. If --counts_no_header, "
                        "expects 3 whitespace-separated columns: count, name, TF_name. "
                        "Otherwise expects a header row with those names.")
    p.add_argument("--inactive_counts", required=True,
                   help="TFBS counts for inactive (not-significant) cCREs (same format).")
    p.add_argument("--counts_no_header", action="store_true",
                   help="Set if count files have no header (raw awk output).")
    p.add_argument("--annotation", required=True,
                   help="TSV with at least: name, normalized_activity_NGN2, sequence")
    p.add_argument("--elastic_net", required=True,
                   help="TSV with elastic-net coefficients: feature, coefficient")
    p.add_argument("--outdir", required=True)
    # Ridge / bootstrap knobs
    p.add_argument("--n_boot", type=int, default=1000)
    p.add_argument("--C", type=float, default=1.0)
    p.add_argument("--min_positive", type=int, default=5)
    p.add_argument("--min_usable_boot", type=int, default=50)
    p.add_argument("--n_jobs", type=int, default=-1)
    p.add_argument("--seed", type=int, default=42)
    # Focused-set filter knobs (must match what is described in the results section)
    p.add_argument("--fisher_fdr", type=float, default=0.05,
                   help="BH FDR cutoff for the Fisher 'differential binding' set (default 0.05)")
    p.add_argument("--or_enrich", type=float, default=2.0,
                   help="Fisher OR cutoff for 'enriched' (default 2.0)")
    p.add_argument("--or_deplete", type=float, default=0.8,
                   help="Fisher OR cutoff for 'depleted' (default 0.8)")
    p.add_argument("--ridge_fdr", type=float, default=0.05,
                   help="BH FDR cutoff for ridge significance flag in supplementary table (default 0.05)")
    p.add_argument("--normalized_activity_col", default="normalized_activity_NGN2", help="Column name in annotation file for the normalized NGN2 activity (default 'normalized_activity_NGN2')")
    return p.parse_args()


# ----------------------------- Data prep ----------------------------- #

def load_long_counts(path, category, no_header):
    if no_header:
        df = pd.read_csv(path, sep=r"\s+", header=None,
                         names=["count", "name", "TF_name"])
    else:
        df = pd.read_csv(path, sep="\t")
    df["activity_category"] = category
    return df


def build_wide_tfbs_table(active_path, inactive_path, no_header):
    df_a = load_long_counts(active_path, "active_element", no_header)
    df_b = load_long_counts(inactive_path, "not_significant_element", no_header)
    combined = pd.concat([df_a, df_b], ignore_index=True)
    wide = pd.pivot_table(
        combined,
        index=["name", "activity_category"],
        columns="TF_name",
        values="count",
        fill_value=0,
    ).reset_index()
    return wide


def compute_gc_content(df, seq_col="sequence"):
    seq = df[seq_col].astype(str).str.upper()
    df["gc_content"] = (seq.str.count("G") + seq.str.count("C")) / seq.str.len()
    return df


# ----------------------------- Fisher (active vs inactive) ----------------------------- #

def run_fisher_active_vs_inactive(wide_df, pseudocount=1):
    """Same logic as test_tf_enrichment, but stripped to what we need downstream."""
    a = wide_df[wide_df["activity_category"] == "active_element"]
    b = wide_df[wide_df["activity_category"] == "not_significant_element"]
    n_a, n_b = len(a), len(b)
    tf_cols = [c for c in wide_df.columns
               if c not in ("name", "activity_category")]
    out = []
    for tf in tf_cols:
        ka = (a[tf] > 0).sum()
        kb = (b[tf] > 0).sum()
        table = [
            [ka + pseudocount, (n_a - ka) + pseudocount],
            [kb + pseudocount, (n_b - kb) + pseudocount],
        ]
        odds, p = fisher_exact(table)
        out.append({
            "TF": tf, "fisher_OR": odds, "fisher_p": p,
            "n_pos_active": int(ka), "n_pos_inactive": int(kb),
        })
    res = pd.DataFrame(out)
    res["fisher_p_adj"] = multipletests(res["fisher_p"], method="fdr_bh")[1]
    return res


# ----------------------------- Ridge bootstrap (activity predictor) ----------------------------- #

def bootstrap_single_tf(activity, density_total, gc, tf_counts,
                        n_boot, C, min_positive, min_usable_boot, seed):
    y_full = (tf_counts > 0).astype(int)
    if y_full.sum() < min_positive:
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0

    n = len(tf_counts)
    rng = np.random.default_rng(seed)
    betas = []

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        y = (tf_counts[idx] > 0).astype(int)
        if len(np.unique(y)) < 2:
            continue
        tf_density_excl = density_total[idx] - tf_counts[idx]
        X = np.column_stack([activity[idx], tf_density_excl, gc[idx]])
        model = LogisticRegression(penalty="l2", C=C, solver="lbfgs", max_iter=3000)
        model.fit(X, y)
        # coef order matches X column order: activity is first
        betas.append(model.coef_[0][0])

    n_usable = len(betas)
    if n_usable < min_usable_boot:
        return np.nan, np.nan, np.nan, np.nan, np.nan, n_usable

    betas = np.array(betas)
    beta_mean = betas.mean()
    OR_mean = np.exp(beta_mean)
    CI_l = np.percentile(betas, 2.5)
    CI_u = np.percentile(betas, 97.5)
    p_emp = 2 * min(np.mean(betas <= 0), np.mean(betas >= 0))
    p_emp = max(p_emp, 1.0 / n_usable)
    return beta_mean, OR_mean, CI_l, CI_u, p_emp, n_usable


def run_ridge_bootstrap(df, tf_cols, n_boot, C, min_positive,
                        min_usable_boot, n_jobs, seed, normalized_activity_col):
    activity = df[normalized_activity_col].values
    density_total = df["tf_density_total"].values
    gc = df["gc_content"].values

    rng = np.random.default_rng(seed)
    tf_seeds = rng.integers(0, 2**31 - 1, size=len(tf_cols))

    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] Bootstrapping {len(tf_cols)} TFs "
          f"x {n_boot} reps on {n_jobs} workers (n_cCRE = {len(df)})...",
          flush=True)

    out = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(bootstrap_single_tf)(
            activity, density_total, gc, df[tf].values,
            n_boot, C, min_positive, min_usable_boot, int(s),
        )
        for tf, s in zip(tf_cols, tf_seeds)
    )

    print(f"[{time.strftime('%H:%M:%S')}] Done in {time.time() - t0:.1f}s",
          flush=True)

    rows = []
    for tf, (beta_mean, OR_mean, CI_l, CI_u, p_emp, n_usable) in zip(tf_cols, out):
        rows.append({
            "TF": tf,
            "beta_activity": beta_mean,
            "OR_per_unit_activity": OR_mean,
            "beta_CI_lower": CI_l,
            "beta_CI_upper": CI_u,
            "OR_CI_lower": np.exp(CI_l) if not np.isnan(CI_l) else np.nan,
            "OR_CI_upper": np.exp(CI_u) if not np.isnan(CI_u) else np.nan,
            "p_empirical": p_emp,
            "n_usable_boot": n_usable,
        })
    res = pd.DataFrame(rows)
    valid = res["p_empirical"].notna()
    res.loc[valid, "p_adj"] = multipletests(
        res.loc[valid, "p_empirical"], method="fdr_bh"
    )[1]
    return res.sort_values("p_adj")


# ----------------------------- Plot ----------------------------- #

def make_concordance_plot(df, out_svg, out_png, ridge_fdr=0.05):
    """
    Replaces the Spearman x-axis with the bootstrap-ridge log-OR-per-unit-activity.
    Closed markers = ridge BH FDR < ridge_fdr; open = not significant on ridge side.
    """
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sig = df["p_adj"] < ridge_fdr
    ax.scatter(df.loc[sig, "beta_activity"], df.loc[sig, "coefficient"],
               s=70, facecolors="#3b6fb6", edgecolors="black", linewidths=0.5,
               label=f"Ridge BH FDR < {ridge_fdr}")
    ax.scatter(df.loc[~sig, "beta_activity"], df.loc[~sig, "coefficient"],
               s=70, facecolors="white", edgecolors="#3b6fb6", linewidths=1.2,
               label=f"Ridge BH FDR ≥ {ridge_fdr}")

    # CI bars (light, only for plotted points)
    ax.errorbar(df["beta_activity"], df["coefficient"],
                xerr=[df["beta_activity"] - df["beta_CI_lower"],
                      df["beta_CI_upper"] - df["beta_activity"]],
                fmt="none", ecolor="grey", alpha=0.4, capsize=0, linewidth=0.8)

    ax.axhline(0, color="grey", linestyle="--", alpha=0.5)
    ax.axvline(0, color="grey", linestyle="--", alpha=0.5)
    ax.set_xlabel("Ridge β (log-OR of TF presence per unit NGN2 activity,\n"
                  "adjusted for GC and TF density)")
    ax.set_ylabel("Elastic-net coefficient (NGN2 activity model)")

    texts = [ax.text(r["beta_activity"], r["coefficient"], r["TF"], fontsize=9)
             for _, r in df.iterrows()]
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_svg, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


# ----------------------------- Main ----------------------------- #

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Build wide TFBS table from the two count files
    print(f"[{time.strftime('%H:%M:%S')}] Loading TFBS counts...", flush=True)
    wide = build_wide_tfbs_table(
        args.active_counts, args.inactive_counts, args.counts_no_header,
    )
    print(f"  cCREs in counts: {len(wide)}")

    # 2. Fisher: active vs inactive (the 'differential binding' set, full table)
    print(f"[{time.strftime('%H:%M:%S')}] Running Fisher...", flush=True)
    fisher = run_fisher_active_vs_inactive(wide)
    fisher_path = os.path.join(args.outdir, "tf_activity_fisher_full.tsv.gz")
    fisher.to_csv(fisher_path, sep="\t", index=False, compression="gzip")
    print(f"  → {fisher_path}")

    # 3. Merge TFBS counts with annotation to get activity + sequence
    print(f"[{time.strftime('%H:%M:%S')}] Loading annotation...", flush=True)
    ann = pd.read_csv(args.annotation, sep="\t")
    needed = ["name", args.normalized_activity_col, "sequence"]
    missing = [c for c in needed if c not in ann.columns]
    if missing:
        raise SystemExit(f"annotation file missing columns: {missing}")
    ann = ann[needed]

    df = wide.merge(ann, on="name", how="inner")
    n_before = len(df)
    df = df.dropna(subset=[args.normalized_activity_col, "sequence"])
    print(f"  merged cCREs: {len(df)} "
          f"(dropped {n_before - len(df)} with missing activity/sequence)")

    df = compute_gc_content(df)

    tf_cols = [c for c in df.columns
               if c not in ("name", "activity_category", "sequence",
                            args.normalized_activity_col, "gc_content")]
    df[tf_cols] = df[tf_cols].astype(float)
    df["tf_density_total"] = df[tf_cols].sum(axis=1)
    print(f"  TFs to test: {len(tf_cols)}")

    # 4. Ridge bootstrap regression on activity (the expensive step)
    ridge = run_ridge_bootstrap(
        df, tf_cols,
        n_boot=args.n_boot, C=args.C,
        min_positive=args.min_positive,
        min_usable_boot=args.min_usable_boot,
        n_jobs=args.n_jobs, seed=args.seed, normalized_activity_col=args.normalized_activity_col,
    )
    ridge_path = os.path.join(args.outdir, "tf_activity_ridge_full.tsv.gz")
    ridge.to_csv(ridge_path, sep="\t", index=False, compression="gzip")
    print(f"  → {ridge_path}")

    # 5. Build the focused supplementary table = Fisher hits ∩ non-zero EN coefs
    print(f"[{time.strftime('%H:%M:%S')}] Building concordance table...", flush=True)
    en = pd.read_csv(args.elastic_net, sep="\t")[["feature", "coefficient"]]\
        .rename(columns={"feature": "TF"})
    en_nz = en[en["coefficient"] != 0]
    print(f"  Elastic-net features with non-zero coefficients: {len(en_nz)}")

    fisher_hit = fisher[
        (fisher["fisher_p_adj"] < args.fisher_fdr)
        & ((fisher["fisher_OR"] > args.or_enrich)
           | (fisher["fisher_OR"] < args.or_deplete))
    ]
    print(f"  Fisher hits (FDR<{args.fisher_fdr}, OR>{args.or_enrich} or <{args.or_deplete}): "
          f"{len(fisher_hit)}")

    focused = (fisher_hit
               .merge(en_nz, on="TF", how="inner")
               .merge(ridge, on="TF", how="left"))
    print(f"  Focused set (intersection): {len(focused)} TFs")

    # Tidy columns / ordering for the supplementary table
    focused = focused[[
        "TF",
        "fisher_OR", "fisher_p_adj", "n_pos_active", "n_pos_inactive",
        "coefficient",
        "beta_activity", "OR_per_unit_activity",
        "beta_CI_lower", "beta_CI_upper", "OR_CI_lower", "OR_CI_upper",
        "p_empirical", "p_adj", "n_usable_boot",
    ]].rename(columns={
        "fisher_OR": "fisher_OR_active_vs_inactive",
        "fisher_p_adj": "fisher_BH_p_adj",
        "coefficient": "elastic_net_coefficient",
        "p_empirical": "ridge_p_empirical",
        "p_adj": "ridge_BH_p_adj",
    })

    sup_path = os.path.join(args.outdir, "supplementary_table_concordance.tsv")
    focused.to_csv(sup_path, sep="\t", index=False)
    print(f"  → {sup_path}")

    # 6. Plot
    out_svg = os.path.join(args.outdir, "figure_concordance.svg")
    out_png = os.path.join(args.outdir, "figure_concordance.png")
    plot_df = focused.rename(columns={
        "elastic_net_coefficient": "coefficient",
        "ridge_BH_p_adj": "p_adj",
    })
    make_concordance_plot(plot_df, out_svg, out_png, ridge_fdr=args.ridge_fdr)
    print(f"  → {out_svg} / {out_png}")

    # 7. Summary
    print("\n=== Summary ===")
    print(f"  TFs in Fisher table:    {len(fisher)}")
    print(f"  TFs in Ridge table:     {ridge['p_empirical'].notna().sum()} tested "
          f"({ridge['p_empirical'].isna().sum()} dropped via min_positive/min_usable_boot)")
    print(f"  Fisher hits:            {len(fisher_hit)}")
    print(f"  Focused (∩ EN non-zero): {len(focused)}")
    n_concordant = ((plot_df["beta_activity"] * plot_df["coefficient"]) > 0).sum()
    print(f"  Same-sign (β_activity & EN coef): {n_concordant} / {len(plot_df)}")
    sig_ridge = (plot_df["p_adj"] < args.ridge_fdr).sum()
    print(f"  Ridge BH-significant at FDR<{args.ridge_fdr}: {sig_ridge} / {len(plot_df)}")


if __name__ == "__main__":
    sys.exit(main())
