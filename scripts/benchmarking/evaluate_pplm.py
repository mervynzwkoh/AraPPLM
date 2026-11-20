#!/usr/bin/env python3
"""
Evaluate PPLM Predictions against DeepAraPPI & ESMAraPPI Benchmarks.

Computes primary metrics matching DeepAraPPI/ESMAraPPI and PPLM's PPI_model.py
(AUPRC, AUC_PR via precision-recall curve trapezoidal auc, AUROC, Precision,
Recall, F1, MCC, Specificity, Top10/20/50), and supports summarizing the
6 pooling combinations (mean, max, min, combine_mean_max, combine_mean_max_min,
combine_mean_min).

Usage:
    # Evaluate multiple DeepAraPPI task predictions at once:
    python scripts/benchmarking/evaluate_pplm.py \
        --c2 results/DeepAraPPI/pretrained/deepara_c2_combine_mean_max_scores.csv \
        --c3 results/DeepAraPPI/pretrained/deepara_c3_combine_mean_max_scores.csv \
        --rice results/DeepAraPPI/pretrained/deepara_rice_combine_mean_max_scores.csv \
        --output_dir results/DeepAraPPI/pretrained/

    # Or summarize all 6 pooling combinations in a directory:
    python scripts/benchmarking/evaluate_pplm.py \
        --pooling_ablation_dir results/DeepAraPPI/pretrained/ \
        --suite DeepAraPPI \
        --output_dir results/DeepAraPPI/pretrained/
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    auc,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
)


POOLING_COMBOS = [
    "mean",
    "max",
    "min",
    "combine_mean_max",
    "combine_mean_max_min",
    "combine_mean_min",
]

# Reference benchmark baselines from DeepAraPPI (Zheng et al., 2023)
DEEPARAPPI_BASELINES = {
    "Task1": {
        "Description": "Random 80/20 Split (Low Difficulty)",
        "RCNN": 0.925,
        "Domain2vec": 0.868,
        "GO2vec": 0.939,
        "DeepAraPPI": 0.965,
    },
    "Task2": {
        "Description": "C2 Test Set: One Unseen Protein (Medium Difficulty)",
        "RCNN": 0.746,
        "Domain2vec": 0.780,
        "GO2vec": 0.871,
        "DeepAraPPI": 0.897,
    },
    "Task3": {
        "Description": "C3 Test Set: Both Unseen Proteins (High Difficulty)",
        "RCNN": 0.481,
        "Domain2vec": 0.681,
        "GO2vec": 0.803,
        "DeepAraPPI": 0.825,
    },
    "Task4_Rice": {
        "Description": "Oryza sativa Cross-Species Generalization",
        "RCNN": 0.248,
        "Domain2vec": 0.279,
        "GO2vec": 0.265,
        "DeepAraPPI": 0.305,
    },
}

ARAPPLM_TASKS = {
    "Task1": {"Description": "AraPPLM C1 Train Partition (theta=0.50, c_min=0.80)"},
    "Task2": {"Description": "AraPPLM C2 Test Set: One Unseen Protein"},
    "Task3": {"Description": "AraPPLM C3 Test Set: Both Unseen Proteins (Strict Zero-Shot)"},
    "Task4_C4": {"Description": "AraPPLM C4 Test Set: Both Seen Proteins, Novel Test Edge"},
}

# Reference benchmark baselines from ESMAraPPI (Zhou et al., 2023)
ESMARAPPI_BASELINES = {
    "ESM_C2": {
        "Description": "ESMAraPPI C2 Test Set (One Unseen Protein, 37,444 pairs)",
        "ESMAraPPI": 0.834,
        "DeepAraPPI": 0.871,
        "TAGPPI": 0.700,
        "PIPR": 0.588,
        "RAPPPID": 0.516,
        "D-SCRIPT": 0.292,
    },
    "ESM_C3": {
        "Description": "ESMAraPPI C3 Test Set (Both Unseen Proteins, 8,866 pairs)",
        "ESMAraPPI": 0.810,
        "DeepAraPPI": 0.785,
        "TAGPPI": 0.554,
        "PIPR": 0.387,
        "RAPPPID": 0.371,
        "D-SCRIPT": 0.291,
    },
}


def compute_metrics(y_true, y_score, threshold=0.5):
    """Compute comprehensive classification, Top-K, and ranking metrics."""
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    # Top-K precision (matching PPI_model.py::evaluate)
    idx_list = np.argsort(-y_score)
    gt_sorted = y_true[idx_list]
    top10 = float(np.sum(gt_sorted[:10]) / min(10, len(gt_sorted))) if len(gt_sorted) else 0.0
    top20 = float(np.sum(gt_sorted[:20]) / min(20, len(gt_sorted))) if len(gt_sorted) else 0.0
    top50 = float(np.sum(gt_sorted[:50]) / min(50, len(gt_sorted))) if len(gt_sorted) else 0.0

    # Ranking metrics (Threshold-independent)
    auprc = float(average_precision_score(y_true, y_score))
    precision_arr, recall_arr, thresh_arr = precision_recall_curve(y_true, y_score)
    auc_pr_trapz = float(auc(recall_arr, precision_arr))

    try:
        auroc = float(roc_auc_score(y_true, y_score))
    except ValueError:
        auroc = float("nan")

    # Threshold-dependent metrics (matching > thresh in PPI_model.py)
    y_pred = (y_score > threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    accuracy = float(accuracy_score(y_true, y_pred))
    precision = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    sensitivity = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    specificity = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    f1 = float(f1_score(y_true, y_pred, zero_division=0))
    mcc = float(matthews_corrcoef(y_true, y_pred))

    # Optimal threshold maximizing F1 score
    f1_scores = np.divide(
        2 * (precision_arr * recall_arr),
        (precision_arr + recall_arr),
        out=np.zeros_like(precision_arr),
        where=(precision_arr + recall_arr) != 0,
    )
    best_idx = int(np.argmax(f1_scores))
    best_f1 = float(f1_scores[best_idx])
    best_thresh = float(thresh_arr[best_idx]) if best_idx < len(thresh_arr) else 0.5

    return {
        "AUPRC": auprc,
        "AUC_PR_Trapz": auc_pr_trapz,
        "AUROC": auroc,
        "Accuracy": accuracy,
        "Precision": precision,
        "Sensitivity_Recall": sensitivity,
        "Specificity": specificity,
        "F1_Score": f1,
        "MCC": mcc,
        "Top10": top10,
        "Top20": top20,
        "Top50": top50,
        "Best_Threshold": best_thresh,
        "Best_F1": best_f1,
        "Total_Samples": len(y_true),
        "Positive_Samples": int(np.sum(y_true)),
        "Negative_Samples": int(len(y_true) - np.sum(y_true)),
    }


def print_task_evaluation(task_name, metrics, baseline_dict=None):
    """Print formatted evaluation report."""
    print("=" * 70)
    print(f"EVALUATION REPORT: {task_name.upper()}")
    if baseline_dict and task_name in baseline_dict:
        print(f"Description: {baseline_dict[task_name].get('Description', '')}")
    print("=" * 70)

    print(
        f"Dataset Size:     {metrics['Total_Samples']:,} pairs "
        f"({metrics['Positive_Samples']:,} positive, {metrics['Negative_Samples']:,} negative | "
        f"Ratio 1:{metrics['Negative_Samples'] / max(1, metrics['Positive_Samples']):.1f})"
    )
    print("-" * 70)
    print(f"  [*] AUPRC (Avg Precision):  {metrics['AUPRC']:.4f}")
    print(f"  [*] AUC_PR (PR-Curve AUC):  {metrics['AUC_PR_Trapz']:.4f}")
    print(f"  [*] AUROC:                  {metrics['AUROC']:.4f}")
    print(f"  F1 Score (thresh=0.5):      {metrics['F1_Score']:.4f}")
    print(f"  Precision:                  {metrics['Precision']:.4f}")
    print(f"  Sensitivity (Recall):       {metrics['Sensitivity_Recall']:.4f}")
    print(f"  Specificity:                {metrics['Specificity']:.4f}")
    print(f"  MCC:                        {metrics['MCC']:.4f}")
    print(f"  Accuracy:                   {metrics['Accuracy']:.4f}")
    print(
        f"  Top-10 / Top-20 / Top-50:   "
        f"{metrics['Top10']:.4f} / {metrics['Top20']:.4f} / {metrics['Top50']:.4f}"
    )
    print(
        f"  Optimal F1:                 {metrics['Best_F1']:.4f} "
        f"(at threshold {metrics['Best_Threshold']:.4f})"
    )

    if baseline_dict and task_name in baseline_dict:
        b = baseline_dict[task_name]
        model_entries = [(k, v) for k, v in b.items() if k != "Description"]
        if model_entries:
            print("-" * 70)
            benchmark_suite = "ESMARAPPI" if "ESM" in task_name else "DEEPARAPPI"
            print(f"COMPARISON WITH {benchmark_suite} BENCHMARKS (AUPRC):")
            print(f"  PPLM:                       {metrics['AUPRC']:.4f}")
            for model_k, val in model_entries:
                print(f"  {model_k:<27} {val:.4f}")

            primary_comp = "ESMAraPPI" if "ESMAraPPI" in b else "DeepAraPPI"
            if primary_comp in b:
                delta = metrics["AUPRC"] - b[primary_comp]
                sign = "+" if delta >= 0 else ""
                print(f"  PPLM vs {primary_comp}:         {sign}{delta:.4f}")
    print("=" * 70 + "\n")


def load_predictions_csv(csv_path):
    """Load prediction CSV supporting both (pred_score, true_label) and (prediction, label)."""
    df = pd.read_csv(csv_path)
    if "pred_score" in df.columns and "true_label" in df.columns:
        return df["true_label"].values, df["pred_score"].values
    if "prediction" in df.columns and "label" in df.columns:
        return df["label"].values, df["prediction"].values
    raise ValueError(
        f"CSV {csv_path} must contain ('pred_score', 'true_label') or ('prediction', 'label') columns."
    )


def evaluate_file(
    csv_path, task_name="Task", output_dir=None, summary_filename=None, suite=None
):
    """Load predictions CSV and compute metrics."""
    if not os.path.exists(csv_path):
        print(f"[!] Error: File not found: {csv_path}")
        return None

    y_true, y_score = load_predictions_csv(csv_path)
    metrics = compute_metrics(y_true, y_score)
    if suite == "AraPPLM":
        baselines = ARAPPLM_TASKS
    else:
        baselines = {**DEEPARAPPI_BASELINES, **ESMARAPPI_BASELINES, **ARAPPLM_TASKS}
    print_task_evaluation(task_name, metrics, baselines)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_name = summary_filename or f"{task_name.lower()}_metrics.txt"
        out_summary_path = os.path.join(output_dir, out_name)
        with open(out_summary_path, "w") as f:
            f.write(f"Task: {task_name}\n")
            for k, v in metrics.items():
                f.write(f"{k}: {v}\n")
        print(f"Saved metrics summary to: {out_summary_path}")

    return metrics


def run_pooling_ablation_summary(ablation_dir, suite, output_dir):
    """Discover and summarize all 6 pooling combinations across tasks."""
    if suite == "DeepAraPPI":
        task_specs = [
            ("c2", "Task2", "deepara"),
            ("c3", "Task3", "deepara"),
            ("rice", "Task4_Rice", "deepara"),
        ]
    elif suite == "AraPPLM":
        task_specs = [
            ("c1", "Task1", "araplm"),
            ("c2", "Task2", "araplm"),
            ("c3", "Task3", "araplm"),
            ("c4", "Task4_C4", "araplm"),
        ]
    else:
        task_specs = [
            ("c2", "ESM_C2", "esmarappi"),
            ("c3", "ESM_C3", "esmarappi"),
        ]

    ablation_rows = []
    for pooling in POOLING_COMBOS:
        for split_tag, task_key, prefix in task_specs:
            csv_path = os.path.join(ablation_dir, f"{prefix}_{split_tag}_{pooling}_scores.csv")
            if not os.path.exists(csv_path):
                continue
            y_true, y_score = load_predictions_csv(csv_path)
            m = compute_metrics(y_true, y_score)
            ablation_rows.append({"Pooling": pooling, "Task": task_key, **m})

    if not ablation_rows:
        print(f"[!] No pooling ablation CSVs found in {ablation_dir}")
        return

    df_ablation = pd.DataFrame(ablation_rows)
    os.makedirs(output_dir, exist_ok=True)
    ablation_csv = os.path.join(output_dir, "pooling_ablation_summary.csv")
    df_ablation.to_csv(ablation_csv, index=False)

    print("\n" + "=" * 80)
    print(f"MULTI-POOLING ABLATION SUMMARY ({suite})")
    print("=" * 80)
    cols = ["Pooling", "Task", "AUPRC", "AUC_PR_Trapz", "AUROC", "F1_Score", "MCC"]
    print(df_ablation[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\n[INFO] Pooling ablation summary saved to: {ablation_csv}")
    print("=" * 80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate PPLM on Plant PPI Benchmarks (DeepAraPPI / ESMAraPPI / AraPPLM)"
    )
    parser.add_argument("--input", help="Single prediction CSV file")
    parser.add_argument(
        "--task",
        default="Task",
        choices=[
            "Task1",
            "Task2",
            "Task3",
            "Task4_Rice",
            "Task4_C4",
            "ESM_C2",
            "ESM_C3",
            "Custom",
        ],
        help="Task name for benchmark lookup",
    )

    # DeepAraPPI / AraPPLM multi-file evaluation
    parser.add_argument("--c1", help="Path to C1 prediction CSV")
    parser.add_argument("--c2", help="Path to C2 test prediction CSV")
    parser.add_argument("--c3", help="Path to C3 test prediction CSV")
    parser.add_argument("--c4", help="Path to AraPPLM C4 test prediction CSV")
    parser.add_argument("--rice", help="Path to DeepAraPPI Rice test prediction CSV")

    # ESMAraPPI multi-file evaluation
    parser.add_argument("--esm_c2", help="Path to ESMAraPPI C2 prediction CSV")
    parser.add_argument("--esm_c3", help="Path to ESMAraPPI C3 prediction CSV")

    # Multi-pooling ablation summary
    parser.add_argument(
        "--pooling_ablation_dir",
        default=None,
        help="Directory containing {prefix}_{split}_{pooling}_scores.csv files across 6 pooling configs",
    )
    parser.add_argument(
        "--suite",
        default="DeepAraPPI",
        choices=["DeepAraPPI", "ESMAraPPI", "AraPPLM"],
        help="Benchmark suite for --pooling_ablation_dir and baseline selection",
    )

    parser.add_argument(
        "--output_dir", default="results", help="Directory to save metric summary reports"
    )
    args = parser.parse_args()

    results_table = []

    if args.input:
        m = evaluate_file(
            args.input, task_name=args.task, output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": args.task, **m})

    # DeepAraPPI / AraPPLM evaluations
    if args.c1:
        m = evaluate_file(
            args.c1, task_name="Task1", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "Task1 (C1)", **m})

    if args.c2:
        m = evaluate_file(
            args.c2, task_name="Task2", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "Task2 (C2)", **m})

    if args.c3:
        m = evaluate_file(
            args.c3, task_name="Task3", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "Task3 (C3)", **m})

    if args.c4:
        m = evaluate_file(
            args.c4, task_name="Task4_C4", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "Task4 (C4)", **m})

    if args.rice:
        m = evaluate_file(
            args.rice, task_name="Task4_Rice", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "Task4 (Rice)", **m})

    # ESMAraPPI evaluations
    if args.esm_c2:
        m = evaluate_file(
            args.esm_c2, task_name="ESM_C2", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "ESMAraPPI C2 (One Unseen)", **m})

    if args.esm_c3:
        m = evaluate_file(
            args.esm_c3, task_name="ESM_C3", output_dir=args.output_dir, suite=args.suite
        )
        if m:
            results_table.append({"Task": "ESMAraPPI C3 (Both Unseen)", **m})

    if len(results_table) > 1:
        df_summary = pd.DataFrame(results_table)
        summary_csv = os.path.join(args.output_dir, "benchmark_summary.csv")
        df_summary.to_csv(summary_csv, index=False)
        print(f"\n[INFO] Consolidated Benchmark Summary saved to: {summary_csv}")

    if args.pooling_ablation_dir:
        run_pooling_ablation_summary(
            args.pooling_ablation_dir, suite=args.suite, output_dir=args.output_dir
        )


if __name__ == "__main__":
    main()
