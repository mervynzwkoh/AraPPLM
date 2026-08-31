#!/usr/bin/env python3
"""
Evaluate PPLM Predictions against DeepAraPPI Benchmarks

Computes primary metrics matching DeepAraPPI (AUPRC, AUROC, Precision, Recall, F1, MCC, Specificity)
and outputs a comparison against the DeepAraPPI paper baselines (RCNN, Domain2vec, GO2vec, DeepAraPPI).

Usage:
    # Evaluate a single task:
    python scripts/evaluate_pplm.py \
        --input results/deepara_c2_scores.csv \
        --task Task2 \
        --output_dir results/

    # Evaluate multiple task predictions at once:
    python scripts/evaluate_pplm.py \
        --c1 results/deepara_c1_scores.csv \
        --c2 results/deepara_c2_scores.csv \
        --c3 results/deepara_c3_scores.csv \
        --rice results/deepara_rice_scores.csv \
        --output_dir results/
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    precision_recall_curve,
    roc_curve,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    accuracy_score
)

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
    """Compute comprehensive classification and ranking metrics."""
    # Ranking metrics (Threshold-independent)
    auprc = average_precision_score(y_true, y_score)
    try:
        auroc = roc_auc_score(y_true, y_score)
    except ValueError:
        auroc = float("nan")

    # Threshold-dependent metrics
    y_pred = (y_score >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    accuracy = accuracy_score(y_true, y_pred)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0  # Recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    f1 = f1_score(y_true, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_true, y_pred)

    # Find optimal threshold that maximizes F1 score
    precision_arr, recall_arr, thresh_arr = precision_recall_curve(y_true, y_score)
    f1_scores = np.divide(
        2 * (precision_arr * recall_arr),
        (precision_arr + recall_arr),
        out=np.zeros_like(precision_arr),
        where=(precision_arr + recall_arr) != 0,
    )
    best_idx = np.argmax(f1_scores)
    best_f1 = f1_scores[best_idx]
    best_thresh = thresh_arr[best_idx] if best_idx < len(thresh_arr) else 0.5

    return {
        "AUPRC": auprc,
        "AUROC": auroc,
        "Accuracy": accuracy,
        "Precision": precision,
        "Sensitivity_Recall": sensitivity,
        "Specificity": specificity,
        "F1_Score": f1,
        "MCC": mcc,
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

    print(f"Dataset Size:     {metrics['Total_Samples']:,} pairs "
          f"({metrics['Positive_Samples']:,} positive, {metrics['Negative_Samples']:,} negative | "
          f"Ratio 1:{metrics['Negative_Samples']/max(1, metrics['Positive_Samples']):.1f})")
    print("-" * 70)
    print(f"  [*] AUPRC (Primary Metric):   {metrics['AUPRC']:.4f}")
    print(f"  [*] AUROC:                    {metrics['AUROC']:.4f}")
    print(f"  F1 Score (thresh=0.5):      {metrics['F1_Score']:.4f}")
    print(f"  Precision:                  {metrics['Precision']:.4f}")
    print(f"  Sensitivity (Recall):       {metrics['Sensitivity_Recall']:.4f}")
    print(f"  Specificity:                {metrics['Specificity']:.4f}")
    print(f"  MCC:                        {metrics['MCC']:.4f}")
    print(f"  Accuracy:                   {metrics['Accuracy']:.4f}")
    print(f"  Optimal F1:                 {metrics['Best_F1']:.4f} (at threshold {metrics['Best_Threshold']:.4f})")

def print_task_evaluation(task_name, metrics, baseline_dict=None):
    """Print formatted evaluation report."""
    print("=" * 70)
    print(f"EVALUATION REPORT: {task_name.upper()}")
    if baseline_dict and task_name in baseline_dict:
        print(f"Description: {baseline_dict[task_name].get('Description', '')}")
    print("=" * 70)

    print(f"Dataset Size:     {metrics['Total_Samples']:,} pairs "
          f"({metrics['Positive_Samples']:,} positive, {metrics['Negative_Samples']:,} negative | "
          f"Ratio 1:{metrics['Negative_Samples']/max(1, metrics['Positive_Samples']):.1f})")
    print("-" * 70)
    print(f"  [*] AUPRC (Primary Metric):   {metrics['AUPRC']:.4f}")
    print(f"  [*] AUROC:                    {metrics['AUROC']:.4f}")
    print(f"  F1 Score (thresh=0.5):      {metrics['F1_Score']:.4f}")
    print(f"  Precision:                  {metrics['Precision']:.4f}")
    print(f"  Sensitivity (Recall):       {metrics['Sensitivity_Recall']:.4f}")
    print(f"  Specificity:                {metrics['Specificity']:.4f}")
    print(f"  MCC:                        {metrics['MCC']:.4f}")
    print(f"  Accuracy:                   {metrics['Accuracy']:.4f}")
    print(f"  Optimal F1:                 {metrics['Best_F1']:.4f} (at threshold {metrics['Best_Threshold']:.4f})")

    # Comparison against baselines
    if baseline_dict and task_name in baseline_dict:
        b = baseline_dict[task_name]
        print("-" * 70)
        benchmark_suite = "ESMARAPPI" if "ESM" in task_name else "DEEPARAPPI"
        print(f"COMPARISON WITH {benchmark_suite} BENCHMARKS (AUPRC):")
        print(f"  PPLM (Zero-Shot):           {metrics['AUPRC']:.4f}")
        for model_k, val in b.items():
            if model_k != "Description":
                print(f"  {model_k:<27} {val:.4f}")

        primary_comp = "ESMAraPPI" if "ESMAraPPI" in b else "DeepAraPPI"
        if primary_comp in b:
            delta = metrics['AUPRC'] - b[primary_comp]
            sign = "+" if delta >= 0 else ""
            print(f"  PPLM vs {primary_comp}:         {sign}{delta:.4f}")
    print("=" * 70 + "\n")

def evaluate_file(csv_path, task_name="Task", output_dir=None):
    """Load predictions CSV and compute metrics."""
    if not os.path.exists(csv_path):
        print(f"[!] Error: File not found: {csv_path}")
        return None

    df = pd.read_csv(csv_path)
    if "pred_score" not in df.columns or "true_label" not in df.columns:
        print(f"[!] Error: CSV must contain 'pred_score' and 'true_label' columns.")
        return None

    y_true = df["true_label"].values
    y_score = df["pred_score"].values

    metrics = compute_metrics(y_true, y_score)
    all_baselines = {**DEEPARAPPI_BASELINES, **ESMARAPPI_BASELINES}
    print_task_evaluation(task_name, metrics, all_baselines)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_summary_path = os.path.join(output_dir, f"{task_name.lower()}_metrics.txt")
        with open(out_summary_path, "w") as f:
            f.write(f"Task: {task_name}\n")
            for k, v in metrics.items():
                f.write(f"{k}: {v}\n")
        print(f"Saved metrics summary to: {out_summary_path}")

    return metrics

def main():
    parser = argparse.ArgumentParser(description="Evaluate PPLM on Plant PPI Benchmarks (DeepAraPPI / ESMAraPPI)")
    # Single file evaluation
    parser.add_argument("--input", help="Single prediction CSV file")
    parser.add_argument("--task", default="Task", choices=["Task1", "Task2", "Task3", "Task4_Rice", "ESM_C2", "ESM_C3", "Custom"], help="Task name for benchmark lookup")
    
    # DeepAraPPI multi-file evaluation
    parser.add_argument("--c1", help="Path to DeepAraPPI C1 test prediction CSV")
    parser.add_argument("--c2", help="Path to DeepAraPPI C2 test prediction CSV")
    parser.add_argument("--c3", help="Path to DeepAraPPI C3 test prediction CSV")
    parser.add_argument("--rice", help="Path to DeepAraPPI Rice test prediction CSV")
    
    # ESMAraPPI multi-file evaluation
    parser.add_argument("--esm_c2", help="Path to ESMAraPPI C2 prediction CSV")
    parser.add_argument("--esm_c3", help="Path to ESMAraPPI C3 prediction CSV")
    
    parser.add_argument("--output_dir", default="results", help="Directory to save metric summary reports")

    args = parser.parse_args()

    results_table = []

    if args.input:
        m = evaluate_file(args.input, task_name=args.task, output_dir=args.output_dir)
        if m:
            results_table.append({"Task": args.task, **m})

    # DeepAraPPI evaluations
    if args.c1:
        m = evaluate_file(args.c1, task_name="Task1", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "Task1 (C1)", **m})

    if args.c2:
        m = evaluate_file(args.c2, task_name="Task2", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "Task2 (C2)", **m})

    if args.c3:
        m = evaluate_file(args.c3, task_name="Task3", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "Task3 (C3)", **m})

    if args.rice:
        m = evaluate_file(args.rice, task_name="Task4_Rice", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "Task4 (Rice)", **m})

    # ESMAraPPI evaluations
    if args.esm_c2:
        m = evaluate_file(args.esm_c2, task_name="ESM_C2", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "ESMAraPPI C2 (One Unseen)", **m})

    if args.esm_c3:
        m = evaluate_file(args.esm_c3, task_name="ESM_C3", output_dir=args.output_dir)
        if m:
            results_table.append({"Task": "ESMAraPPI C3 (Both Unseen)", **m})

    if len(results_table) > 1:
        df_summary = pd.DataFrame(results_table)
        summary_csv = os.path.join(args.output_dir, "benchmark_summary.csv")
        df_summary.to_csv(summary_csv, index=False)
        print(f"\n[INFO] Consolidated Benchmark Summary saved to: {summary_csv}")

if __name__ == "__main__":
    main()
