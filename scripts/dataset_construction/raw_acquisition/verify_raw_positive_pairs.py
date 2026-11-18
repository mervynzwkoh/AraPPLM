"""
Legacy wrapper for scripts/dataset_construction/01_verify_positive_pairs.py
"""
import os
import sys

target = os.path.join(os.path.dirname(__file__), "..", "01_positive_pairs", "01_verify_positive_pairs.py")
with open(target, "r", encoding="utf-8") as f:
    code = compile(f.read(), target, "exec")
exec(code, globals())
