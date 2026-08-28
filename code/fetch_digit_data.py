#!/usr/bin/env python3
"""Fetch UCI optdigits (8x8) and MNIST (28x28) from OpenML and write the raw
binaries that V2_digit_simulation.R reads. Float32 row-major (C-order), int32
labels, little-endian. Run once before the R digit analyses."""
import sys
import numpy as np
from sklearn.datasets import fetch_openml

OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp"

# UCI optical recognition of handwritten digits: 5620 x 64, values 0..16
opt = fetch_openml("optdigits", version=1, as_frame=False)
Xo = np.ascontiguousarray(np.asarray(opt.data, dtype=np.float32))
yo = np.asarray(opt.target, dtype=np.int32)
print("optdigits:", Xo.shape, "range", float(Xo.min()), float(Xo.max()))
Xo.tofile(f"{OUT}/digits8_X.bin")
yo.tofile(f"{OUT}/digits8_y.bin")

# MNIST: 70000 x 784, values 0..255
mn = fetch_openml("mnist_784", version=1, as_frame=False)
Xm = np.ascontiguousarray(np.asarray(mn.data, dtype=np.float32))
ym = np.asarray(mn.target, dtype=np.int32)
print("mnist:", Xm.shape, "range", float(Xm.min()), float(Xm.max()))
Xm.tofile(f"{OUT}/mnist_X.bin")
ym.tofile(f"{OUT}/mnist_y.bin")
print("wrote binaries to", OUT)
