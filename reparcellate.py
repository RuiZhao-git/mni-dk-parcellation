"""
reparcellate.py — Re-parcellate fMRI BOLD into 84-region timeseries
using the Desikan-Killiany atlas in MNI152NLin6Asym 1mm space.

Usage:
    python reparcellate.py --bold <bold.nii.gz> --out <output_dir>

The default output is a TSV with the same format (column ordering + naming)
as the sample TSVs shared in the collaboration: 16 subcortical columns
followed by 68 cortical columns, uppercase + underscore naming. Preprocessing
state by default:
  - Linear detrending applied during voxel-to-region averaging
  - First 4 TRs discarded (T1 saturation transient)
  - NO bandpass filtering
  - NO confound regression

If bandpass filtering should be applied during parcellation, pass
--bandpass 0.01,0.1 (or any low,high band in Hz).

Notes on atlas:
  - The bundled atlas is the smoothed maximum-probability variant of the
    Desikan-Killiany aparc+aseg parcellation in MNI152NLin6Asym 1mm space,
    from Lotter et al. (g-node DOI 10.12751/g-node.2mnxpm, CC BY-SA 4.0).
  - The "tight" variant (strict cortical ribbon, no smoothing) is also
    bundled and can be selected via --atlas tight.
  - The aparc+aseg combination includes Ventral DC (Ventral Diencephalon)
    in the subcortical set rather than cerebellum. The output column names
    reflect this: LEFT_VENTRAL_DC and RIGHT_VENTRAL_DC are at the column
    positions that other DK pipelines may use for cerebellum.
"""
import os
import sys
import time
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.maskers import NiftiLabelsMasker
from scipy.signal import butter, filtfilt


# ---------------------------------------------------------------------------
# Output column ordering + naming (matches the collaboration's TSV format)
# ---------------------------------------------------------------------------
# Column order: 16 subcortical (LEFT 8, RIGHT 8) then 68 cortical (CTX_LH 34, CTX_RH 34).
# Naming: uppercase + underscore; '_PROPER' suffix for THALAMUS; VentralDC -> VENTRAL_DC.

GNODE_TO_OUTPUT_NAME = {
    # subcortical
    "Left-Thalamus":       "LEFT_THALAMUS_PROPER",
    "Left-Caudate":        "LEFT_CAUDATE",
    "Left-Putamen":        "LEFT_PUTAMEN",
    "Left-Pallidum":       "LEFT_PALLIDUM",
    "Left-Hippocampus":    "LEFT_HIPPOCAMPUS",
    "Left-Amygdala":       "LEFT_AMYGDALA",
    "Left-Accumbens-area": "LEFT_ACCUMBENS_AREA",
    "Left-VentralDC":      "LEFT_VENTRAL_DC",
    "Right-Thalamus":       "RIGHT_THALAMUS_PROPER",
    "Right-Caudate":        "RIGHT_CAUDATE",
    "Right-Putamen":        "RIGHT_PUTAMEN",
    "Right-Pallidum":       "RIGHT_PALLIDUM",
    "Right-Hippocampus":    "RIGHT_HIPPOCAMPUS",
    "Right-Amygdala":       "RIGHT_AMYGDALA",
    "Right-Accumbens-area": "RIGHT_ACCUMBENS_AREA",
    "Right-VentralDC":      "RIGHT_VENTRAL_DC",
}


def rename_gnode_label(label: str) -> str:
    """Convert a g-node atlas label to the output naming convention."""
    if label in GNODE_TO_OUTPUT_NAME:
        return GNODE_TO_OUTPUT_NAME[label]
    # cortical: 'ctx-lh-bankssts' -> 'CTX_LH_BANKSSTS'
    if label.startswith("ctx-"):
        return label.upper().replace("-", "_")
    raise ValueError(f"Unrecognized atlas label: {label!r}")


def build_output_order(gnode_labels: list[str]) -> list[int]:
    """
    Return a list of g-node 0-based indices in the desired output order:
      slots 0-7:  LEFT subcortical (VentralDC, Thalamus, Caudate, Putamen,
                                    Pallidum, Hippocampus, Amygdala, Accumbens)
      slots 8-15: RIGHT subcortical (same order)
      slots 16-49: CTX_LH cortical (alphabetical)
      slots 50-83: CTX_RH cortical (alphabetical)
    """
    target_order = [
        # Slot 0-7: LEFT subcortical (VentralDC first, matching position of
        # CEREBELLUM in other DK pipelines)
        "Left-VentralDC", "Left-Thalamus", "Left-Caudate", "Left-Putamen",
        "Left-Pallidum", "Left-Hippocampus", "Left-Amygdala", "Left-Accumbens-area",
        # Slot 8-15: RIGHT subcortical
        "Right-VentralDC", "Right-Thalamus", "Right-Caudate", "Right-Putamen",
        "Right-Pallidum", "Right-Hippocampus", "Right-Amygdala", "Right-Accumbens-area",
    ]
    # Slot 16-49: cortical LH (already alphabetical in g-node labels 0-33)
    target_order += [lbl for lbl in gnode_labels if lbl.startswith("ctx-lh-")]
    # Slot 50-83: cortical RH (g-node labels 34-67)
    target_order += [lbl for lbl in gnode_labels if lbl.startswith("ctx-rh-")]

    label_to_idx = {lbl: i for i, lbl in enumerate(gnode_labels)}
    return [label_to_idx[lbl] for lbl in target_order]


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def derive_output_filename(bold_path: Path, out_dir: Path) -> Path:
    """
    Derive output TSV filename from BOLD filename.
    BIDS pattern: replace '_desc-preproc_bold.nii.gz' with
                  '_desc-timeseries_desikan_killiany.tsv'.
    Fallback: use BOLD stem + '_desikan_killiany.tsv'.
    """
    name = bold_path.name
    if "_desc-preproc_bold.nii.gz" in name:
        tsv_name = name.replace(
            "_desc-preproc_bold.nii.gz",
            "_desc-timeseries_desikan_killiany.tsv",
        )
    else:
        stem = re.sub(r"\.nii(\.gz)?$", "", name)
        tsv_name = f"{stem}_desikan_killiany.tsv"
    return out_dir / tsv_name


def verify_alignment(bold_img: nib.Nifti1Image, atlas_img: nib.Nifti1Image) -> None:
    if bold_img.shape[:3] != atlas_img.shape:
        raise ValueError(
            f"Shape mismatch: BOLD spatial {bold_img.shape[:3]} vs atlas {atlas_img.shape}"
        )
    if not np.allclose(bold_img.affine, atlas_img.affine, atol=1e-3):
        raise ValueError(
            "Affine mismatch between BOLD and atlas. They must be in the same template space."
        )


def bandpass_filter(ts: np.ndarray, low: float, high: float, tr: float) -> np.ndarray:
    """Second-order Butterworth bandpass with zero-phase filtfilt, along axis 0."""
    fs = 1.0 / tr
    nyq = fs / 2.0
    b, a = butter(2, [low / nyq, high / nyq], btype="bandpass")
    return filtfilt(b, a, ts, axis=0)


def parse_band(s: str) -> tuple[float, float]:
    """Parse '0.01,0.1' into (0.01, 0.1)."""
    parts = s.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Bandpass must be 'low,high' in Hz, e.g. 0.01,0.1")
    low, high = float(parts[0]), float(parts[1])
    if not (0 < low < high):
        raise argparse.ArgumentTypeError(f"Invalid band: low={low}, high={high}")
    return low, high


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Re-parcellate fMRI BOLD using the Desikan-Killiany atlas in "
            "MNI152NLin6Asym 1mm space, outputting an 84-region TSV."
        ),
    )
    parser.add_argument("--bold", required=True, type=Path,
                        help="Path to BOLD NIfTI (must be in MNI152NLin6Asym 1mm space).")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output directory (will be created if missing).")
    parser.add_argument("--atlas", default="smoothed", choices=["smoothed", "tight"],
                        help="Atlas variant (default: smoothed maximum-probability).")
    parser.add_argument("--atlas-path", type=Path, default=None,
                        help="Override atlas variant with an explicit NIfTI path.")
    parser.add_argument("--n-cut", type=int, default=4,
                        help="Number of TRs to discard from the start (default: 4 for T1 "
                             "saturation transient). Use 0 to skip.")
    parser.add_argument("--bandpass", type=parse_band, default=None,
                        help="Optional bandpass band as 'low,high' in Hz (e.g. 0.01,0.1). "
                             "If not given, no bandpass is applied.")
    parser.add_argument("--tr", type=float, default=None,
                        help="TR in seconds, used for bandpass. If not given, taken from "
                             "BOLD header pixdim[4].")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite the output TSV if it already exists. By default, "
                             "the script skips subjects whose output already exists so that "
                             "interrupted batch runs can be safely resumed by re-running.")
    args = parser.parse_args()

    # ----- Resolve atlas path -----
    if args.atlas_path is not None:
        atlas_path = args.atlas_path
    else:
        repo_root = Path(__file__).resolve().parent
        if args.atlas == "smoothed":
            atlas_path = repo_root / "atlas" / "seg-aparc.aseg_space-MNI152NLin6Asym_desc-smoothed.nii.gz"
        else:
            atlas_path = repo_root / "atlas" / "seg-aparc.aseg_space-MNI152NLin6Asym.nii.gz"
    atlas_tsv = atlas_path.parent / "seg-aparc.aseg_space-MNI152NLin6Asym.tsv"
    if not atlas_path.exists():
        sys.exit(f"Atlas file not found: {atlas_path}")
    if not atlas_tsv.exists():
        sys.exit(f"Atlas label TSV not found: {atlas_tsv}")

    args.out.mkdir(parents=True, exist_ok=True)

    # ----- Resolve output filename early + skip if it exists (resume support) -----
    out_path = derive_output_filename(args.bold, args.out)
    if out_path.exists() and not args.overwrite:
        print(f"Output already exists, skipping: {out_path}")
        print("(use --overwrite to force re-processing)")
        return 0

    t_start = time.time()
    print(f"BOLD:     {args.bold}")
    print(f"Atlas:    {atlas_path.name}")
    print(f"Output:   {out_path}")
    if out_path.exists() and args.overwrite:
        print("WARNING: --overwrite is set; existing output will be replaced.")

    # ----- Load + verify alignment -----
    bold_img = nib.load(str(args.bold))
    atlas_img = nib.load(str(atlas_path))
    print(f"BOLD shape:  {bold_img.shape}, voxel size: {bold_img.header.get_zooms()}")
    print(f"Atlas shape: {atlas_img.shape}, voxel size: {atlas_img.header.get_zooms()}")
    verify_alignment(bold_img, atlas_img)
    print("Alignment check: OK")

    # ----- Resolve TR -----
    if args.tr is not None:
        tr = args.tr
    else:
        tr = float(bold_img.header["pixdim"][4])
    print(f"TR: {tr:.3f}s (from {'argument' if args.tr is not None else 'header'})")

    # ----- Load g-node labels (84 in atlas file order) -----
    atlas_labels_df = pd.read_csv(atlas_tsv, sep="\t", header=None, names=["idx", "name"])
    gnode_labels = atlas_labels_df["name"].tolist()
    if len(gnode_labels) != 84:
        sys.exit(f"Expected 84 atlas labels, got {len(gnode_labels)}")

    # ----- Extract region timeseries -----
    print("Extracting region timeseries (nilearn NiftiLabelsMasker, detrend=True)...")
    masker = NiftiLabelsMasker(
        labels_img=str(atlas_path),
        strategy="mean",
        detrend=True,
        standardize=False,
        verbose=0,
    )
    region_ts = masker.fit_transform(str(args.bold))  # (T, n_regions_with_voxels)
    T = region_ts.shape[0]
    print(f"Extracted shape: ({T}, {region_ts.shape[1]})")

    # ----- Map extracted columns to g-node label index -----
    # nilearn 0.13: masker.region_ids_ is dict with 'background' first, then labels
    region_ids_attr = getattr(masker, "region_ids_", None)
    if region_ids_attr is not None:
        region_ids = [v for k, v in region_ids_attr.items() if k != "background"]
    else:
        region_ids = list(range(1, 85))
    # Build full (T, 84) array indexed by g-node label position (label - 1)
    full_ts = np.full((T, 84), np.nan)
    for col, lbl in zip(range(region_ts.shape[1]), region_ids):
        if 1 <= int(lbl) <= 84:
            full_ts[:, int(lbl) - 1] = region_ts[:, col]
    n_missing = int(np.isnan(full_ts).any(axis=0).sum())
    if n_missing:
        print(f"WARNING: {n_missing} region(s) had no voxels in the mask and are NaN.")

    # ----- Reorder + rename to output convention -----
    order = build_output_order(gnode_labels)
    output_labels = [rename_gnode_label(gnode_labels[i]) for i in order]
    output_ts = full_ts[:, order]
    assert len(output_labels) == 84

    # ----- Discard first N TRs -----
    if args.n_cut > 0:
        if args.n_cut >= T:
            sys.exit(f"--n-cut ({args.n_cut}) >= number of TRs ({T})")
        output_ts = output_ts[args.n_cut:]
        print(f"Discarded first {args.n_cut} TRs. Remaining: {output_ts.shape[0]}")

    # ----- Optional bandpass -----
    if args.bandpass is not None:
        low, high = args.bandpass
        nyq = 1.0 / tr / 2.0
        if high >= nyq:
            sys.exit(f"Bandpass high cutoff ({high} Hz) >= Nyquist ({nyq:.3f} Hz at TR={tr}s)")
        print(f"Applying bandpass {low}-{high} Hz (Butterworth 2nd order, filtfilt)...")
        output_ts = bandpass_filter(output_ts, low, high, tr)
    else:
        print("Bandpass: not applied (use --bandpass low,high to enable).")

    # ----- Save TSV -----
    out_df = pd.DataFrame(output_ts, columns=output_labels)
    out_df.to_csv(out_path, sep="\t", index=False, float_format="%.10g")
    elapsed = time.time() - t_start
    print(f"Saved: {out_path}")
    print(f"Final shape: {output_ts.shape}, mean={output_ts.mean():.4f}, std={output_ts.std():.4f}")
    print(f"Elapsed: {elapsed:.1f} seconds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
