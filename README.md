# MNI-DK parcellation

Parcellation pipeline for fMRIPrep-preprocessed BOLD data into 84-region
timeseries using the **Desikan-Killiany (DK)** cortical and subcortical
atlas in **MNI152NLin6Asym 1mm** space.

The pipeline uses `nilearn.maskers.NiftiLabelsMasker` (`strategy='mean'`,
`detrend=True`) for voxel-to-region averaging, discards the first 4 TRs
to remove the T1 saturation transient, and saves an 84-column TSV per
subject. Bandpass filtering is **not** applied by default (see "Output"
below).

---

## Requirements

- Python 3.11 or newer
- Approximately 6 GB RAM per subject (one at a time; no GPU needed)
- No HPC / cluster / specialized hardware required — runs locally on
  laptop or workstation

Dependencies (`requirements.txt`):

```
nilearn>=0.13.0
nibabel>=5.0
scipy>=1.10
numpy>=1.24
pandas>=2.2,<3
```

Install:

```bash
pip install -r requirements.txt
```

---

## Usage

### Single subject

```bash
python reparcellate.py \
    --bold /path/to/sub-XXX_..._space-MNI152NLin6Asym_res-01_desc-preproc_bold.nii.gz \
    --out  /path/to/output_dir/
```

The output TSV is named by replacing `_desc-preproc_bold.nii.gz` with
`_desc-timeseries_desikan_killiany.tsv`, and placed in `--out`.

### Multiple subjects

For batch processing of a directory of BOLD files, see
`examples/process_directory.sh`. This is a **template**: open it in a
text editor, set `BOLD_DIR` and `OUTPUT_DIR` to your real paths, then
run it from a terminal:

```bash
bash examples/process_directory.sh
```

The script will loop over every `*_desc-preproc_bold.nii.gz` in
`BOLD_DIR` and write each output TSV to `OUTPUT_DIR`.

### Options

| Flag | Default | Description |
|---|---|---|
| `--bold` | (required) | Path to BOLD NIfTI in MNI152NLin6Asym 1mm space |
| `--out` | (required) | Output directory |
| `--atlas` | `smoothed` | Atlas variant: `smoothed` (recommended for fMRI) or `tight` |
| `--atlas-path` | (none) | Override atlas with an explicit NIfTI path |
| `--n-cut` | `4` | Number of TRs to discard from the start (0 = none) |
| `--bandpass` | (none) | Optional bandpass band as `low,high` in Hz, e.g. `0.01,0.1` |
| `--tr` | (from header) | TR in seconds (used for bandpass) |

---

## Output

The output TSV has **84 columns** in this fixed order:

| Column index (1-based) | Region |
|---|---|
| 1–8 | LEFT subcortical: VENTRAL_DC, THALAMUS_PROPER, CAUDATE, PUTAMEN, PALLIDUM, HIPPOCAMPUS, AMYGDALA, ACCUMBENS_AREA |
| 9–16 | RIGHT subcortical (same order) |
| 17–50 | CTX_LH_* cortical (34 regions, alphabetical) |
| 51–84 | CTX_RH_* cortical (34 regions, alphabetical) |

### Default preprocessing state

By default, the output TSV has:

- **Spatial:** mean voxel-to-region averaging
- **Temporal — applied:**
  - Linear detrending (during extraction, removes per-region linear drift)
  - First 4 TRs discarded (T1 saturation transient)
- **Temporal — NOT applied:**
  - Bandpass filtering
  - Confound regression (motion / WM / CSF / GSR)

This matches the preprocessing state expected by downstream FC pipelines
that apply their own bandpass at the FC computation step. **Applying
bandpass twice (here and downstream) would over-filter the signal**, so
the default leaves bandpass to the downstream code.

If your downstream code does **not** apply bandpass and you would like
it applied during parcellation, use `--bandpass 0.01,0.1` (or any
appropriate band).

### Verifying output preprocessing state

The output state is directly verifiable from the TSV itself:

| Check | How |
|---|---|
| Number of TRs (cut applied?) | `pd.read_csv(...).shape[0]` |
| Detrend applied? | Each column's mean should be ≈ 0 |
| Bandpass applied? | Power spectral density of each column should fall outside the filter band |

---

## Methodology rationale

| Choice | Rationale |
|---|---|
| Atlas: smoothed maximum-probability variant | The tight (strict cortical ribbon) variant excludes voxels in the partial-volume periphery of the cortex, which are functionally informative for fMRI BOLD. The smoothed variant (6 mm FWHM Gaussian + Schaefer 2018 gray-matter mask + max-probability label) covers these voxels. Per the atlas author's README, smoothed is recommended for fMRI use cases. |
| Strategy: `mean` (NiftiLabelsMasker) | Standard voxel-to-region aggregation for fMRI BOLD. |
| Detrend: `True` (during extraction) | Removes per-region linear drift accumulated during the scan. Standard preprocessing step; baseline assumption for rs-fMRI analyses. |
| Standardize: `False` | Preserves original signal amplitude. Downstream metrics (Pearson FC, etc.) are scale-invariant, so standardization is not needed at the parcellation step. |
| Discard first 4 TRs | T1 saturation transient: spin equilibrium is not yet reached in the first ~4 TRs, signal baseline is unstable. Universally discarded in modern fMRI pipelines. |
| Bandpass: off by default | See "Default preprocessing state" above. |
| Confound regression: not applied | Out of scope for parcellation. Apply downstream via your preferred confound model (e.g., aCompCor, motion regression from fMRIPrep's `desc-confounds_timeseries.tsv`). |

---

## Atlas attribution

The bundled atlas files in `atlas/` are derived from:

> Lotter, L. D., Saberi, A., Hansen, J. Y., Misic, B., Paquola, C.,
> Barker, G. J., Bokde, A. L. W., Desrivières, S., Flor, H., Grigis, A.,
> Garavan, H., Gowland, P., Heinz, A., Brühl, R., Martinot, J.-L.,
> Paillère Martinot, M.-L., Artiges, E., Nees, F., Orfanos, D. P.,
> Lemaitre, H., Paus, T., Poustka, L., Hohmann, S., Holz, N., Fröhner,
> J. H., Smolka, M. N., Vaidya, N., Walter, H., Whelan, R., Schumann,
> G., the IMAGEN consortium, Eickhoff, S. B., Bauer, T., & Dukart, J.
> "Maximum-probability atlas in MNI space" via `aparc+aseg` from
> FreeSurfer. g-node DOI 10.12751/g-node.2mnxpm, CC BY-SA 4.0.

The two variants bundled:

- `seg-aparc.aseg_space-MNI152NLin6Asym_desc-smoothed.nii.gz` — smoothed
  maximum-probability variant (default; recommended for fMRI).
- `seg-aparc.aseg_space-MNI152NLin6Asym.nii.gz` — tight cortical-ribbon
  variant (alternative).
- `seg-aparc.aseg_space-MNI152NLin6Asym.tsv` — label list (84 labels).

The full archive (including additional atlases like Destrieux and DKT) is
available at <https://doi.gin.g-node.org/10.12751/g-node.2mnxpm/>.

### Region set note: VentralDC vs. Cerebellum

The `aparc+aseg` combination from FreeSurfer includes **Ventral DC**
(Ventral Diencephalon) in the subcortical set rather than the cerebellum
cortex. The output TSV reflects this: columns `LEFT_VENTRAL_DC` and
`RIGHT_VENTRAL_DC` occupy slot positions 1 and 9, where other DK
pipelines (using `aseg` with cerebellum) may have `LEFT_CEREBELLUM_CORTEX`
and `RIGHT_CEREBELLUM_CORTEX`. The remaining 14 subcortical regions
(7 per hemisphere: Thalamus, Caudate, Putamen, Pallidum, Hippocampus,
Amygdala, Accumbens) match exactly.

---

## License

MIT — see `LICENSE`.

The bundled atlas files are licensed under CC BY-SA 4.0 (see atlas
attribution above).
