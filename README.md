# Figures for the SharpClaw polynomial-reconstruction docs page

This branch holds the code that draws the seven figures on the Clawpack documentation page `pyclaw/sharpclaw_poly.rst` (clawpack/doc), which documents the polynomial reconstructions added to SharpClaw in clawpack/pyclaw#492 (`solver.lim_type` 4, 5 and 6). The docs themselves are static: they contain only the PNGs, and building them runs none of this code. This branch exists so the figures can be reproduced and extended later.

It is a standalone branch with no shared history with the docs or with PyClaw. It is meant to stay as it is; the docs page cites the commit.

| file | purpose |
|---|---|
| `make_figures.py` | redraws `figures/*.png` |
| `poly_weights.py` | exact weights, stencils and Fourier symbols of the reconstructions |
| `sine_reconstruction.py` | the sine example and the observed-order table |
| `figures/*.png` | the figures as committed to the docs |
| `requirements.txt` | the NumPy and matplotlib versions the figures were made with |

## Reproducing the figures

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python make_figures.py              # rewrites figures/*.png
python sine_reconstruction.py       # observed-order table on the page
python poly_weights.py --table      # weight tables on the page
```

The PNGs are written without a Software tag, so with the versions in `requirements.txt` the output is byte-for-byte identical to the files in `figures/` and in the docs.

## Where the code comes from

`poly_weights.py` and `sine_reconstruction.py` are unmodified copies of the files in PyClaw's `development/poly_reconstruction/` directory at commit `3dd5e747` (clawpack/pyclaw#492). PyClaw's copies are canonical: there `poly_weights.py` also generates the Fortran weights `src/pyclaw/sharpclaw/poly.f90`, and its `--check` and `--write` options only work inside a PyClaw checkout. `make_figures.py` comes from the same directory at commit `02f32a68` of the pull request branch, where it lived before the figures moved to the docs.

If the stencils change in PyClaw, copy its `poly_weights.py` and `sine_reconstruction.py` here, rerun `make_figures.py`, and update the figures and tables in the docs.
