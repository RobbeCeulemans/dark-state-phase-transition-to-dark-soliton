# Dark-state phase transition to a dark soliton in a dissipative Bose–Hubbard chain

Code to reproduce the figures of:

> R. Ceulemans, S. E. Begg, M. J. Davis, and M. Wouters,
> *"Dark-state phase transition to a dark soliton in a dissipative Bose–Hubbard chain"*,
> accepted (2026). Preprint: [arXiv:2509.25707](https://doi.org/10.48550/arXiv.2509.25707)

## Data

The datasets required to run the plotting scripts are hosted separately on Zenodo:

> **[Zenodo DOI here]**

Download the archive and unzip it so the resulting `darkstatephase_data/` folder sits
directly next to the `Fig_1`–`Fig_4` folders, i.e. at the repository root:

```
DarkStatePhase/
├── darkstatephase_data/   <- from Zenodo
├── Fig_1/
├── Fig_2/
├── Fig_3/
├── Fig_4/
├── Project.toml
└── Manifest.toml
```

`darkstatephase_data/` is git-ignored and not part of this repository.

## Repository structure

| Folder | Script | Produces |
|---|---|---|
| `Fig_1/` | *(none — hand-drawn schematic)* | Figure 1(a) |
| `Fig_2/` | `switchingtimes.jl` | Figure 2(b) and 2(c) |
| `Fig_2/` | `soliton.jl` | Figure B1 (Appendix B) |
| `Fig_3/` | `soliton_instab.jl` | Figures 3(a)/3(b), 3(c)/3(d), and B2 (Appendix B) |
| `Fig_4/` | `phasediagram.jl` | Figure 4(a)/4(b) and E1 (Appendix E) |

Each script saves its output PDF(s) alongside itself in the same `Fig_X/` folder.

## Requirements

- [Julia](https://julialang.org/) (developed and tested with 1.12.6)
- A Python installation with `matplotlib`, accessible to [PyCall.jl](https://github.com/JuliaPy/PyCall.jl)/[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl)
- A working LaTeX installation (the plots render text via `text.usetex=true`)

All Julia package dependencies are pinned in `Project.toml`/`Manifest.toml`.

## Running

From the repository root, instantiate the environment once:

```julia
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then run any figure script, e.g.:

```
julia --project=. Fig_4/phasediagram.jl
```

Scripts resolve their own data/output paths relative to their own location, so this
works regardless of the directory you launch `julia` from — as long as `--project=.`
points at the repository root (where `Project.toml` lives) and `darkstatephase_data/`
has been placed as described above.

## License

This code is released under the MIT License — see [LICENSE](LICENSE).
