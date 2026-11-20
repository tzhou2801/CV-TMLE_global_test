# CV-TMLE global test: simulation code

R code for the two simulation studies in

> Zhou T, Gruber S, Lee H, Lee W, Nie L, van der Laan MJ (2026). A CV-TMLE global test approach to improve power in rare disease clinical studies with multiple-component endpoints. arXiv:2605.02851. https://arxiv.org/abs/2605.02851

| Script | Simulation study |
|---|---|
| [`sim_study_1.R`](sim_study_1.R) | Study 1 (Section 6.1) |
| [`sim_study_2.R`](sim_study_2.R) | Study 2 (Section 6.2 and Appendix A) |

## Requirements

R and the following CRAN packages:

```r
install.packages(c("tmle", "SuperLearner", "dplyr", "foreach", "doParallel"))
```

`MASS` and `parallel` ship with R. The results in the paper were produced with R 4.2.1, tmle 2.0.1.1,
SuperLearner 2.0-28, MASS 7.3-60, dplyr 1.1.2, foreach 1.5.2 and doParallel 1.0.17.

## Running the simulations

Each script is self-contained. Running it executes the full study (1,000 replications per scenario, seed 202701)
in parallel and writes timestamped CSV files with the raw and summary results to the working directory:

```sh
Rscript sim_study_1.R
Rscript sim_study_2.R
```

`sim_study_1.R` uses 7 worker processes (`n_cores = 7` at the bottom of the script); `sim_study_2.R` uses all but
one of the available cores. Each study takes about 2 hours on an 8-core laptop. Every replication is seeded
separately, so the results do not depend on the number of cores.

## Citation

```bibtex
@misc{zhou2026cvtmle,
  author        = {Zhou, Tianyue and Gruber, Susan and Lee, Hana and Lee, Wonyul and Nie, Lei and van der Laan, Mark J.},
  title         = {A {CV-TMLE} global test approach to improve power in rare disease clinical studies with multiple-component endpoints},
  year          = {2026},
  eprint        = {2605.02851},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME}
}
```

## License

MIT; see [LICENSE](LICENSE).
