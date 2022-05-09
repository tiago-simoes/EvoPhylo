# EvoPhylo

A package to perform automated morphological character partitioning for phylogenetic analyses and analyze macroevolutionary parameter outputs from clock (time-calibrated) Bayesian inference analyses.

The ideas and rationale behind the functionality and objectives of the analyses available in this package were first presented and can be referenced to [Simões, T. R. & Pierce, S. E. 2021. Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. *Nature Ecology & Evolution*.](https://doi.org/10.1038/s41559-021-01532-x)

EvoPhylo is currently designed to pre- and postprocess morphological data for relaxed clock Bayesian phylogenetic analyses using the software [Mr. Bayes](https://nbisweden.github.io/MrBayes/). Subsequent releases will also implement functions to postprocess output data produced by the software package [BEAST2](http://www.beast2.org/).


## Installing package **EvoPhylo**

Install the release version directly from CRAN 
```{r}
install.packages("EvoPhylo")
```

or from Github (check this for the latest updates)
```{r}
# install.packages("devtools")
devtools::install_github("tiago-simoes/EvoPhylo")
```

## Tutorials

Please check the [EvoPhylo website](https://tiago-simoes.github.io/EvoPhylo/) for step by step tutorial vignettes.


## Authors

* [**Tiago R. Simões**](https://tiago-simoes.com/): *Conceptualization and code development* - [Github Profile](https://github.com/tiago-simoes)
* [**Noah Greifer**](https://github.com/ngreifer): *Code development* 
* [**Stephanie E. Pierce**](https://projects.iq.harvard.edu/spierce/home): *Project supervision* 

## License

This project is licensed under General Public License, version 2.

