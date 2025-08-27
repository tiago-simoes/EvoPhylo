
---
`EvoPhylo` News and Updates
======
## Development: EvoPhylo 0.3.5
(August 26th, 2025)
* Function `combine_log` now also reads and combine log files from BEAST2 and MCMCTREE
* Added Z-transformation for functions `plot_back_rates` and `plot_treerates_sgn`
* New function `write_partitioned_alignments2` write alignment partitions as separate alignment files for various data types
* New function `clade_membership`designates clade membership for each tip for downstream analyses summarizing rates for each clade

## Development: EvoPhylo 0.3.3
(July 21st, 2023)

* Function `plot_back_rates` corrected for MrBayes.
* New warnings and import of MrBayes single {all} clock partition for `plot_treerates_sgn`.

## Development: EvoPhylo 0.3.2
(October 31st, 2022)

* New function `plot_back_rates` and updates on vignettes.

## Development: EvoPhylo 0.3.1
(OCtober 20th, 2022)

* Bug fixes and adjustment of selection functions for BEAST2: old `get_pwt_rates` function exists now in two versions`get_pwt_rates_MrBayes` and `get_pwt_rates_BEAST2`, and updates on `plot_treerates_sgn` to handle BEAST2.

## Development: EvoPhylo 0.3
(August 1st, 2022)

Implementation of functions for BEAST2 data now complete! 

 * Added function (`write_partitioned_alignments`) to write partitioned alignments based on the inferred clusters.
 * Added BEAST2 support for clock rate table import - old `get_clockrate_table` function exists now in two versions, `get_clockrate_table_MrBayes` and `get_clockrate_table_BEAST2`
  
## Development: EvoPhylo 0.2.1
(July 8th, 2022)

 * Added function (`drop.dummy.mb`) to remove "dummy" tip from Mr. Bayes summary trees and export tree accounting for metadata on the tips (for fully extinct clades).

## Development: EvoPhylo 0.2.0

(June 27th, 2022)

 * Added functions (`drop.dummy.beast`, `offset.to.dummy`, `write.beast.treedata`, `offset.to.dummy.metadata`) to handle trees with offsets as produced by the BEAST2 SA package for fully extinct clades
 * FBD functions can now handle BEAST2 log files produced by several packages (SA, BDSKY)
 * `plot_treerates_sgn()` function update: outputs chosen thresholds and estimated rate values for each threshold. 

## EvoPhylo 0.1.0

(May 12th, 2022)

* First version!
