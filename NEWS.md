
---
`EvoPhylo` News and Updates
======
## Development: EvoPhylo 0.2.1
(July 8th, 2022)

 * Added function (`drop.dummy.mb`) to remove "dummy" tip from Mr. Bayes summary trees and export tree accounting for metadata on the tips (for fully extinct clades).

## Development: EvoPhylo 0.2.0

(June 27th, 2022)

 * Added functions (`drop.dummy.beast`, `offset.to.dummy`, `write.beast.treedata`, `offset.to.dummy.metadata`) to handle trees with offsets as produced by the BEAST2 SA package for fully extinct clades
 * FBD functions can now handle BEAST2 log files produced by several packages (SA, BDSKY)
 * Added BEAST2 support for clock rate table import - function exists now in two versions, `get_clockrate_table_MrBayes` and `get_clockrate_table_BEAST2`
 * `plot_treerates_sgn()` function update: outputs chosen thresholds and estimated rate values for each threshold. 

## EvoPhylo 0.1.0

(May 12th, 2022)

* First version!
