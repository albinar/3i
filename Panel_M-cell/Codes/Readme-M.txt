Order of Code:

Seperate_Files_Folders-M.R
GlobalFrame657-M.R 
lgl411-M.R
main3iTcellClean-M.R
main3iTcell-M.R (gating)
gthresOutlierCheck-M.R
main3iTcellFT-M.R
main3iTcellsecondPart-M.R   -> Phenotypes.Rdata, Matric3iTcell.csv
main3iTcellthirdPart-M.R    -> T_Test and WilcoxTest Results
main3iTcellthirdPart_B-M.R  -> Creates RchyOptimyx plots, AdjPvalues plots
ComparisonOfProportions-M.R -> FinalTableOfProportions.csv

Figures:
FancyFigures-M.R (Creates Box Plots - requires FinalTableOfProportions.csv)


Dec. 29th 2016
The above is outdated.
New Order of code:

preProcessingMain.R (SPLEEN) AND preProcessingMain_MLN.R (MLN)
SangerCleanMain.R (SPLEEN & MLN)
main3iTcell-M_parallel_SangerClean.R (SPLEEN) AND main3iTcell-M_parallel_SangerClean_MLN.R (MLN)

Aug. 2 2017
REPLACE main3iTcell-M_parallel_SangerClean.R (SPLEEN) in the above order with main3iTcell-M_parallel_SangerClean2.R (SPLEEN)