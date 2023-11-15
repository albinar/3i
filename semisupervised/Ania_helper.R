###3i semi-supervised pipeline parts
##Written by : Ania Lorenc, Aug 2018
#

#### Functions

get_subset_wt_ranks <- function(wt_ranks_all, ko_ranks){
  #removing WT ranks corresponding to KO ranks. Deals with duplicated ranks by removing closest rank (first one in the case of tie)
  wt_ranks_after_removal <- setdiff(wt_ranks_all, ko_ranks)
  for(i in ko_ranks[duplicated(ko_ranks)]){
    to_remove <- which.min(abs(wt_ranks_after_removal - i))
    wt_ranks_after_removal <- wt_ranks_after_removal[-to_remove]
  }
  return(wt_ranks_after_removal)
}