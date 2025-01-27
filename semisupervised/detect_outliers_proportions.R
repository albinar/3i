require(data.table)
require(ggplot2)
require(dplyr)
require(readr)








source("3i-Helper.R")
source("Ania_helper.R")
#Input.path needs to be changed based on user's directory
input.path <- "/home/rstudio/results/3i/"
##This is for BM:


######################INPUT FILES' locations
file_with_hits <- paste(input.path,panel,"/semisupervised/",organ,"/",panel,"-",organ,"_hits_IMPC_ML_080618.csv",sep="")
file_with_Events_Proportions <- paste(input.path,panel,"semisupervised",organ,"allProportions_Table.csv",sep="/")
file_with_LMnumbers_map <-paste(input.path,panel,"semisupervised",organ,"Mnumber_Lnumber.csv" ,sep="/")


######################FUNCTIONS 

add_zygosity <- function(dt, Genotype_col="Genotype", new_format=FALSE){
  if(new_format){
    dt[, zygosity:=sub(".*:","",get(Genotype_col))]
    dt[zygosity=="Hemi", zygosity:="Hom"]
  }
  else{
    dt[, zygosity:=case_when(get(Genotype_col)%in%c("+/+","+/Y")~"WT",
                             grepl("/\\+",get(Genotype_col)) ~"Het",
                             grepl("\\+/",get(Genotype_col)) ~"Het",
                             is.character(get(Genotype_col) )~"Hom")]
  }
  return(dt)
  
}

take_only_WTs <- function(all_data_dt){
  #for given data dt only WT individuals, requires the column "Genotype"
  return(all_data_dt[Genotype%in%c("+/+","+/Y","WT")])
  
}

get_a_balanced_subset <- function(one_day_data) {
  #takes only one animal per sex and strain
  one_day_data_collapsed <-
    one_day_data[, .SD[1], by = .(Genotype, Gender), .SDcols = colnames(one_day_data)]
  how_many_one_sex <-
    one_day_data_collapsed [, .SD[1], by = .(Genotype, Gender), .SDcols =
                              colnames(one_day_data)][, .N, by = Gender][, min(N)]
  if (length(unique(one_day_data_collapsed$Gender)) > 1) {
    which_sex <-
      one_day_data_collapsed [, .SD[1], by = .(Genotype, Gender), .SDcols =
                                colnames(one_day_data)][, .N, by = Gender][order(Gender)][, c("Female", "Male")[which.min(N)]]
    res = rbind(one_day_data_collapsed[Gender == which_sex], one_day_data_collapsed[!Gender ==
                                                                                      which_sex][sample(1:nrow(one_day_data_collapsed[!Gender == which_sex]), how_many_one_sex)])
    return(res)
  } else{
    return(NULL)
  }
}


######################
library(readr)
  BM_hits <- data.table(read_csv(file_with_hits))
 # BM_hits<-BM_hits[BM_hits$zygosity!="WT",]
  LMnumbers_map <- data.table(read_csv(file_with_LMnumbers_map))
  LMnumbers_map[,Genotype:=NULL]
  LMnumbers_map[,Gender:=NULL]
  LMnumbers_map[,Assay_Date:=NULL]
  LMnumbers_map[,zygosity:=NULL]
  BM_for_M <- read.csv(file_with_Events_Proportions,check.names = F)
  to.bo.removed <- c(grep(x = colnames(BM_for_M),pattern = "All Events-%Parent"),grep(x = colnames(BM_for_M),pattern = "Singlets-%Parent" ),
                     grep(x = colnames(BM_for_M),pattern = "Live-%Parent"),grep(x = colnames(BM_for_M),pattern = "Lymphocytes-%Parent" ),
                     grep(x = colnames(BM_for_M),pattern = "CD45-%Parent"),grep(x = colnames(BM_for_M),pattern = "NOT*") )
  BM_for_M<-BM_for_M[,-to.bo.removed]
  BM_for_M <- data.table(BM_for_M)
  add_zygosity(BM_for_M,new_format=TRUE)
  BM_for_M[,`Assay Date`:=as.Date(`Assay Date`, format="%d-%b-%y")]
  #detach("package:RchyOptimyx", unload=TRUE)
  BM_M <- merge(x=BM_for_M, y=LMnumbers_map, by.x="Barcodes", by.y="Label_Barcode")
  

#Here I assume you are interested in outlier detection based on "fraction of parental gate"- if you would prefer % of CD45, you would need to prepare apropriate columns and preselect for them.
  #long format, column selection
  BM_for_M_m <- melt(BM_M[,c("Genotype","FCS files","Barcodes","Assay Date","Gender","zygosity","Mouse", grep("%", colnames(BM_for_M), val=T)), with=F],
                    id.vars=c("Genotype","FCS files","Barcodes","Assay Date","Gender","zygosity","Mouse"))[order(Genotype)]

  
  
  #WT-based overall and daily medians (!not sex-specific, assuming that days off would shift both sexes i na  similar way)  
  overall_medians_BM <- take_only_WTs(BM_for_M_m)[, .(medval=median(value, na.rm=T), madval=mad(value, na.rm=T), N = sum(!is.na(value))), by=.( variable)]
  daily_medians_BM <- take_only_WTs(BM_for_M_m)[, .(med_day=median(value, na.rm=T), mad_day=mad(value, na.rm=T), N = sum(!is.na(value)) ), by=.( variable, `Assay Date`)]
  

  all_samples_without_hits_BM <- BM_M[!Mouse%in%BM_hits$Mouse] 
  
  #cast them into broad format
  all_samples_without_hits_BM_d <- BM_M[,c("Genotype","FCS files","Barcodes","Assay Date","Gender","zygosity","Mouse", grep("%", colnames(BM_for_M), val=T)), with=F]
  
  days_to_use_balanced_ko <- c(unique(BM_M[!`Assay Date`%in%daily_medians_BM$`Assay Date`]$`Assay Date`),unique(daily_medians_BM[N<=2,`Assay Date`]))
  
  
  #obtain balanced medians
  library(dplyr)
  balanced_ko_medians_BM <- bind_rows(lapply( days_to_use_balanced_ko, function(x){
    bind_rows(lapply(unique(BM_for_M_m $variable), function(par){
      #  print(x)
      #  print(par)
      bs <- get_a_balanced_subset(all_samples_without_hits_BM_d[`Assay Date`==x, c("Mouse","Genotype", "Gender", "Assay Date", "zygosity", as.character(par)), with=F])
      if(!is.null(bs)){
        bsm <-melt(bs[,c("Mouse","Genotype", "Gender", "Assay Date", "zygosity", as.character(par)), with=F],
                   id.vars = c("Mouse","Genotype", "Gender", "Assay Date", "zygosity"))
        bsm[,.(`Assay Date`=x, variable=par, med_day_ko=median(value, na.rm=TRUE), mad_day_ko=mad(value, na.rm=TRUE), N= sum(!is.na(value)))]
      }
      else{bsm=NULL}
    }))
  }))
  
  
  
  together_BM <- 
    data.table:::merge.data.table(balanced_ko_medians_BM, daily_medians_BM, by=c("variable", "Assay Date"), suffixes = c(".ko", ".wt"), all=T) %>%
    data.table:::merge.data.table(., overall_medians_BM, by=c("variable"), all=T)
  
  
  together_BM[,is_med_ko_within_1_medval:=(medval-madval <= med_day_ko )&(medval+madval >=med_day_ko )]
  together_BM[,is_med_wt_within_1_medval:=(medval-madval <= med_day )&(medval+madval >=med_day )]
  together_BM[,is_med_ko_within_2_medval:=(medval-2*madval <= med_day_ko )&(medval+2*madval >=med_day_ko )]
  together_BM[,is_med_wt_within_2_medval:=(medval-2*madval <= med_day )&(medval+2*madval >=med_day )]
  together_BM[,is_med_ko_within_3_medval:=(medval-3*madval <= med_day_ko )&(medval+3*madval >=med_day_ko )]
  together_BM[,is_med_wt_within_3_medval:=(medval-3*madval <= med_day )&(medval+3*madval >=med_day )]
  
  #Dates when:
  ## wt present, some params outside of the 2 medians
  ##without wt/2 or less wt, balanced ko outside of the 2 medians
  ##
  suspected_dates_BM <- data.table:::merge.data.table(together_BM[N.wt>2][is_med_wt_within_2_medval==FALSE, .(params_wt=length(variable)),by="Assay Date"],
                                                      together_BM[N.ko>2][(is.na(is_med_wt_within_2_medval)&is_med_ko_within_2_medval==FALSE), .(params_ko=length(variable)),by=`Assay Date`], all=T, by="Assay Date", suffixes=c(".wt", ".ko")) %>%
    data.table:::merge.data.table(., together_BM[(is.na(N.wt)&N.ko<=2), .(Ko_param_with_less_than_3= length(variable)), by="Assay Date"], all=T) %>%
    data.table:::merge.data.table(., together_BM[N.wt<=2, .(wt_param_with_less_than_3 = length(variable)),by="Assay Date"], all=T)
  
  #How many strains we would loose by removing susp.days?
  data.table(dcast(BM_M[Genotype!="WT"][,.N, by=.(`Assay Date`%in%suspected_dates_BM$`Assay Date`, Genotype)],Genotype~`Assay Date`))[,.(good_samples=`FALSE`,susp_samples=`TRUE`,all=sum(c(`FALSE`,`TRUE`), na.rm=TRUE)),by=Genotype][all>=4][good_samples<4][order(all, decreasing = TRUE)]
  
  #Susp.days versus number of suspected parameters
  together_BM[is_med_wt_within_2_medval==FALSE|(is.na(is_med_wt_within_2_medval)&is_med_ko_within_2_medval==FALSE),.(n_params=.N),by=.(`Assay Date`)][order(`Assay Date`)]
