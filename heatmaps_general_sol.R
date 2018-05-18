#####
# scooperdock@whrc.org
# Modified from Kylen's scripts
# Make heatmaps with time window on the x-axis and mortality threshold on the Y.
#####


#### Clear environment and load packages ####
rm(list=ls(all=TRUE))
require(doParallel)
require(gplots)
require(earlywarnings)
require(plyr)

#### Set parameters ####
# Run on lcomp?
NASAcomp_run = T

#
outpath_parent_dir = "Heatmapv7_20"

# Add labels or not?
label_heatmap = T

#Using cutoff data? What year is cutoff?
use_year_cutoff=F
year_cutoff=1980

#Set longest length of time to use between tree surveys
max_interval=10

# Mortality thresholds 
mort_thresh_list =  seq(0,25,1)


# List to define whether to make percentage heatmaps in addition to continuous figs
perc_list = list(
  "trend" = T,
  "pulse" = T,
  "ddj" = F,
  "ews" = F
)

# Metric thresholds for percentage heatmaps. Not currently making perc heatmaps for ddj, ews
metric_thresh_list = list(
  "trend" = c(0),
  "pulse" = c(-1.5,-1.75,-2),
  "ews" = c(0),
  "ddj" = c(0)
)

# Which VI data sets to include in analysis
data_set_vec=c("GIMMS")#,"MODIS","LANDSAT","LANDSATMAX","MODISTERRA")

# For each dataset, which to analysis to include in LM?
analysis_list = list(
    "GIMMS" = c("trend","pulse"))#,"ddj","ews"))#,
  #   "MODIS" = c("trend","pulse","ddj"),
  #   "MODISTERRA" = c("trend","pulse","ddj"),
  #   "LANDSAT" = c("trend","pulse","ddj"),
  #   "LANDSATMAX" = c("trend","pulse","ddj")
  # )

# Sub metrics. This is only really used for ddj and ews, so it's a bit inefficient to have it for both but oh well...
sub_metrics = list(
  "trend" = c("trend"),
  "pulse" = c("pulse"),
  "ews" = c("ar1","sd","kurt","densratio"),
  "ddj" = c("S2.t")
)

# Color ranges for continuous heatmaps. The third element of the list is T if you want the colorbar to have red be lower than blue
col_range_list = list(
  "trend" =c(-.5,.5,T),
  "pulse" =c(-3,-1,T),
  "ar1" = c(-1,1,F),
  "sd" = c(-1,1,F),
  "kurt" = c(-1,1,F),
  "densratio" = c(-1,1,F),
  "S2.t" = c(-1,1,F)
)

# Set figure dimensions (differs a little for ews and ddj so they can fit more closely)
fig_dims = list(
  "trend" = c(4.9,4.9),
  "pulse" = c(4.9,4.9),
  "ews" = c(6.125,4.9),
  "ddj" = c(6.125,4.9)
  )
# For cafi only
shift_list = list(
  "trend" = .5,
  "pulse"= 0,
  "ews" = 0,
  "ddj" = 0
)

# Set timelag start
tl_start = 0

# Time lag lengths for VI trend, pulse, ews, and ddj, separated by data set
tl_list = list(
  # GIMMS
  "GIMMS" = list(
    "trend" = 5:20,
    "pulse" = 5:20,
    "ews" = 5:20,
    "ddj"  = 5:20
  ),
  # MODIS
  "MODIS" = list(
    "trend" = 5:20,
    "pulse" = 5:20,
    "ews" = 5:20,
    "ddj"  = 5:20
  ),
  "MODISTERRA" = list(
    "trend" = 5:20,
    "pulse" = 5:20,
    "ews" = 5:20,
    "ddj"  = 5:20
  ),
  "LANDSAT" = list(
    "trend" = 5:20,
    "pulse" = 5:20,
    "ews" = 5:20,
    "ddj"  = 5:20
  ),
  "LANDSATMAX" = list(
    "trend" = 5:20,
    "pulse" = 5:20,
    "ews" = 5:20,
    "ddj"  = 5:20
  )
)
# 3x3 or 1pixel selection?
pix_select_list = list(
  "GIMMS" = "1pixa",
  "MODIS" = "1pixa",
  "MODISTERRA" = "1pixa",
  "LANDSAT" = "3x3a",
  "LANDSATMAX" = "3x3a"
)

# Parkland, Boreal, or All CIPHA plots? Species select for CAFI plots?
reg_species_list = "All"

reg_select = "All"

# For cafi plots, species select
spec_select = 

# Mortality type?
mtype = "biomass"

# Set exclusion plots, presumably some judgment call
exc_list = list(
  "cipha" = c("HAR1","HAR2","HAR3"),
  "cafi" = c()
)
# Null value
null_val = -32768

#### Some quick variable definition based on the inputs ####
#Set path lead based on whether we're on lcomp or not:
if(NASAcomp_run) {
  leadpath = "/att/nobackup/bmroger1/"
} else {
  leadpath = "/Volumes/"
}

# Source heatmap functions
#source("/home/bmroger1/r_scripts/sol/code/analysis/functions/create_heatmap_sol.R")
source("H:/bmroger1/r_scripts/sol/code/analysis/functions/create_heatmap_sol.R")
# Set full y_range
y_range = 1981:2015

# Set fire threshold as % of plot burned with set timewindow before observation
burn_thresh = 1 #Percent of pixel burned
burn_window = 50 #Years since burn

# Set rule for missing years. 
miss_rule = c(25,2)

#### Functions ####

# Get list of sites to eliminate
pick_sites = function(site_names,reg_names) {
  elim_names = c()
  if(reg_select != "All") {
    elim_names = c(elim_names,site_names[reg_names!=reg_select])
  }
  
  elim_names = c(elim_names, exc_list[[mort_set]])
  keep_names = site_names[!(site_names %in% elim_names)] 
  
  return(keep_names) 
}

# Clean up data frames
clean_df = function(in_df) {
  # Eliminate any exclusion sites
  in_df = in_df[(rownames(in_df) %in% keep_site_list),]
  
  # Eliminate null vals
  in_df[in_df == null_val] = NA
  
  # #Eliminate columns region and ref
  # elim_cols = which(colnames(in_df) %in% c("region","ref"))
  # if(length(elim_cols>0)) {
  #   in_df = in_df[ ,-elim_cols]
  # }
  
  # Match up all column names
  add_years = y_range[!(paste0("X",y_range) %in% colnames(in_df))]
  if(length(add_years)>0) {
    in_df[,paste0("X",add_years)] = NA
  }
  
  # Ensure column order is set properly
  in_df = in_df[,paste0("X",y_range)]
  
  return(in_df)
}

# Process burn df
#Complicated, but essentially returns 1s for all sites that were not burned within 50 years and NAs for any that were
process_burn = function(in_df) {
  burn_site_list = as.factor(in_df[,"Plot_ID"])
  burn_windows = lapply(y_range,function(x) ((x-burn_window):x))
  unburned_df = sapply(burn_windows,function(x) rowSums(in_df[,paste0("X",as.character(x)[paste0("X",as.character(x)) %in% colnames(in_df)])]>burn_thresh))==0
  ### Notes on above                                                       
  unburned_df<-data.frame(in_df[,"Plot_ID"],unburned_df*1)
  colnames(unburned_df) = c("Plot_ID",paste0("X",y_range))
  unburned_df[unburned_df==0]=NA
  #unburned_df = data.matrix(ddply(data.frame(unburned_df),.(Plot_ID),colwise(mean)))
  unburned_df = unburned_df[,-1]
  colnames(unburned_df) = paste0("X",y_range)
  rownames(unburned_df) = burn_site_list
  return(unburned_df)
}


# Shift vi results for cafi
shift_vi = function(in_df,shift_fraction) {
  temp_shift_df = round(shift_df*shift_fraction,0)
  new_df = in_df
  for(r in 1:dim(in_df)[1]) {
    for(y in y_range) {
      if(!is.na(temp_shift_df[r,paste0("X",y)])) {
        new_y = y-temp_shift_df[r,paste0("X",y)]
        new_df[r,paste0("X",y)] = in_df[r,paste0("X",new_y)]
      }
    }
  }
  return(new_df)
}

#### Mortality read ####
# Read in mortality data for each
if (use_year_cutoff){
  raw_mort_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/",year_cutoff,"_as_cutoff/psp/mort_biomass_",year_cutoff,".csv"),row.names = "Plot_ID")
  raw_damaged_df<-read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/",year_cutoff,"_as_cutoff/psp/psp_flags_burns_",year_cutoff,".csv"),stringsAsFactors = FALSE,row.names = "Plot_ID")
}else{
  raw_mort_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/psp/mort_biomass.csv"),row.names = "Plot_ID")
  raw_damaged_df<-read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/psp/psp_flags_burns.csv"),stringsAsFactors = FALSE,row.names = "Plot_ID")
}
# Clean up dataframe
# keep_site_list = pick_sites(rownames(raw_mort_df),raw_mort_df["region"])
keep_site_list = rownames(raw_mort_df)[!(rownames(raw_mort_df) %in% exc_list)] #Creates list that removes some sites
mort_df = clean_df(raw_mort_df)

#Read in damage dataframe from flags in field data

raw_damaged_df<-raw_damaged_df[,-1]
raw_damaged_df[is.na(raw_damaged_df)]=0
burn_windows = lapply(y_range,function(x) ((x-burn_window):x))
undamaged_df<- sapply(burn_windows,function(x) rowSums(raw_damaged_df[,paste0("X",as.character(x)[paste0("X",as.character(x)) %in% colnames(raw_damaged_df)])]>0))==0
undamaged_df=undamaged_df*1
colnames(undamaged_df)=paste0("X",y_range)
undamaged_df[undamaged_df==0]=NA
undamaged_df<-clean_df(undamaged_df)


# #### Select topmort only ####
# topmort = 3
# for(r in 1:(dim(mort_df)[1])) {
#   mort_df[r,head(order(mort_df[r,],na.last=F),-topmort)] = NA
# }

#Read in lengths of time since last survey
raw_shift_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/psp/survey_interval.csv"),row.names = "Plot_ID")
shift_df = clean_df(raw_shift_df)

#Remove data that has invervals longer than max_interval
shift_df[shift_df>max_interval]=NA

#remove sites from mort_df that were measured more than max_interval years after the previous measurement
keep_max_interval=shift_df
keep_max_interval[!is.na(keep_max_interval)]=1
mort_df=mort_df*keep_max_interval #At max interval of 5, reduces number of sites from 6550 to 2324

#### Create DF and Populate ####
mort_vi_df = data.frame("mort" = unlist(mort_df))

#### Read and Process VI Results ####
# Populate with premade csvs
for(data_set in data_set_vec) {#Loops through data sets
  # Read burned data
  if (use_year_cutoff){
    raw_burn_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/",year_cutoff,"_as_cutoff/fire/",data_set,"_fire",year_cutoff,".csv"))
  }else{
     raw_burn_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/fire/",data_set,"_fire.csv"))
  }
  ##Creates dataframe called raw_burned_df which has the fire data
  unburned_df = clean_df(process_burn(raw_burn_df))
  #Complicated, but essentially returns 1s for all sites that were not burned within 50 years and NAs for any that were
  for(an_type in analysis_list[[data_set]]) {#Loop through each analysis
    for(metric in sub_metrics[[an_type]]) {#Loop through each metric
      for(tl_length in tl_list[[data_set]][[an_type]]) {#Loop through each timelag
        temp_df = read.csv(paste0(leadpath,"scooperdock/GIMMS_Tests/data/ndvi_metrics/",data_set,"/",an_type,"/",metric,"_tlength",tl_length,".csv"),row.names=1)
        temp_df = clean_df(temp_df)
        # Remove burned sites
        temp_df = temp_df*unburned_df #Multiplies by 1 or 0, thus removing burned sites
        temp_df = temp_df*undamaged_df #Removes damaged sites
        
        mort_vi_df[paste(data_set,an_type,metric,tl_length,sep="_")] = unlist(temp_df)
      }
    }
  }
}


#### Plot Creation ####

# Create heatmaps
for(data_set in data_set_vec) {#loops through data sets
  for(an_type in analysis_list[[data_set]]) {#loops through analyses
    for(metric in sub_metrics[[an_type]]) {#loops through metrics
      for(metric_thresh in metric_thresh_list[[an_type]]) {#loops through metric thresholds
        # Create outpaths for pngs
        outpath_sub_dir = paste0(leadpath,"scooperdock/GIMMS_Tests/data/plots/pretty/",data_set,"/AllObs/",an_type,"/",outpath_parent_dir,"/",mtype,"/")
        dir.create(outpath_sub_dir,showWarnings=F,recursive=T)
        outpath_pngs = outpath_func_list[[an_type]](outpath_sub_dir,mtype,metric,"site",reg_species_list,
                                                    burn_thresh,burn_window,metric_thresh,pix_select_list[[data_set]])
        print(outpath_pngs[[1]])
        # Continuous
        create_heatmap_cont(mort_vi_df,data_set,an_type,metric,
                            tl_list[[data_set]][[an_type]],mort_thresh_list,metric_thresh,
                            outpath_pngs[["cont"]],label_heatmap,fig_dims[[an_type]][1],fig_dims[[an_type]][2],col_range_list[[metric]])
        # Percentage
        if(perc_list[[an_type]]) {
          create_heatmap_perc(mort_vi_df,data_set,an_type,metric,
                              tl_list[[data_set]][[an_type]],mort_thresh_list,metric_thresh,
                              outpath_pngs[["perc"]],label_heatmap,fig_dims[[an_type]][1],fig_dims[[an_type]][2])
        }
      }
    }
  }
}
