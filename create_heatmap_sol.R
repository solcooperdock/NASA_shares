#####
# scooperdock@whrc.org
# Modified from Kylen's script
# Functions for creating heatmaps (both continuous and percent based)
#####

# Create outpaths
#Functions which tell where to write the .png files of heatmaps and what to call them
create_outpaths_trend = function(outdir,m,met,plot_lev,reg_spec_text,bthresh,bwin,met_thresh,pix_select) {
  return( 
    list( 
      "perc" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","ViThresh=",met_thresh,"_bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_",met,".png"),
      "cont" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","ViThresh=",met_thresh,"_bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_",met,"Cont.png")
    )
  )
}
create_outpaths_pulse = function(outdir,m,met,plot_lev,reg_spec_text,bthresh,bwin,met_thresh,pix_select) {
  return( 
      list( 
        "perc" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","ZvalThresh=",met_thresh,"_bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_",met,".png"),
        "cont" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_",met,"Cont.png")
      )
    )
}
create_outpaths_ews = function(outdir,m,met,plot_lev,reg_spec_text,bthresh,bwin,met_thresh,pix_select) {
  return( 
    list( 
      "perc" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","metThresh=",met_thresh,"_bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_ews_",met,".png"),
      "cont" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_ews_",met,"Cont.png")
    )
  )
}
create_outpaths_ddj = function(outdir,m,met,plot_lev,reg_spec_text,bthresh,bwin,met_thresh,pix_select) {
  return( 
    list( 
      "perc" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","metThresh=",met_thresh,"_bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_ddj_",met,".png"),
      "cont" = paste0(outdir,m,"_",plot_lev,"_",reg_spec_text,"_tstart0_","bWindowThresh=",bwin,"-",bthresh,"_",pix_select,"_ddj_",met,"Cont.png")
    )
  )
}

#Creates a list to loop through to use the above functions
outpath_func_list = list(
  "trend" = create_outpaths_trend,
  "pulse" = create_outpaths_pulse,
  "ews" = create_outpaths_ews,
  "ddj" = create_outpaths_ddj
)


# Create percentage heatmap
create_heatmap_perc = function(in_df,ds,an,metric,tl_range,mort_range,met_thresh,out_png,add_labels,fig_w,fig_h) {
  #Function with the following inputs:(mort_vi_df,data_set,an_type,metric,tl_list[[data_set]][[an_type]],
  #mort_thresh_list[[mort_set]],metric_thresh,outpath_pngs[["perc"]],label_heatmap,fig_dims[[an_type]][1],
  #fig_dims[[an_type]][2],col_range_list[[metric]])
  #Color palette
  my_palette = colorRampPalette(c("blue", "white","red"))(n = 299)
  # col_breaks = c(seq(0,33,length=100), # for blue
  #                seq(33.1,66,length=100),  # for white
  #                seq(66.1,100,length=100)) # for red
  
  #For subtracting control
  col_breaks = c(seq(-50,-16.67,length=100), # for blue
                 seq(-16.68,16.67,length=100),  # for white
                 seq(16.68,50,length=100)) # for red
  
  heatmap_mat = matrix(NA,nrow=length(tl_range),ncol=length(mort_range),dimnames=list(rev(tl_range),mort_range))
  #creates an empty matrix 16x26
  cellnote_array = heatmap_mat

  # Populate heatmap array
  for(mort in mort_range) {#Loops through mort_thresh_list
    # Remove NAs
    mort_list = in_df[!is.na(in_df[,"mort"]),"mort"] #This removes NAs from column "mort" in df, confusing because its also looping through a variable called "mort"
    mort_pass_list = mort_list>=mort #this removes any values that are below the threshold in mort_thresh_list
    mort_Nopass_list = mort_list<mort
    
    
    temp_vi_df = in_df[,paste(ds,an,metric,tl_range,sep="_")] #Takes columns from tl_range for specific data_set, analysis combo
    temp_vi_df = sapply(colnames(temp_vi_df),function(x) temp_vi_df[!is.na(in_df[,"mort"]),x]) #removes all sites when no mort obs existed
    colnames(temp_vi_df) = tl_range #sets column names to timelags 5:20
    for(tl in tl_range) {#Loops for lags in timelag range
      # Remove more NAs
      temp_vi_list = temp_vi_df[,as.character(tl)]#All data for each time lag
      #heatmap_mat[as.character(tl),as.character(mort)] = 100*sum(temp_vi_list[mort_pass_list]<met_thresh,na.rm=T)/sum(!is.na(temp_vi_list[mort_pass_list])) #Original
      heatmap_mat[as.character(tl),as.character(mort)] = 100*(sum(temp_vi_list[mort_pass_list]<met_thresh,na.rm=T)/sum(!is.na(temp_vi_list[mort_pass_list])) - sum(temp_vi_list[mort_Nopass_list]<met_thresh,na.rm=T)/sum(!is.na(temp_vi_list[mort_Nopass_list]))) #subtract control
      #Populates each cell of heatmap_mat with value 100*the sum of the number of temp_vi_list values that are above mort threshold and less than the metric threshold
      #divided by the number of data points that have values, ie non NAs (output is % of values below met threshold)
      cellnote_array[as.character(tl),as.character(mort)] = paste0(as.character(round(100*(sum(temp_vi_list[mort_pass_list]<met_thresh,na.rm=T)/sum(!is.na(temp_vi_list[mort_pass_list]))
                                                                    - sum(temp_vi_list[mort_Nopass_list]<met_thresh,na.rm=T)/sum(!is.na(temp_vi_list[mort_Nopass_list])))),0),"\n", 
                                                                    as.character(sum(!is.na(temp_vi_list[mort_pass_list]))))
      #Fills in each cell of cellnote_array with "rounded result of heatmap_mat equation" space "number of valid values"
      }
  }
  
  # Reorder heatmap array and change column names
  colnames(heatmap_mat) = c("All",tail(mort_range,-1)) #Changes column 0 to All
  rownames(heatmap_mat) = paste("0",rev(tl_range),sep="-") #Changes row names to 0-# instead of just #

  # Write out to png
  png(filename=out_png,units="cm",width = fig_w, height = fig_h,res=600)       
  
  if(add_labels) {
    heatmap.2(heatmap_mat,
              Rowv=NA, Colv=NA,density.info="none",dendrogram=c("none"),trace="none",
              cellnote=cellnote_array,notecol="black",notecex=.15,
              col = my_palette,breaks=col_breaks,
              cexRow=.5,cexCol=.5,
              key=F, lwid=c(0.1,4), lhei=c(0.1,4),margins=c(2,2),
              colsep=c(0,ncol(heatmap_mat)),rowsep=c(0,nrow(heatmap_mat)),
              sepcolor="black", sepwidth=c(0,0),
              na.color="grey"
    )
  } else {
    heatmap.2(heatmap_mat,
              Rowv=NA, Colv=NA,density.info="none",dendrogram=c("none"),trace="none",
              col = my_palette,breaks=col_breaks,
              cexRow=.5,cexCol=.5,
              key=F, lwid=c(0.1,4), lhei=c(0.1,4),margins=c(2,2),
              colsep=c(0,ncol(heatmap_mat)),rowsep=c(0,nrow(heatmap_mat)),
              sepcolor="black", sepwidth=c(0,0),
              na.color="grey"
              )
  }
  dev.off()
}

# Create continuous heatmap
create_heatmap_cont = function(in_df,ds,an,metric,tl_range,mort_range,met_thresh,out_png,add_labels,fig_w,fig_h,col_range) {
  
  heatmap_mat = matrix(NA,nrow=length(tl_range),ncol=length(mort_range),dimnames=list(rev(tl_range),mort_range))
  cellnote_array = heatmap_mat
  
  # Populate heatmap array
  for(mort in mort_range) {
    # Remove NAs
    mort_list = in_df[!is.na(in_df[,"mort"]),"mort"]
    mort_pass_list = mort_list>=mort
    mort_Nopass_list = mort_list<mort
    
    temp_vi_df = in_df[,paste(ds,an,metric,tl_range,sep="_")] 
    temp_vi_df = sapply(colnames(temp_vi_df),function(x) temp_vi_df[!is.na(in_df[,"mort"]),x])
    colnames(temp_vi_df) = tl_range
    for(tl in tl_range) {
      # Remove more NAs
      temp_vi_list = temp_vi_df[,as.character(tl)]#All data for each time lag
      #heatmap_mat[as.character(tl),as.character(mort)] = mean(temp_vi_list[mort_pass_list],na.rm=T) #Original, Average of that time lag
      heatmap_mat[as.character(tl),as.character(mort)] = mean(temp_vi_list[mort_pass_list],na.rm=T) - mean(temp_vi_list[mort_Nopass_list],na.rm=T) #modification, subtract out control plots
      if(an!="ddj") {
      cellnote_array[as.character(tl),as.character(mort)] = paste0(as.character(round(mean(temp_vi_list[mort_pass_list],na.rm=T) - mean(temp_vi_list[mort_Nopass_list],na.rm=T),2)),
                                                                   "\n","(",as.character(sum(!is.na(temp_vi_list[mort_pass_list]))),")")
      
      } else {
        cellnote_array[as.character(tl),as.character(mort)] = paste0(format(mean(temp_vi_list[mort_pass_list],na.rm=T) - mean(temp_vi_list[mort_Nopass_list],na.rm=T),digits=2,scientific=T),
                                                                     "\n","(",as.character(sum(!is.na(temp_vi_list[mort_pass_list]))),")")
      }
    }
  }
  
  # Reorder heatmap array and change column names
  colnames(heatmap_mat) = c("All",tail(mort_range,-1))
  rownames(heatmap_mat) = paste("0",rev(tl_range),sep="-")
  
  #Color palette
  if(an == "ews" | an == "ddj") {
    col_breaks = seq(min(heatmap_mat,na.rm=T),max(heatmap_mat,na.rm=T),length=300)
    my_palette = colorRampPalette(c("blue", "white","red"))(n = 299)
  } else {
    col_breaks = seq(col_range[1],col_range[2],length=300)
    if(col_range[3]) { # reverse palette
      my_palette = colorRampPalette(c("red", "white","blue"))(n = 299)
    } else {
      my_palette = colorRampPalette(c("blue", "white","red"))(n = 299)
    }
  }
  
  # Make plot and write out to png
  png(filename=out_png,units="cm",width = fig_w, height = fig_h,res=600)       
  
  if(add_labels) {
    heatmap.2(heatmap_mat,
              Rowv=NA, Colv=NA,density.info="none",dendrogram=c("none"),trace="none",
              cellnote=cellnote_array,notecol="black",notecex=.15,
              col = my_palette,breaks=col_breaks,
              cexRow=.5,cexCol=.5,
              key=F, lwid=c(0.1,4), lhei=c(0.1,4),margins=c(2,2),
              colsep=c(0,ncol(heatmap_mat)),rowsep=c(0,nrow(heatmap_mat)),
              sepcolor="black", sepwidth=c(0,0),
              na.color="grey"
    )
  } else {
    heatmap.2(heatmap_mat,
              Rowv=NA, Colv=NA,density.info="none",dendrogram=c("none"),trace="none",
              col = my_palette,breaks=col_breaks,
              cexRow=.5,cexCol=.5,
              key=F, lwid=c(0.1,4), lhei=c(0.1,4),margins=c(2,2),
              colsep=c(0,ncol(heatmap_mat)),rowsep=c(0,nrow(heatmap_mat)),
              sepcolor="black", sepwidth=c(0,0),
              na.color="grey"
    )
  }
  dev.off()
}

