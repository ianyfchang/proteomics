# generate chromatogame plots
library(scales)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(doSNOW)
cl <- makeSOCKcluster(16) # number of CPU cores
registerDoSNOW(cl)

outputDir="output_dir"
expand_boundary=0.7
legend_text_size=12
pdf_ncol = 6
pdf_nrow = 6
pdf_width = 23.4 # A1 paper width in inch
pdf_height = 33.1 # A1 paper height in inch
pdf_title_font_size = 24
relative_line_thicks = 1.2 # 110%
chromatogram_from_skyline="chromatogram_exported_from_skyline.tsv"
boundary_from_skyline="boundary_exported_from_skyline.csv"
uniprot_id=TRUE # use uniport protein name?
number_of_time_axis_breaks=30
empty_step_of_time_axis_breaks=5
number_of_intensity_axis_breaks=30
empty_step_of_intensity_axis_breaks=5
use_fix_intensity_scale=T

# color for fragment ion
colors<-c("#0000FF","#A52A2A","#8A2BE2","##e31a1c","#33a02c","#ff7f00" )


# read chromatogram & boundary file
chroma<-fread(chromatogram_from_skyline,sep="\t",header=T,stringsAsFactors = F,na.string="#N/A")
# clean isotopelabeltype
chroma<-subset(chroma, IsotopeLabelType %in% c('light','heavy'))

boundary<-read.delim(boundary_from_skyline,sep=",",header=T,stringsAsFactors = F,na.strings = "#N/A")
# get peptide:protein name mapping
peptide_protein<-unique(boundary[,c("Peptide.Sequence", "Protein.Name")])

if(uniprot_id){
  peptide_protein$Protein.Name<-gsub("^.+\\|","",peptide_protein$Protein.Name)
}
dir.create(paste(outputDir,"/single_chromatogram",sep=""))

# create transition id column
chroma$Transition.id<-paste(chroma$PrecursorCharge,chroma$FragmentIon,chroma$ProductCharge,sep=".")
boundary$Transition.id<-paste(boundary$Precursor.Charge,boundary$Fragment.Ion,boundary$Product.Charge,sep=".")


# identify unique filename, peptide, transition.id
tmp_peptide<-unique(chroma[,c('FileName','PeptideModifiedSequence','Transition.id','FragmentIon')])
tmp_peptide<-tmp_peptide[order(tmp_peptide$FileName,tmp_peptide$PeptideModifiedSequence, tmp_peptide$Transition.id),]
tmp_fragmention    <- aggregate( FragmentIon ~ FileName+PeptideModifiedSequence, tmp_peptide, paste0, collapse=",")
unique_wif_peptide <- aggregate( Transition.id ~ FileName+PeptideModifiedSequence, tmp_peptide, paste0, collapse=",")
unique_wif_peptide <- merge(unique_wif_peptide, tmp_fragmention, by=c("FileName","PeptideModifiedSequence"))
unique_wif_peptide <- merge(unique_wif_peptide, peptide_protein, by.x="PeptideModifiedSequence",by.y="Peptide.Sequence",all.x=T)
unique_wif_peptide <- unique_wif_peptide[order(unique_wif_peptide$Protein.Name, unique_wif_peptide$FileName),]
row.names(unique_wif_peptide) <- 1:dim(unique_wif_peptide)[1]
protein_fragmention_for_paging <- unique(unique_wif_peptide[,c("Protein.Name","PeptideModifiedSequence","Transition.id","FragmentIon")])



# SavitzkyGolaySmooth
# https://sourceforge.net/p/proteowizard/code/HEAD/tree/trunk/pwiz/pwiz_tools/Skyline/
# skyline-source-code/proteowizard-code/pwiz/pwiz_tools/Skyline/Model/Results/ChromHeaderInfo.cs
# rewirte sg smooth to rcpp later
SavitzkyGolaySmooth<-function(intensities){
  # if (length(intensities) < 9)
  #   return(intensities)
  intRaw = intensities
  intSmooth = intRaw
  # 1:4 is original data
  tmp<-sapply(5:(length(intRaw) - 4),function(i){
    sum = 59 * intRaw[i] +
      54 * (intRaw[i - 1] + intRaw[i + 1]) +
      39 * (intRaw[i - 2] + intRaw[i + 2]) +
      14 * (intRaw[i - 3] + intRaw[i + 3]) -
      21 * (intRaw[i - 4] + intRaw[i + 4])
      return(sum / 231)
  })
  # last 4 data are the same
  #Array.Copy(intRaw, intRaw.Length - 4, intSmooth, intSmooth.Length - 4, 4);
  intSmooth[5:(length(intSmooth)-4)]<-tmp
  return(intSmooth);
}

get_peptide_max_light_heavy_intensity<-function(pro_frag){
  profrag<-as.list(pro_frag)
  # get boundary first
  tmp_transitionid<-unlist(strsplit(profrag$Transition.id,","))
  tmp_boundary<-subset(boundary,Peptide.Sequence==profrag$PeptideModifiedSequence & Transition.id %in% tmp_transitionid)
  tmp_melted_chroma_light<-subset(melted_chroma,PeptideModifiedSequence==profrag$PeptideModifiedSequence & Transition.id %in% tmp_transitionid & IsotopeLabelType=="light")
  tmp_melted_chroma_heavy<-subset(melted_chroma,PeptideModifiedSequence==profrag$PeptideModifiedSequence & Transition.id %in% tmp_transitionid & IsotopeLabelType=="heavy")
  
  tmp_merged_chroma_boundary_light<-merge(tmp_melted_chroma_light, tmp_boundary[,c('File.Name','Peptide.Sequence','Transition.id','light.Start.Time','light.End.Time')], 
                                          by.x=c("FileName","PeptideModifiedSequence","Transition.id"),by.y=c("File.Name","Peptide.Sequence","Transition.id"),all.x=T)
  tmp_merged_chroma_boundary_heavy<-merge(tmp_melted_chroma_heavy, tmp_boundary[,c('File.Name','Peptide.Sequence','Transition.id','heavy.Start.Time','heavy.End.Time')], 
                                          by.x=c("FileName","PeptideModifiedSequence","Transition.id"),by.y=c("File.Name","Peptide.Sequence","Transition.id"),all.x=T)
  tmp_light_intensity<-subset(tmp_merged_chroma_boundary_light, Times>= light.Start.Time & Times<=light.End.Time,select="IntensitiesSmoothed")
  tmp_heavy_intensity<-subset(tmp_merged_chroma_boundary_heavy, Times>= heavy.Start.Time & Times<=heavy.End.Time,select="IntensitiesSmoothed")
  if(dim(tmp_light_intensity)[1]!=0){
    profrag$light_max<- max(tmp_light_intensity,na.rm=T)
  }else{
    profrag$light_max=NA
  }
  if(dim(tmp_heavy_intensity)[1]!=0){
    profrag$heavy_max<- max(tmp_heavy_intensity,na.rm=T)
  }else{
    profrag$heavy_max=NA
  }
  # get intensities between boundary for every chromatogram
  return(data.frame(profrag,stringsAsFactors = F))
}

every_nth <- function(x, nth, empty = TRUE, inverse = T) # claning spcing axis label
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

plot_light_or_heavy<-function(data,light=T,
                              peak_time=0,
                              #peak_intensity=0,
                              start.time=0,
                              end.time=0,
                              global_max_intensity){
  if(light==T){
    title="light"
  }else{
    title="heavy"
  }
  # max intensities between boundary
  if(length(start.time)!=0){
    # test ks.test of every fragment ion
    
    peak_intensity_for_arrow<-max(subset(data, Times>=start.time & Times<=end.time, select=IntensitiesSmoothed),na.rm = T)
    if( use_fix_intensity_scale==F){
      y_max<-peak_intensity_for_arrow
    }
  }else{
    peak_intensity_for_arrow<-numeric(0)
  }

  x_min<-min(data$Times,na.rm=T)
  x_max<-max(data$Times,na.rm=T)
  
  x_axis_break<-pretty_breaks(number_of_time_axis_breaks)(floor(x_min):ceiling(x_max))
  x_axis_label<-every_nth(x_axis_break,empty_step_of_time_axis_breaks)
  
  y_min<-0
  y_max <- max(data$IntensitiesSmoothed) # if data have no boundary, use max intensity
  
  # if data have boundary
  if(length(start.time)!=0 && use_fix_intensity_scale==F){
      y_max<-peak_intensity_for_arrow
  }
  
  if(use_fix_intensity_scale==T && is.na(global_max_intensity)==F){
    y_max<-global_max_intensity
  }
  
  if (y_max>=1000){
    power<-floor(log10(y_max))- 1
    ylabel<-paste("Intensity (10^",power,")", sep="")
    data$IntensitiesSmoothed<-data$IntensitiesSmoothed/(10^power)
    peak_intensity_for_arrow<-peak_intensity_for_arrow/(10^power)
    if(use_fix_intensity_scale || is.na(global_max_intensity)==F){
      y_max<-y_max/(10^power)
    }else{
      y_max<-max(data$IntensitiesSmoothed)
    }
  }else{
    ylabel<-"Intensity"
  }
  
  # create y axis breaks & labels
  y_max<-y_max*1.2
  y_axis_break<-pretty_breaks(number_of_intensity_axis_breaks)(y_min:y_max)
  y_axis_label<-every_nth(y_axis_break,empty_step_of_intensity_axis_breaks)
    
  p<-ggplot(data,aes(x=Times,y=IntensitiesSmoothed,color=TransitionId_ProductMz))+
    geom_line(size=rel(relative_line_thicks))+
    scale_x_continuous("Retention Time",limits = c(x_min,x_max),expand=c(0,0), breaks=x_axis_break, labels = x_axis_label)+
    scale_y_continuous(ylabel,          limits = c(y_min,y_max),expand=c(0,0), breaks=y_axis_break, labels = y_axis_label)+
    scale_color_manual(values = colors)+
    guides(color=guide_legend(nrow=2))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.box.background = element_rect(color="black"),
          #legend.spacing = unit(1.5, 'lines'),
          #legend.box.spacing = unit(1.5, 'lines'),
          legend.position = "top",
          legend.justification=c(0, 0),
          legend.text = element_text(size=rel(0.7)),
          legend.title = element_blank())
  
    intensities_between_start_end<-subset(data, Times>=start.time & Times<=end.time, select=IntensitiesSmoothed,drop=T)
    
    if(length(intensities_between_start_end)!=0){
      y_max_between_start_end<-max(intensities_between_start_end,na.rm=T)
    }else{
      y_max_between_start_end<-numeric(0)
    }
    
    y_for_vline<-rep(0,length(start.time)+length(end.time))
    yend_for_vline<-rep(y_max_between_start_end,length(start.time)+length(end.time))
    x_for_vline<-c(start.time,end.time)
    xend_for_vline<-x_for_vline
    
    if(length(start.time)!=0){
      p <- p+ 
        annotate('segment' ,x=x_for_vline,   y=y_for_vline, xend=xend_for_vline, yend=yend_for_vline, color="black",linetype="dashed")+
        annotate('segment', x=peak_time-0.02,y=(peak_intensity_for_arrow*1.05),xend=peak_time-0.01,yend=(peak_intensity_for_arrow*1.03),arrow = arrow(angle=15,type="closed"))+
        annotate('text',x=peak_time,y=(peak_intensity_for_arrow*1.03),label=peak_time, vjust=0,hjust=0)
    }
    
  return(p)
}

generate_chromatogram<-function(wif_peptide,max.heavy.intensity, max.light.intensity){
  wif_peptide<-as.list(wif_peptide)
  
  # get boundary data from boundary file
  peptide_boundary<-subset(boundary,File.Name==wif_peptide$FileName & Peptide.Sequence==wif_peptide$PeptideModifiedSequence)
  
  productMz<-as.double(c(peptide_boundary$light.Product.Mz,peptide_boundary$heavy.Product.Mz))
  
  
  peptide_chromatogram<-subset(melted_chroma, FileName==wif_peptide$FileName & PeptideModifiedSequence == wif_peptide$PeptideModifiedSeq )
  peptide_chromatogram$TransitionId_ProductMz <- paste(peptide_chromatogram$Transition.id,peptide_chromatogram$ProductMz, sep="-")
  #peptide_chromatogram$Times<-round(peptide_chromatogram$Times,digits = 2)
  
  # append space after legend
  heavy_idx<-which(peptide_chromatogram$IsotopeLabelType=="heavy")
  peptide_chromatogram[heavy_idx,'TransitionId_ProductMz']    <-paste(peptide_chromatogram[heavy_idx,'TransitionId_ProductMz'],"(heavy) ",sep=" ")
  light_index<-which(peptide_chromatogram$IsotopeLabelType=="light")
  peptide_chromatogram[light_index,'TransitionId_ProductMz']<-paste(peptide_chromatogram[light_index,'TransitionId_ProductMz'],"(light) ",sep=" ")
  
  
  # find rank 1 heavy.rentention.time
  heavy.retention.time<-peptide_boundary[which(peptide_boundary$heavy.Peak.Rank==1),"heavy.Retention.Time"]
  light.retention.time<-peptide_boundary[which(peptide_boundary$light.Peak.Rank==1),"light.Retention.Time"]
  
  # get start, end times
  light.start.time<-peptide_boundary[which(peptide_boundary$light.Peak.Rank==1),"light.Start.Time"]
  light.end.time<-peptide_boundary[which(peptide_boundary$light.Peak.Rank==1),"light.End.Time"]
  
  heavy.start.time<-peptide_boundary[which(peptide_boundary$heavy.Peak.Rank==1),"heavy.Start.Time"]
  heavy.end.time<-peptide_boundary[which(peptide_boundary$heavy.Peak.Rank==1),"heavy.End.Time"]
  
  mid_point<-heavy.start.time+(heavy.end.time - heavy.start.time)/2
  
  seg<-subset(peptide_chromatogram,
              Times >=  mid_point - expand_boundary & 
              Times <=  mid_point + expand_boundary )

  if(length(heavy.start.time)==0){ # if No heavy boundary, plot all chromatogram
    seg<-peptide_chromatogram
  }
  
  lightSeg<-subset(seg,IsotopeLabelType=="light")
  
  heavySeg<-subset(seg,IsotopeLabelType=="heavy")
  
  lightP<-plot_light_or_heavy(lightSeg,light=T,
                              peak_time=light.retention.time,
                              # peak_intensity=light_peak_intensity,
                              start.time=light.start.time,
                              end.time=light.end.time,
                              global_max_intensity=max.light.intensity)
  heavyP<-plot_light_or_heavy(heavySeg,light=F,
                              peak_time=heavy.retention.time,
                              # peak_intensity=heavy_peak_intensity,
                              start.time=heavy.start.time,
                              end.time=heavy.end.time,
                              global_max_intensity=max.heavy.intensity)
  
  p1<-arrangeGrob(lightP,heavyP,ncol=1, top=textGrob(unique(peptide_boundary$Replicate.Name),gp=gpar(fontface="bold")))
  p2<-rectGrob(gp = gpar(color="#555555"))
  pp<-gTree(children = gList(p1,p2))
  return(pp)
  # add retangel
}

# generate PDFs

generate_pdfs<-function(protein_fragmention,pdf_ncol, pdf_nrow, pdf_width, pdf_height, pdf_title_font_size){
  pro_frag <- as.list(protein_fragmention)
  transition.id.list<-unlist(strsplit(pro_frag$Transition.id,","))
  # subseting chromatogram list from unique_wif_peptide
  tmp_wif_peptide <- subset(unique_wif_peptide,PeptideModifiedSequence == pro_frag$PeptideModifiedSequence & Protein.Name == pro_frag$Protein.Name & Transition.id == pro_frag$Transition.id)
  
  # generate plots
  plot_list<-apply(tmp_wif_peptide,1,generate_chromatogram, max.light.intensity=pro_frag$light_max, max.heavy.intensity=pro_frag$heavy_max)
  
  totalPlot<-dim(tmp_wif_peptide)[1] 
  #plot_per_page<-pdf_ncol*pdf_nrow
  #pages<-ceiling(totalPlot/plot_per_page)
  nlay <- max(matrix(seq_len(pdf_nrow * pdf_ncol)), na.rm = TRUE)
  npages <- totalPlot%/%nlay + as.logical(totalPlot%%nlay)
  groups <- split(plot_list, rep(seq_len(npages), each = nlay, length.out = totalPlot))
  pl <- vector(mode = "list", length = npages)
  nfrag<-length(transition.id.list)
  
  pdf(paste(outputDir,"/single_chromatogram/chromatogram_",pro_frag$Protein.Name,"_",pro_frag$PeptideModifiedSequence,".pdf",sep=""), width = pdf_width, height = pdf_height)
  
  for (g in seq_along(groups)) {
    grobwidth<-NULL
    grobs <- grobTree(
      gp = gpar(fontsize = pdf_title_font_size, fontface = "bold"), 
      textGrob(label = paste(pro_frag$Protein.Name,"_",pro_frag$PeptideModifiedSequence,sep=""), 
               name = "title1", 
               x = unit(0.2, "lines"), 
               y = unit(1, "lines"), 
               hjust = 0, 
               vjust = 0)
    )
    grobwidth<-grobWidth("title1")
    i=1
    for (i in 1:nfrag){
      grobs<-addGrob(grobs, textGrob(label = paste("_", transition.id.list[i],sep=""), 
                                     name = paste("title",(i+1),sep=""), 
                                     x = grobwidth + unit(0.2, "lines"), 
                                     y = unit(1, "lines"), 
                                     hjust = 0, 
                                     vjust = 0, 
                                     gp = gpar(col = colors[i])))
      grobwidth<-grobwidth + grobWidth(paste("title",(i+1),sep=""))
    }
    grobs<-addGrob(grobs, textGrob(label = paste("Page ",g," of ",npages,sep=""), 
                                   name = "titleP", 
                                   x = unit(0.99, "npc"),
                                   y = unit(1, "lines"),
                                   just = "right", 
                                   gp = gpar(col = "#333333")))
    gg <- arrangeGrob(grobs=groups[[g]],ncol=pdf_ncol,nrow=pdf_nrow,top=grobs,padding = unit(5, "line"))
    grid.newpage()
    grid.draw(gg)
  }
  dev.off()
  write(paste(pro_frag$Protein.Name,"chromatogram done!",sep=" "),stdout())
}


##########################################################################################################
write("Melting chromatogram times and intensities",stdout())
start.time<-Sys.time()
tmp_intensity<-lapply(strsplit(chroma$Intensities,","),as.numeric)
tmp_smoothed_intensity<-lapply(tmp_intensity,SavitzkyGolaySmooth)
tmp_times<-lapply(strsplit(chroma$Times,","),as.numeric)
tmp_times<-lapply(tmp_times,round, digits=2)
melted_chroma<-data.frame(FileName=rep(chroma$FileName,sapply(tmp_intensity,length)),
                          PeptideModifiedSequence=rep(chroma$PeptideModifiedSequence,sapply(tmp_intensity,length)),
                          PrecursorCharge=rep(chroma$PrecursorCharge,sapply(tmp_intensity,length)),
                          ProductMz=rep(chroma$ProductMz,sapply(tmp_intensity,length)),
                          FragmentIon=rep(chroma$FragmentIon,sapply(tmp_intensity,length)),
                          ProductCharge=rep(chroma$ProductCharge,sapply(tmp_intensity,length)),
                          IsotopeLabelType=rep(chroma$IsotopeLabelType,sapply(tmp_intensity,length)),
                          TotalArea=rep(chroma$TotalArea,sapply(tmp_intensity,length)),
                          Transition.id=rep(chroma$Transition.id,sapply(tmp_intensity,length)),
                          Times=unlist(tmp_times),
                          Intensities=unlist(tmp_intensity),
                          IntensitiesSmoothed=unlist(tmp_smoothed_intensity))
melted_chroma[which(melted_chroma$IntensitiesSmoothed < 0 ),'IntensitiesSmoothed'] <- 0
Sys.time()-start.time

write("Get max intensities of light and heavy chromatogram between identified boundary",stdout())
start.time<-Sys.time()
protein_fragmention_for_paging <- foreach(i = 1:dim(protein_fragmention_for_paging)[1],.combine='rbind') %dopar% {
  return(get_peptide_max_light_heavy_intensity(protein_fragmention_for_paging[i,]))
}
Sys.time()-start.time

write("Generate PDFs",stdout())
start.time<-Sys.time()
# generate PDFS
foreach(i=c(1:dim(protein_fragmention_for_paging)[1]), .packages=c("scales","ggplot2","data.table","gridExtra","grid")) %dopar% {
  write(paste0(i, protein_fragmention_for_paging[i,], collapse = ", "),stdout())
  generate_pdfs(protein_fragmention=protein_fragmention_for_paging[i,],
                pdf_ncol = pdf_ncol,
                pdf_nrow = pdf_nrow,
                pdf_width = pdf_width,
                pdf_height = pdf_height,
                pdf_title_font_size = pdf_title_font_size)
}
Sys.time()-start.time


##
# unite pdfs
old_merged_pdf<-list.files(path=outputDir,pattern = "merged_chromatogram.pdf",recursive = T )
if(length(old_merged_pdf)>0){
  file.remove(paste(outputDir,"/",old_merged_pdf,sep=""))
}
pdf_files <- paste0(outputDir,"/single_chromatogram/",list.files(path=paste(outputDir,"/single_chromatogram",sep=""),pattern = "*.pdf",recursive = T ),collapse =" ")
try(
  if(system("pdftk",ignore.stderr=T)==99){
    write(paste("merge pdfs by pdftk",sep=""),stdout())
    system(paste("pdftk ",pdf_files," output ",outputDir,"/merged_chromatogram.pdf",sep=""))
  }else if (system("pdfunite",ignore.stderr=T)==99){
    write(paste("merge pdfs by pdfunite",sep=""),stdout())    
    system(paste("pdfunite ",pdf_files," ",outputDir,"/merged_chromatogram.pdf",sep=""))
  }else{
    print("no pdf utilities available: pdfunite (centos 7) or pdftk (windows + linux)")
  })
