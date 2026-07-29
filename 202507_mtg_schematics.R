rm(list=ls())
setwd("C:/Users/dita5057/OneDrive - UCB-O365/Gryffindor/Fraginae_mitochondria/Manuscript/mtg_annotation")
library(MetBrewer)
library(scales)
library(extrafont)
windowsFonts(A=windowsFont("Seaford"))
windowsFonts(B=windowsFont("Gill Sans MT"))


# FFRA ----------------------------------------------------------------------------
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Ffra_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( "#00906aaa", "#ff6167ee","#8f1400ee","#440154FF","#2d71fcaa","#00000088","#285d59bb","#ff9f00", "#689cd6aa", "black", "#011640dd", "#f09500aa" )

# warna <- c( "#29AF7Faa","#DD5E66FF", "#DD5E66FF","#00000044","#285d59bb","#ff9f00bb","#00000088","#2d71fcaa", "black", "#011640dd", "#f09500aa", "#2d71fcaa",  "#689cd6aa" )


warna <- c( "#00000000","#285d59bb","#29AF7Faa","#00000044", "#2d71fcaa", "black","#f09707bb","#9a75d1aa","#f09707bb", "#011640dd","#f09500aa", "#2d71fcaa",  "#689cd6aa" )

# c("#27187e", "#758bfd","#AEB8FE")

frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*3))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,9, 
       expression(paste(italic('Fragum fragum'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,8,
         unique(dat$Type),
         fill=alpha(warna, 1), 
         horiz=T, cex=0.75,
         text.col ="black"
  )
}


# dev.off()


# FWHI------------------------------------------------------------------------------
rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Fwhi_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088",
            # "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "black", "#011640dd", "#ffaf00bb" )


warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044",  "#2d71fcaa", "black", "#011640dd", "#78f0e2aa" ,"#9AB082ff","#9a75d1aa","#f09707bb", "#f09707bb","#f09707bb", "#f0950a99" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,9, 
       expression(paste(italic('Fragum whitleyi'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,8,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.45,
         text.col ="black"
  )
}


# dev.off()


# TGIG ------------------------------------------------------------------------------

### TGIG

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Tgig_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
#   "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044", "#2d71fcaa", "#011640cc","#9AB08fcc","#9a75d1aa","#f09707bb","#f09707bb","#f09707bb", "#78f0e2aa" , "#f0950a99" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Tridacna gigas'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.65,
         text.col ="black"
  )
}


# dev.off()



# TCRO ------------------------------------------------------------------------------

### TCRO

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Tcro_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
  "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044","#2d71fcaa", "#011640cc","#21bdd988","#9AB08fcc", "#f09707bb","#9a75d1aa","#f09707bb","#f09707bb","#78f0e2aa" , "#f0950a99" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Tridacna crocea'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.6,
         text.col ="black"
  )
}


# dev.off()




# FSUE ------------------------------------------------------------------------------

### FSUE

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Fsue_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
#   "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044","#00599499", "#689cd6aa","#011640cc","#fd8d00aa","#9fc2cc","#2d71fcaa", "#2712fedd","#9AB08fcc","#f09707aa", "#fcdc83aa", "#9a75d1aa","#f0950a99" ,"#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Fragum sueziense'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.55,
         text.col ="black"
  )
}


# dev.off()




# HHIP------------------------------------------------------------------------------

### HHIP

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/Hhip_schematics.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
#   "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044","#2d71fcaa","#f09707bb","#f09707dd","#f29600bb","#f09707dd", "#9a75d1aa","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Hippopus hippopus'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.7,
         text.col ="black"
  )
}


 # dev.off()



# FMUT ------------------------------------------------------------------------------

### FMUT

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Fmut_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values

warna <- c( "#00000000","#285d59bb","#29AF7Faa", "#9AB08fcc","black", "#00000044","#2d71fcaa","#011640cc","#9a75d1aa","#f09707bb","#f09707bb","#f09707bb","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*3))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Fulvia mutica'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.6,
         text.col ="black"
  )
}


# dev.off()



# CEDU ------------------------------------------------------------------------------

### CEDU

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Cedu_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
  "#2d71fcaa", "#00000088","#285d59bb","#000000dd","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "#000000dd","#9AB08fcc","#00000044", "#9a75d1aa", "#f09707bb","#f09707bb","#f09707bb","#f09707bb","#2d71fcaa","#011640cc","#f29600bb","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*4))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Cerastoderma edule'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.7,
         text.col ="black"
  )
}


# dev.off()




# ATUB ------------------------------------------------------------------------------

### ATUB

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Atub_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
  "#2d71fcaa", "#00000088","#285d59bb","#000000dd","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "#000000dd","#9AB08fcc","#00000044","#2d71fcaa", "#9a75d1aa","#f09707bb","#f09707bb","#f09707bb","#011640cc","#f29600bb","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000


{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*4))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.7
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Acanthocardia tuberculata'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.6,
         text.col ="black"
  )
}


# dev.off()



# TDER------------------------------------------------------------------------------

### TDER

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/20240831_Tder_schematic.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values
# warna <- c( #"#00906aaa", "#ff6167ee","#8f1400ee","#00000088","#00000088",
#   "#2d71fcaa", "#00000088","#285d59bb","black","#f29600", "#689cd6aa", "#011640dd", "black", "#ffaf00bb" )



warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044", "#2d71fcaa","#011640cc","#9AB08fcc","#9a75d1aa","#f09707bb","#f09707dd","#f29600bb","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*2))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Tridacna derasa'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.7,
         text.col ="black"
  )
}


# dev.off()





# CGLA------------------------------------------------------------------------------

### CGLA


rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/Cgla_schematics.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values


warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#9AB08fcc","#00000044", "#9a75d1aa","#f09707bb","#f09707dd","#f29600bb","#f09707dd","#9fc2cc", "#2712fedd", "#78f0e2aa" ,"#AEB8FEFF" ,"#f0950a99","#f09707dd" )


frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000

{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*4))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.8
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Cerastoderma glaucum'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.7,
         text.col ="black"
  )
}



# dev.off()



# AMED ------------------------------------------------------------------------------

### AMED

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/Amed_schematics.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values

warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044",  "#011640dd", "#2d71fcaa", "#689cd6aa","#f09707bb", "#f09707bb","#f09707bb", "#f0950a99","#9a75d1aa" )



frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000


{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*4))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.7
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Americardia media'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.6,
         text.col ="black"
  )
}


# dev.off()




# LCRA ------------------------------------------------------------------------------

### LCRA

rm(list=ls())
dat <- read.csv("./ORF_TE_search/20240831_schematics_updated/Lcra_schematics.csv")

head(dat)

types <- unique(dat$Type)
typeid <- c()

for (jj in 1:nrow(dat)){
  tmp <- dat[jj,]
  for (ii in 1:length(types)){
    if (tmp$Type == types[ii]){
      typeid <- c(typeid,ii)
    }
  }
}

dat <- data.frame(dat,typeid)

Yax <- dat$Y1 ##Input vector for y-axis values

warna <- c( "#00000000","#285d59bb","#29AF7Faa", "black","#00000044",  "#2d71fcaa", "#f09707bb", "#f09707bb","#f09707bb", "#9a75d1aa" )



frame_start = 0
frame_stop = max(dat$Stop)+2000
frame_interval = 2000


{ 
  plot.new()
  par(oma=c(2,1,2,1))
  par(mar=c(2,1,2,1) + 0.2)
  par(family="Seaford")
  par(bg="white")
  plot.window(xlim = c(frame_start,frame_stop), ylim = c(0,length(unique(Yax))*4))
  
  
  axis(1, pos = 0, at = seq(frame_start,frame_stop,by=frame_interval),
       #labels=F,
       las=2, tck = 0.02, cex.axis=0.7
       # col = "white",        # axis line colour
       # col.ticks = "white",  # tick colour
       # col.axis = "white"    #label colour
  )
  axis(2, pos = 0,at = seq(0,length(unique(Yax)),by=1),labels=F,tick=F)
  
  # we will modify
  width_rec = 0.25
  
  for (i in 1:nrow(dat)) {
    start = dat$Start[i]
    stop = dat$Stop[i]
    length_rec = abs(start - stop)
    length_tri = 100
    
    if (length_tri > length_rec) {
      length_tri = length_rec
      flag = "small"
    } else {
      length_rec = abs(length_rec - length_tri)
      flag = "large"
    }
    
    if (start < stop) { # flagging right
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                border=NA, col = warna[dat$typeid[i]])
      } else {
        # first our rectangle
        polygon(x = c(start, start, start + length_rec, stop, start + length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start + length_rec, start + length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
    
    if (start > stop) { #flagging left
      
      if (flag == "small") {
        polygon(x = c(start, start, stop), y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
                col = warna[dat$typeid[i]], border = NA)
      } else {
        # first our rectangle
        polygon(x = c(start, start, start - length_rec, stop, start - length_rec),
                y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i] + width_rec, Yax[i], Yax[i] - width_rec),
                col = warna[dat$typeid[i]], border = NA)
        # now our triangle
        # polygon(x = c(start - length_rec, start - length_rec, stop),
        #         y = c(Yax[i] - width_rec, Yax[i] + width_rec, Yax[i]),
        #         col = alpha(warna[dat$typeid[i]], 0.55), border = NA)
      }
    }
  }
  
  
  # using segments instead:
  # segments(x0 = dat$Start, y0 = Yax,
  #          x1 = dat$Stop, y1 = Yax,
  #         lwd=9,
  #         lend=1,
  #         col=alpha(warna[dat$typeid], 0.5))
  text(frame_start,10, 
       expression(paste(italic('Laevicardium crassum'),' Tandem Repeat and Sequencing Technology Comparison')), 
       adj=0, cex=1.2 
       # col = "white"
  )
  # text(seq(1500,11000,by = 1000), par("usr")[3],
  #      srt = 45, adj = 1, xpd = TRUE,
  #      labels = seq(1500,11000,by = 1000), 
  #      cex = 1)
  
  legend(frame_start,9,
         unique(dat$Type),
         fill=warna, 
         horiz=T, cex=0.6,
         text.col ="black"
  )
}


# dev.off()


