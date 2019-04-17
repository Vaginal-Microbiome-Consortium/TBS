## ------------------------------------------------------------------------
mybarplot <- function(mypropdata, plottype="ggplot2", addlegend=FALSE, dxbv=NULL, colScheme="classic", alphaDiversity=NULL, outpath=NULL) { # pass proportion data, samples on rows, taxa on columns.

  mypropdata <- mypropdata[,apply(mypropdata,2,sum) > 0]

  # assign vagitype
  mytypes <- apply(mypropdata,1,which.max)
  maxprop <- mypropdata[matrix(c(1:nrow(mypropdata),mytypes), ncol=2)]
  mytypes <- colnames(mypropdata)[mytypes]
  mytypes[maxprop < 0.3] <- "No Type"
  
  # set type order
  uniqueTypes <- unique(mytypes)
  lactoTypes <- c(
                  grep("crispatus", uniqueTypes, value=TRUE),
                  grep("delbrueckii", uniqueTypes, value=TRUE),  
                  grep("jensenii", uniqueTypes, value=TRUE), 
                  grep("gasseri", uniqueTypes, value=TRUE), 
                  grep("iners", uniqueTypes, value=TRUE)
                    ) 
  otherLacto <- setdiff(grep("Lactobacillus", uniqueTypes, value=TRUE), lactoTypes)
  nonLactoTypes <- c(
                     grep("Gardnerella", uniqueTypes, value=TRUE), 
                     grep("BVAB", uniqueTypes, value=TRUE), 
                     grep("Sneathia amnii", uniqueTypes, value=TRUE)
                     )
  myTypeOrder <- c(lactoTypes, otherLacto, nonLactoTypes)
  myTypeOrder <- c(myTypeOrder, setdiff(uniqueTypes,myTypeOrder))
  if (length(grep("No Type", myTypeOrder)) > 0) {
    myTypeOrder <- c(myTypeOrder[-grep("No Type", myTypeOrder)], "No Type")
  }
  myTypeOrder <- data.frame(typeOrder=c(1:length(myTypeOrder)),row.names=myTypeOrder)
  
  mypropdata <- 100*mypropdata
 
  # order data by type then proportion
  mypropdata <- mypropdata[order(myTypeOrder[mytypes,],-maxprop),]
  barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]


  

  # Hardik PTB
  if (colScheme == "hardik") {
    mycolors <- rep("#faf0e6", ncol(mypropdata)+1)
    mycolors[grep("Clostridiales.BVAB2", colnames(mypropdata))] <- "#dee0e5"
    mycolors[grep("Coriobacteriaceae.OTU27", colnames(mypropdata))] <- "#bdc2cc"
    mycolors[grep("Megasphaera.OTU71.type2", colnames(mypropdata))] <- "#60636a"
    mycolors[grep("Aerococcus.christensenii", colnames(mypropdata))] <- "#a5acaf"
    mycolors[grep("Streptococcus.cluster29", colnames(mypropdata))] <- "#9edae5"
    mycolors[grep("Dialister.cluster51", colnames(mypropdata))] <- "#dbdb8d"
    mycolors[grep("Dialister", colnames(mypropdata))] <- "#dbdb8d"
    mycolors[grep("TM7.OTU.H1", colnames(mypropdata))] <- "#e377c2"
    mycolors[grep("Mycoplasma.girerdii", colnames(mypropdata))] <- "#7f7f7f"
    mycolors[grep("Prevotella.amnii", colnames(mypropdata))] <- "#c7c7c7"
    mycolors[grep("Megasphaera.OTU70.type1", colnames(mypropdata))] <- "#f7b6d2"
    mycolors[grep("Sneathia.amnii", colnames(mypropdata))] <- "#9467bd"
    mycolors[grep("Lactobacillus.gasseri_cluster", colnames(mypropdata))] <- "#c5b0d5"
    mycolors[grep("Lactobacillus.jensenii", colnames(mypropdata))] <- "#ffbb78"
    #mycolors[grep("Lactobacillus.jensenii/fornicalis/psittaci", colnames(mypropdata))] <- "#ffbb78"
    mycolors[grep("Prevotella.cluster2", colnames(mypropdata))] <- "#1f77b4"
    mycolors[grep("Atopobium.vaginae", colnames(mypropdata))] <- "#c49c94"
    mycolors[grep("Gardnerella.vaginalis", colnames(mypropdata))] <- "#d62728"
    mycolors[grep("Lachnospiraceae.BVAB1", colnames(mypropdata))] <- "#ff7f0e"
    #mycolors[grep("Lactobacillus.crispatus.cluster", colnames(mypropdata))] <- "#fff5aa"
    mycolors[grep("Lactobacillus.crispatus", colnames(mypropdata))] <- "#fff5aa"
    mycolors[grep("Lactobacillus.iners", colnames(mypropdata))] <- "#aec7e8"
    mycolors[grep("Prevotella.bivia", colnames(mypropdata))] <- "#17becf"
    mycolors[grep("Sneathia.sanguinegens", colnames(mypropdata))] <- "#ff9896"
    mycolors[grep("Parvimonas.OTU142", colnames(mypropdata))] <- "#bcbd22"
    mycolors[grep("Ureaplasma.cluster23", colnames(mypropdata))] <- "#2ca02c"
    mycolors[grep("Dialister.micraerophilus", colnames(mypropdata))] <- "#8c564b"
    mycolors[grep("Other", colnames(mypropdata))] <- "#faf0e6" 
  } else if (colScheme == "classic") {
    if (plottype=="ggplot2") {
      #mycolors <- c(rainbow(ncol(mypropdata), start=0.2, end=0.90))
      mycolors <- rep("gray", ncol(mypropdata))
    } else {
      mycolors <- c(rainbow(ncol(mypropdata)+1, start=0.2, end=0.90))
    }
    # for genus-level
    mycolors[grep("Lactobacillus",   colnames(mypropdata))]      <- "yellow"
    mycolors[grep("Sneathia",   colnames(mypropdata))]           <- "purple"
    mycolors[grep("Gardnerella",     colnames(mypropdata))]      <- "red"
    mycolors[grep("Lachnospiraceae", colnames(mypropdata))]      <- "orange"
    mycolors[grep("BVAB", colnames(mypropdata))]                 <- "orange"
    mycolors[grep("Prevotella",        colnames(mypropdata))]    <- "blue"
    mycolors[grep("Atopobium",        colnames(mypropdata))]     <- "brown"
    
    
    mycolors[grep("crispatus",   colnames(mypropdata))]           <- "yellow"
    mycolors[grep("iners",   colnames(mypropdata))]               <- "lightblue"
    mycolors[grep("gasseri",   colnames(mypropdata))]             <- "pink"
    mycolors[grep("delbrueckii",   colnames(mypropdata))]         <- "black"
    mycolors[grep("jensenii",   colnames(mypropdata))]            <- "green"
    mycolors[grep("Sneathia.amnii",        colnames(mypropdata))] <- "purple"
    mycolors[grep("Streptococcus",     colnames(mypropdata))]     <- "orange"
  } else {
    mycolors <- c(rainbow(ncol(mypropdata)+1, start=0.2, end=0.90))
    mycolors <- rep("gray", ncol(mypropdata)+1)
    # for genus-level
    mycolors[grep("Lactobacillus",   colnames(mypropdata))]      <- "yellow"
    mycolors[grep("Sneathia",   colnames(mypropdata))]           <- "purple"
    mycolors[grep("Gardnerella",     colnames(mypropdata))]      <- "red"
    mycolors[grep("Lachnospiraceae", colnames(mypropdata))]      <- "orange"
    mycolors[grep("BVAB", colnames(mypropdata))]                 <- "orange"
    mycolors[grep("Prevotella",        colnames(mypropdata))]    <- "blue"
    mycolors[grep("Atopobium",        colnames(mypropdata))]     <- "brown"
    
    
    mycolors[grep("crispatus",   colnames(mypropdata))]           <- "yellow"
    mycolors[grep("iners",   colnames(mypropdata))]               <- "lightblue"
    mycolors[grep("gasseri",   colnames(mypropdata))]             <- "pink"
    mycolors[grep("delbrueckii",   colnames(mypropdata))]         <- "black"
    mycolors[grep("jensenii",   colnames(mypropdata))]            <- "green"
    mycolors[grep("Sneathia.amnii",        colnames(mypropdata))] <- "purple"
    mycolors[grep("Streptococcus",     colnames(mypropdata))]     <- "orange"
  }


  mylabels <- rep("", nrow(mypropdata))
  mylabels[!duplicated(barplottypes)] <- barplottypes[!duplicated(barplottypes)]

  if (length(dxbv) > 0) {
    mylabels[dxbv == "Yes"] <- paste(mylabels[dxbv=="Yes"],"---",sep=" ")
  }

  if (length(outpath) > 0) {
    png(paste(outpath, "barplot.png", sep=""), width=2560, height=1440) 
  }

  if (plottype == "ggplot2") {
    library(ggplot2)
    mycolors <- factor(mycolors)
    mycolordf <- data.frame(mycolors, Taxa=colnames(mypropdata))
    sampleorder <- factor(rownames(mypropdata))
    
    ggplotdata <- melt(as.matrix(mypropdata),  
                   id.vars=names(mypropdata), 
                   varnames=c("SampleID", "Taxa"), 
                   value.name="ATprop")
    ggplotdata$SampleID <- factor(ggplotdata$SampleID, levels=sampleorder)
    ggplotdata <- merge(ggplotdata, mycolordf, by="Taxa")
    ggplotdata <- ggplotdata[ggplotdata$ATprop != 0,]

    p <- ggplot(ggplotdata, aes(SampleID, ATprop, fill=mycolors, group=ATprop)) + 
      geom_bar(stat="identity", position="stack", width=1) +
      scale_fill_manual(values = levels(mycolors)) +
      labs(x="Sample", y="Relative Abundance") + 
      theme(axis.text.x=element_blank(), 
	    legend.position="none", 
	    axis.ticks=element_blank()) + 
      ggtitle("All Vaginal Samples")
    print(p)
  } else {
    barplot(as.matrix(t(mypropdata)), ylim=c(0,100), col=mycolors, 
            las=2, names.arg = mylabels, axisnames=TRUE, 
            xlab="samples", ylab="Percentage of Reads (%)", 
            cex.names=0.5, border="NA", space=0)
  }
  if (length(outpath) > 0) {
    dev.off()
  }
 

  if (addlegend==TRUE) {
    print("Legend")
    mynrow <- ceiling(ncol(mypropdata)/5)

    x <- (1:4)-1
    y <- (1:mynrow)-1
    x <- rep(x, each=mynrow)
    y <- rep(y,4)
    
    if (length(outpath) > 0) {
      png(paste(outpath, "barplotKey.png", sep=""), width=2560, height=1440) 
    }
    par(cex=0.5)
    plot(0,0,type="n", xlab="", ylab="", xaxs="i", yaxs="i", 
             xlim=c(0,max(x)+2), ylim=c(max(y)+0.5,-0.5), xaxt="n", yaxt="n")
    rect(x+0.2/4, y-0.4, x+0.2, y+0.4, border="black", col=mycolors[1:(4*mynrow)])
    text(x+0.2*1.2, y, substr(names(mypropdata[1:(4*mynrow)]),0,25), cex=0.5, 
           adj=c(0,0.5))

    lastnrow <- ncol(mypropdata) - mynrow*4
    x <- rep(4, mynrow)
    y <- 1:lastnrow-1 
    rect(x+0.2/4, y-0.4, x+0.2, y+0.4, border="black",  
           col=mycolors[(4*mynrow+1):ncol(mypropdata)])
    text(x+0.2*1.2, y, substr(names(mypropdata)[(4*mynrow+1):ncol(mypropdata)],0,25),
         cex=0.5, adj=c(0,0.5))
    if (length(outpath) > 0) {
      dev.off()
    }
  }
  if (length(alphaDiversity) > 0) {
    alphaDiversity <- alphaDiversity[order(myTypeOrder[mytypes,],-maxprop)]
    if (length(outpath) > 0) {
      png(paste(outpath, "alphabar.png", sep=""), width=2560, height=1440) 
    }
    barplot(alphaDiversity, names.arg=rep("",length(alphaDiversity)))
    if (length(outpath) > 0) {
      dev.off()
    }
  }
}

