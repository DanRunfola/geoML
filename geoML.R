#Generalized Implementation of Propensity Matching, CT, RF, and MatchIt Linear Models
#Inputs: 
#Required:
#Dataset
#List of Control Variables
#Single Outcome Variable
#Single Treatment Variable
#Folder for output
#Prefix for file outputs
#Optional: 
#four key variables for balance plotting (list)
#geog.fields = alternative columns for latitude and longitude.

#Outputs:
#Stargazer table of descriptive statistics (prefix_desc.html)
#Map of all locations (prefix_map.png)
#Pre-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_prebalance.html)
#Stargazer table for propensity linear model (prefix_propensityModel.html)
#Post-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_postbalance.html)
#Stargazer table of balance statistics (prefix_balancestats.html)
#Causal Tree results (prefix_ct.png)
#Stargazer of linear model using matched cases and interactions identified in Causal Tree (prefix_geoMath.html)
#Map of predicted results from Causal Tree (prefix_prediction.png)
#CSV of predicted results from Causal Tree (prefix_prediction.csv)
#Random Forest results (prefix_rf.csv)
#Figure of most important variables in RF (prefix_purity.png)
#Map of uncertainty ("% of observations within 1 standard deviation of the mean") from RF (prefix_rfUnc.png)

#Libraries
library(sp)
library(ggplot2)
library(rgdal)
#library(spatstat)
#library(maptools)
#library(rgeos)
#library(geojsonio)
library(MatchIt)
library(rpart)
library(Rcpp)
library(stargazer)
#library(doBy)
sourceCpp("/home/aiddata/Desktop/Github/GEF_MFA/Analyses/split.cpp")

#============================================================
#============================================================
#Helper Functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#============================================================
#============================================================
#Main Function
geoML <- function(dta, 
                  trt, 
                  ctrl, 
                  outcome, 
                  pth, 
                  file.prefix, 
                  kvar, 
                  geog.fields = c("latitude", "longitude"), 
                  caliper=0.5, 
                  counterfactual.name="Control")
{
  #Truncate to relevant variables
  cvar.lim <- ceiling(length(ctrl) / 2)
  ctrl.vars <- ctrl[1:cvar.lim]
  ctrl.names <- ctrl[(ceiling(length(ctrl) / 2)+1):length(ctrl)]
  
  sub.dta <- dta[c(ctrl.vars, trt[1], outcome[1])]
  
  labels <- c(ctrl.names, trt[2], outcome[2])
  
  #============================================================
  #============================================================
  #Stargazer table of descriptive statistics (prefix_desc.html)
  stargazer(sub.dta, 
            type="html", 
            median=TRUE, 
            digits=2, 
            title=paste("Descriptive Statistics: ",trt[2]," (Treated), ",counterfactual.name, "(Control)", sep=""),
            covariate.labels = labels,
            out = paste(pth,file.prefix,"_desc.html",sep="")
  )
  
  
  #============================================================
  #============================================================
  #Map of all locations (prefix_map.png)
  lonlat <- sub.dta[,c(geog.fields[2], geog.fields[1])]
  spdf <- SpatialPointsDataFrame(coords = lonlat, data = sub.dta,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
  land.mask <- readOGR("data/countries.geojson", "OGRGeoJSON")
  p <- spplot(spdf, zcol=trt[1], cex=0.25,
              cuts=2,
              at=0.5,
              col.regions=c("black", "whitesmoke"),
              xlim=c(-180,180),
              ylim=c(-60,90),
              #key=list(lines=TRUE),
              pch=c(2,5),
              main=list(label=paste("Treatment (",trt[2],") and Control (",counterfactual.name,") Locations.")
                        ,cex=0.5),
              sp.layout = list(list(land.mask, fill="grey", first=TRUE)))
  p$legend$bottom$args$key$text[[1]][1] <- counterfactual.name
  p$legend$bottom$args$key$text[[1]][2] <- trt[2]
  png(paste(out_path,file.prefix,"_map.png",sep=""),
      width = 6, 
      height = 4, 
      units = 'in', 
      res = 300)
  print(p)
  dev.off()
  
  #============================================================
  #============================================================
  #Run MatchIt

  exec_str = paste("matchit(",trt[1], "~", paste(ctrl.vars, collapse="+"),
                   ",data=na.omit(sub.dta),caliper=",caliper,")",sep="")
  
  print(exec_str)
  
  m.ret <- eval(parse(text=exec_str))
  
  #============================================================
  #============================================================
  #Pre-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_prebalance.png)
  plot.dta <- sub.dta
  plot.dta$propensity <- predict(m.ret$model, newdata=plot.dta, type="response")
  plot.dta$Type <- counterfactual.name
  plot.dta[plot.dta[trt[1]] == 1,]["Type"] <- trt[2]
  
  a0 <- ggplot(data=plot.dta, aes(x=propensity,fill=Type)) + 
    geom_histogram(binwidth=.05, alpha=.5, position="identity")+ 
    ggtitle("Pre-Balance Metrics") + 
    theme(legend.position="bottom") + xlab("Propensity Scores")+
    labs(fill="")
  
  a1 <- ggplot(data=plot.dta, aes(x=plot.dta[kvar[1]],fill=Type)) + 
    geom_histogram(binwidth=NULL, alpha=.5, position="identity")+ 
    ylab("") + theme(legend.position="none") + xlab(ctrl.names[match(kvar[1],ctrl.vars)])
  
  a2 <- ggplot(data=plot.dta, aes(x=plot.dta[kvar[2]],fill=Type)) +  
    geom_histogram(binwidth=NULL, alpha=.5, position="identity")+ 
    ylab("") + theme(legend.position="none") + xlab(ctrl.names[match(kvar[2],ctrl.vars)])
  
  a3 <- ggplot(data=plot.dta, aes(x=plot.dta[kvar[3]],fill=Type))+ 
    geom_histogram(binwidth=NULL, alpha=.5, position="identity")+ 
    ylab("") + theme(legend.position="none") + xlab(ctrl.names[match(kvar[3],ctrl.vars)])
  
  a4 <- ggplot(data=plot.dta, aes(x=plot.dta[kvar[4]],fill=Type)) + 
    geom_histogram(binwidth=NULL, alpha=.5, position="identity") + 
    ylab("") + theme(legend.position="none") + xlab(ctrl.names[match(kvar[4],ctrl.vars)])
  
  png(paste(out_path,file.prefix,"_prebalance.png",sep=""),
      width = 7, 
      height = 6, 
      units = 'in', 
      res = 300)
  multiplot(a0,a1, a2, a3, a4, layout=matrix(c(1,1,2,3,4,5), nrow=3, byrow=TRUE))
  dev.off()
  
  #============================================================
  #============================================================ 
  #Stargazer table for propensity linear model (prefix_propensityModel.html)
  stargazer(sub.dta, 
            type="html", 
            median=TRUE, 
            digits=2, 
            title=paste("Descriptive Statistics: ",trt[2]," (Treated), ",counterfactual.name, "(Control)", sep=""),
            covariate.labels = labels,
            out = paste(pth,file.prefix,"_desc.html",sep="")
  )
  
  #============================================================
  #============================================================
  #Post-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_postbalance.html)
  #============================================================
  #============================================================
  #Stargazer table of balance statistics (prefix_balancestats.html)
  #============================================================
  #============================================================
  #Causal Tree results (prefix_ct.png)
  #============================================================
  #============================================================
  #Stargazer of linear model using matched cases and interactions identified in Causal Tree (prefix_geoMath.html)
  #============================================================
  #============================================================
  #Map of predicted results from Causal Tree (prefix_prediction.png)
  #============================================================
  #============================================================
  #CSV of predicted results from Causal Tree (prefix_prediction.csv)
  #============================================================
  #============================================================
  #Random Forest results (prefix_rf.csv)
  #============================================================
  #============================================================
  #Figure of most important variables in RF (prefix_purity.png)  
  #============================================================
  #============================================================
  #Map of uncertainty ("% of observations within 1 standard deviation of the mean") from RF (prefix_rfUnc.png)
}
