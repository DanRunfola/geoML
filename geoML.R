#
#
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
library(rpart.plot)
library(Rcpp)
library(stargazer)
library(matrixStats)
#library(doBy)

if (Sys.getenv("USER") == "vagrant" | Sys.getenv("USER") == "ubuntu") {
  sourceCpp(path.expand("~/geoML/split.cpp"))
} else {
  sourceCpp("/home/aiddata/Desktop/Github/geoML/split.cpp")
}



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

ctev <- function(y, wt,parms) {
  out = node_evaluate(y)
  list(label= out[1], deviance=out[2])
}

ctsplit <- function(y, wt, x, parms, continuous) {
  if (continuous) {
    n = nrow(y)
    res = splitc(y)
    list(goodness=res[1:(n-1)], direction=res[n:(2*(n-1))])
  }
  else{
    res = splitnc(y,x)
    n=(length(res)+1)/2
    list(goodness=res[1:(n-1)], direction=res[n:(2*n-1)])
  }
}


ctinit <- function(y, offset, parms, wt) {
  sfun <- function(yval, dev, wt, ylevel, digits ) {
    print(yval)
    paste("events=", round(yval[,1]),
          ", coef= ", format(signif(yval[,2], digits)),
          ", deviance=" , format(signif(dev, digits)),
          sep = '')}
  environment(sfun) <- .GlobalEnv
  list(y =y, parms = 0, numresp = 1, numy = 4,
       summary = sfun)
}

#============================================================
#============================================================
#Main Function
geoML <- function(dta,
                  trt,
                  ctrl,
                  outcome,
                  out_path,
                  file.prefix,
                  kvar,
                  top.rep,
                  geog.fields = c("latitude", "longitude"),
                  caliper=0.5,
                  counterfactual.name="Control",
                  tree.ctrl = c(5,10),
                  tree.cnt = 1000,
                  col.invert=FALSE,
		  tree.cex=1.0)
{

  #Ensure no names have spaces
  names(dta) <- make.names(names(dta))

  #Truncate to relevant variables
  cvar.lim <- ceiling(length(ctrl) / 2)
  ctrl.vars <- ctrl[1:cvar.lim]
  ctrl.names <- ctrl[(ceiling(length(ctrl) / 2)+1):length(ctrl)]

  sub.dta <- dta[c(ctrl.vars, trt[1], outcome[1])]

  ret.vec.dta <- dta[top.rep]

  lab.ids <- sapply(sub.dta, is.numeric)

  sub.dta.desc <- sub.dta[ , lab.ids]

  for(k in 1:length(names(sub.dta.desc)))
  {

    id <- ctrl.names[match(names(sub.dta.desc)[k],ctrl.vars)]
    if(names(sub.dta.desc)[k] == trt[1])
    {
      id = trt[2]
    }
    if(names(sub.dta.desc)[k] == outcome[1])
    {
      id = outcome[2]
    }

    if(k == 1)
    {
      if(!is.na(id))
      {
        labels <- id
      }
      else
      {
        labels <- names(sub.dta.desc)[k]
      }
    }
    else
    {
      if(!is.na(id))
      {
        labels <- c(labels, id)
      }
      else
      {
        labels <- c(labels,names(sub.dta.desc)[k])
      }
    }
  }




  #============================================================
  #============================================================
  #Stargazer table of descriptive statistics (prefix_desc.html)
  desc.out = paste(out_path,file.prefix,"_desc.html",sep="")
  len_o <- 1:length(names(sub.dta.desc))
  stargazer(sub.dta.desc[sub.dta.desc[trt[1]] == 1,],
            type="html",
            median=TRUE,
            digits=2,
            title=paste("Descriptive Statistics for ",trt[2]," (Treated)", sep=""),
            covariate.labels = labels,
            omit.summary.stat=c("n"),
            order=len_o,
            out = paste(out_path,file.prefix,"_desc.html",sep="")
  )



  #============================================================
  #============================================================
  #Map of all locations (prefix_map.png)
  map.all.out <- paste(out_path,file.prefix,"_map.png",sep="")
  lonlat <- sub.dta[,c(geog.fields[2], geog.fields[1])]
  spdf <- SpatialPointsDataFrame(coords = lonlat, data = sub.dta,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))


  if (Sys.getenv("USER") == "vagrant" | Sys.getenv("USER") == "ubuntu") {
    land.mask <- readOGR(path.expand("~/geoML/data/countries.geojson"), "OGRGeoJSON")
  } else {
    land.mask <- readOGR("data/countries.geojson", "OGRGeoJSON")
  }


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
  sub.dta$propensity <- predict(m.ret$model, newdata=sub.dta, type="response")

  #============================================================
  #============================================================
  #Pre-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_prebalance.png)
  pre.balance.path <- paste(out_path,file.prefix,"_prebalance.png",sep="")
  plot.dta <- sub.dta
  plot.dta$Type <- counterfactual.name
  plot.dta[plot.dta[trt[1]] == 1,]["Type"] <- trt[2]

  #plot.dta <- lapply(plot.dta, function(x) as.numeric(
  #			as.character(x)))

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
  prop.model.out <- paste(out_path,file.prefix,"_propensityModel.html",sep="")

  for(k in 2:length(m.ret$model$coefficients))
  {

        id <- ctrl.names[match(names(m.ret$model$coefficients[k]),ctrl.vars)]

  if(k == 2)
      {
        if(!is.na(id))
          {
            prop.labels <- id
          }
        else
          {
            prop.labels <- names(m.ret$model$coefficients[k])
          }
      }
  else
      {
        if(!is.na(id))
        {
            prop.labels <- c(prop.labels, id)
        }
          else
        {
          prop.labels <- c(prop.labels, names(m.ret$model$coefficients[k]))
        }
      }
  }


  stargazer(m.ret$model,
            type="html",
            median=TRUE,
            digits=2,
            title=paste("Propensity Model: ",trt[2]," (Treated), ",counterfactual.name, " (Control)", sep=""),
            covariate.labels = prop.labels,
            font.size="small",
            no.space=TRUE,
            omit.stat=c("LL","ser","f"),
            single.row=TRUE,
            ci=TRUE, ci.level=0.95,
            out = paste(out_path,file.prefix,"_propensityModel.html",sep="")
  )


  #============================================================
  #============================================================
  #Post-balance plot - propensity score and four key variables (first four in list taken otherwise; prefix_postbalance.png)
  post.balance.outpath <- paste(out_path,file.prefix,"_postbalance.png",sep="")
  plot.dta <- match.data(m.ret)
  plot.dta$Type <- counterfactual.name
  plot.dta[plot.dta[trt[1]] == 1,]["Type"] <- trt[2]

   a0 <- ggplot(data=plot.dta, aes(x=distance,fill=Type)) +
    geom_histogram(binwidth=.05, alpha=.5, position="identity")+
    ggtitle("Post-Balance Metrics") +
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

  png(paste(out_path,file.prefix,"_postbalance.png",sep=""),
      width = 7,
      height = 6,
      units = 'in',
      res = 300)
  multiplot(a0,a1, a2, a3, a4, layout=matrix(c(1,1,2,3,4,5), nrow=3, byrow=TRUE))
  dev.off()



  #============================================================
  #============================================================
  #Stargazer table of balance statistics (prefix_balancestats.html)
  balance.stats.out <- paste(out_path,file.prefix,"_balancestats.html",sep="")
  s.rm <- summary(m.ret)$reduction
  row.names(s.rm)[1] <- "Propensity Score"
  for(k in 2:length(row.names(s.rm)))
  {

    id <- ctrl.names[match(row.names(s.rm)[k],ctrl.vars)]
      if(!is.na(id))
      {
        row.names(s.rm)[k] <- id
      }
      else
      {
        row.names(s.rm)[k] <- row.names(s.rm)[k]
      }
  }


  stargazer(s.rm,
            type="html",
            title="Percent Balance Improvement:",
            dep.var.labels.include = TRUE,
            summary=FALSE,
            out = paste(out_path,file.prefix,"_balancestats.html",sep="")
  )


  #============================================================
  #============================================================
  #Causal Tree results (prefix_ct.png)
  #Edit 0 and 1 propensity cases
  sub.dta$ML.prop <- sub.dta$propensity

  for(i in 1:length(sub.dta['ML.prop']))
  {

    if(sub.dta['ML.prop'][i] >= 0.99)
    {
      sub.dta['ML.prop'][i] = 0.99
    }
    if(sub.dta['ML.prop'][i] <= 0.01)
    {
      sub.dta['ML.prop'][i] = .01
    }
  }

  #Tree propensity calculations
  transOutcome <- list(rep(0,nrow(sub.dta)))
  trans.rf.prop <- list(rep(0,nrow(sub.dta)))

  for(i in 1:nrow(sub.dta))
  {
    if(sub.dta[trt[1]][[1]][i] == 1)
    {
      #Treated
      transOutcome[i] = sub.dta[outcome[1]][[1]][i] *
        (1 / sub.dta$ML.prop[i])
      trans.rf.prop[i] = sub.dta["ML.prop"][[1]][i]
    }
    else
    {
      #Untreated
      transOutcome[i] = -1 * (sub.dta[outcome[1]][[1]][i] *
                                ((1-0) / (1 - sub.dta$ML.prop[i])))
      trans.rf.prop[i] = sub.dta["ML.prop"][[1]][i] * -1
    }
  }
  sub.dta$transOutcome <- unlist(transOutcome)
  sub.dta$transProp <- unlist(trans.rf.prop)

  #------------------
  #------------------
  #CT
  #------------------
  #------------------
  tree.dta <- sub.dta[!is.na(sub.dta$transOutcome),]
  match.ids <- m.ret$match.matrix[!is.na(m.ret$match.matrix)]
  match.ids <- c(match.ids, names(m.ret$match.matrix[!is.na(m.ret$match.matrix),]))

  tree.dta <- tree.dta[match.ids,]

  alist <- list(eval=ctev, split=ctsplit, init=ctinit)
  dbb = tree.dta
  k = tree.ctrl[2]
  n = dim(dbb)[1]
  crxvdata = dbb
  crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
  list = 1:k
  m.split = tree.ctrl[1]

  errset = list()

  for (i in 1:k){
    errset[[i]] = list()
    trainingset <- subset(crxvdata, id %in% list[-i])

    tree.exec <- paste("sub.fit = rpart(cbind(",
                       outcome[1],",",
                       trt[1],",ML.prop,transOutcome)~",
                       paste(ctrl.vars, collapse="+"),
                       ",trainingset, control=rpart.control(cp = 0,minsplit = m.split),method=alist)", sep="")
    sub.fit <- eval(parse(text=tree.exec))
    sub.fit.dm = data.matrix(sub.fit$frame)
    index = as.numeric(rownames(sub.fit$frame))
    removed_nodes = 0
    removed_nodes = cross_validate(sub.fit.dm, index,removed_nodes)
    removed_nodes = removed_nodes[-1]

    for(l in 1:length(removed_nodes)){
      error = 0
      sub.fit.pred = snip.rpart(sub.fit, removed_nodes[1:l])

      #Subset Fit
      testset <- subset(crxvdata, id %in% c(i))
      pt = predict(sub.fit.pred,testset,type = "matrix")
      y = data.frame(pt)
      val = data.matrix(y)
      idx = as.numeric(rownames(testset))
      dbidx = as.numeric(rownames(dbb))


      for(pid in 1:(dim(y)[1])){
        id = match(idx[pid],dbidx)
        error = error + (dbb$transOutcome[id] - val[pid])^2
      }

	if(length(pt) == 0){error = 0}

      if(error == 0){
        errset[[i]][l] = 1000000
      }
      else{
        errset[[i]][l] = error/k
      }
    }
  }

  #Identify the average error to depth ratio across all cross-validations
  avg.index <- vector()
  for(e in 1:length(errset))
  {
    avg.index[e] <- which.min(errset[[e]])
  }

  #---------------
  #Build Final Tree
  #---------------
  ct.png.out <- paste(out_path,file.prefix,"_ct.png",sep="")
  final.tree.exec <- paste("sub.fit = rpart(cbind(",
                     outcome[1],",",
                     trt[1],",ML.prop,transOutcome)~",
                     paste(ctrl.vars, collapse="+"),
                     ",crxvdata, control=rpart.control(cp = 0,minsplit = m.split),method=alist)", sep="")

  fit1 <- sub.fit <- eval(parse(text=final.tree.exec))

  fit = data.matrix(fit1$frame)
  index = as.numeric(rownames(fit1$frame))


  removed_nodes = 0
  removed_nodes = cross_validate(fit, index,removed_nodes)
  removed_nodes = removed_nodes[-1]
  pruned_nodes = removed_nodes[1:round(mean(avg.index))]
  final.tree <- snip.rpart(fit1, pruned_nodes)

  print.tree <- final.tree
  var.rec <- ""
  het.lab <- ""
  for(i in 1:length(levels(print.tree$frame$var)))
  {
    if(levels(print.tree$frame$var)[i] != "<leaf>")
    {
      var.rec <- c(var.rec, ctrl.vars[match(levels(print.tree$frame$var)[i],ctrl.vars)])
      het.lab <- c(het.lab, paste(ctrl.names[match(levels(print.tree$frame$var)[i],ctrl.vars)],"*Treatment"))
      levels(print.tree$frame$var)[i] <- ctrl.names[match(levels(print.tree$frame$var)[i],ctrl.vars)]
    }
  }

  png(paste(out_path,file.prefix,"_ct.png",sep=""),
      width = 7,
      height = 4,
      units = 'in',
      res = 600)
  if(col.invert == TRUE)
  {
    rpart.plot(print.tree, extra=1, branch=1, type=4, tweak=1, cex=tree.cex, clip.right.labs=FALSE,
               box.col=c("palegreen3", "pink")[findInterval(print.tree $frame$yval, v = c(-999999999999999999,0))],
               faclen=0,
               varlen=0,fallen.leaves=FALSE)
  }
  else
  {
    rpart.plot(print.tree, extra=1, branch=1, type=4, tweak=1, cex=tree.cex, clip.right.labs=FALSE,
               box.col=c("pink", "palegreen3")[findInterval(print.tree $frame$yval, v = c(-999999989999999999,0))],
               faclen=0,
               varlen=0,fallen.leaves=FALSE)
  }

  #title(paste("Causal Tree: ",trt[2]," (Treated), ",counterfactual.name, " (Control)", sep=""), cex.main=0.75)
  dev.off()

  #============================================================
  #============================================================
  #Stargazer of linear model using matched cases and interactions identified in Causal Tree (prefix_linearMatch.html)
  linear.het.out <- paste(out_path,file.prefix,"_linearMatch.html",sep="")
  lm.exec <- paste("lm(",outcome[1],"~",trt[1],"+",paste(ctrl.vars, collapse="+"), sep="")

  print(match.data(m.ret))
  het.dta <- lapply(match.data(m.ret), function(x) scale(as.numeric(x)))
  print(het.dta)
  for(i in 2:length(var.rec))
  {
    lm.exec <- paste(lm.exec, "+", var.rec[i],"*",trt[1])
  }


  lm.exec <- paste(lm.exec, ",data=het.dta)")


  het.model <- eval(parse(text=lm.exec))
  v.rec.add <- paste("treatment:", var.rec, sep="")
  for(k in 2:length(het.model$coefficients))
  {

    id <- ctrl.names[match(names(het.model$coefficients[k]),ctrl.vars)]

    if(is.na(id))
    {
     id <- het.lab[match(names(het.model$coefficients[k]),v.rec.add)]
    }

    if(k == 2)
    {
      if(!is.na(id))
      {
        all.lab <- id
      }
      else
      {
        all.lab <- names(het.model$coefficients[k])
      }
    }
    else
    {
      if(!is.na(id))
      {
        all.lab <- c(all.lab, id)
      }
      else
      {
        all.lab <- c(all.lab, names(het.model$coefficients[k]))
      }
    }
  }

  print(all.lab)
  #all.lab <- c("Treatment", ctrl.names, het.lab[2:length(het.lab)])
  print(het.model)
  print(lm.exec)
  stargazer(het.model,
            type="html",
            median=TRUE,
            digits=2,
	    digits.extra = 100,
            title=paste("Matched Model: ",trt[2]," (Treated), ",counterfactual.name, " (Control)", sep=""),
            covariate.labels = all.lab,
            font.size="small",
            no.space=TRUE,
            omit.stat=c("LL","ser","f"),
            single.row=TRUE,
            ci=TRUE, ci.level=0.95,
            dep.var.labels = outcome[2],
            out = paste(out_path,file.prefix,"_linearMatch.html",sep="")
  )
  #============================================================
  #============================================================
  #Map of predicted results from Causal Tree (prefix_map_estimate.png)
  sub.dta$tree.pred <- predict(final.tree, newdata=sub.dta)

  out.sum = paste(out_path,file.prefix,"_summary_estimates.txt",sep="")

  print("SUMMARY STATISTICS OF TREE PREDICTION:")
  print(summary(sub.dta$tree.pred))

  write("Treatment Cases:", file=out.sum, append=TRUE)
  write(paste(dim(sub.dta.desc[sub.dta.desc[trt[1]] == 1,])[1]), file = out.sum,
        append = TRUE)
  write("Control Cases:", file=out.sum, append=TRUE)
  write(paste(dim(sub.dta.desc[sub.dta.desc[trt[1]] == 0,])[1]), file = out.sum,
        append = TRUE)

  ct.mean.est <- summary(sub.dta$tree.pred)[4]

  sub.dta <- cbind(sub.dta, ret.vec.dta)



  trt.dta <- sub.dta[sub.dta[trt[1]] == 1,]


  lonlat <- trt.dta[c(geog.fields[2], geog.fields[1])]

  spdf <- SpatialPointsDataFrame(coords = lonlat, data = trt.dta,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

  if(col.invert == TRUE)
  {
  p <- spplot(spdf, zcol="tree.pred", cex=0.4,
              cuts=5,
              col.regions=terrain.colors(5, alpha=1),
              xlim=c(-180,180),
              ylim=c(-60,90),
              main=list(label=paste("Impact of Treatment (",trt[2],")", sep="")
                        ,cex=0.5),
              sp.layout = list(list(land.mask, fill="grey", first=TRUE)))
  }
  else
  {
    p <- spplot(spdf, zcol="tree.pred", cex=0.4,
                cuts=5,
                col.regions=rev(terrain.colors(5, alpha=1)),
                xlim=c(-180,180),
                ylim=c(-60,90),
                main=list(label=paste("Impact of Treatment (",trt[2],")", sep="")
                          ,cex=0.5),
                sp.layout = list(list(land.mask, fill="grey", first=TRUE)))
  }
  map.est.out <- paste(out_path,file.prefix,"_map_estimate.png",sep="")
  png(paste(out_path,file.prefix,"_map_estimate.png",sep=""),
      width = 6,
      height = 4,
      units = 'in',
      res = 300)
  print(p)
  dev.off()


  #============================================================
  #============================================================
  #CSV of predicted results from Causal Tree and relevant covariates (prefix_prediction.csv)
  #write trt only
  # write.csv(trt.dta, paste(out_path,file.prefix,"_prediction.csv",sep=""))
  # write trt and cnt
  write.csv(sub.dta, paste(out_path,file.prefix,"_prediction.csv",sep=""))

  #Biggset and smallest
  ret.trt.dta <- trt.dta[trt.dta$tree.pred == max(trt.dta$tree.pred) |
                           trt.dta$tree.pred == min(trt.dta$tree.pred),]

  write.csv(ret.trt.dta, paste(out_path, file.prefix, "_BW.csv", sep=""))
  best.worst.out <- paste(out_path, file.prefix, "_BW.csv", sep="")

  #============================================================
  #============================================================
  #Random Forest results (prefix_rf.csv)
  if (Sys.getenv("USER") == "vagrant" | Sys.getenv("USER") == "ubuntu") {
    python.path <- path.expand("~/CausalForest/CF.py")
  } else {
    python.path <- "/home/aiddata/Desktop/Github/CausalForest/CF.py"
  }


  csv.str <- tempfile(fileext=".csv")
  nums <- sapply(tree.dta, is.numeric)

  rf.tree.dta <- tree.dta[ , nums]
  rf.tree.dta <- rf.tree.dta[complete.cases(rf.tree.dta),]
  print(summary(rf.tree.dta))
  write.csv(rf.tree.dta, csv.str)

  c.vars.pass = paste(ctrl.vars[ctrl.vars %in% names(rf.tree.dta)], collapse=",")
  sys.call.str <- paste("python", python.path, csv.str, c.vars.pass, outcome[1], "transProp",
                        paste(out_path,file.prefix,"_rf.csv",sep=""), tree.cnt)

  tree.resp <- system(sys.call.str, intern=TRUE)
  #CSV is written in the python script in the above line.

  #============================================================
  #============================================================
  #Figure of most important variables in RF (prefix_purity.png)
  purity <- tail(tree.resp, n=1)
  p.A <- gsub("\\[|\\]","", purity)
  p.B <- gsub("\\(", "", p.A)
  p.C <- gsub("'","", p.B)
  p.D <- gsub("\\)", "", p.C)
  pur.out <- strsplit(gsub(" ","",p.D), ",")
  df.len <- length(pur.out[[1]]) / 2
  purity.df <- data.frame(var=1:df.len, purity=1:df.len, col=1:df.len)


  for(i in 1:df.len)
  {

    purity.df[1][i,] <- ctrl.names[match(as.character(pur.out[[1]][(2 * i - 1)]), ctrl.vars)]
    purity.df[2][i,] <- as.character(pur.out[[1]][2 * i])

    #Calculate vector of bar colors based on the initial CT.
    het.id <- match(paste("treatment:",as.character(pur.out[[1]][(2 * i - 1)]),sep=""),
          names(het.model$coefficients))

    if(!is.na(het.id))
    {
      if(col.invert == FALSE)
      {
        if(het.model$coefficients[het.id][[1]] > 0)
        {
          purity.df[3][i,] <- "green"
        }
        if(het.model$coefficients[het.id][[1]] < 0)
        {
          purity.df[3][i,] <- "red"
        }
      }
      else
      {
        if(het.model$coefficients[het.id][[1]] > 0)
        {
          purity.df[3][i,] <- "red"
        }
        if(het.model$coefficients[het.id][[1]] < 0)
        {
          purity.df[3][i,] <- "green"
        }
      }
      if(het.model$coefficients[het.id][[1]] == 0)
      {
        purity.df[3][i,] <- "black"
      }
    }
    else
    {
      purity.df[3][i,] <- colors()[309]
    }

  }



  ##==========================
  ##Purity calculations are temporarily disabled
  #Jianing is looking into SciKit errors.
  # png(paste(out_path,file.prefix,"_purity.png",sep=""),
  #     width = 6,
  #     height = 4,
  #     units = 'in',
  #     res = 300)
  # par(mai=c(1,2,1,1))
  # barplot(height=as.numeric(purity.df["purity"][[1]]),
  #         names.arg=purity.df["var"][[1]], horiz=TRUE,
  #         cex.names=0.5, las=2, xlab="Tree Purity Contribution",
  #         cex.axis = 0.6,
  #         xlim=c(0,max(as.numeric(purity.df["purity"][[1]]))),
  #         col=purity.df["col"][[1]])
  # dev.off()

  #Map of uncertainty ("% of observations within 1 standard deviation of the mean") from RF (prefix_rf_unc.png)
  rf.res <- read.csv(paste(out_path,file.prefix,"_rf.csv",sep=""), header=FALSE)
  rf.tree.dta$unc <- colSds(as.matrix(rf.res)) * 1.96
  rf.tree.dtaB <- rf.tree.dta[rf.tree.dta[trt[1]] == 1,]
  lonlat <- rf.tree.dtaB[,c(geog.fields[2], geog.fields[1])]
  spdf <- SpatialPointsDataFrame(coords = lonlat, data = rf.tree.dtaB,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

  #write("\nUncertainty Mean (+/- 95%):",file=out.sum, append=TRUE)
  #write(summary(rf.tree.dta$unc)[4], file=out.sum, append=TRUE)

  p <- spplot(spdf, "unc", cex=0.4,
              cuts=5,
              col.regions=rev(heat.colors(6, alpha=1)),
              xlim=c(-180,180),
              ylim=c(-60,90),
              main=list(label="Uncertainty in Estimates (+/- @ 95% Confidence Interval)"
                        ,cex=0.5),
              sp.layout = list(list(land.mask, fill="grey", first=TRUE)))

  uncertainty.map.out <- paste(out_path,file.prefix,"_map_uncertainty.png",sep="")
  png(paste(out_path,file.prefix,"_map_uncertainty.png",sep=""),
      width = 6,
      height = 4,
      units = 'in',
      res = 300)
  print(p)
  dev.off()

##Produce histogram of uncertainty, with a red line marking the CT result.
  hist.out <- paste(out_path,file.prefix,"_uncertHist.png", sep="")
  rf.means <- rowMeans(rf.res)
  mean.rf.est <- mean(rf.means)
L <-  hist(rf.means, main="Model Uncertainty",
       xlab=paste("Estimated Mean Impact",outcome[2]),
       ylab="Number of Simulations"
       )

if(col.invert == FALSE)
{
clr <- ifelse(L$breaks < 0, "#F8766D", "#00BA38")[-length(L$breaks)]
}
else
{
clr <- ifelse(L$breaks < 0, "#00BA38", "#F8766D")[-length(L$breaks)]
}


  png(paste(out_path,file.prefix,"_uncertHist.png", sep=""),
	width = 6,
	height = 4,
	units='in',
	res=300)

  hist(rf.means, main="Model Uncertainty",
       xlab=paste("Estimated Mean Impact",outcome[2]),
       ylab="Number of Simulations",
       col=clr)
  abline(v=mean.rf.est, col="blue", lwd=2)
  #abline(v=ct.mean.est, col="red", lwd=1)
  legend("topleft", c("Random Forest Mean"),
fill=c("blue"), bty="n", cex=0.5)
dev.off()


  # ----------------------------------

  # out.sum              # txt  - Text file containing N for Control and Treatment cases
  # hist.out             # png  - Histogram of estimate / overall summary
  # best.worst.out       # csv  - CSV of best and worst projects
  # desc.out             # html - Descriptive Stats
  # map.all.out          # png  - Map of all locations
  # prop.model.out       # html - Propensity model
  # pre.balance.path     # png  - Pre balance figure
  # post.balance.outpath # png  - Post-balance figure
  # balance.stats.out    # html - Balance summary table
  # ct.png.out           # png  - Causal Tree Figure
  # linear.het.out       # html - Linear model with het effects table
  # map.est.out          # png  - Map of Estimates
  # uncertainty.map.out  # png  - Map of uncertainty


  htmler <- setRefClass(
    "htmler",
    fields = list(
      html = "character"
    ),
    methods = list(
      initialize = function() {
        "Init html"
        html <<- "<!DOCTYPE html>
                   <html>
                   <head>
                   <meta charset='UTF-8'>
                   <title>geoML Output</title>
                   </head>
                   <body>"
      },
      complete = function() {
        "Complete html"
        html <<- paste(html, "</body></html>", sep="")
        return(html)
      },
      view = function() {
        print(html)
      },
      add = function(text, newLine=TRUE) {
        n <- ''
        if (newLine == TRUE) {
          n <- '\n'
        }

        html <<- paste(html, text, n, sep="")
      }
    )
  )

  html <- htmler$new()
  html$add(paste("<h2> Contrast: ",trt[2]," (Treated), ",counterfactual.name, " (Control)</h2>", sep=""))
  for (i in readLines(out.sum)) {
    html$add(paste('<b>',i,'</b><br \\>'))
  }
  html$add('<br \\>')

  html$add(paste('<img src="', basename(hist.out), '" width="600px">', sep=""))

  #Description of histogram
  html$add("<p>")
  html$add("Figure 1. Mean Estimated Impacts from Causal Tree and Random Forest Estimation Strategies.")
  html$add("</p>")
  html$add("<p>")
  #html$add("Figure 1 illustrates the overall mean impact estimate from the Random Forest.")
  html$add("In this illustration, the blue vertical line is the global estimate of impact from the random forest.")
  html$add("The range of uncertainty is illustrated using a histogram, indicating the proportion of simulations which result in positive or negative estimates.")
  html$add("Specific confidence intervals can be calculated for this global estimate, or for individual project locations.  A map of project location uncertainties is presented below.")
  html$add("</p>")

  for (i in readLines(desc.out)) {
    html$add(i, newLine=FALSE)
  }
  html$add('\n')


  html$add(paste('<img src="', basename(map.all.out), '" width="600px">', sep=""))

  for (i in readLines(prop.model.out)) {
    html$add(i, newLine=FALSE)
  }

  html$add(paste('<img src="', basename(pre.balance.path), '" width="600px">', sep=""))
  html$add(paste('<img src="', basename(post.balance.outpath), '" width="600px">', sep=""))

  for (i in readLines(balance.stats.out)) {
    html$add(i, newLine=FALSE)
  }
  html$add('\n')

  html$add(paste('<img src="', basename(ct.png.out), '" width="1200px">', sep=""))

  for (i in readLines(linear.het.out)) {
    html$add(i, newLine=FALSE)
  }
  html$add("<p>Standardized Beta Coefficients for Linear Model with Heterogeneous Effects</p>")
  html$add('\n')

  html$add(paste('<img src="', basename(map.est.out), '" width="600px">', sep=""))

  html$add(paste('<img src="', basename(uncertainty.map.out), '" width="600px">', sep=""))

  html_out <- html$complete()

  print(html_out)
  write(html_out, file=paste(out_path,file.prefix,"_complete.html",sep=""))

  # ----------------------------------


  return(trt.dta)


  }
