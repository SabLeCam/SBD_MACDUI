library(hzar)
library(doMC)
library(foreach)
library(iterators)
library(parallel)
library(extrafont)
library(coda)

################################################################################# 

### prepare data
setwd("path_to_datafile")
read.csv("../hapgp_sab.txt",header=T, sep="\t") -> macMolecular
macMolecular


##################################################
##### Cline analyses #############################
##################################################
#based on supplementary information R script examples from Derryberry et al 2013

## A typical chain length. This value is the default setting in the package.
chainLength=1e5


## Make each model run off a separate seed
mainSeed=
  list(A=c(as.integer(runif(6, 1, 999))), 
       B=c(as.integer(runif(6, 1, 999))), 
       C=c(as.integer(runif(6, 1, 999))),
       D=c(as.integer(runif(6, 1, 999))))


if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC()
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}

##############################################################################
## Molecular Analysis
## Load example Molecular data from the data table.
macMolecular

## ## Picking an allele for a locus

locnames<-"b1a_f"    #  do that for all haplogroups
hap_data<-macMolecular$b1a_f
hap_nsamples<-macMolecular$f_nSamples


## Blank out space in memory to hold molecular analysis
if(length(apropos("^mac$",ignore.case=FALSE)) == 0 ||
   !is.list(mac) ) mac <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
mac$hap <- list();
## Space to hold the observed data
mac$hap$obs <- list();
## Space to hold the models to fit
mac$hap$models <- list();
## Space to hold the compiled fit requests
mac$hap$fitRs <- list();
## Space to hold the output data chains
mac$hap$runs <- list();
## Space to hold the analysed data
mac$hap$analysis <- list();

## Locus Ada, Allele A from Brumfield et al 2001
mac$hap$obs <-hzar.doMolecularData1DPops(macMolecular$distance, hap_data, hap_nsamples)

## Look at a graph of the observed data
hzar.plot.obsData(mac$hap$obs);


## Make a helper function
mac.loadmodel <- function(scaling,tails,id=paste(scaling,tails,sep=".")){
  mac$hap$models[[id]] <<- hzar.makeCline1DFreq(mac$hap$obs, scaling, tails)
}

mac.loadmodel("none","none","model1");
mac.loadmodel("fixed" ,"none","model6");


## Observations were between 0 and 1511 km
mac$hap$models <- sapply(mac$hap$models,
                         hzar.model.addBoxReq,
                         0 , 1600,
                         simplify=FALSE)

## Check the updated settings
print(mac$hap$models)

## Check the updated settings
print(mac$hap$models)

#save output
#checkupdatedsettings<-capture.output(print(mac$hap$models))

#cat(checkupdatedsettings,file=paste0(locnames,"_checkupdatedsettings.txt"),sep="",append=TRUE)

## Compile each of the models to prepare for fitting
mac$hap$fitRs$init <- sapply(mac$hap$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=mac$hap$obs,
                             verbose=TRUE,
                             simplify=FALSE)

mac$hap$fitRs$init$model1$mcmcParam$chainLength<-chainLength; #1e6

mac$hap$fitRs$init$model1$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

mac$hap$fitRs$init$model1$mcmcParam$seed[[1]] <-
  mainSeed$A

mac$hap$fitRs$init$model6$mcmcParam$chainLength <-
  chainLength; #1e5

mac$hap$fitRs$init$model6$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

mac$hap$fitRs$init$model6$mcmcParam$seed[[1]] <-
  mainSeed$B


## Check fit request settings
print(mac$hap$fitRs$init)


## Run just one of the models for an initial chain
mac$hap$runs$init <- 
  list()

mac$hap$runs$init$model1 <-
  hzar.doFit(mac$hap$fitRs$init$model1)

## Run another model for an initial chain
mac$hap$runs$init$model6 <-
  hzar.doFit(mac$hap$fitRs$init$model6)

## Compile a new set of fit requests using the initial chains
mac$hap$fitRs$chains <-
  lapply(mac$hap$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
mac$hap$fitRs$chains <-
  hzar.multiFitRequest(mac$hap$fitRs$chains,
                       each=3,
                       baseSeed=NULL)


# Get parameters in each model.

## center
random_center <- runif(6, 0, 1600) # center for all models
for(i in seq(1, 6)){
  mac$hap$fitRs$chains[[i]]$modelParam$init["center"] <-  random_center[i]
}

## width
random_width <- runif(6, 0, 1600)
for(i in seq(1, 6)){
  mac$hap$fitRs$chains[[i]]$modelParam$init["width"] <- random_width[i]
}



## Go ahead and run a chain of 3 runs for every fit request
mac$hap$runs$chains <- hzar.doChain.multi(mac$hap$fitRs$chains,
                                          doPar=TRUE,
                                          inOrder=FALSE,
                                          count=3)




## Did model1 converge?
summary(do.call(mcmc.list,
                lapply(mac$hap$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

plot(do.call(mcmc.list,
             lapply(mac$hap$runs$chains[1:3],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_1_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(mac$hap$runs$chains[1:3],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_1_convergence,file=paste0(locnames,"_Check_model_1_convergence.txt"),sep="",append=TRUE)


## Did model6 converge?
summary(do.call(mcmc.list,
                lapply(mac$hap$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )


par("mar")
par(mar=c(1,1,1,1))
plot(do.call(mcmc.list,
             lapply(mac$hap$runs$chains[4:6],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )


#save output
Checkmodel_6_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(mac$hap$runs$chains[4:6],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_6_convergence,file=paste0(locnames,"_Check_model_6_convergence.txt"),sep="",append=TRUE)


## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
mac$hap$analysis$initDGs <- list(
  nullModel = hzar.dataGroup.null(mac$hap$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
mac$hap$analysis$initDGs$model1 <-
  hzar.dataGroup.add(mac$hap$runs$init$model1)

mac$hap$analysis$initDGs$model6 <-
  hzar.dataGroup.add(mac$hap$runs$init$model6)



##Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, model1......to model15).
mac$hap$analysis$oDG <-
  hzar.make.obsDataGroup(mac$hap$analysis$initDGs)

mac$hap$analysis$oDG <-
  hzar.copyModelLabels(mac$hap$analysis$initDGs,
                       mac$hap$analysis$oDG)


#save output
Checkdatagroupobjs<-capture.output(print(summary(mac$hap$analysis$oDG$data.groups)))

cat(Checkdatagroupobjs,file=paste0(locnames,"_Check_dataGroup_objs.txt"),sep="",append=TRUE)

## Compare the 2 cline models to the null model graphically
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_comparing2models",".png"),pointsize=8)
hzar.plot.cline(mac$hap$analysis$oDG);
dev.off()

## Do model selection based on the AICc scores
print(mac$hap$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(mac$hap$analysis$oDG));


#save output
AICctableforallmodels<-capture.output(print(mac$hap$analysis$AICcTable <-
                                              hzar.AICc.hzar.obsDataGroup(mac$hap$analysis$oDG)))

cat(AICctableforallmodels,file=paste0(locnames,"_AICc_table_for_all_models.txt"),sep="",append=TRUE)

## Print out the model with the minimum AICc score
print(mac$hap$analysis$model.name <-
        rownames(mac$hap$analysis$AICcTable
        )[[ which.min(mac$hap$analysis$AICcTable$AICc )]])

#save output
SelectedModel<-capture.output(print(mac$hap$analysis$model.name <-
                                      rownames(mac$hap$analysis$AICcTable
                                      )[[ which.min(mac$hap$analysis$AICcTable$AICc )]]))

cat(SelectedModel,file=paste0(locnames,"_selected_model.txt"),sep="",append=TRUE)

## Extract the hzar.dataGroup object for the selected model
mac$hap$analysis$model.selected <-
  mac$hap$analysis$oDG$data.groups[[mac$hap$analysis$model.name]]

#save the oupt
Var_params<-capture.output(print(hzar.getLLCutParam(mac$hap$analysis$model.selected,
                                                    names(mac$hap$analysis$model.selected$data.param))))
cat(Var_params,file=paste0(locnames,"_MaxLL_var_params_for_selected_model.txt"),sep="",append=TRUE)

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(mac$hap$analysis$model.selected))

#save the oupt
MaxLL_params<-capture.output(print(hzar.get.ML.cline(mac$hap$analysis$model.selected)))
cat(MaxLL_params,file=paste0(locnames,"_MaxLL_params_for_selected_model.txt"),sep="",append=TRUE)

## Plot the maximum likelihood cline for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_maxLL_selectedmodel",".png"),pointsize=8)
hzar.plot.cline(mac$hap$analysis$model.selected);
dev.off()

## Plot the 95% credible cline region for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_fuzzycline_selectedmodel",".png"),pointsize=8)
hzar.plot.fzCline(mac$hap$analysis$model.selected);
dev.off()

#Rename loc <- replace "cur" with loc name
#names(gd)[[j]]=as.character(locnames)

model.selected <- mac$hap$analysis$model.selected

# Save selected.model object to .rds file for later analysis
loc_obj_rds <- paste0(locnames, ".rds")
saveRDS(model.selected, loc_obj_rds)

####plots clines####

#load selected models

hapgp_b1a_f<-readRDS(file = "b1a_f.rds")
hapgp_b1b_f<-readRDS(file = "b1b_f.rds")
hapgp_b2_3_f<-readRDS(file = "b2_b3_f.rds")


library("scales") 

hzar.plot.obsData(hapgp_b2_3_f, main="cline bestfit models for cox1f haplogroups in females")
hzar.plot.fzCline(hapgp_b1a_f, lty=2, fzCol=alpha("#b3b3b3ff",0.5),pch=1,add=TRUE)
hzar.plot.fzCline(hapgp_b1b_f, lty=2, fzCol=alpha("#ff6600ff",0.5),pch=2,add=TRUE)
legend('topright', legend=c("b2/b3", "b1a", "b1b"),pch=c(3,1,2), cex=0.4)




