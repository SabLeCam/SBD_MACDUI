setwd("~/Documents/MAC_DUI/SBD_MS/analyses_eric/IBD/")

# geographic data

read.table('sites.txt', sep='\t', header=T) -> sites

library(marmap)

getNOAA.bathy(lon1=-5,lon2=3,lat1=45,lat2=51, resolution=1, keep=T) -> a

trans <- trans.mat(a, min.depth = 0)
out.path <- lc.dist(trans, sites[,3:2], res = "path")
geo_dist <- as.matrix(lc.dist(trans, sites[,3:2], res = "dist"))
rownames(geo_dist) <- c("CRO","BRI","MAH","BRE","AYT","FOU")
colnames(geo_dist) <- c("CRO","BRI","MAH","BRE","AYT","FOU")

plot(a, image=F, deep=0, shallow=0, step=1, lwd=0.8 )
points(sites[,3], sites[,2], cex=1, pch=21, col="black", bg="lightblue") 
text(sites[,3], sites[,2], sites[,1], adj=c(-0.2,0.2), cex=1, pch=21, col="black", bg="blue") 
text(-3.5, 46, "Bay of\nBiscay", font=3)
text(-1.528188, 50.20050, "English Channel", font=3)



########## slope comparison ###############
fst_f<-"0.00000
0.10564 0.00000				
0.21094	0.13146 0.00000			
0.31973	0.16998	0.00198	0.00000	
0.20347	0.07103	0.00008	0.00684	0.00000
0.26211	0.16145	-0.02257	-0.00543	0.00833 0.00000"

mat_fst_f <- data.matrix( read.table(text=fst_f, fill=TRUE, col.names=paste("V", 1:6))  )
mat_fst_f[upper.tri(mat_fst_f)] <- t(mat_fst_f)[upper.tri(mat_fst_f)]
rownames(mat_fst_f)<-c("CRO","BRI","MAH","BRE","AYT","FOU")
colnames(mat_fst_f)<-c("CRO","BRI","MAH","BRE","AYT","FOU")
mat_fst_f
mat_fst_f<-mat_fst_f/(1-mat_fst_f)

fst_m<-"0.00000
0.061350 0.00000
0.48989	0.41417 0.00000
0.28952	0.24094	0.05498 0.00000
0.39586	0.33159	0.02117	-0.01044 0.00000
0.41733	0.35154	-0.01035 0.0057	-0.0148 0.00000"

mat_fst_m <- data.matrix( read.table(text=fst_m, fill=TRUE, col.names=paste("V", 1:6))  )
mat_fst_m[upper.tri(mat_fst_f)] <- t(mat_fst_m)[upper.tri(mat_fst_m)]
rownames(mat_fst_m)<-c("CRO","BRI","MAH","BRE","AYT","FOU")
colnames(mat_fst_m)<-c("CRO","BRI","MAH","BRE","AYT","FOU")
mat_fst_m
mat_fst_m<-mat_fst_m/(1-mat_fst_m)

log_geo_dist<-log(geo_dist)

library(bootstrap)
library(prabclus)
library(reshape2)

m_geo<-melt(log_geo_dist)
geo <- m_geo[m_geo$Var1 != m_geo$Var2,]


mfst<-melt(mat_fst_m)
mfst <- mfst[mfst$Var1 != mfst$Var2,]

ffst<-melt(mat_fst_f)
ffst <- ffst[ffst$Var1 != ffst$Var2,]

fst<-rbind(ffst,mfst)
#1=F, 2=M
#fst$sex<-c(rep("1",30),rep("1",30))
geo<-rbind(geo,geo)
#geo$sex<-c(rep("1",30),rep("2",30))
 
sex<-c(rep("F",30),rep("M",30))

log_geo_dist <- log(geo_dist+quantile(as.vector(as.dist(geo_dist)),0.25))

 
H01_reg <- regeqdist(geo$value, fst$value, grouping=sex, groups=c("F","M"))


diffreg <- function(dmx1 , dmx2 , dmy1, dmy2){
  # dmx1, dmx2 = matrices of geographic distances 1 & 2
  # dmy1, dmy2 = matrices of genetic distances 1 & 2
  
  dmxc <- dmyc <- jr <- lmfit <- xvi <- yvi <- list()
  nc <- sediff <- coefdiff <- pval <- condition <- numeric(0)
  dmxc[[1]] <- dmx1
  dmxc[[2]] <- dmx2
  dmyc[[1]] <- dmy1
  dmyc[[2]] <- dmy2
  
  groups <- c(1,2)
  
  for (i in 1:2) {
    jr[[i]] <- list()
    # N ind in matrix i
    nc[i] <- sum(dim(dmxc[[i]])[1])
    # sous matrice geo en vecteur
    xvi[[i]] <- as.vector(as.dist(dmxc[[i]]))
    # sous matrice genet en vecteur
    yvi[[i]] <- as.vector(as.dist(dmyc[[i]]))
  }
  
  xall <- c(xvi[[1]], xvi[[2]])
  xcenter <- mean(xall)
  
  clm <- jackpseudo <- jackestcl <- jackvarcl <- list()
  jackse <- jackest <- tstat <- tdf <- numeric(0)
  for (i in 1:2) {
    jackpseudo[[i]] <- list()
    jackestcl[[i]] <- jackvarcl[[i]] <- numeric(0)
    for (j in 1:2) jackpseudo[[i]][[j]] <- numeric(0)
  }
  for (i in 1:2) {
    #sous matrice geo pondérée
    xvi[[i]] <- xvi[[i]] - xcenter
    lmfit[[i]] <- lm(yvi[[i]] ~ xvi[[i]])
    mm <- model.matrix(~xvi[[i]])
    condition[i] <- kappa(mm)
    clm[[i]] <- coef(lmfit[[i]])
  }
  for (i in 1:2) for (j in 1:2) jr[[i]][[j]] <- bootstrap::jackknife(1:sum(dim(dmxc[[i]])[1]), regdist, dmx = dmxc[[i]], dmy = dmyc[[i]], xcenter = xcenter, param = j)
  
  for (j in 1:2) {
    for (i in 1:2) if (is.na(jr[[i]][[j]]$jack.se)) 
      computable <- FALSE
    coefdiff[j] <- clm[[1]][j] - clm[[2]][j]
    for (i in 1:2) {
      for (k in 1:nc[i]) jackpseudo[[i]][[j]][k] <- nc[i] * 
          clm[[i]][j] - (nc[i] - 1) * jr[[i]][[j]]$jack.values[k]
      jackestcl[[i]][j] <- mean(jackpseudo[[i]][[j]])
      jackvarcl[[i]][j] <- var(jackpseudo[[i]][[j]])
    }
    jackest[j] <- jackestcl[[1]][j] - jackestcl[[2]][j]
    jackse[j] <- sqrt(jackvarcl[[1]][j]/nc[1] + jackvarcl[[2]][j]/nc[2])
    tstat[j] <- jackest[j]/jackse[j]
    tdf[j] <- jackse[j]^4/((jackvarcl[[1]][j]/nc[1])^2/(nc[1] - 
                                                          1) + (jackvarcl[[2]][j]/nc[2])^2/(nc[2] - 1))
    if (tstat[j] > 0) 
      pval[j] <- 2 * (pt(tstat[j], tdf[j], lower.tail = FALSE))
    else pval[j] <- 2 * (pt(tstat[j], tdf[j]))
  }
  
  out <- list(pval = pval, coefdiff = coefdiff, condition = condition, 
              lmfit = lmfit, jr = jr, xcenter = xcenter, tstat = tstat, 
              tdf = tdf, jackest = jackest, jackse = jackse, jackpseudo = jackpseudo, 
              groups = groups)
  class(out) <- "regeqdist"
  out
}

H_reg <- diffreg(log_geo_dist,log_geo_dist,mat_fst_f,mat_fst_m)

#library(emmeans)

#mod <- lm(mat_fst_f[!is.na(mat_fst_f)] ~ mat_fst_m[!is.na(mat_fst_m)])
#summary(mod)


###convert arlequin to genepop (via genind..)
