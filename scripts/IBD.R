setwd("~/ownCloud/Stages/2022_CDD_Kreckelbergh_Eugenie/manuscript/analyses_eric/IBD")

# genetic data

as.matrix(read.table("pairwise_fst.txt")) -> fst
fst/(1-fst) -> gen_dist

# geographic data

read.table('sites.txt', sep='\t', header=T) -> sites

library(marmap)

getNOAA.bathy(lon1=-5,lon2=3,lat1=45,lat2=51, resolution=1, keep=T) -> a

trans <- trans.mat(a, min.depth = 0)
out.path <- lc.dist(trans, sites[,3:2], res = "path")
geo_dist <- as.matrix(lc.dist(trans, sites[,3:2], res = "dist"))
rownames(geo_dist) <- rownames(gen_dist)
colnames(geo_dist) <- colnames(gen_dist)

plot(a, image=F, deep=0, shallow=0, step=1, lwd=0.8 )

points(sites[,3], sites[,2], cex=1, pch=21, col="black", bg="lightblue") 
text(sites[,3], sites[,2], sites[,1], adj=c(-0.2,0.2), cex=1, pch=21, col="black", bg="blue") 
text(-3.5, 46, "Bay of\nBiscay", font=3)
text(-1.528188, 50.20050, "English Channel", font=3)

lapply(out.path, lines, col = "orange", lwd = 5, lty = 1) -> dummy

# test of Isolation By Distance (IBD)

geo_dist_F <- geo_dist ; gen_dist_F <- gen_dist
geo_dist_M <- geo_dist ; gen_dist_M <- gen_dist


############# all points #############
# males

geo_dist_F[lower.tri(geo_dist_F)] <- NA
gen_dist_F[lower.tri(geo_dist_F)] <- NA

# females

geo_dist_M[upper.tri(geo_dist_M)] <- NA
gen_dist_M[upper.tri(gen_dist_M)] <- NA

# males vs. females

plot(gen_dist_F ~ geo_dist_F, 
	pch=21, bg=rgb(1,0,0, alpha=0.5), 
	ylim=c(0,1), xlab="geographic distance (km)", ylab="genetic distance (Fst/(1-Fst))")
points(geo_dist_M, gen_dist_M, pch=22, bg=rgb(0,0,1, alpha=0.5))

abline(lm(gen_dist_F[!is.na(gen_dist_F)] ~ geo_dist_F[!is.na(geo_dist_F)]), col=rgb(1,0,0, alpha=0.5))
abline(lm(gen_dist_M[!is.na(gen_dist_M)] ~ geo_dist_M[!is.na(geo_dist_M)]), col=rgb(0,0,1, alpha=0.5))

legend("topleft", legend=c("males", "females"), pch=c(22, 21), 
						pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5))
						)

############# without LeCro et StBrieuc #############

geo_dist_F[ -c(1,2), -c(1,2)] -> geo_dist_F2
gen_dist_F[ -c(1,2), -c(1,2)] -> gen_dist_F2

geo_dist_M[-c(1,2),-c(1,2)] -> geo_dist_M2
gen_dist_M[-c(1,2),-c(1,2)] -> gen_dist_M2

plot(gen_dist_F2 ~ geo_dist_F2, 
	pch=21, bg=rgb(1,0,0, alpha=0.5), #ylim=c(0,1), 
	xlab="geographic distance (km)", ylab="genetic distance (Fst/(1-Fst))")
points(geo_dist_M2, gen_dist_M2, pch=22, bg=rgb(0,0,1, alpha=0.5))

abline(lm(gen_dist_F2[!is.na(gen_dist_F2)] ~ geo_dist_F2[!is.na(geo_dist_F2)]), col=rgb(1,0,0, alpha=0.5))
abline(lm(gen_dist_M2[!is.na(gen_dist_M2)] ~ geo_dist_M2[!is.na(geo_dist_M2)]), col=rgb(0,0,1, alpha=0.5))

legend("topright", legend=c("males", "females"), pch=c(22, 21), 
						pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5))
						)
										
read.table("dist_sex_fst.txt", head=T) -> tab

summary(lm(formula = fst ~ sex + dist + sex*dist, data = tab))

				

read.table("pairwise_fst_astable.txt", head=T) -> tab
plot(male~female, data=tab)
abline(a=0,b=1)
abline(lm(male~female, data=tab), col=rgb(0,0,1, alpha=0.5))

read.table("fstMF.txt", head=T) -> tab
boxplot(tab$fstMF~tab$site)
