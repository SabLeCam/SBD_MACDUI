setwd("Path_to_your_datafiles")

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


#############################
# Sex-specific IBD analysis
#############################
##### Phist values estimated using Arlequin
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



# Example: mat_fst_m, mat_fst_f, geo_dist already defined
# Transform FST/(1-FST)
mat_fst_m <- mat_fst_m / (1 - mat_fst_m)
mat_fst_f <- mat_fst_f / (1 - mat_fst_f)

# Mirror lower triangle to upper
mat_fst_m[upper.tri(mat_fst_m)] <- t(mat_fst_m)[upper.tri(mat_fst_m)]
mat_fst_f[upper.tri(mat_fst_f)] <- t(mat_fst_f)[upper.tri(mat_fst_f)]

# Function to convert pairwise matrices to long dataframe
pairwise_df <- function(mat_fst, geo_dist, sex_label){
  idx <- lower.tri(mat_fst)
  data.frame(
    pop1 = rownames(mat_fst)[row(mat_fst)[idx]],
    pop2 = colnames(mat_fst)[col(mat_fst)[idx]],
    distance = geo_dist[idx],
    phist = mat_fst[idx],
    sex = sex_label
  )
}

# Create dataframe for males and females
df_m <- pairwise_df(mat_fst_m, geo_dist, "Male")
df_f <- pairwise_df(mat_fst_f, geo_dist, "Female")

# Combine
df <- rbind(df_m, df_f)
df$sex_num <- ifelse(df$sex == "Female", 1, 0)
df$log_distance <- log(df$distance)

# -------------------
# 1. Sex-specific slopes
# -------------------
lm_male <- lm(phist ~ distance, data = df[df$sex=="Male",])
lm_female <- lm(phist ~ distance, data = df[df$sex=="Female",])

slope_m <- coef(lm_male)["distance"]
slope_f <- coef(lm_female)["distance"]

cat("Male slope:", slope_m, "\n")
#Male slope: 0.000901103 
cat("Female slope:", slope_f, "\n")
#Female slope: 0.0004252056 

# -------------------
# 2. Male:Female dispersal ratio
# sigma_m / sigma_f = sqrt(slope_f / slope_m)
# -------------------
disp_ratio <- sqrt(slope_f / slope_m)
cat("Male:Female dispersal ratio:", disp_ratio, "\n")
#Male:Female dispersal ratio: 0.6869297 

# -------------------
# 3. Combined regression with interaction
# -------------------
lm_ibd <- lm(phist ~ distance * sex, data = df)
summary(lm_ibd)

#Coefficients:
#.                Estimate Std. Error t value Pr(>|t|)   
#(Intercept)      -0.0537623  0.0651509  -0.825  0.41677   
#distance          0.0004252  0.0001195   3.557  0.00147 **
#  sexMale          -0.0263514  0.0921373  -0.286  0.77714   
#distance:sexMale  0.0004759  0.0001690   2.815  0.00917 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.1405 on 26 degrees of freedom
#Multiple R-squared:  0.7616,	Adjusted R-squared:  0.7341 
#F-statistic: 27.69 on 3 and 26 DF,  p-value: 2.963e-08




# -------------------
# 4. Permutation test for slope difference
# -------------------
set.seed(123)
nperm <- 9999
perm_diff <- numeric(nperm)

for(i in 1:nperm){
  df$sex_perm <- sample(df$sex)
  
  lm_m_perm <- lm(phist ~ distance, data = df[df$sex_perm=="Male",])
  lm_f_perm <- lm(phist ~ distance, data = df[df$sex_perm=="Female",])
  
  perm_diff[i] <- coef(lm_m_perm)["distance"] - coef(lm_f_perm)["distance"]
}

obs_diff <- slope_m - slope_f
p_perm <- mean(abs(perm_diff) >= abs(obs_diff))

cat("Observed slope difference:", obs_diff, "\n")
#Observed slope difference: 0.0004758973 
cat("Permutation p-value:", p_perm, "\n")
#Permutation p-value: 0.01880188
# -------------------
# 5. Optional: visualize
# -------------------

library(ggplot2)
ggplot(df, aes(distance, phist, color=sex)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=TRUE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(Phi[ST] / (1 - Phi[ST])),
    color = "Sex"
  )

