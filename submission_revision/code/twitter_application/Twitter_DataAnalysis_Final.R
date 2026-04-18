##############################################################
# Twitter Application: Cluster Curve Visualization
#
# Loads pre-computed estimates and cluster labels, then produces
# four output figures corresponding to paper figures:
#   - Paper Figure 2 (right panel): scatter plot of the first two FPC scores
#   - Paper Figure 3 (top row):    estimated curves for Cluster 1 (vec1, large)
#   - Paper Figure 3 (middle row): estimated curves for Cluster 2 (vec2, medium)
#   - Paper Figure 3 (bottom row): estimated curves for Cluster 0 (vec0, small)
# Each curve plot has 5 panels (left to right): p1, p2, p3, Z1, Z2
#
# Input:  ../../data/processed/Twiiter_figure_logit_mul_final_Dec.RData
# Output: ../../outputs/figures/
##############################################################
library(dbscan)
library(elbow)
library(ggplot2)

# Load pre-computed cluster assignments and estimated probability/intensity curves
load("../../data/processed/Twiiter_figure_logit_mul_final_Dec.RData")

# Optionally load the raw W matrix to plot example individual trajectories
# (uncomment the block below and library calls to reproduce trajectory plots)
# load("../../data/processed/W_matrix_final.RData")
# 
# library(cfda)
# start_number <- 1670
# end_number <- 1730
# user_index_sub <- c(sample(vec0,2),sample(vec1,2))
# graph_data_set <- data.frame(W_matrix_final[c(user_index_sub[c(1,2)],
#                                               user_index_sub[c(3,4)]),start_number:end_number])
# graph_data_set$id <- c(user_index_sub )
# 
# library(reshape2)
# df_long <- melt(graph_data_set,
#                 id.vars = "id",
#                 variable.name = "time",
#                 value.name = "state")
# df_long$time <- as.character(df_long$time)
# df_long$time <- substr(df_long$time, 2, nchar(df_long$time))
# df_long$time <- as.numeric(df_long$time)
# library(cfda)
# plotData(df_long, col = c("red",  "green", "blue")) +
#   labs(title = "")+
#   scale_x_continuous(
#     breaks = c(min(timestamps01[start_number:end_number]), max(timestamps01[start_number:end_number])),
#     labels = c("March 11, 2017","March 12, 2017")
#   ) +
#   theme(legend.position = "none",
#         text=element_text(size = 15))
# 
# plotData(df_long, col = c("red",  "green", "blue")) +
#   labs(title = "")+
#   scale_x_continuous(
#     breaks = c(min(timestamps01[start_number:end_number]), max(timestamps01[start_number:end_number])),
#     labels = c("March 11, 2017","March 12, 2017")
#   ) +
#   theme(legend.position = "none",
#         text=element_text(size = 12))


# --- Paper Figure 2 (right panel): Scatter plot of FPC scores coloured by cluster assignment ---
tps2 <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + 
  geom_point(aes(shape=Cluster),size=3)+
  ggtitle(paste0("Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2]))) + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))+
  scale_x_continuous(limits = c(-5500, 500)) +  # Set x-axis limits
  scale_y_continuous(limits = c(-1600, 1300))

out_dir <- "../../outputs/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "twitter_cluster_scatter.png"), tps2, width = 10, height = 8, dpi = 300)

t <- timestamps01

# Shared y-axis limits for p2 are computed from clusters 1 and 2 only;
# cluster 0 uses its own wider limits (see below)
p_min22 <- min(c(p2_figure_final[, c(vec2, vec1)]))
p_max22 <- max(c(p2_figure_final[, c(vec2, vec1)]))

# --- Paper Figure 3 (top row): Estimated curves for Cluster 1 (large cluster, vec1) ---
# Panels (left to right): p1, p2, p3, Z1, Z2
# Individual trajectories shown in grey; cluster median overlaid in red dashed
png(file.path(out_dir, "twitter_cluster1_curves.png"), width = 3000, height = 600, res = 150)
par(mfrow=c(1,5))

meanp=apply(t(p1_figure_final)[vec1,],2,median)
matplot(t, t(t(p1_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp2=apply(t(p2_figure_final)[vec1,],2,median)
matplot(t, t(t(p2_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min22,p_max22))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp3=apply(t(p3_figure_final)[vec1,],2,median)
matplot(t,t(t(p3_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min3,p_max3))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz=apply(t(Z1_figure)[vec1,],2,median)
matplot(t, t(t(Z1_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz2=apply(t(Z2_figure)[vec1,],2,median)
matplot(t, t(t(Z2_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
dev.off()

# --- Paper Figure 3 (middle row): Estimated curves for Cluster 2 (medium cluster, vec2) ---
png(file.path(out_dir, "twitter_cluster2_curves.png"), width = 3000, height = 600, res = 150)
par(mfrow=c(1,5))
meanp=apply(t(p1_figure_final)[vec2,],2,median)
matplot(t, t(t(p1_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp2=apply(t(p2_figure_final)[vec2,],2,median)
matplot(t, t(t(p2_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min22,p_max22))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp3=apply(t(p3_figure_final)[vec2,],2,median)
matplot(t, t(t(p3_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min3,p_max3))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
meanz=apply(t(Z1_figure)[vec2,],2,median)
matplot(t, t(t(Z1_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz2=apply(t(Z2_figure)[vec2,],2,median)
matplot(t, t(t(Z2_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
dev.off()

# --- Paper Figure 3 (bottom row): Estimated curves for Cluster 0 (small outlier cluster, vec0) ---
# p2 and Z1 use wider y-axis limits to accommodate more extreme values in this cluster
png(file.path(out_dir, "twitter_cluster0_curves.png"), width = 3000, height = 600, res = 150)
par(mfrow=c(1,5))
meanp=apply(t(p1_figure_final)[vec0,],2,median)
matplot(t, t(t(p1_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp2=apply(t(p2_figure_final)[vec0,],2,median)
matplot(t, t(t(p2_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(-600,p_max22))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp3=apply(t(p3_figure_final)[vec0,],2,median)
matplot(t, t(t(p3_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min3,p_max3))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
meanz=apply(t(Z1_figure)[vec0,],2,median)
matplot(t, t(t(Z1_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(-1200,z_max1))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz2=apply(t(Z2_figure)[vec0,],2,median)
matplot(t, t(t(Z2_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
dev.off()

message("Figures saved to: ", normalizePath(out_dir))
