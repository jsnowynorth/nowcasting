library(R.utils)
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)

# josh box plots ----------------------------------------------------------
rigid <- readRDS("data/Qmse_reviewer.Rdata")
rigid <- rigid[,,-c(59:70)]
warp <- read.csv("data/mses/Warping_All_MSE.csv")

warp <- warp %>% filter(start_step >= 4)

# warp <- Rigid_Warping_All_RMSE %>% filter(method == 'warping')

rigid_df <- NULL

for(i in 4:58){
  df_tmp <- data.frame(start_step = rep(i, 10), 
                       end_step = rep(i, 10) + seq(1,10),
                       forecast_horizon = seq(5, 50, by = 5),
                       forecast_MSE = rigid[,1,i],
                       persistence_MSE = rigid[,2,i],
                       percent_improvement = rigid[,4,i],
                       method = 'rigid')
  
  rigid_df <- rbind(rigid_df, df_tmp)
  
}

df_all <- rbind(rigid_df, warp)

plot3 <- ggplot(df_all, aes(x=factor(forecast_horizon), y=percent_improvement, fill = factor(method))) +
  geom_boxplot(outlier.colour=NA) +
  coord_cartesian(ylim = c(-70, 70)) +
  scale_fill_manual(values=c("white", "darkgrey")) +
  theme_classic() +
  theme(legend.position="none") +
  geom_hline(yintercept=0)+
  labs(x="Forecast Horizon (min)", y="Percent Improvement")

png("/Users/joshuanorth/Desktop/rigid_warping_withoutlier_MSE.png", units="in", width=7, height=5, res=500)
plot3
dev.off()


# percent rigid is better
warp <- warp %>% filter(start_step >= 4)
sum(if_else(rigid_df$percent_improvement - warp$percent_improvement > 0, 1, 0))/length(rigid_df$percent_improvement)

rigid_df %>% group_by(forecast_horizon) %>% select(percent_improvement) %>% summarise(mean(percent_improvement))
rigid_df %>% group_by(forecast_horizon) %>% select(percent_improvement) %>% summarise(median(percent_improvement))

warp %>% group_by(forecast_horizon) %>% select(percent_improvement) %>% summarise(mean(percent_improvement))
warp %>% group_by(forecast_horizon) %>% select(percent_improvement) %>% summarise(median(percent_improvement))











