

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  ## fonction qui permet de faire des images avec plusieurs graphiques avec ggplot2
  library(grid)
  
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

make_plot = function(data, family, rate_cens, n, h_0.01, title = ""){
  if (h_0.01){
    sub_data = data[(data$family == family) &
                      (data$rate_cens == rate_cens) &
                      (data$n == n), ]
  } else {
    sub_data = data[(data$family == family) &
                      (data$rate_cens == rate_cens) &
                      (data$n == n) &
                      (data$h != 0.01), ]
  }
  
  to_plot = melt(data = sub_data,
                 id.vars = "h",
                 measure.vars = c("est_par_var",
                                  "est_par_bias_sq",
                                  "est_par_error_sq_mean"))
  
  plot = ggplot(data = to_plot, aes(x = factor(h), y =  value, group = variable)) +
    geom_line(size = 1, aes(linetype = variable)) +
    geom_point(aes(shape = variable), size = 6) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 24),
          axis.title.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          title = element_text(size = 27),
          legend.title = element_text(size = 27),
          legend.text = element_text(size=24),
          legend.key.size = unit(1,"cm")) +
    ylab("error") +
    xlab("h") +
    # scale_y_continuous(#"", 
    #                    limits=c(0,max(to_plot$value))
    #                    #labels = c("variance", "bias", "tot. error")
    #                    ) +
    scale_shape_discrete("",
                         labels  = c("variance", "bias", "tot. error")) +
    scale_linetype_discrete("", labels = c("variance", "bias", "tot. error")) +
    ggtitle(title)
  
  return(list(plot = plot, to_plot = to_plot))
}


library(ggplot2)
library(reshape2)
library(dplyr)


#setwd("~/Google Drive/GitHub/simulations copules/")
setwd("~/Dropbox/Yohann/copules/application - simulated data")

results = read.csv("output_2018-04-19/result_simulated_data.csv")



results[(results$family == 4 &
           results$rate_cens == 0.3 &
           results$n == 1000 &
           results$h == 0.15), ]


results_summary1 = results %>%
  group_by(family, rate_cens, n, h, x) %>%
  mutate(est_par_mean = mean(est_par, na.rm = T),
         est_par_var = sd(est_par, na.rm = T)^2) %>%
  mutate(est_par_bias_sq = (exact_par - est_par_mean)^2,
         est_par_error_sq = (est_par - exact_par)^2)

head(data.frame(results_summary1))
head(data.frame(results_summary1)[results_summary1$h == 0.05,])


results_summary2 = results_summary1 %>%
  group_by(family, rate_cens, n, h, x) %>%
  summarise(
    est_par_mean = unique(est_par_mean),
    exact_par = unique(exact_par),
    est_par_var = unique(est_par_var),
    est_par_bias_sq = unique(est_par_bias_sq),
    est_par_error_sq_mean = mean(est_par_error_sq, na.rm = T))

head(data.frame(results_summary2)[(results_summary2$h == 0.05 & results_summary2$n == 500),])


results_summary_x = results_summary2 %>%
  group_by(family, rate_cens, n, h) %>%
  summarise(est_par_var = mean(est_par_var),
            est_par_bias_sq = mean(est_par_bias_sq),
            est_par_error_sq_mean = mean(est_par_error_sq_mean))

data.frame(results_summary_x)[(results_summary_x$h == 0.05 & results_summary_x$n == 1000),]

data.frame(results_summary_x)[(results_summary_x$family == 4 & 
                                 #results_summary_x$rate_cens == 0.05 &
                                 results_summary_x$n == 1000),]


# --------------------------------------------------
#   make graphics
# --------------------------------------------------

### n = 500

plot_1_03_500 = make_plot(data = results_summary_x,
                          family = 1,
                          rate_cens = 0.3,
                          n = 500,
                          h_0.01 = F,
                          title = "Gaussian, q = 0.3, n = 500")

plot_3_03_500 = make_plot(data = results_summary_x,
                          family = 3,
                          rate_cens = 0.3,
                          n = 500,
                          h_0.01 = F,
                          title = "Clayton, q = 0.3, n = 500")

plot_4_03_500 = make_plot(data = results_summary_x,
                          family = 4,
                          rate_cens = 0.3,
                          n = 500,
                          h_0.01 = F,
                          title = "Gumbel, q = 0.3, n = 500")

plot_5_03_500 = make_plot(data = results_summary_x,
                          family = 5,
                          rate_cens = 0.3,
                          n = 500,
                          h_0.01 = F,
                          title = "Frank, q = 0.3, n = 500")


plot_1_05_500 = make_plot(data = results_summary_x,
                          family = 1,
                          rate_cens = 0.5,
                          n = 500,
                          h_0.01 = F,
                          title = "Gaussian, q = 0.5, n = 500")

plot_3_05_500 = make_plot(data = results_summary_x,
                          family = 3,
                          rate_cens = 0.5,
                          n = 500,
                          h_0.01 = F,
                          title = "Clayton, q = 0.5, n = 500")

plot_4_05_500 = make_plot(data = results_summary_x,
                          family = 4,
                          rate_cens = 0.5,
                          n = 500,
                          h_0.01 = F,
                          title = "Gumbel, q = 0.5, n = 500")

plot_5_05_500 = make_plot(data = results_summary_x,
                          family = 5,
                          rate_cens = 0.5,
                          n = 500,
                          h_0.01 = F,
                          title = "Frank, q = 0.5, n = 500")


multiplot(plot_1_03_500$plot + theme(legend.position = "None") + ylim(c(0, 0.11)),
          plot_3_03_500$plot + theme(legend.position = "None") + ylim(c(0, 1.1)),
          plot_4_03_500$plot + theme(legend.position = "None") + ylim(c(0, 0.5)),
          plot_5_03_500$plot + theme(legend.position = "None") + ylim(c(0, 10)),
          plot_1_05_500$plot + theme(legend.position = "None") + ylim(c(0, 0.11)),
          plot_3_05_500$plot + theme(legend.position = "None") + ylim(c(0, 1.1)),
          plot_4_05_500$plot + theme(legend.position = "None") + ylim(c(0, 0.5)),
          plot_5_05_500$plot + theme(legend.position = "None") + ylim(c(0, 10)),
          cols = 2
)

#####################

plot_1_0_500 = make_plot(data = results_summary_x,
                         family = 1,
                         rate_cens = 0,
                         n = 500,
                         h_0.01 = F,
                         title = "Gaussian, q = 0, n = 500")

plot_3_0_500 = make_plot(data = results_summary_x,
                         family = 3,
                         rate_cens = 0,
                         n = 500,
                         h_0.01 = F,
                         title = "Clayton, q = 0, n = 500")

plot_4_0_500 = make_plot(data = results_summary_x,
                         family = 4,
                         rate_cens = 0,
                         n = 500,
                         h_0.01 = F,
                         title = "Gumbel, q = 0, n = 500")


plot_1_01_500 = make_plot(data = results_summary_x,
                          family = 1,
                          rate_cens = 0.1,
                          n = 500,
                          h_0.01 = F,
                          title = "Gaussian, q = 0.1, n = 500")

plot_3_01_500 = make_plot(data = results_summary_x,
                          family = 3,
                          rate_cens = 0.1,
                          n = 500,
                          h_0.01 = F,
                          title = "Clayton, q = 0.1, n = 500")

plot_4_01_500 = make_plot(data = results_summary_x,
                          family = 4,
                          rate_cens = 0.1,
                          n = 500,
                          h_0.01 = F,
                          title = "Gumbel, q = 0.1, n = 500")

multiplot(plot_1_0_500$plot + theme(legend.position = "None") + ylim(c(0, 0.042)),
          plot_3_0_500$plot + theme(legend.position = "None") + ylim(c(0, 0.44)),
          plot_4_0_500$plot + theme(legend.position = "None") + ylim(c(0, 0.22)),
          plot_1_01_500$plot + theme(legend.position = "None") + ylim(c(0, 0.042)),
          plot_3_01_500$plot + theme(legend.position = "None") + ylim(c(0, 0.44)),
          plot_4_01_500$plot + ylim(c(0, 0.22)),
          cols = 2
)

### n = 1000

plot_1_03_1000 = make_plot(data = results_summary_x,
                           family = 1,
                           rate_cens = 0.3,
                           n = 1000, 
                           h_0.01 = F,
                           title = "Gaussian, q = 0.3, n = 1000")

plot_3_03_1000 = make_plot(data = results_summary_x,
                           family = 3,
                           rate_cens = 0.3,
                           n = 1000, 
                           h_0.01 = F,
                           title = "Clayton, q = 0.3, n = 1000")

plot_4_03_1000 = make_plot(data = results_summary_x,
                           family = 4,
                           rate_cens = 0.3,
                           n = 1000, 
                           h_0.01 = F,
                           title = "Gumbel, q = 0.3, n = 1000")

plot_5_03_1000 = make_plot(data = results_summary_x,
                           family = 5,
                           rate_cens = 0.3,
                           n = 1000, 
                           h_0.01 = F,
                           title = "Frank, q = 0.3, n = 1000")

plot_1_05_1000 = make_plot(data = results_summary_x,
                           family = 1,
                           rate_cens = 0.5,
                           n = 1000,
                           h_0.01 = F,
                           title = "Gaussian, q = 0.5, n = 1000")


plot_3_05_1000 = make_plot(data = results_summary_x,
                           family = 3,
                           rate_cens = 0.5,
                           n = 1000,
                           h_0.01 = F,
                           title = "Clayton, q = 0.5, n = 1000")

plot_4_05_1000 = make_plot(data = results_summary_x,
                           family = 4,
                           rate_cens = 0.5,
                           n = 1000,
                           h_0.01 = F,
                           title = "Gumbel, q = 0.5, n = 1000")

plot_5_05_1000 = make_plot(data = results_summary_x,
                           family = 5,
                           rate_cens = 0.5,
                           n = 1000,
                           h_0.01 = F,
                           title = "Frank, q = 0.5, n = 1000")

multiplot(plot_1_03_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.07)),
          plot_3_03_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.5)),
          plot_4_03_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.3)),
          plot_5_03_1000$plot + theme(legend.position = "None") + ylim(c(0, 4)),
          plot_1_05_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.07)),
          plot_3_05_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.5)),
          plot_4_05_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.3)),
          plot_5_05_1000$plot + theme(legend.position = "None") + ylim(c(0, 4)),
          cols = 2
)

####################

plot_1_0_1000 = make_plot(data = results_summary_x,
                          family = 1,
                          rate_cens = 0,
                          n = 1000,
                          h_0.01 = F,
                          title = "Gaussian, q = 0, n = 1000")

plot_3_0_1000 = make_plot(data = results_summary_x,
                          family = 3,
                          rate_cens = 0,
                          n = 1000,
                          h_0.01 = F,
                          title = "Clayton, q = 0, n = 1000")

plot_4_0_1000 = make_plot(data = results_summary_x,
                          family = 4,
                          rate_cens = 0,
                          n = 1000,
                          h_0.01 = F,
                          title = "Gumbel, q = 0, n = 1000")


plot_1_01_1000 = make_plot(data = results_summary_x,
                           family = 1,
                           rate_cens = 0.1,
                           n = 1000,
                           h_0.01 = F,
                           title = "Gaussian, q = 0.1, n = 1000")

plot_3_01_1000 = make_plot(data = results_summary_x,
                           family = 3,
                           rate_cens = 0.1,
                           n = 1000,
                           h_0.01 = F,
                           title = "Clayton, q = 0.1, n = 1000")

plot_4_01_1000 = make_plot(data = results_summary_x,
                           family = 4,
                           rate_cens = 0.1,
                           n = 1000,
                           h_0.01 = F,
                           title = "Gumbel, q = 0.1, n = 1000")

multiplot(plot_1_0_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.02)),
          plot_3_0_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.25)),
          plot_4_0_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.08)),
          plot_1_01_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.02)),
          plot_3_01_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.25)),
          plot_4_01_1000$plot + theme(legend.position = "None") + ylim(c(0, 0.08)),
          cols = 2)


