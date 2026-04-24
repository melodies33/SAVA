library(ggpubr)
library(Cairo)
library(tidyverse)
library(latex2exp)
# parameters to control the size of figures
fd = 1.5
pointsize = 4*fd
lnwidth = 0.8*fd
axissizett = 24*fd
axissizetext = 20*fd
legendsize = 20*fd
titlesize = 26*fd
totaltitlesize = 30*fd
#### plot the TSR for different methods #####
  dat.f = read.csv('amazon-fashion-comprison.csv')
  dat.b = read.csv('amazon-allbeauty-comprison.csv')
  dat.l = read.csv('amazon-luxurybeauty-comprison.csv')
  ggplot()+
    geom_point(aes(x = decision.time, y = Totalselection, shape = Method, color = Method), size = pointsize, data = dat.f)+
    geom_line(aes(x = decision.time, y = Totalselection, color = Method),
              linewidth = lnwidth, data = dat.f)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    scale_x_continuous(breaks = seq(100,1200,200), limits = c(100,1200))+
    labs(x = 'Indices of decision times')+
    ggtitle('Amazon Fashion')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = decision.time, y = Totalselection, shape = Method, color = Method), size = pointsize, data = dat.b)+
    geom_line(aes(x = decision.time, y = Totalselection, color = Method),
              linewidth = lnwidth, data = dat.b)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,400)+
    labs(x = 'Indices of decision times')+
    ggtitle('All Beauty')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  
  ggplot()+
    geom_point(aes(x = decision.time, y = Totalselection, shape = Method, color = Method), size = pointsize, data = dat.l)+
    geom_line(aes(x = decision.time, y = Totalselection, color = Method),
              linewidth = lnwidth, data = dat.l)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1000)+
    labs(x = 'Indices of decision times')+
    ggtitle('Luxury Beauty')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p3 
  
  ggarrange(p1,p2,p3, ncol = 3, common.legend = T, legend = 'bottom')-> p
  annotate_figure(p, top = text_grob(TeX("Selection Number Comparison"), size = totaltitlesize)) -> p
  
  ggsave('amazon-sava.pdf', plot = p, width = 24, height = 9)

#### plot test levels for selected decision times ####


total = read_csv('threshold-compare.csv')
total = total[-which(total$Index == 500),]
total$Index = as.factor(total$Index)
write.csv(total, 'plotdata.csv')

sava_data <- filter(total, Method == "SAVA")

other_methods <- filter(total, Method != "SAVA")
sava_color <- "#E69F00" 
other_colors <- c("LORD++" = "#56B4E9", "SAFFRON" = "#009E73", "ADDIS" = "#CC79A7")
colors <- c("SAVA" = sava_color, other_colors)

shapes <- c("ADDIS" = NA, "LORD++" = NA, "SAFFRON" = NA , "SAVA" = 16)
line_types <- c("ADDIS" = "dotdash", "LORD++" = "dashed","SAFFRON" = "dotted","SAVA" = "solid")

indexlabel = c('Item index : 100', 'Item index : 200', 'Item index : 300', 'Item index : 400')
names(indexlabel) = c('100','200','300','400')

sava_times <- unique(sava_data$Decision_time)

other_methods_ext <- other_methods %>%
  group_by(Index, Arm, Method) %>%
  expand(Decision_time = sava_times) %>%
  
  mutate(Test_level = first(other_methods$Test_level[other_methods$Method == unique(Method)]))

plot_data_all <- bind_rows(
  sava_data,
  other_methods_ext
)

ggplot(filter(plot_data_all, Arm == "A"), 
       aes(x = Decision_time, y = Test_level, 
           color = Method, linetype = Method, shape = Method,
           group = interaction(Method, Index))) +
  geom_line(linewidth = lnwidth) +
  geom_point(data = filter(plot_data_all, Arm == "A" & Method == "SAVA"), 
             size = pointsize, show.legend = TRUE) +
  facet_wrap(~Index, ncol = 2, labeller = labeller(Index = indexlabel)) +
  scale_color_manual(values = colors, name = "Method") +
  scale_linetype_manual(values = line_types, name = "Method") +
  scale_shape_manual(values = shapes, name = "Method") +
  scale_y_continuous(name = "Test Level", limits = c(0, 0.1)) +
  labs(x = "Decision Time", title = "Arm A") +
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext)) +
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize)) +
  theme(strip.text.x = element_text(size = legendsize)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(plot.title = element_text(size = titlesize, hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(
    linetype = line_types,
    shape = shapes
  ))) -> plot_A



ggplot(filter(plot_data_all, Arm == "B"), 
       aes(x = Decision_time, y = Test_level, 
           color = Method, linetype = Method, shape = Method,
           group = interaction(Method, Index))) +
  geom_line(linewidth = lnwidth) +
  geom_point(data = filter(plot_data_all, Arm == "B" & Method == "SAVA"), 
             size = pointsize, show.legend = TRUE) +
  facet_wrap(~Index, ncol = 2, labeller = labeller(Index = indexlabel)) +
  scale_color_manual(values = colors, name = "Method") +
  scale_linetype_manual(values = line_types, name = "Method") +
  scale_shape_manual(values = shapes, name = "Method") +
  scale_y_continuous(name = "Test Level", limits = c(0, 0.1)) +
  labs(x = "Decision Time", title = "Arm B") +
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext)) +
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize)) +
  theme(strip.text.x = element_text(size = legendsize)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(plot.title = element_text(size = titlesize, hjust = 0.5)) +
  
  guides(color = guide_legend(override.aes = list(
    linetype = line_types,
    shape = shapes
  ))) -> plot_B


ggarrange(plot_A,plot_B, ncol = 2, common.legend = T, legend = 'bottom')-> pp
annotate_figure(pp, top = text_grob(TeX(" Comparison of Test Levels"), size = totaltitlesize)) -> pp

ggsave('test level compare.pdf', plot = pp, width = 26, height = 13)
