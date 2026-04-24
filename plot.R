library(ggpubr)
library(Cairo)
library(tidyverse)
library(latex2exp)
library(RColorBrewer)
# total parameters setting for controlling sizes of figures
fd = 1.5
pointsize = 4*fd
lnwidth = 0.8*fd
axissizett = 24*fd
axissizetext = 20*fd
legendsize = 20*fd
titlesize = 26*fd
totaltitlesize = 30*fd
#### counterexample 1 ####
counterex = read_csv('sava_counterexample_1.csv')
pointindex = seq(10,length(counterex$time),10)
ggplot()+
  geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = counterex[pointindex,])+
  geom_line(aes(x = time, y = FSR, color = method),
            linewidth = 0.8, data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits = c(0,0.4))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p1

ggplot()+
  geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = counterex[pointindex,])+
  geom_line(aes(x = time, y = TSR, color = method),
            linewidth = 0.8, data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  ylim(0,0.2)+
  ggtitle('Intermediate TSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p2

ggplot()+
  geom_point(aes(x = time, y = mFSR, shape = method, color = method), size = pointsize, data = counterex[pointindex,])+
  geom_line(aes(x = time, y = mFSR, color = method),
            linewidth = 0.8, data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits = c(0,0.4))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p3

ggarrange(p3,p1,p2, ncol = 3, common.legend = T, legend = 'bottom')-> p
annotate_figure(p, top = text_grob("Methods Comparison", size = 30)) -> p

CairoPNG(filename = 'counter example 1.png', width = 1800,height = 600)
p
dev.off()

ggsave('counter example 1.pdf', plot = p, width = 21, height = 7)

#### counterexample 2 ####
counterex = read_csv('sava_counterexample_2.csv')

pointindex = seq(10,length(counterex$time),10)
ggplot()+
  geom_point(aes(x = time, y = FSR, color = method, shape = method),data = counterex[pointindex,], size = pointsize)+
  geom_line(aes(x = time, y = FSR, color = method),
            linewidth = 0.8,data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  scale_y_continuous(breaks = seq(0,0.5,0.1), limits = c(0,0.5))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p1

ggplot(counterex)+
  geom_point(aes(x = time, y = TSR, color = method, shape = method),data = counterex[pointindex,], size = pointsize)+
  geom_line(aes(x = time, y = TSR, color = method),
            linewidth = 0.8, data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  ylim(0,1)+
  ggtitle('Intermediate TSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p2

ggplot()+
  geom_point(aes(x = time, y = mFSR, color = method, shape = method),data = counterex[pointindex,], size = pointsize)+
  geom_line(aes(x = time, y = mFSR, color = method),
            linewidth = 0.8, data = counterex)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(strip.text.x = element_text(size = 20))+
  scale_y_continuous(breaks = seq(0,0.5,0.1), limits = c(0,0.5))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = 26, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p3


ggarrange(p3,p1,p2, ncol = 3, common.legend = T, legend = 'bottom')-> p
annotate_figure(p, top = text_grob("Methods Comparison", size = 30)) -> p

CairoPNG(filename = 'counter example 2.png', width = 1800,height = 600)
p
dev.off()

ggsave('counter example 2.pdf', plot = p, width = 21, height = 7)

#### results for different k #####
re1 = read_csv('sava oracle k curve p0.02.csv')
re2 = read_csv('sava oracle k curve p0.33.csv')
re3 = read_csv('sava oracle k curve p0.66.csv')

re1 = re1[which(re1$time == 3000),]
re2 = re2[which(re2$time == 3000),]
re3 = re3[which(re3$time == 3000),]
l = length(re1$time)
fsr = c(re1$FSR, re2$FSR, re3$FSR)
tsr = c(re1$TSR, re2$TSR, re3$TSR)
pi = rep(c('0.05','0.33','0.67'), each = l)
k = c(re1$k, re2$k, re3$k)
re = tibble(k = k, FSR = fsr, TSR = tsr, p = pi)
pointindex = seq(5,length(re$k),5)

ggplot()+
  geom_point(aes(x = k, y = FSR, shape = p, color = p), size = pointsize, data = re[pointindex,])+
  geom_line(aes(x = k, y = FSR, color = p),
            linewidth = lnwidth, data = re)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  theme(strip.text.x = element_text(size = legendsize))+
  labs(x = 'k')+
  scale_y_continuous(breaks = seq(0,0.05,0.01), limits = c(0,0.05))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p1

ggplot()+
  geom_point(aes(x = k, y = TSR, shape = p, color = p), size = pointsize, data = re[pointindex,])+
  geom_line(aes(x = k, y = TSR, color = p),
            linewidth = lnwidth, data = re)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  theme(strip.text.x = element_text(size = legendsize))+
  scale_y_continuous(breaks = seq(0.8,1,0.05), limits = c(0.8,1))+
  labs(x = 'k')+
  ggtitle('Intermediate TSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p2

ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom')-> p
annotate_figure(p, top = text_grob(TeX("Effects of different k"), size = totaltitlesize)) -> p

CairoPNG(filename = 'oracle k.png', width = 1600,height = 800)
p
dev.off()
ggsave('oracle k.pdf', plot = p, width = 20, height = 10)
#### Gaussian case ####

#### pi+ ####
{
  sava = read_csv('sava_gaussian_pi.csv')
onlinelord = read_csv('lordpp_gaussian_pi.csv')
onlinesaff = read_csv('saffron_gaussian_pi.csv')
onlineaddis = read_csv('addis_gaussian_pi.csv')
selectindex = c(1:800,1001:1800,2001:2800,3001:3800)

sava = sava[selectindex,]
onlinelord = onlinelord[selectindex,]
onlinesaff = onlinesaff[selectindex,]
onlineaddis = onlineaddis[selectindex,]
l = length(selectindex)
fsr = c(sava$FSR, onlinelord$FSR, onlinesaff$FSR, onlineaddis$FSR)
pi = c(sava$pi, onlinelord$pi, onlinesaff$pi, onlineaddis$pi)
tsr = c(sava$TSR, onlinelord$TSR, onlinesaff$TSR, onlineaddis$TSR)
method = c('SAVA','LORD++', 'SAFFRON', 'ADDIS')
time = c(sava$time, onlinelord$time, onlinesaff$time, onlineaddis$time)
methodvec = rep(method, each = l)
totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
names(pilabel) = c('0.2','0.4','0.6','0.8')

pointindex = seq(80,length(totalpi$time),80)

ggplot()+
  geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
  geom_line(aes(x = time, y = FSR, color = method),
            linewidth = lnwidth, data = totalpi)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
  theme(strip.text.x = element_text(size = legendsize))+
  theme(panel.spacing = unit(2, "lines"))+
  labs(x = 'Indices of decision times')+
  scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p1

ggplot()+
  geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
  geom_line(aes(x = time, y = TSR, color = method),
            linewidth = lnwidth, data = totalpi)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
  theme(strip.text.x = element_text(size = legendsize))+
  ylim(0,1)+
  theme(panel.spacing = unit(2, "lines"))+
  labs(x = 'Indices of decision times')+
  ggtitle('Intermediate TSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p2

ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom')-> p
annotate_figure(p, top = text_grob(TeX("Methods Comparison (constant $\\mu=0.1$)"), size = totaltitlesize)) -> p

CairoPNG(filename = 'gaussian_pi.png', width = 1800,height = 900)
  p
dev.off()
  
ggsave('gaussian_pi.pdf', plot = p, width = 26, height = 13)
}
#### mu #####
{
  sava = read_csv('sava_gaussian_mu.csv')
  onlinelord = read_csv('lordpp_gaussian_mu.csv')
  onlinesaff = read_csv('saffron_gaussian_mu.csv')
  onlineaddis = read_csv('addis_gaussian_mu.csv')
  selectindex = c(1:800,1001:1800,2001:2800,3001:3800)
  sava = sava[selectindex,]
  onlinelord = onlinelord[selectindex,]
  onlinesaff = onlinesaff[selectindex,]
  onlineaddis = onlineaddis[selectindex,]
  l = length(selectindex)
  
  fsr = c(sava$FSR, onlinelord$FSR, onlinesaff$FSR, onlineaddis$FSR)
  mu = c(sava$mu,   onlinelord$mu, onlinesaff$mu, onlineaddis$mu)
  tsr = c(sava$TSR,   onlinelord$TSR, onlinesaff$TSR, onlineaddis$TSR)
  time = c(sava$time, onlinelord$time, onlinesaff$time, onlineaddis$time)
  totalmu = tibble( time = time, FSR = fsr, TSR = tsr, mu = as.factor(mu), method = methodvec)
  
  mulabel = c('mu : 0.05', 'mu : 0.10', 'mu : 0.15', 'mu : 0.20')
  names(mulabel) = c('0.05','0.1','0.15','0.2')
  pointindex = seq(80,length(totalmu$time),80)
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalmu[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totalmu)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~mu, labeller = labeller(mu = mulabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalmu[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totalmu)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~mu, labeller = labeller(mu = mulabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom')-> pp
  annotate_figure(pp, top = text_grob(TeX("Methods Comparison (constant $\\pi^+=0.5$)"), size = totaltitlesize)) -> pp
  
  
  CairoPNG(filename = 'gaussian_mu.png', width = 1800,height = 900)
  pp
  dev.off()
  
  
  ggsave('gaussian_mu.pdf', plot = pp, width = 26, height = 13)
}

#### truncated Gaussian case ####

#### pi+ ####
{
  sava = read_csv('sava_general_pi.csv')
  onlinelord = read_csv('lordpp_general_pi.csv')
  onlinesaff = read_csv('saffron_general_pi.csv')
  onlineaddis = read_csv('addis_general_pi.csv')
  sava = sava[selectindex,]
  onlinelord = onlinelord[selectindex,]
  onlinesaff = onlinesaff[selectindex,]
  onlineaddis = onlineaddis[selectindex,]
  
  fsr = c(sava$FSR,  onlinelord$FSR, onlinesaff$FSR, onlineaddis$FSR)
  pi = c(sava$pi, onlinelord$pi, onlinesaff$pi, onlineaddis$pi)
  tsr = c(sava$TSR, onlinelord$TSR, onlinesaff$TSR, onlineaddis$TSR)
  method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
  
  time = c(sava$time, onlinelord$time, onlinesaff$time, onlineaddis$time)
  methodvec = rep(method, each = l)
  totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
  
  pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
  names(pilabel) = c('0.2','0.4','0.6','0.8')
  
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totalpi)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    labs(x = 'Indices of decision times')+
    scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totalpi)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom')-> p
  annotate_figure(p, top = text_grob(TeX("Methods Comparison (constant $\\mu = 1$)"), size = totaltitlesize)) -> p
  
  CairoPNG(filename = 'general_pi.png', width = 1800,height = 900)
  p
  dev.off()
  
  
  ggsave('general_pi.pdf', plot = p, width = 26, height = 13)
}
#### mu #####
{
  sava = read_csv('sava_general_mu.csv')
  onlinelord = read_csv('lordpp_general_mu.csv')
  onlinesaff = read_csv('saffron_general_mu.csv')
  onlineaddis = read_csv('addis_general_mu.csv')
  sava = sava[selectindex,]
  onlinelord = onlinelord[selectindex,]
  onlinesaff = onlinesaff[selectindex,]
  onlineaddis = onlineaddis[selectindex,]
  
  fsr = c(sava$FSR, onlinelord$FSR, onlinesaff$FSR, onlineaddis$FSR)
  mu = c(sava$mu, onlinelord$mu, onlinesaff$mu, onlineaddis$mu)
  tsr = c(sava$TSR,  onlinelord$TSR, onlinesaff$TSR, onlineaddis$TSR)
  totalmu = tibble( time = time, FSR = fsr, TSR = tsr, mu = as.factor(mu), method = methodvec)
  
  mulabel = c('mu : 0.8', 'mu : 1.0', 'mu : 1.2', 'mu : 1.4')
  names(mulabel) = c('0.8','1','1.2','1.4')
  
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalmu[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totalmu)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~mu, labeller = labeller(mu = mulabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalmu[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totalmu)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~mu, labeller = labeller(mu = mulabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom')-> pp
  annotate_figure(pp, top = text_grob(TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"), size = totaltitlesize)) -> pp
  
  ggsave('general_mu.pdf', plot = pp, width = 26, height = 13)
  
  CairoPNG(filename = 'general_mu.png', width = 1800,height = 900)
  pp
  dev.off()
  
  
}

#### Gamma #####

# change pi+
{
  sava = read_csv('sava_gamma_pi.csv')
  lord = read_csv('lordpp_gamma_pi.csv')
  saff = read_csv('saffron_gamma_pi.csv')
  addi = read_csv('addis_gamma_pi.csv')
  total_rows <- nrow(sava)
  block_size <- 1000
  keep_per_block <- 800
  full_blocks <- total_rows %/% block_size
  indices <- unlist(lapply(1:full_blocks, function(b) {
    start <- (b - 1) * block_size + 1
    end   <- (b - 1) * block_size + keep_per_block
    start:end
  }))
  selectindex = indices
  l = length(selectindex)
  sava = sava[selectindex,]
  lord = lord[selectindex,]
  saff = saff[selectindex,]
  addi = addi[selectindex,]
  
  fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
  pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
  tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
  method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
  
  time = c(sava$time, lord$time, saff$time, addi$time)
  methodvec = rep(method, each = l)
  totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
  

pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
names(pilabel) = c('0.2','0.4','0.6','0.8')

pointindex = seq(80, length(totalpi$time), 80)

ggplot()+
  geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
  geom_line(aes(x = time, y = FSR, color = method),
            linewidth = lnwidth, data = totalpi)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
  theme(strip.text.x = element_text(size = legendsize))+
  theme(panel.spacing = unit(2, "lines"))+
  labs(x = 'Indices of decision times')+
  scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
  ggtitle('Intermediate FSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p1

ggplot()+
  geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
  geom_line(aes(x = time, y = TSR, color = method),
            linewidth = lnwidth, data = totalpi)+
  theme(axis.title = element_text(size = axissizett),
        axis.text = element_text(size = axissizetext))+
  theme(legend.text = element_text(size = legendsize),
        legend.title = element_text(size = legendsize))+
  facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
  theme(strip.text.x = element_text(size = legendsize))+
  theme(panel.spacing = unit(2, "lines"))+
  ylim(0,1)+
  labs(x = 'Indices of decision times')+
  ggtitle('Intermediate TSR level Comparison')+
  theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
  scale_color_brewer(palette = 'Set1') -> p2

ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p

annotate_figure(
  p,
  top = text_grob(
    TeX("Methods Comparison (constant $\\mu_\\delta = 1.3$)"),
    size = totaltitlesize
  )
) -> p

ggsave('total_gamma_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_gamma_multiply.csv')
  lord = read_csv('lordpp_gamma_mumulti.csv')
  saff = read_csv('saffron_gamma_mumulti.csv')
  addi = read_csv('addis_gamma_mumulti.csv')
  
  sava = sava[selectindex,]
  lord = lord[selectindex,]
  saff = saff[selectindex,]
  addi = addi[selectindex,]
  
  fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
  delta = c(sava$multiply, lord$mu, saff$mu, addi$mu)
  tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
  method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
  
  time = c(sava$time, lord$time, saff$time, addi$time)
  methodvec = rep(method, each = l)
  totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
  
  
  deltalabel = c('mu_delta: 1.2', 'mu_delta: 1.3', 'mu_delta: 1.4', 'mu_delta: 1.5')
  names(deltalabel) = c('1.2','1.3','1.4','1.5')
  
  pointindex = seq(80, length(totaldelta$time), 80)
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totaldelta)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    labs(x = 'Indices of decision times')+
    scale_y_continuous(breaks = seq(0,0.01,0.002), limits = c(0,0.01))+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totaldelta)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
  
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_gamma_delta.pdf', plot = p, width = 26, height = 13)
}

#### Beta #####
# change pi
{
  sava = read_csv('sava_beta_pi.csv')
  lord = read_csv('lordpp_beta_pi.csv')
  saff = read_csv('saffron_beta_pi.csv')
  addi = read_csv('addis_beta_pi.csv')
  
  total_rows <- nrow(sava)
  block_size <- 1000
  keep_per_block <- 800
  full_blocks <- total_rows %/% block_size
  indices <- unlist(lapply(1:full_blocks, function(b) {
    start <- (b - 1) * block_size + 1
    end   <- (b - 1) * block_size + keep_per_block
    start:end
  }))
  selectindex = indices
  l = length(selectindex)
  sava = sava[selectindex,]
  lord = lord[selectindex,]
  saff = saff[selectindex,]
  addi = addi[selectindex,]
  
  fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
  pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
  tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
  method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
  
  time = c(sava$time, lord$time, saff$time, addi$time)
  methodvec = rep(method, each = l)
  totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
  
  
  pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
  names(pilabel) = c('0.2','0.4','0.6','0.8')
  
  pointindex = seq(80, length(totalpi$time), 80)
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totalpi)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    labs(x = 'Indices of decision times')+
    scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totalpi)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
  
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\mu_\\delta = 0.14$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_beta_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_beta_delta.csv')
  lord = read_csv('lordpp_beta_mu.csv')
  saff = read_csv('saffron_beta_mu.csv')
  addi = read_csv('addis_beta_mu.csv')
  
  sava = sava[selectindex,]
  lord = lord[selectindex,]
  saff = saff[selectindex,]
  addi = addi[selectindex,]
  
  fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
  delta = c(sava$delta, lord$mu, saff$mu, addi$mu)
  tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
  method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
  
  time = c(sava$time, lord$time, saff$time, addi$time)
  methodvec = rep(method, each = l)
  totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
  
  deltalabel = c('mu_delta: 0.12', 'mu_delta: 0.13', 'mu_delta: 0.14', 'mu_delta: 0.15')
  names(deltalabel) = c('0.12','0.13','0.14','0.15')
  
  
  pointindex = seq(80, length(totaldelta$time), 80)
  
  ggplot()+
    geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
    geom_line(aes(x = time, y = FSR, color = method),
              linewidth = lnwidth, data = totaldelta)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    labs(x = 'Indices of decision times')+
    scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
    ggtitle('Intermediate FSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p1
  
  ggplot()+
    geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
    geom_line(aes(x = time, y = TSR, color = method),
              linewidth = lnwidth, data = totaldelta)+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext))+
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
    theme(strip.text.x = element_text(size = legendsize))+
    theme(panel.spacing = unit(2, "lines"))+
    ylim(0,1)+
    labs(x = 'Indices of decision times')+
    ggtitle('Intermediate TSR level Comparison')+
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
    scale_color_brewer(palette = 'Set1') -> p2
  
  ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
  
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_beta_delta.pdf', plot = p, width = 26, height = 13)
}

#### Beta distribution with coin betting method #####
# change pi+
{
    sava = read_csv('sava_beta_pi_coin.csv')
    lord = read_csv('lordpp_beta_coin_pi.csv')
    saff = read_csv('saffron_beta_coin_pi.csv')
    addi = read_csv('addis_beta_coin_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\mu_\\delta = 0.05$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_beta_coin_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
    sava = read_csv('sava_beta_delta_coin.csv')
  lord = read_csv('lordpp_beta_coin_mu.csv')
  saff = read_csv('saffron_beta_coin_mu.csv')
  addi = read_csv('addis_beta_coin_mu.csv')
  
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$delta, lord$mu, saff$mu, addi$mu)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    
    
    deltalabel = c('mu_delta: 0.04', 'mu_delta: 0.05', 'mu_delta: 0.06', 'mu_delta: 0.07')
    names(deltalabel) = c('0.04','0.05','0.06','0.07')
    
    
    pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
  
  ggsave('total_beta_coin_delta.pdf', plot = p, width = 26, height = 13)
}

#### sub-Gaussian: Gaussian model #####
# change pi+
{
  sava = read_csv('sava_subgau_gauss_pi.csv')
   lord = read_csv('lordpp_sub_gauss_pi.csv')
    saff = read_csv('saffron_sub_gauss_pi.csv')
    addi = read_csv('addis_sub_gauss_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\mu_\\delta = 0.5$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_subGauss_Gauss_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_subgau_gauss_delta.csv')
    lord = read_csv('lordpp_sub_gauss_mu.csv')
    saff = read_csv('saffron_sub_gauss_mu.csv')
    addi = read_csv('addis_sub_gauss_mu.csv')
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$delta, lord$mu, saff$mu, addi$mu)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    
    deltalabel = c('mu_delta: 0.4', 'mu_delta: 0.45', 'mu_delta: 0.5', 'mu_delta: 0.55')
    names(deltalabel) = c('0.4','0.45', '0.5','0.55')    
    
    pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
  ggsave('total_subGauss_Gauss_delta.pdf', plot = p, width = 26, height = 13)
}

#### sub-Gaussian: Uniform #####
# change pi+
{
  sava = read_csv('sava_subgau_unif_pi.csv')
    lord = read_csv('lordpp_sub_unif_pi.csv')
    saff = read_csv('saffron_sub_unif_pi.csv')
    addi = read_csv('addis_sub_unif_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\mu_\\delta = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    

  ggsave('total_subGauss_unif_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_subgau_unif_delta.csv')
    lord = read_csv('lordpp_sub_unif_mu.csv')
    saff = read_csv('saffron_sub_unif_mu.csv')
    addi = read_csv('addis_sub_unif_mu.csv')
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$delta, lord$mu, saff$mu, addi$mu)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    deltalabel = c('mu_delta: 0.45', 'mu_delta: 0.5', 'mu_delta: 0.55', 'mu_delta: 0.6')
    names(deltalabel) = c('0.45','0.5','0.55','0.6')
    
    
    pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
  
  ggsave('total_subGauss_unif_delta.pdf', plot = p, width = 26, height = 13)
}

#### sub-Gaussian: Bernoulli #####
# change pi+
{
  sava = read_csv('sava_subgau_ber_pi.csv')
   lord = read_csv('lordpp_sub_ber_pi.csv')
    saff = read_csv('saffron_sub_ber_pi.csv')
    addi = read_csv('addis_sub_ber_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\mu_\\delta = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
  ggsave('total_subGauss_ber_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_subgau_ber_delta.csv')
    
    lord = read_csv('lordpp_sub_ber_mu.csv')
    saff = read_csv('saffron_sub_ber_mu.csv')
    addi = read_csv('addis_sub_ber_mu.csv')
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$delta, lord$mu, saff$mu, addi$mu)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    deltalabel = c('mu_delta: 0.5', 'mu_delta: 0.6', 'mu_delta: 0.7', 'mu_delta: 0.8')
    names(deltalabel) = c('0.5','0.6','0.7','0.8')
        pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
  ggsave('total_subGauss_ber_delta.pdf', plot = p, width = 26, height = 13)
}

#### Misspecified models: Uniform case: change a #####
# change a
{
  sava1 = read_csv('sava_subgau_unif_a.csv')
  sava1$mu = rep(0.5, length(sava1$time))
  sava2 = read_csv('sava_subgau_unif_a_mu1.csv')
  sava2$mu = rep(1, length(sava1$time))
  sava = rbind(sava1, sava2)
  total_rows <- nrow(sava)     
  block_size <- 1000           
  keep_per_block <- 800       
  full_blocks <- total_rows %/% block_size  
  indices <- unlist(lapply(1:full_blocks, function(b) {
    start <- (b - 1) * block_size + 1
    end   <- (b - 1) * block_size + keep_per_block
    start:end
  }))
  sava <- sava[indices, ]
  sava$mu = as.factor(sava$mu)
  sava$a = as.factor(sava$a)
  alabel = c('a: 0.5', 'a: 1', 'a: 1.25', 'a: 1.5', 'a: 2', 'a: 4')
  names(alabel) = c('0.5','1','1.25','1.5','2','4')
  
  pointindex = seq(80, length(sava$time), 80)
  
  ggplot() +
    geom_point(
      aes(x = time, y = FSR, shape = mu, color = mu),
      size = pointsize,
      data = sava[pointindex,]
    ) +
    geom_line(
      aes(x = time, y = FSR, color = mu),
      linewidth = lnwidth,
      data = sava
    ) +
    labs(color = TeX("$\\mu_{\\delta}$"), 
         shape = TeX("$\\mu_{\\delta}$"))+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext)) +
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~a, labeller = labeller(a = alabel), ncol = 3) +
    theme(strip.text.x = element_text(size = legendsize)) +
    theme(panel.spacing = unit(2, "lines")) +
    labs(x = 'Indices of decision times') +
    scale_y_continuous(breaks = seq(0,0.15,0.05), limits = c(0,0.15)) +
    ggtitle('Intermediate FSR level') +
    theme(plot.title = element_text(size = titlesize, hjust = 0.5)) -> p1
  
  ggplot() +
    geom_point(
      aes(x = time, y = TSR, shape = mu, color = mu),
      size = pointsize,
      data = sava[pointindex,]
    ) +
    geom_line(
      aes(x = time, y = TSR, color = mu),
      linewidth = lnwidth,
      data = sava
    ) +
    labs(color = TeX("$\\mu_{\\delta}$"), 
         shape = TeX("$\\mu_{\\delta}$"))+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext)) +
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~a, labeller = labeller(a = alabel), ncol = 3) +
    theme(strip.text.x = element_text(size = legendsize)) +
    ylim(0,1) +
    theme(panel.spacing = unit(2, "lines")) +
    labs(x = 'Indices of decision times') +
    ggtitle('Intermediate TSR level') +
    theme(plot.title = element_text(size = titlesize, hjust = 0.5)) -> p2
  
  ggarrange(p1, p2, ncol = 2,
            common.legend = T, legend = 'bottom') -> p
  
  annotate_figure(
    p,
    top = text_grob(
      TeX("SAVA Performance"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('sava_unif_a.pdf', plot = p, width = 35, height = 12)
}

####  Misspecified models: Truncated Gaussian case: change K #####
# change K
{
  sava1 = read_csv('sava_truncgauss_K_mu1.5.csv')
  sava1$mu = rep(1.5, length(sava1$time))
  sava2 = read_csv('sava_truncgauss_K.csv')
  sava2$mu = rep(1, length(sava1$time))
  sava = rbind(sava1, sava2)
  total_rows <- nrow(sava)     
  block_size <- 1000           
  keep_per_block <- 800       
  full_blocks <- total_rows %/% block_size  
  indices <- unlist(lapply(1:full_blocks, function(b) {
    start <- (b - 1) * block_size + 1
    end   <- (b - 1) * block_size + keep_per_block
    start:end
  }))
  sava <- sava[indices, ]
  sava$mu = as.factor(sava$mu)
  sava$K = as.factor(sava$K)
  Klabel = c('K: 2', 'K: 3','K: 4','K: 5','K: 6','K: 10')
  names(Klabel) = c('2','3','4','5','6','10')
  
  pointindex = seq(80, length(sava$time), 80)
  
  ggplot() +
    geom_point(
      aes(x = time, y = FSR, shape = mu, color = mu),
      size = pointsize,
      data = sava[pointindex,]
    ) +
    geom_line(
      aes(x = time, y = FSR, color = mu),
      linewidth = lnwidth,
      data = sava
    ) +
    labs(color = TeX("$\\mu_{\\delta}$"), 
         shape = TeX("$\\mu_{\\delta}$"))+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext)) +
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~K, labeller = labeller(K = Klabel), ncol = 3) +
    theme(strip.text.x = element_text(size = legendsize)) +
    theme(panel.spacing = unit(2, "lines")) +
    labs(x = 'Indices of decision times') +
    scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005)) +
    ggtitle('Intermediate FSR level') +
    theme(plot.title = element_text(size = titlesize, hjust = 0.5)) -> p1
  
  ggplot() +
    geom_point(
      aes(x = time, y = TSR, shape = mu, color = mu),
      size = pointsize,
      data = sava[pointindex,]
    ) +
    geom_line(
      aes(x = time, y = TSR, color = mu),
      linewidth = lnwidth,
      data = sava
    ) +
    labs(color = TeX("$\\mu_{\\delta}$"), 
         shape = TeX("$\\mu_{\\delta}$"))+
    theme(axis.title = element_text(size = axissizett),
          axis.text = element_text(size = axissizetext)) +
    theme(legend.text = element_text(size = legendsize),
          legend.title = element_text(size = legendsize))+
    facet_wrap(~K, labeller = labeller(K = Klabel), ncol = 3) +
    theme(strip.text.x = element_text(size = legendsize)) +
    ylim(0,1) +
    theme(panel.spacing = unit(2, "lines")) +
    labs(x = 'Indices of decision times') +
    ggtitle('Intermediate TSR level') +
    theme(plot.title = element_text(size = titlesize, hjust = 0.5)) -> p2
  
  ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
  
  annotate_figure(
    p,
    top = text_grob(
      TeX("SAVA Performance"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('sava_truncgauss_k.pdf', plot = p, width = 35, height = 12)
}


#### dependent case: Gaussian model #####
# change pi+
{
  sava = read_csv('sava_depend_subgau_gauss_pi.csv')
   lord = read_csv('lordpp_depend_sub_gauss_pi.csv')
    saff = read_csv('saffron_depend_sub_gauss_pi.csv')
    addi = read_csv('addis_depend_sub_gauss_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
  annotate_figure(
    p,
    top = text_grob(
      TeX("Methods Comparison (constant $\\rho = 0.5$)"),
      size = totaltitlesize
    )
  ) -> p
  
  ggsave('total_depend_subGauss_Gauss_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_depend_subgau_gauss_rho.csv')
 lord = read_csv('lordpp_depend_sub_gauss_mu.csv')
    saff = read_csv('saffron_depend_sub_gauss_mu.csv')
    addi = read_csv('addis_depend_sub_gauss_mu.csv')
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$rho, lord$rho, saff$rho, addi$rho)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    deltalabel = c('rho: 0.2', 'rho: 0.3', 'rho: 0.4', 'rho: 0.5')
    names(deltalabel) = c('0.2','0.3','0.4','0.5')
    pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
    
  ggsave('total_depend_subGauss_Gauss_delta.pdf', plot = p, width = 26, height = 13)
}


#### dependent case: Gamma model #####
# change pi+
{
  sava = read_csv('sava_depend_gamma_pi.csv')
  lord = read_csv('lordpp_depend_gamma_pi.csv')
  saff = read_csv('saffron_depend_gamma_pi.csv')
  addi = read_csv('addis_depend_gamma_pi.csv')
    
    total_rows <- nrow(sava)
    block_size <- 1000
    keep_per_block <- 800
    full_blocks <- total_rows %/% block_size
    indices <- unlist(lapply(1:full_blocks, function(b) {
      start <- (b - 1) * block_size + 1
      end   <- (b - 1) * block_size + keep_per_block
      start:end
    }))
    selectindex = indices
    l = length(selectindex)
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    pi = c(sava$pi, lord$pi, saff$pi, addi$pi)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totalpi = tibble( time = time, FSR = fsr, TSR = tsr, pi = as.factor(pi), method = methodvec)
    
    pilabel = c('pi: 0.2', 'pi: 0.4', 'pi: 0.6', 'pi: 0.8')
    names(pilabel) = c('0.2','0.4','0.6','0.8')
    pointindex = seq(80, length(totalpi$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totalpi[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totalpi)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~pi, labeller = labeller(pi = pilabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\rho = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
   
  ggsave('total_depend_gamma_pi.pdf', plot = p, width = 26, height = 13)
}
# change pdelta
{
  sava = read_csv('sava_depend_gamma_multiply.csv')
  lord = read_csv('lordpp_depend_gamma_mu.csv')
    saff = read_csv('saffron_depend_gamma_mu.csv')
    addi = read_csv('addis_depend_gamma_mu.csv')
    
    sava = sava[selectindex,]
    lord = lord[selectindex,]
    saff = saff[selectindex,]
    addi = addi[selectindex,]
    
    fsr = c(sava$FSR,  lord$FSR, saff$FSR, addi$FSR)
    delta = c(sava$multiply, lord$rho, saff$rho, addi$rho)
    tsr = c(sava$TSR, lord$TSR, saff$TSR, addi$TSR)
    method = c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS')
    
    time = c(sava$time, lord$time, saff$time, addi$time)
    methodvec = rep(method, each = l)
    totaldelta = tibble( time = time, FSR = fsr, TSR = tsr, delta = as.factor(delta), method = methodvec)
    
    deltalabel = c('rho: 0.2', 'rho: 0.4', 'rho: 0.6', 'rho: 0.8')
    names(deltalabel) = c('0.2','0.4','0.6','0.8')
    
    pointindex = seq(80, length(totaldelta$time), 80)
    
    ggplot()+
      geom_point(aes(x = time, y = FSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = FSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      labs(x = 'Indices of decision times')+
      scale_y_continuous(breaks = seq(0,0.005,0.001), limits = c(0,0.005))+
      ggtitle('Intermediate FSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p1
    
    ggplot()+
      geom_point(aes(x = time, y = TSR, shape = method, color = method), size = pointsize, data = totaldelta[pointindex,])+
      geom_line(aes(x = time, y = TSR, color = method),
                linewidth = lnwidth, data = totaldelta)+
      theme(axis.title = element_text(size = axissizett),
            axis.text = element_text(size = axissizetext))+
      theme(legend.text = element_text(size = legendsize),
            legend.title = element_text(size = legendsize))+
      facet_wrap(~delta, labeller = labeller(delta = deltalabel), ncol = 2)+
      theme(strip.text.x = element_text(size = legendsize))+
      theme(panel.spacing = unit(2, "lines"))+
      ylim(0,1)+
      labs(x = 'Indices of decision times')+
      ggtitle('Intermediate TSR level Comparison')+
      theme(plot.title = element_text(size = titlesize, hjust = 0.5))+
      scale_color_brewer(palette = 'Set1') -> p2
    
    ggarrange(p1,p2, ncol = 2, common.legend = T, legend = 'bottom') -> p
    
    annotate_figure(
      p,
      top = text_grob(
        TeX("Methods Comparison (constant $\\pi^+ = 0.5$)"),
        size = totaltitlesize
      )
    ) -> p
    
    
  ggsave('total_depend_gamma_multiply.pdf', plot = p, width = 26, height = 13)
}

