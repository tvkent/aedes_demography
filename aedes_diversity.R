#library(edgeR)
library(data.table)
library(tidyverse)
library(matrixStats)
library(viridis)
library(boot)
library(cowplot)
library(gridExtra)
library(tanagR)
library(rjson)

setwd('~/aedes/Results/')
pal = tanagr_palette("chlorochrysa_nitidissima", n = 7,discrete = F)


#sweep pi
brazil_0fold_sweep_pi <- read_tsv('Brazil_0fold_sweeps.01_pi.txt') %>% filter(pop=='Brazil')
brazil_4fold_sweep_pi <- read_tsv('Brazil_4fold_sweeps.01_pi.txt') %>% filter(pop=='Brazil')
brazil_04fold_pi <- merge(brazil_0fold_sweep_pi,brazil_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1')) %>% na.omit()

colombia_0fold_sweep_pi <- read_tsv('Colombia_0fold_sweeps.01_pi.txt') %>% filter(pop=='Colombia')%>%mutate(pop = replace(pop, pop == "Colombia", "Rio_Claro"))
colombia_4fold_sweep_pi <- read_tsv('Colombia_4fold_sweeps.01_pi.txt') %>% filter(pop=='Colombia')%>%mutate(pop = replace(pop, pop == "Colombia", "Rio_Claro"))
colombia_04fold_pi <- merge(colombia_0fold_sweep_pi,colombia_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1')) %>% na.omit()

kenya_0fold_sweep_pi <- read_tsv('Kenya_0fold_sweeps.01_pi.txt') %>% filter(pop=='Kenya')
kenya_4fold_sweep_pi <- read_tsv('Kenya_4fold_sweeps.01_pi.txt') %>% filter(pop=='Kenya')
kenya_04fold_pi <- merge(kenya_0fold_sweep_pi,kenya_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1'))%>% na.omit()

senegal_0fold_sweep_pi <- read_tsv('Senegal_0fold_sweeps.01_pi.txt') %>% filter(pop=='Senegal')
senegal_4fold_sweep_pi <- read_tsv('Senegal_4fold_sweeps.01_pi.txt') %>% filter(pop=='Senegal')
senegal_04fold_pi <- merge(senegal_0fold_sweep_pi,senegal_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1')) %>% na.omit()

gabon_0fold_sweep_pi <- read_tsv('Gabon_0fold_sweeps.01_pi.txt') %>% filter(pop=='Gabon')
gabon_4fold_sweep_pi <- read_tsv('Gabon_4fold_sweeps.01_pi.txt') %>% filter(pop=='Gabon')
gabon_04fold_pi <- merge(gabon_0fold_sweep_pi,gabon_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1')) %>% na.omit()

usa_0fold_sweep_pi <- read_tsv('USA_0fold_sweeps.01_pi.txt') %>% filter(pop=='USA')
usa_4fold_sweep_pi <- read_tsv('USA_4fold_sweeps.01_pi.txt') %>% filter(pop=='USA')
usa_04fold_pi <- merge(usa_0fold_sweep_pi,usa_4fold_sweep_pi,by=c('pop','chromosome','window_pos_1')) %>% na.omit()

all_04_sweep_pi <- rbind(brazil_04fold_pi,colombia_04fold_pi,kenya_04fold_pi,senegal_04fold_pi,gabon_04fold_pi,usa_04fold_pi)

all_04_sweep_pi %>% 
  mutate(pinpis=(count_diffs.x/count_comparisons.x)/(count_diffs.y/count_comparisons.y)) %>% 
  na.omit() %>% 
  filter_all(all_vars(!is.infinite(.))) %>%
  group_by(pop) %>% 
  summarise(mean=((sum(count_diffs.x)/sum(count_comparisons.x))/(sum(count_diffs.y)/sum(count_comparisons.y))),se=sd(pinpis)/n()) #%>%
  #spread(pop,mean,se)

boot_mean_04_sweep_pi <- boot(all_04_sweep_pi,boot_04_pi,R=1000)
boot_mean_04_sweep_pi_ci <- t(sapply(1:length(boot_mean_04_sweep_pi$t0), function(x)
  boot.ci(boot_mean_04_sweep_pi,index=x,type='perc')$perc[4:5])) %>%
  data.frame(.)%>%
  data.table::setnames(c('low','high')) %>%
  dplyr::mutate(set=names(boot_mean_04_sweep_pi$t0),estimate=boot_mean_04_sweep_pi$t0)

boot_mean_04_sweep_pi_ci <- boot_mean_04_sweep_pi_ci %>% mutate(region='sweeps') %>%
  mutate(pop=factor(set,levels=c('Brazil','Rio_Claro','Gabon','Kenya','Senegal','USA')))
boot_mean_04_pi_ci <- boot_mean_04_pi_ci %>% mutate(region='genome wide') %>%
  mutate(pop=factor(set,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

boot_all <- rbind(boot_mean_04_pi_ci,boot_mean_04_sweep_pi_ci)

pi04_sweep_plot <- boot_all %>%
  mutate(pop=factor(set,levels=c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya'))) %>%
  ggplot(aes(y=estimate,x=pop,col=pop)) +
  geom_point(size=2,position = position_dodge(width = 0.9),aes(shape=region)) +
  geom_errorbar(aes(ymin=low,ymax=high,linetype=region),width=0.25,linewidth=0.9,position = position_dodge(width = 0.9)) +
  ylab(expression(pi[0]/pi[4]))+
  labs(col='pop')+
  xlab(element_blank())+
  scale_x_discrete(labels=c('Brazil','Cali (Colombia)','Río Claro (Colombia)','USA','Senegal','Gabon','Kenya'))+
  theme_classic() +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=20))+
  scale_colour_manual(values = pal,guide='none')+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=10),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

#pi

lg1_4fold_pi <- read_tsv('AaegL5_1_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz_500kb_pi.txt') %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia') %>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))


lg1_4fold_pi_plot <- lg1_4fold_pi %>% group_by('pop') %>%ggplot(aes(x=midpoint/1000000,y=avg_pi,colour=pop)) +
  facet_grid(.~chromosome,space = 'free_x',scales='free') +
  geom_smooth(method='loess',span=0.1,se = F) +
  ylim(c(0,0.05))+
  scale_x_continuous()+
  theme_classic() +
  scale_color_manual(values = pal)+
  theme(axis.text.x=element_text(angle=30,hjust=1),axis.title.y=element_text(size=14))+
  xlab('Window midpoint (Mb)')+
  ylab(expression(paste('Average 4-fold',pi)))

lg2_4fold_pi <- read_tsv('AaegL5_2_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

lg2_4fold_pi_plot <- lg2_4fold_pi %>% group_by('pop') %>%ggplot(aes(x=midpoint/1000000,y=avg_pi,colour=pop)) +
  facet_grid(.~chromosome,space = 'free_x',scales='free') +
  geom_smooth(method='loess',span=0.1,se = F) +
  ylim(c(0,0.05))+
  scale_x_continuous()+
  theme_classic() +
  scale_color_manual(values = pal)+
  theme(axis.text.x=element_text(angle=30,hjust=1),axis.title.y=element_text(size=14))+
  xlab('Window midpoint (Mb)')+
  ylab(expression(paste('Average 4-fold',pi)))

lg3_4fold_pi <- read_tsv('AaegL5_3_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

lg3_4fold_pi_plot <- lg3_4fold_pi %>% group_by('pop') %>%ggplot(aes(x=midpoint/1000000,y=avg_pi,colour=pop)) +
  facet_grid(.~chromosome,space = 'free_x',scales='free') +
  geom_smooth(method='loess',span=0.1,se = F) +
  ylim(c(0,0.05))+
  scale_x_continuous()+
  theme_classic() +
  scale_color_manual(values = pal)+
  theme(axis.text.x=element_text(angle=30,hjust=1),axis.title.y=element_text(size=14))+
  xlab('Window midpoint (Mb)')+
  ylab(expression(paste('Average 4-fold',pi)))

plot_grid(lg1_4fold_pi_plot,lg2_4fold_pi_plot,lg3_4fold_pi_plot,nrow = 3)

all_4fold_pi <- rbind(lg1_4fold_pi,lg2_4fold_pi,lg3_4fold_pi)

#bootstrap total pi levels
boot_pi <- function(d,i){
  mean <- d[i,] %>% group_by(pop) %>% summarise(mean_pi=(sum(count_diffs)/sum(count_comparisons))) %>% spread(pop,mean_pi) %>% unlist()
  return(mean)
}
#lowexpcutoff
boot_mean_pi <- boot(all_4fold_pi,boot_pi,R=1000)
boot_mean_pi_ci <- t(sapply(1:length(boot_mean_pi$t0), function(x)
  boot.ci(boot_mean_pi,index=x,type='perc')$perc[4:5])) %>%
  data.frame(.)%>%
  data.table::setnames(c('low','high')) %>%
  dplyr::mutate(set=names(boot_mean_pi$t0),estimate=boot_mean_pi$t0)

pi4_plot <- boot_mean_pi_ci %>% 
  mutate(pop=factor(set,levels=c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya'))) %>%
  ggplot(aes(y=estimate,x=pop,col=pop)) +
  geom_point(size=2,position = position_dodge2(width=0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high),width=0.25,linewidth=.8) +
  ylab(expression(paste("Average ",pi[4])))+
  labs(col='population')+
  xlab(element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = 15,linetype=0)))+
  theme_classic() +
  scale_x_discrete(labels=c('Brazil','Cali (Colombia)','Río Claro (Colombia)','USA','Senegal','Gabon','Kenya'))+
  theme(axis.text.x=element_text(angle=30,hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=20),legend.position = 'none')+
  scale_colour_manual(values = pal)+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))


##############
# 0-fold vs 4-fold
##############
lg1_0fold_pi <- read_tsv('AaegL5_1_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

lg2_0fold_pi <- read_tsv('AaegL5_2_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Río_Claro','Gabon','Kenya','Senegal','USA')))

lg3_0fold_pi <- read_tsv('AaegL5_3_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

all_0fold_pi <- rbind(lg1_0fold_pi,lg2_0fold_pi,lg3_0fold_pi)

all_04fold_pi <- merge(all_0fold_pi,all_4fold_pi,by=c('pop','chr','midpoint'))

all_04fold_pi %>% group_by(pop) %>% summarise(mean_pi=sum(count_diffs.y)/sum(count_comparisons.y)) %>% spread(pop,mean_pi) %>% unlist()

#bootstrap total pi levels
boot_04_pi <- function(d,i){
  mean <- d[i,] %>% group_by(pop) %>% summarise(mean_pi=((sum(count_diffs.x)/sum(count_comparisons.x))/(sum(count_diffs.y)/sum(count_comparisons.y)))) %>% spread(pop,mean_pi) %>% unlist()
  return(mean)
}
#lowexpcutoff
boot_mean_04_pi <- boot(all_04fold_pi,boot_04_pi,R=10000)
boot_mean_04_pi_ci <- t(sapply(1:length(boot_mean_04_pi$t0), function(x)
  boot.ci(boot_mean_04_pi,index=x,type='perc')$perc[4:5])) %>%
  data.frame(.)%>%
  data.table::setnames(c('low','high')) %>%
  dplyr::mutate(set=names(boot_mean_04_pi$t0),estimate=boot_mean_04_pi$t0)
pal=c("#D4940E","#976407","#36260A","#F1ED7F","#286F29","#062C3B","#00504F")
pi04_plot <- boot_mean_04_pi_ci %>% 
  mutate(pop=factor(set,levels=c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya')))%>% 
  ggplot(aes(y=estimate,x=pop,col=pop)) +
  geom_point(size=2,position = position_dodge2(width=0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high),width=0.25,linewidth=0.8) +
  ylab(expression(pi[0]/pi[4]))+
  labs(col='population')+
  xlab(element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = 15,linetype=0)))+
  theme_classic() +
  scale_x_discrete(labels=c('Brazil','Cali (Colombia)','Río Claro (Colombia)','USA','Senegal','Gabon','Kenya'))+
  theme(axis.text.x=element_text(angle=30,hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=20),legend.position = 'none')+
  scale_colour_manual(values = pal)+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

##############################
# FIGURE 5
##############################
top_row <- plot_grid(pi4_plot,pi04_plot,labels = c('A','B'),nrow=1)
bottom_row <- plot_grid(NULL,pi04_sweep_plot,NULL,labels=c('','C',''),nrow=1,rel_widths = c(0.3,1,0.2))
plot_grid(top_row,bottom_row,nrow=2)

##############
# 0-fold vs intergenic
##############
lg1_intergenic_pi <- read_tsv('AaegL5_1_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

lg2_intergenic_pi <- read_tsv('AaegL5_2_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

lg3_intergenic_pi <- read_tsv('AaegL5_3_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz_500kb_pi.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>%
  filter(pop!='Colombia')%>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

all_intergenic_pi <- rbind(lg1_intergenic_pi,lg2_intergenic_pi,lg3_intergenic_pi)

all_intergenic0fold_pi <- merge(all_0fold_pi,all_intergenic_pi,by=c('pop','chr','midpoint'))

#bootstrap total pi levels
boot_intergenic_pi <- function(d,i){
  mean <- d[i,] %>% group_by(pop) %>% summarise(mean_pi=((sum(count_diffs.x)/sum(count_comparisons.x))/(sum(count_diffs.y)/sum(count_comparisons.y)))) %>% spread(pop,mean_pi) %>% unlist()
  return(mean)
}
#lowexpcutoff
boot_mean_intergenic_pi <- boot(all_intergenic0fold_pi,boot_intergenic_pi,R=1000)
boot_mean_intergenic_pi_ci <- t(sapply(1:length(boot_mean_intergenic_pi$t0), function(x)
  boot.ci(boot_mean_intergenic_pi,index=x,type='perc')$perc[4:5])) %>%
  data.frame(.)%>%
  data.table::setnames(c('low','high')) %>%
  dplyr::mutate(set=names(boot_mean_intergenic_pi$t0),estimate=boot_mean_intergenic_pi$t0)

piintergenic4_plot <- boot_mean_intergenic_pi_ci %>% 
  mutate(pop=factor(set,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))%>% 
  ggplot(aes(y=estimate,x=pop,col=pop)) +
  geom_point(size=5,position = position_dodge2(width=0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high),width=0.25,linewidth=1) +
  ylab(expression(pi[0]/pi[intergenic]))+
  labs(col='population')+
  xlab(element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = 15,linetype=0)))+
  theme_classic() +
  scale_colour_manual(values = pal)


###########################
# FST
###########################

lg3_0fold_fst <- read_tsv('AaegL5_3_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz_500kb_fst.txt')  %>% 
  mutate(midpoint=window_pos_2-250000,chr=chromosome) %>% unite(pops,c(pop1,pop2),sep='-',remove=F)

lg3_0fold_fst_plot <- lg3_0fold_fst %>% 
  subset(pops=="Cali-Brazil" | pops=="Cali-USA" | pops=='Rio_Claro-Cali' | pops=='Rio_Claro-Brazil' | pops=='Rio_Claro-USA' | pops=='Brazil-USA')%>%group_by('pops') %>%
  ggplot(aes(x=midpoint/1000000,y=avg_wc_fst,colour=pops)) +
  geom_smooth(method='loess',span=0.2,se = F) +
  scale_x_continuous()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=30,hjust=1),axis.title.y=element_text(size=14))+
  xlab('Window midpoint (Mb)')+
  ylab(expression('Average 0-fold F'['ST']))


###########################
# DXY
###########################
setwd('~/aedes/Results/pixy/')
lg1_dxy <- read_tsv('AaegL5_1_invariant_variant_rep_map_depth_masked_fix.vcf.gz_100kb_dxy.txt')  %>% 
  mutate(midpoint=window_pos_2-50000,chr=chromosome) %>% unite(pops,c(pop1,pop2),sep='-',remove=F)
lg2_dxy <- read_tsv('AaegL5_2_invariant_variant_rep_map_depth_masked_fix.vcf.gz_100kb_dxy.txt')  %>% 
  mutate(midpoint=window_pos_2-50000,chr=chromosome) %>% unite(pops,c(pop1,pop2),sep='-',remove=F)
lg3_dxy <- read_tsv('AaegL5_3_invariant_variant_rep_map_depth_masked_fix.vcf.gz_100kb_dxy.txt')  %>% 
  mutate(midpoint=window_pos_2-50000,chr=chromosome) %>% unite(pops,c(pop1,pop2),sep='-',remove=F)

full_dxy <- rbind(lg1_dxy,lg2_dxy,lg3_dxy) %>% na.omit()
full_dxy_genomewide <- full_dxy %>% group_by(pops) %>% na.omit() %>% summarise(mean_dxy=sum(count_diffs)/sum(count_comparisons))


#bootstrap total pi levels
boot_dxy <- function(d,i){
  mean <- d[i,] %>% group_by(pops) %>% summarise(mean_dxy=sum(count_diffs)/sum(count_comparisons)) %>% spread(pops,mean_dxy) %>% unlist()
  return(mean)
}
#lowexpcutoff
boot_mean_dxy <- boot(full_dxy,boot_dxy,R=100)
boot_mean_dxy_ci <- t(sapply(1:length(boot_mean_dxy$t0), function(x)
  boot.ci(boot_mean_dxy,index=x,type='perc')$perc[4:5])) %>%
  data.frame(.)%>%
  data.table::setnames(c('low','high')) %>%
  dplyr::mutate(set=names(boot_mean_dxy$t0),estimate=boot_mean_dxy$t0)

boot_mean_dxy_ci <- boot_mean_dxy_ci %>% separate(set, c("pop1", "pop2"), "-", remove=F) %>%
  mutate(region1=if_else(str_detect(pop1, 'Brazil|Cali|Rio_Claro|Medellin'),"SA",
                         ifelse(str_detect(pop1, 'Senegal|Kenya|Gabon'), "Africa",
                         ifelse(str_detect(pop1, 'SoCal|NorCal|Exeter|Clovis|Florida'), "NA","none")))) %>%
  mutate(region2=if_else(str_detect(pop2, 'Brazil|Cali|Rio_Claro|Medellin'),"SA",
                         ifelse(str_detect(pop2, 'Senegal|Kenya|Gabon'), "Africa",
                                ifelse(str_detect(pop2, 'SoCal|NorCal|Exeter|Clovis|Florida'), "NA","none")))) %>%
  mutate(contrast=if_else((region1=='SA' & region2=='NA') | (region2=='SA' & region1=="NA"), "SA-NA",
                          if_else((region1=='SA' & region2=='Africa') | (region2=='SA' & region1=="Africa"),'SA-Africa',
                                  if_else((region1=='NA' & region2=='Africa') | (region2=='NA' & region1=="Africa"),'NA-Africa',
                                          if_else((region1=='Africa' & region2=='Africa') | (region2=='Africa' & region1=="Africa"),'Africa',
                                                  if_else((region1=='SA' & region2=='SA') | (region2=='SA' & region1=="SA"),'SA',
                                                          if_else((region1=='NA' & region2=='NA') | (region2=='NA' & region1=="NA"),'NA','none')))))))

boot_mean_dxy_ci$contrast <- factor(boot_mean_dxy_ci$set, levels=unique(boot_mean_dxy_ci$set))
boot_mean_dxy_ci <- transform(boot_mean_dxy_ci, set=reorder(set, contrast) ) 

boot_mean_dxy_ci$set <- boot_mean_dxy_ci$set %>% str_replace("Rio_Claro","Río Claro")

boot_mean_dxy_ci %>% 
  group_by(contrast) %>%
  arrange(desc(estimate),.by_group = TRUE)%>%
  #mutate(set=fct_reorder(set,contrast)) %>% 
  filter(!str_detect(set, 'Medellin'))%>%
  ggplot(aes(y=estimate,x=fct_inorder(set),col=contrast, label=contrast)) +
  facet_grid( .~factor(contrast,levels=c("Africa","SA-Africa","NA-Africa","SA-NA","NA","SA")), space='free_x', switch = 'x', scales = "free")+
  geom_point(size=1,position = position_dodge2(width=0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high),width=0.25,linewidth=1) +
  geom_vline(xintercept = c(1,3,5,7,9,11,13,15,17,19,21),col='#D9D9D9')+
  geom_point(size=1,position = position_dodge2(width=0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high),width=0.25,linewidth=1) +
  labs(col='population pair')+
  xlab(element_blank())+
  ylab(expression("d"["xy"]))+
  #ylab(expression("F"["ST"]))+
  guides(colour = guide_legend(override.aes = list(shape = 15,linetype=0)))+
  theme_minimal() +
  theme(axis.text.x=element_text(angle=40,hjust=1,size=10),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=16),legend.position = 'none',
        strip.placement = "outside") 


lg2_dxy %>%
  #filter(str_detect(pops, 'SoCal')) %>%
  #filter(str_detect(pops, 'Kenya|Gabon|Senegal')) %>%
  filter(str_detect(pops,'NorCal|Florida|Clovis|SoCal|Exeter')) %>%
  filter(str_detect(pops,'Rio_Claro|Cali|Brazil')) %>%
  group_by('pops') %>%ggplot(aes(x=midpoint/1000000,y=avg_dxy,colour=pops)) +
  geom_smooth(method='loess',span=0.1,se = F) +
  scale_x_continuous()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=30,hjust=1),axis.title.y=element_text(size=14))+
  xlab('Window midpoint (Mb)')+
  ylab(expression('Average dxy'))

###########################
# SFS
###########################

brazil_sfs<- transpose(fread('~/aedes/brazil_new_sfs.txt')) %>%
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Brazil',rel_af=(rowid/max(rowid))/2)

cali_sfs <- transpose(fread('~/aedes/cali_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Cali',rel_af=(rowid/max(rowid))/2)

clovis_sfs <- transpose(fread('~/aedes/clovis_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Clovis',rel_af=(rowid/max(rowid))/2)

gabon_sfs <- transpose(fread('~/aedes/gabon_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Gabon',rel_af=(rowid/max(rowid))/2)

kenya_sfs <- transpose(fread('~/aedes/kenya_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Kenya',rel_af=(rowid/max(rowid))/2)

rio_claro_sfs<- transpose(fread('~/aedes/rio_claro_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Rio_Claro',rel_af=(rowid/max(rowid))/2)

senegal_sfs <- transpose(fread('~/aedes/senegal_new_sfs.txt')) %>% 
  tibble::rowid_to_column() %>% 
  mutate(rel_freq=V1/sum(V1),pop='Senegal',rel_af=(rowid/max(rowid))/2)

full_sfs <- rbind(brazil_sfs,cali_sfs,clovis_sfs,gabon_sfs,kenya_sfs,rio_claro_sfs,senegal_sfs) %>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','Clovis')))

full_sfs %>% ggplot(aes(x=rel_af,y=rel_freq,col=pop))+
  geom_line() +
  ylab('Relative count')+
  labs(col='population')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))+
  scale_colour_manual(values = pal)


####################################
# DFE
####################################

#4fold
############
brazil_4fold_nes <- fread('Brazil_4fold.nes') %>% mutate(pop='Brazil')
cali_4fold_nes <- fread('Cali_4fold.nes') %>% mutate(pop='Cali')
Rio_Claro_4fold_nes <- fread('Rio_Claro_4fold.nes') %>% mutate(pop='Rio_Claro')
Gabon_4fold_nes <- fread('Gabon_4fold.nes') %>% mutate(pop='Gabon')
Kenya_4fold_nes <- fread('Kenya_4fold.nes') %>% mutate(pop='Kenya')
Senegal_4fold_nes <- fread('Senegal_4fold.nes') %>% mutate(pop='Senegal')
USA_4fold_nes <- fread('USA_4fold.nes') %>% mutate(pop='USA')
nes_4fold <- rbind(brazil_4fold_nes,cali_4fold_nes,Rio_Claro_4fold_nes,Gabon_4fold_nes,Kenya_4fold_nes,Senegal_4fold_nes,USA_4fold_nes) %>%
  group_by(pop) %>%
  mutate(Strong=sum(`10-100`,`100-inf`)) %>%
  ungroup() %>%
  melt() %>% 
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA'))) %>%
  group_by(pop,variable) %>%
  rename(NEs=variable,Prop_muts=value)


brazil_4fold_nes_boot <- fread('Brazil_4fold_bootstrapped1000.nes') %>% mutate(pop='Brazil')
cali_4fold_nes_boot <- fread('Cali_4fold_bootstrapped1000.nes') %>% mutate(pop='Cali')
Rio_Claro_4fold_nes_boot <- fread('Rio_Claro_4fold_bootstrapped1000.nes') %>% mutate(pop='Rio_Claro')
Gabon_4fold_nes_boot <- fread('Gabon_4fold_bootstrapped1000.nes') %>% mutate(pop='Gabon')
Kenya_4fold_nes_boot <- fread('Kenya_4fold_bootstrapped1000.nes') %>% mutate(pop='Kenya')
Senegal_4fold_nes_boot <- fread('Senegal_4fold_bootstrapped1000.nes') %>% mutate(pop='Senegal')
USA_4fold_nes_boot <- fread('USA_4fold_bootstrapped1000.nes') %>% mutate(pop='USA')
fullbootnes_4fold <- rbind(brazil_4fold_nes_boot,
                          cali_4fold_nes_boot,
                          Rio_Claro_4fold_nes_boot,
                          Gabon_4fold_nes_boot,
                          Kenya_4fold_nes_boot,
                          Senegal_4fold_nes_boot,
                          USA_4fold_nes_boot)
fullbootnes_4fold$Strong<-rowSums(fullbootnes_4fold [,c("10-100","100-inf")])
fullbootnes_4fold_melt <- fullbootnes_4fold %>% 
  melt() %>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

#95% percentile CI for bootstrapped dfe
fullbootnes_CI <- fullbootnes_4fold_melt %>%
  group_by(pop,variable) %>% 
  summarise(CIlow = quantile(value,probs = c(0.025)),CIhigh=quantile(value,probs=c(0.975))) %>% 
  rename(NEs=variable)  %>%filter(NEs != "10-100") %>% filter(NEs != "100-inf") %>% #filter(NEs != "Strong")
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya'))) #%>% filter(NEs != "10-100") %>% filter(NEs != "100-inf")
fullbootnes_CI$pop <- factor(fullbootnes_CI$pop, levels = c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya'))
                        
pal=c("#D4940E","#976407","#36260A","#F1ED7F","#286F29","#062C3B","#00504F")

nes_4fold %>% filter(NEs != "10-100") %>% filter(NEs != "100-inf") %>% #filter(NEs != "Strong")
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','USA','Senegal','Gabon','Kenya'))) %>% 
  ggplot(aes(y=Prop_muts,x=NEs,fill=pop)) +
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(data=fullbootnes_CI,aes(ymin=CIlow,ymax=CIhigh,x=NEs,fill=pop),inherit.aes=F,position = position_dodge2(),col='black')+
  scale_fill_manual(values = pal,labels=c('Brazil','Cali (Colombia)','Río Claro (Colombia)','USA','Senegal','Gabon','Kenya'),name="Region") +
  theme_classic() +
  ylab("Proportion of mutations") +
  xlab(expression(italic("N"["e"]*"s"))) +
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=11),axis.title.y = element_text(size=11),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))


#intergenic
############
brazil_intergenic_nes <- fread('Brazil_intergenic.nes') %>% mutate(pop='Brazil')
cali_intergenic_nes <- fread('Cali_intergenic.nes') %>% mutate(pop='Cali')
Rio_Claro_intergenic_nes <- fread('Rio_Claro_intergenic.nes') %>% mutate(pop='Rio_Claro')
Gabon_intergenic_nes <- fread('Gabon_intergenic.nes') %>% mutate(pop='Gabon')
Kenya_intergenic_nes <- fread('Kenya_intergenic.nes') %>% mutate(pop='Kenya')
Senegal_intergenic_nes <- fread('Senegal_intergenic.nes') %>% mutate(pop='Senegal')
USA_intergenic_nes <- fread('USA_intergenic.nes') %>% mutate(pop='USA')
nes_intergenic <- rbind(brazil_intergenic_nes,cali_intergenic_nes,Rio_Claro_intergenic_nes,Gabon_intergenic_nes,Kenya_intergenic_nes,Senegal_intergenic_nes,USA_intergenic_nes) %>%
  melt() %>% 
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA'))) %>%
  group_by(pop,variable) %>%
  rename(NEs=variable,Prop_muts=value)


brazil_intergenic_nes_boot <- fread('Brazil_intergenic_bootstrapped1000.nes') %>% mutate(pop='Brazil')
cali_intergenic_nes_boot <- fread('Cali_intergenic_bootstrapped1000.nes') %>% mutate(pop='Cali')
Rio_Claro_intergenic_nes_boot <- fread('Rio_Claro_intergenic_bootstrapped1000.nes') %>% mutate(pop='Rio_Claro')
Gabon_intergenic_nes_boot <- fread('Gabon_intergenic_bootstrapped1000.nes') %>% mutate(pop='Gabon')
Kenya_intergenic_nes_boot <- fread('Kenya_intergenic_bootstrapped1000.nes') %>% mutate(pop='Kenya')
Senegal_intergenic_nes_boot <- fread('Senegal_intergenic_bootstrapped1000.nes') %>% mutate(pop='Senegal')
USA_intergenic_nes_boot <- fread('USA_intergenic_bootstrapped1000.nes') %>% mutate(pop='USA')
fullbootnes_intergenic <- rbind(brazil_intergenic_nes_boot,
                           cali_intergenic_nes_boot,
                           Rio_Claro_intergenic_nes_boot,
                           Gabon_intergenic_nes_boot,
                           Kenya_intergenic_nes_boot,
                           Senegal_intergenic_nes_boot,
                           USA_intergenic_nes_boot)
fullbootnes_intergenic_melt <- fullbootnes_intergenic %>% 
  melt() %>%
  mutate(pop=factor(pop,levels=c('Brazil','Cali','Rio_Claro','Gabon','Kenya','Senegal','USA')))

#95% percentile CI for bootstrapped dfe
fullbootnes_CI <- fullbootnes_intergenic_melt %>%
  group_by(pop,variable) %>% 
  summarise(CIlow = quantile(value,probs = c(0.025)),CIhigh=quantile(value,probs=c(0.975))) %>% 
  rename(NEs=variable)
nes_intergenic %>% ggplot(aes(y=Prop_muts,x=NEs,fill=pop)) +
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(data=fullbootnes_CI,aes(ymin=CIlow,ymax=CIhigh,x=NEs),inherit.aes=F,position = position_dodge2(),col='black')+
  scale_fill_manual(values = pal) +
  theme_classic()
