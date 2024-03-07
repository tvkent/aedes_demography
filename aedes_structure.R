library(tidyverse)
library(cowplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(pophelper)
library(tanagR)
library(scales)
library(rjson)
library(calecopal)
library(MetBrewer)
names(tanagr_palettes)
setwd('~/aedes/Results/')


#pal <- lacroix_palette(type = 'paired',n = 9)
pal = tanagr_palette("chlorochrysa_nitidissima", n = 20,discrete = F)
pal_americas_countries=tanagr_palette("chlorochrysa_nitidissima")[c(1,2,3,6)]
pal_americas_countries=c("#D4940E","#976407","#36260A","#F1ED7F")
pal_americas_locations=pal[c(1:4,8:20)]


#-------------------------------------#
# MAP OF SAMPLES
#-------------------------------------#

library(ggspatial)
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

sites_rc <- data.frame(longitude = c(-74.86), latitude = c(5.9))
sites_rc <- st_as_sf(sites_rc, coords = c("longitude", "latitude"), crs = 4326, 
                     agr = "constant")

sites_cali <- data.frame(longitude = c(-76.55), latitude = c(3.38))
sites_cali <- st_as_sf(sites_cali, coords = c("longitude", "latitude"), crs = 4326, 
                       agr = "constant")

sites_brazil <- data.frame(longitude = c(-54.72), latitude = c(2.44))
sites_brazil <- st_as_sf(sites_brazil, coords = c("longitude", "latitude"), crs = 4326, 
                         agr = "constant")

sites_gabon <- data.frame(longitude = c(13.58), latitude = c(-1.64))
sites_gabon <- st_as_sf(sites_gabon, coords = c("longitude", "latitude"), crs = 4326, 
                        agr = "constant")

sites_kenya <- data.frame(longitude = c(39.6), latitude = c(-3.93))
sites_kenya <- st_as_sf(sites_kenya, coords = c("longitude", "latitude"), crs = 4326, 
                        agr = "constant")

sites_senegal <- data.frame(longitude = c(-16.43), latitude = c(14.64))
sites_senegal <- st_as_sf(sites_senegal, coords = c("longitude", "latitude"), crs = 4326, 
                          agr = "constant")

sites_usa <- data.frame(longitude = c(-122.2,-120.11,-119.86,-119.67,-119.51,-119.16,-117.92,-117.67,-115.53,-117.05,-80.38,-81.78), latitude = c(37.44,36.96,36.83,36.81,36.73,36.3,33.76,33.61,32.97,32.56,29.65,24.56))
sites_usa <- st_as_sf(sites_usa, coords = c("longitude", "latitude"), crs = 4326, 
                      agr = "constant")

map_fig<- (gworld <- ggplot(data = world) +
             geom_sf(fill = "lightgrey") +
             geom_sf(data = sites_rc, size = 3, shape = 23, fill = "#36260A") +
             geom_sf(data = sites_cali, size = 3, shape = 23, fill = "#976407") +
             geom_sf(data = sites_brazil, size = 3, shape = 23, fill = "#D4940E") +
             geom_sf(data = sites_kenya, size = 3, shape = 23, fill = "#00504F") +
             geom_sf(data = sites_gabon, size = 3, shape = 23, fill = "#062C3B") +
             geom_sf(data = sites_senegal, size = 3, shape = 23, fill = "#286F29") +
             geom_sf(data = sites_usa, size = 3, shape = 23, fill = "#F1ED7F") +
             coord_sf(xlim = c(-128, 50), ylim = c(-10, 47), expand = FALSE) +
             theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), 
                   panel.border = element_rect(fill = NA)))

#all samples
pca <- fread("./admixture/all/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.eigenvec",header=F,data.table=F,stringsAsFactors = F)
pca_named<-pca[,-1]
names(pca_named)<-c("sample",paste0("PC",seq(1:20)))

eigenvals <- fread("./admixture/all/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.eigenval",header=F,data.table=F,stringsAsFactors = F)
eigenvals$frac<-round( (eigenvals$V1 / sum(eigenvals$V1)*100), digits=2)

samples <- fread('../whole_sample_sorted_country.031522.csv',data.table=F,stringsAsFactors = FALSE) %>%
  rename(sample=sample_id) %>% mutate()
pca_named <- merge(samples,pca_named,by='sample')

#scree plot by country
eigenvals %>% mutate(pc=row_number()) %>% ggplot(aes(x=pc,y=V1)) +
  geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+
  ylab('Percent variation explained') +
  xlab('Principal component')
  ggtitle("scree plot for country-wise PCA")

#plot by country
country_pca_fig <- pca_named %>%
  ggplot(aes(PC1,PC2, color=country)) +
  geom_point(size=2,alpha=0.6) +
  scale_color_tanagr("chlorochrysa_nitidissima")+
  xlab("PC1 (12.35%)") +
  ylab("PC2 (7.25%)")+
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

# plot by location
pca_named %>%
  ggplot(aes(PC1,PC2, color=location)) +
  geom_point(size=3,alpha=0.8) +
  scale_color_manual(values = pal)+
  theme_cowplot()

#grid 
pc12 <- pca_named %>%
  ggplot(aes(PC1,PC2, color=country)) +
  geom_point(size=3,alpha=0.8) +
  xlab(paste0("PC1 (", eigenvals$frac[1], "%)")) +
  ylab(paste0("PC2 (", eigenvals$frac[2], "%)")) +
  scale_color_tanagr("chlorochrysa_nitidissima")+
  theme_cowplot()

pc23 <- pca_named %>%
  ggplot(aes(PC2,PC3, color=country)) +
  geom_point(size=3,alpha=0.8) +
  xlab(paste0("PC2 (", eigenvals$frac[2], "%)")) +
  ylab(paste0("PC3 (", eigenvals$frac[3], "%)")) +
  scale_color_tanagr("chlorochrysa_nitidissima")+
  theme_cowplot()

pc34 <- pca_named %>%
  ggplot(aes(PC3,PC4, color=country)) +
  geom_point(size=3,alpha=0.8) +
  xlab(paste0("PC3 (", eigenvals$frac[3], "%)")) +
  ylab(paste0("PC4 (", eigenvals$frac[4], "%)")) +
  scale_color_tanagr("chlorochrysa_nitidissima")+
  theme_cowplot()

pc45 <- pca_named %>%
  ggplot(aes(PC4,PC5, color=country)) +
  geom_point(size=3,alpha=0.6) +
  xlab(paste0("PC4 (", eigenvals$frac[4], "%)")) +
  ylab(paste0("PC5 (", eigenvals$frac[5], "%)")) +
  scale_color_tanagr("chlorochrysa_nitidissima")+
  theme_cowplot()

ggarrange(pc12,pc23,pc34,pc45,common.legend = T)


#americas only
pca_americas <- fread("./admixture/americas/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas.eigenvec",header=F,data.table=F,stringsAsFactors = F)
pca_americas_named<-pca_americas[,-1]
names(pca_americas_named)<-c("sample",paste0("PC",seq(1:20)))

eigenvals_americas <- fread("./admixture/americas/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas.eigenval",header=F,data.table=F,stringsAsFactors = F)
eigenvals_americas$frac<-round( (eigenvals_americas$V1 / sum(eigenvals_americas$V1)*100), digits=2)

pca_americas_named <- merge(samples,pca_americas_named,by='sample')

# scree plot americas
eigenvals_americas %>% mutate(pc=row_number()) %>% ggplot(aes(x=pc,y=V1)) +
  geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+
  ylab('Percent variation explained') +
  xlab('Principal component')
  ggtitle("scree plot for Americas PCA")

#plot by country
pca_americas_named %>%
  ggplot(aes(PC1,PC2, color=country)) +
  geom_point(size=3,alpha=0.8) +
  scale_color_manual(values=pal_americas_countries)+
  geom_text_repel(aes(label=location_label))+
  xlab("PC1 (10.86%)") +
  ylab("PC2 (7.72%)") +
  theme_cowplot()

# plot by location
#labels for pc12 plot
label_ix <- c(1,19,29,36,54:71,75,76,79)
pca_americas_named$location_label <- ""
pca_americas_named$location_label[label_ix] <- pca_americas_named$location[label_ix]
#subsetted labels for smaller plots
#label_ix <- c(1,19,29,36,54,55,58,59,60,61,62,63,65,67,68,71,75,76)
label_ix <- c(1,19,36,54,58,59,60,61,63,65,68,75,76)
pca_americas_named$location_label_reduced <- ""
pca_americas_named$location_label_reduced[label_ix] <- pca_americas_named$location[label_ix]

#labels for country level plot
pca_americas_named$location_label_coarse <- ""
pca_americas_named$location_label_rc <- ""
pca_americas_named$location_label_rc[1] <- "Colombia:
Río Claro"
pca_americas_named$location_label_coarse[19] <- "Colombia: Cali"
pca_americas_named$location_label_coarse[64] <- "USA: SoCal"
pca_americas_named$location_label_norcal <- ""
pca_americas_named$location_label_norcal[54] <- "USA: NorCal"
pca_americas_named$location_label_brazil <- ""
pca_americas_named$location_label_brazil[39] <- "Brazil"

#coarse countries for colour
pca_americas_named$country_label_coarse <- pca_americas_named$location
pca_americas_named$country_label_coarse[1:18] <- "Río Claro (CO)"
pca_americas_named$country_label_coarse[19:28] <- "Cali (CO)"
pca_americas_named$country_label_coarse[29:34] <- "Río Claro (CO)"
pca_americas_named$country_label_coarse[53:79] <- "USA"
pca_americas_named$country_label_coarse[35:52] <- "Brazil"

NA_pca_fig <- pca_americas_named %>%
  ggplot(aes(PC1,PC2, color=country_label_coarse)) +
  geom_point(size=2,alpha=0.6) +
  geom_text_repel(aes(label=location_label_coarse),size=2.5,min.segment.length = 10,nudge_x = -.04,show.legend = F)+
  geom_text_repel(aes(label=location_label_brazil),size=2.5,min.segment.length = 10,nudge_y = -0.05,show.legend = F)+
  geom_text_repel(aes(label=location_label_norcal),size=2.5,min.segment.length = 10,nudge_x = 0.07,nudge_y=0.05,show.legend = F)+
  geom_text_repel(aes(label=location_label_rc),size=2.5,min.segment.length = 10,nudge_x = -0.00,nudge_y=0.1,show.legend = F)+
  scale_color_manual(values=pal_americas_countries,name='country')+
  xlab("PC1 (10.86%)") +
  ylab("PC2 (7.72%)") +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+
  theme(legend.position="none")

##############################
# FIGURE 1
##############################
bottom_row <- plot_grid(country_pca_fig, NA_pca_fig, labels = c('B', 'C'), label_size = 14, ncol=2,rel_widths = c(1,0.8))
plot_grid(map_fig, bottom_row, labels = c('A', ''), label_size = 14, ncol = 1,rel_heights = c(0.6,1))




#grid 
pc12 <- pca_americas_named %>%
  ggplot(aes(PC1,PC2, color=location)) +
  geom_point(size=3,alpha=0.8) +
  geom_text_repel(aes(label=location_label_reduced))+
  scale_color_manual(values=pal_americas_locations)+
  xlab(paste0("PC1 (", eigenvals_americas$frac[1], "%)")) +
  ylab(paste0("PC2 (", eigenvals_americas$frac[2], "%)")) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

pc23 <- pca_americas_named %>%
  ggplot(aes(PC2,PC3, color=location)) +
  geom_point(size=3,alpha=0.8) +
  geom_text_repel(aes(label=location_label_reduced))+
  scale_color_manual(values=pal_americas_locations)+
  xlab(paste0("PC2 (", eigenvals_americas$frac[2], "%)")) +
  ylab(paste0("PC3 (", eigenvals_americas$frac[3], "%)")) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

pc34 <- pca_americas_named %>%
  ggplot(aes(PC3,PC4, color=location)) +
  geom_point(size=3,alpha=0.8) +
  geom_text_repel(aes(label=location_label_reduced))+
  scale_color_manual(values=pal_americas_locations)+
  xlab(paste0("PC3 (", eigenvals_americas$frac[3], "%)")) +
  ylab(paste0("PC4 (", eigenvals_americas$frac[4], "%)")) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

pc45 <- pca_americas_named %>%
  ggplot(aes(PC4,PC5, color=location)) +
  geom_point(size=3,alpha=0.8) +
  geom_text_repel(aes(label=location_label_reduced))+
  scale_color_manual(values=pal_americas_locations)+
  xlab(paste0("PC4 (", eigenvals_americas$frac[4], "%)")) +
  ylab(paste0("PC5 (", eigenvals_americas$frac[5], "%)")) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=11),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

ggarrange(pc12,pc23,pc34,pc45,common.legend = F)


#----------------------------------------------------#
# ADMIXTURE
#----------------------------------------------------#

#everything

all_samples_loc <- samples %>% select(location,country) %>% filter(location != "Colombia") %>% 
  mutate(across(everything(), gsub, pattern = ".*:", replacement =""))

pal = tanagr_palette("chlorochrysa_nitidissima", n = 10,discrete = FALSE)
pal = tanagr_palette("dacnis_berlepschi", n=5, discrete = F)

#setwd('~/aedes/Results/admixture/all/plotting/')
#sfiles <- list.files(pattern = "AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked_numericChr_LDpruned_100ld.*Q")
setwd('~/aedes/Results/admixture/all/plotting')
pattern = 'AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.*Q'
sfiles <- list.files(pattern = pattern)
all <- readQ(files=sfiles)
p2 <- plotQ(all,imgoutput="join",
            returnplot=TRUE,
            exportplot=FALSE,
            basesize=11,
            clustercol = pal,
            splab = c('K=2','K=3','K=4','K=5'),
            grplab = all_samples_loc,
            grplabsize=1.9,
            panelratio = c(3,2),
            grplabangle = 60,
            grplabpos = 0.65,
            linesize=0.8,
            pointsize=4,
            selgrp = 'country',
            sortind = 'all',
            ordergrp = TRUE,
            sharedindlab = FALSE
)
admix1 <- p2$plot[[1]]

#americas
setwd('~/aedes/Results/admixture/americas/plotting/')
pal = tanagr_palette("buthraupis_montana",n=4, discrete = F)
americas_samples <- fread('../AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked_numericChr_LDpruned_100ld_americas.fam',data.table = FALSE,header=FALSE,stringsAsFactors=FALSE) %>%
  rename(sample=V1) %>% select(sample)
americas_samples_loc <- merge(americas_samples,samples,by='sample',sort = FALSE) %>% 
  select(location,country) %>%
  mutate(across(everything(), gsub, pattern = ".*:", replacement ="")) %>%
  filter(location != "Colombia")

sfiles <- list.files(pattern = 'AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas.*.Q')
all <- readQ(files=sfiles)
p2 <- plotQ(all,imgoutput="join",
            returnplot=T,
            exportplot=F,
            basesize=11,
            clustercol = pal,
            splab = c('K=2','K=3','K=4'),
            grplab = americas_samples_loc,
            grplabsize=1.9,
            panelratio = c(3,2),
            grplabangle = 60,
            grplabpos = 0.65,
            
            linesize=0.8,
            pointsize=4,
            selgrp = 'country',
            ordergrp = T,
            sortind = 'all',
            sharedindlab = F
            )
admix_2 <- p2$plot[[1]]

plot_grid(NULL,admix1,NULL,admix_2,labels = c('','A','','B'),nrow=2,ncol=2,rel_widths = c(0.03,.97,0.03,.97),hjust = 1)

#----------------------------------------------------#
# STAIRWAYPLOT
#----------------------------------------------------#
setwd('~/aedes/Results/')
pal=tanagr_palette("chlorochrysa_nitidissima")


brazil_stairway <- fread('brazil_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Brazil')
cali_stairway <- fread('cali_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Cali (Colombia)')
clovis_stairway <- fread('clovis_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Clovis')
gabon_stairway <- fread('gabon_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Gabon')
kenya_stairway <- fread('kenya_100k_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Kenya')
rio_claro_stairway <- fread('rio_claro_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Rio Claro (Colombia)')
senegal_stairway <- fread('senegal_masked.final.summary',data.table = F,header=T) %>% mutate(pop='Senegal')

all_stairway <- rbind(brazil_stairway,cali_stairway,rio_claro_stairway,gabon_stairway,kenya_stairway,senegal_stairway)
all_stairway$pop <- factor(all_stairway$pop,levels=c('Brazil','Cali (Colombia)','Rio Claro (Colombia)','Gabon','Kenya','Senegal'))

all_stairway$pop <- all_stairway$pop %>% str_replace("Rio Claro","Río Claro")
all_stairway$pop <- factor(all_stairway$pop,levels=c('Brazil','Cali (Colombia)','Río Claro (Colombia)','Gabon','Kenya','Senegal'))

all_stairway %>% ggplot(aes(y=Ne_median/1000,x=year/1000,col=pop)) + 
  geom_line(size=1.2) + 
  scale_x_log10(breaks = breaks_log(n=15),labels=label_comma(drop0trailing = TRUE),limits=c(0.001,110)) +
  scale_y_log10(breaks = breaks_log(n=15), labels = scales::comma) +
  scale_color_manual(values=c("#D4940E","#976407","#36260A","#062C3B","#00504F","#286F29"),name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  theme(axis.text.x=element_text(angle=40))+
  theme_bw() 

#----------------------------------------------------#
# SMC++
#----------------------------------------------------#
setwd("~/aedes/Results/smc++/")
pal=tanagr_palette("chlorochrysa_nitidissima",n = 7,discrete = F)

brazil_cubic_smc <- fread('Brazil_cubic_smc.csv',data.table = F)
cali_cubic_smc <- fread('Cali_cubic_smc.csv',data.table = F)
usa_cubic_smc <- fread('USA_cubic_smc.csv',data.table = F)
gabon_cubic_smc <- fread('Gabon_cubic_smc.csv',data.table = F)
kenya_cubic_smc <- fread('Kenya_cubic_smc.csv',data.table = F)
rio_claro_cubic_smc <- fread('Rio_Claro_cubic_smc.csv',data.table = F)
senegal_cubic_smc <- fread('Senegal_cubic_smc.csv',data.table = F)

full_cubic <- rbind(brazil_cubic_smc,cali_cubic_smc,usa_cubic_smc,gabon_cubic_smc,kenya_cubic_smc,rio_claro_cubic_smc,senegal_cubic_smc)

full_cubic %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Cubic spline Ne')+
  theme_classic() 

brazil_smc <- fread('Brazil_smc.csv',data.table = F)
cali_smc <- fread('Cali_smc.csv',data.table = F)
usa_smc <- fread('USA_smc.csv',data.table = F)
gabon_smc <- fread('Gabon_smc.csv',data.table = F)
kenya_smc <- fread('Kenya_smc.csv',data.table = F)
rio_claro_smc <- fread('Rio_Claro_smc.csv',data.table = F)
senegal_smc <- fread('Senegal_smc.csv',data.table = F)

full_piecewise <- rbind(brazil_smc,cali_smc,usa_smc,gabon_smc,kenya_smc,rio_claro_smc,senegal_smc)

full_piecewise %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Piecewise Ne 10 knots')+
  theme_classic()

#rp4 low regularization to smooth oscillations
brazil_smc <- fread('Brazil_rp4_smc.csv',data.table = F)
cali_smc <- fread('Cali_rp4_smc.csv',data.table = F)
usa_smc <- fread('USA_rp4_smc.csv',data.table = F)
gabon_smc <- fread('Gabon_rp4_smc.csv',data.table = F)
kenya_smc <- fread('Kenya_rp4_smc.csv',data.table = F)
rio_claro_smc <- fread('Rio_Claro_rp4_smc.csv',data.table = F)
senegal_smc <- fread('Senegal_rp4_smc.csv',data.table = F)

full_rp4 <- rbind(brazil_smc,cali_smc,usa_smc,gabon_smc,kenya_smc,rio_claro_smc,senegal_smc)

full_rp4 %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Piecewise Ne 10 knots, regularization penalty 4')+
  theme_classic()

#rp5 low regularization to smooth oscillations
brazil_smc <- fread('Brazil_100ksites_smc.csv',data.table = F)
cali_smc <- fread('Cali_rp5_smc.csv',data.table = F)
usa_smc <- fread('USA_100ksites_smc.csv',data.table = F)
gabon_smc <- fread('Gabon_rp5_smc.csv',data.table = F)
kenya_smc <- fread('Kenya_rp5_smc.csv',data.table = F)
rio_claro_smc <- fread('Rio_Claro_100k_sites_smc.csv',data.table = F)
senegal_smc <- fread('Senegal_rp5_smc.csv',data.table = F)

full_rp5_africa <- rbind(gabon_smc,kenya_smc,senegal_smc)
full_rp5_americas <- rbind(brazil_smc,cali_smc,rio_claro_smc,usa_smc)
full_rp5 <- rbind(brazil_smc,cali_smc,rio_claro_smc,usa_smc,gabon_smc,kenya_smc,senegal_smc)

full_smc <- full_rp5 %>% mutate(label=ifelse(label=="RC","Río Claro",label)) %>% mutate(label=factor(label,levels=c('Brazil','Cali','Río Claro','USA','Gabon','Kenya','Senegal'))) %>%
  ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=c("#D4940E","#976407","#36260A","#F1ED7F","#00504F","#062C3B","#286F29"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

smc_legend <- get_legend(
  # create some space to the left of the legend
 full_smc + theme(legend.box.margin = margin(0, 0, 0, 12))
)

africas_smc <- full_rp5_africa %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=c("#00504F","#062C3B","#286F29"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+ 
  theme(legend.position="none")

americas_smc <- full_rp5_africa %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=0,col='white',show.legend = F)+
  geom_line(data=full_rp5_americas,size=1) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=c("#D4940E","#976407","#36260A","#F1ED7F"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+ 
  theme(legend.position="none")

##########
# Senegal boot
##########
setwd('~/aedes/Results/smc++/Senegal_bootstraps/')
senegal_boot1 <- fread('Senegal_bootstrap_1_100ksites_rp5_smc.csv',data.table = F)
senegal_boot2 <- fread('Senegal_bootstrap_2_100ksites_rp5_smc.csv',data.table = F)
senegal_boot3 <- fread('Senegal_bootstrap_3_100ksites_rp5_smc.csv',data.table = F)
senegal_boot4 <- fread('Senegal_bootstrap_4_100ksites_rp5_smc.csv',data.table = F)
senegal_boot5 <- fread('Senegal_bootstrap_5_100ksites_rp5_smc.csv',data.table = F)
senegal_boot6 <- fread('Senegal_bootstrap_6_100ksites_rp5_smc.csv',data.table = F)
senegal_boot7 <- fread('Senegal_bootstrap_7_100ksites_rp5_smc.csv',data.table = F)
senegal_boot8 <- fread('Senegal_bootstrap_8_100ksites_rp5_smc.csv',data.table = F)
senegal_boot9 <- fread('Senegal_bootstrap_9_100ksites_rp5_smc.csv',data.table = F)
senegal_boot10 <- fread('Senegal_bootstrap_10_100ksites_rp5_smc.csv',data.table = F)
senegal_boot11 <- fread('Senegal_bootstrap_11_100ksites_rp5_smc.csv',data.table = F)
senegal_boot12 <- fread('Senegal_bootstrap_12_100ksites_rp5_smc.csv',data.table = F)
senegal_boot13 <- fread('Senegal_bootstrap_13_100ksites_rp5_smc.csv',data.table = F)
senegal_boot14 <- fread('Senegal_bootstrap_14_100ksites_rp5_smc.csv',data.table = F)
senegal_boot15 <- fread('Senegal_bootstrap_15_100ksites_rp5_smc.csv',data.table = F)
senegal_boot16 <- fread('Senegal_bootstrap_16_100ksites_rp5_smc.csv',data.table = F)
senegal_boot17 <- fread('Senegal_bootstrap_17_100ksites_rp5_smc.csv',data.table = F)
senegal_boot18 <- fread('Senegal_bootstrap_18_100ksites_rp5_smc.csv',data.table = F)
senegal_boot19 <- fread('Senegal_bootstrap_19_100ksites_rp5_smc.csv',data.table = F)
senegal_boot20 <- fread('Senegal_bootstrap_20_100ksites_rp5_smc.csv',data.table = F)

senegal_smc %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(linewidth=1.3,col="#286F29") + 
  geom_line(data=senegal_boot1,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot1,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot2,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot3,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot4,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot5,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot6,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot7,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot8,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot9,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot10,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot11,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot12,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot13,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot14,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot15,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot16,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot17,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot18,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot19,linewidth=1,col='lightgreen',alpha=0.3)+
  geom_line(data=senegal_boot20,linewidth=1,col='lightgreen',alpha=0.3)+
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  ylab('Ne') +
  xlab('Years ago') +
  #ggtitle('Piecewise Ne 10 knots, regularization penalty 5, masked')+
  theme_classic()+
  theme(legend.justification = c(0,0.5))

setwd('~/aedes/Results/smc++/Gabon_bootstraps/')
gabon_boot1 <- fread('gabon_bootstrap_1_100ksites_rp5_smc.csv',data.table = F)
gabon_boot2 <- fread('gabon_bootstrap_2_100ksites_rp5_smc.csv',data.table = F)
gabon_boot3 <- fread('gabon_bootstrap_3_100ksites_rp5_smc.csv',data.table = F)
gabon_boot4 <- fread('gabon_bootstrap_4_100ksites_rp5_smc.csv',data.table = F)
gabon_boot5 <- fread('gabon_bootstrap_5_100ksites_rp5_smc.csv',data.table = F)
gabon_boot6 <- fread('gabon_bootstrap_6_100ksites_rp5_smc.csv',data.table = F)
gabon_boot7 <- fread('gabon_bootstrap_7_100ksites_rp5_smc.csv',data.table = F)
gabon_boot8 <- fread('gabon_bootstrap_8_100ksites_rp5_smc.csv',data.table = F)
gabon_boot9 <- fread('gabon_bootstrap_9_100ksites_rp5_smc.csv',data.table = F)
gabon_boot10 <- fread('gabon_bootstrap_10_100ksites_rp5_smc.csv',data.table = F)
gabon_boot11 <- fread('gabon_bootstrap_11_100ksites_rp5_smc.csv',data.table = F)
gabon_boot12 <- fread('gabon_bootstrap_12_100ksites_rp5_smc.csv',data.table = F)
gabon_boot13 <- fread('gabon_bootstrap_13_100ksites_rp5_smc.csv',data.table = F)
gabon_boot14 <- fread('gabon_bootstrap_14_100ksites_rp5_smc.csv',data.table = F)
gabon_boot15 <- fread('gabon_bootstrap_15_100ksites_rp5_smc.csv',data.table = F)
gabon_boot16 <- fread('gabon_bootstrap_16_100ksites_rp5_smc.csv',data.table = F)
gabon_boot17 <- fread('gabon_bootstrap_17_100ksites_rp5_smc.csv',data.table = F)
gabon_boot18 <- fread('gabon_bootstrap_18_100ksites_rp5_smc.csv',data.table = F)
gabon_boot19 <- fread('gabon_bootstrap_19_100ksites_rp5_smc.csv',data.table = F)
gabon_boot20 <- fread('gabon_bootstrap_20_100ksites_rp5_smc.csv',data.table = F)

gabon_smc %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(linewidth=1.3,col="#062C3B") + 
  geom_line(data=gabon_boot1,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot1,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot2,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot3,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot4,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot5,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot6,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot7,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot8,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot9,linewidth=1,col='#0b5470',alpha=0.3)+
  geom_line(data=gabon_boot10,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot11,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot12,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot13,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot14,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot15,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot16,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot17,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot18,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot19,linewidth=1,col="#0b5470",alpha=0.3)+
  geom_line(data=gabon_boot20,linewidth=1,col="#0b5470",alpha=0.3)+
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  ylab('Ne') +
  xlab('Years ago') +
  #ggtitle('Piecewise Ne 10 knots, regularization penalty 5, masked')+
  theme_classic()+
  theme(legend.justification = c(0,0.5))

# full centromere and 100k sites around genes mask
# rp5
brazil_100k_smc <- fread('Brazil_100ksites_smc.csv',data.table = F)
cali_100k_smc <- fread('Cali_100ksites_smc.csv',data.table = F)
usa_100k_smc <- fread('USA_100ksites_smc.csv',data.table = F)
gabon_100k_smc <- fread('Gabon_100ksites_smc.csv',data.table = F)
kenya_100k_smc <- fread('Kenya_100ksites_smc.csv',data.table = F)
rio_claro_100k_smc <- fread('Rio_Claro_100k_sites_smc.csv',data.table = F)
senegal_100k_smc <- fread('Senegal_100ksites_smc.csv',data.table = F)

full_100k <- rbind(brazil_100k_smc,cali_100k_smc,usa_100k_smc,gabon_100k_smc,kenya_100k_smc,rio_claro_100k_smc,senegal_100k_smc)

full_100k_plot <- full_100k %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne') +
  xlab('Years ago') +
  #ggtitle('Piecewise Ne 10 knots, regularization penalty 5, masked')+
  theme_classic()+
  theme(legend.justification = c(0,0.5))

full_rp5_africa <- rbind(gabon_100k_smc,kenya_100k_smc,senegal_100k_smc)
full_rp5_americas <- rbind(brazil_100k_smc,cali_100k_smc,rio_claro_100k_smc,usa_100k_smc)
full_rp5 <- rbind(brazil_100k_smc,cali_100k_smc,rio_claro_100k_smc,usa_100k_smc,gabon_100k_smc,kenya_100k_smc,senegal_100k_smc)

full_smc <- full_rp5 %>% mutate(label=ifelse(label=="RC","Río Claro",label)) %>% mutate(label=factor(label,levels=c('Brazil','Cali','Río Claro','USA','Gabon','Kenya','Senegal'))) %>%
  ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=c("#D4940E","#976407","#36260A","#F1ED7F","#062C3B","#00504F","#286F29"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

smc_legend <- get_legend(
  # create some space to the left of the legend
  full_smc + theme(legend.box.margin = margin(0, 0, 0, 12))
)

africas_smc <- full_rp5_africa %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=5),labels = label_log()) +
  scale_color_manual(values=c("#062C3B","#00504F","#286F29"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+ 
  theme(legend.position="none")

americas_smc <- full_rp5_africa %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=0,col='white',show.legend = F)+
  geom_line(data=full_rp5_americas,size=1) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=c("#D4940E","#976407","#36260A","#F1ED7F"),name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+ 
  theme(legend.position="none")


# full centromere and 100k sites around genes mask
# rp4
brazil_100k_smc_rp4 <- fread('Brazil_100ksites_rp4_smc.csv',data.table = F)
cali_100k_smc_rp4 <- fread('Cali_100ksites_rp4_smc.csv',data.table = F)
usa_100k_smc_rp4 <- fread('USA_100ksites_rp4_smc.csv',data.table = F)
gabon_100k_smc_rp4 <- fread('Gabon_100ksites_rp4_smc.csv',data.table = F)
kenya_100k_smc_rp4 <- fread('Kenya_100ksites_rp4_smc.csv',data.table = F)
rio_claro_100k_smc_rp4 <- fread('Rio_Claro_100ksites_rp4_smc.csv',data.table = F)
senegal_100k_smc_rp4 <- fread('Senegal_100ksites_rp4_smc.csv',data.table = F)

full_100k_rp4 <- rbind(brazil_100k_smc,cali_100k_smc_rp4,usa_100k_smc_rp4,gabon_100k_smc_rp4,kenya_100k_smc_rp4,rio_claro_100k_smc_rp4,senegal_100k_smc_rp4)

full_100k_rp4 %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Piecewise Ne 10 knots, regularization penalty 4, masked')+
  theme_classic()

# full centromere and 100k sites around genes mask
# rp6
brazil_100k_smc_rp6 <- fread('Brazil_100ksites_rp6_smc.csv',data.table = F)
cali_100k_smc_rp6 <- fread('Cali_100ksites_rp6_smc.csv',data.table = F)
usa_100k_smc_rp6 <- fread('USA_100ksites_rp6_smc.csv',data.table = F)
gabon_100k_smc_rp6 <- fread('Gabon_100ksites_rp6_smc.csv',data.table = F)
kenya_100k_smc_rp6 <- fread('Kenya_100ksites_rp6_smc.csv',data.table = F)
rio_claro_100k_smc_rp6 <- fread('Rio_Claro_100ksites_rp6_smc.csv',data.table = F)
senegal_100k_smc_rp6 <- fread('Senegal_100ksites_rp6_smc.csv',data.table = F)

full_100k_rp6 <- rbind(brazil_100k_smc_rp6,cali_100k_smc_rp6,usa_100k_smc_rp6,gabon_100k_smc_rp6,kenya_100k_smc_rp6,rio_claro_100k_smc_rp6,senegal_100k_smc_rp6)

full_100k_rp6 %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Piecewise Ne 10 knots, regularization penalty 6, masked')+
  theme_classic()

# RC different params
rio_claro_rp4_smc <- fread('Rio_Claro_rp4_smc.csv',data.table = F) %>% mutate(param='rp4')
rio_claro_smc <- fread('Rio_Claro_smc.csv',data.table = F) %>% mutate(param='default')
rio_claro_100k_smc <- fread('Rio_Claro_100k_sites_smc.csv',data.table = F) %>% mutate(param='100k_cent')

full_rc <- rbind(rio_claro_100k_smc,rio_claro_rp4_smc,rio_claro_smc)

full_rc %>% ggplot(aes(x=x,y=y,lty=param)) + 
  geom_line(size=1.3,col="#00504F") + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Rio Claro parameter space')+
  theme_classic()



full_rc %>% ggplot(aes(x=x,y=y,lty=param)) + 
  geom_line(size=1.3,col="#00504F") + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  ggtitle('Rio Claro parameter space')+
  theme_classic()

# USA populations
# full centromere and 100k sites around genes mask
pal=met.brewer("OKeeffe1", 6, type="continuous")
pal[6]="#F1ED7F"
# rp5
NorCal_100k_smc_rp5 <- fread('NorCal_100ksites_rp5_smc.csv',data.table = F)
SoCal_100k_smc_rp5 <- fread('SoCal_100ksites_rp5_smc.csv',data.table = F)
Florida_100k_smc_rp5 <- fread('Florida_100ksites_rp5_smc.csv',data.table = F)
SoCalExeter_100k_smc_rp5 <- fread('SoCalExeter_100ksites_rp5_smc.csv',data.table = F)
FloridaExeter_100k_smc_rp5 <- fread('FloridaExeter_100ksites_rp5_smc.csv',data.table = F)
Exeter_100k_smc_rp5 <- fread('Exeter_100ksites_rp5_smc.csv',data.table = F)
Clovis_100k_smc_rp5 <- fread('Clovis_100ksites_rp5_smc.csv',data.table = F)
USA_100k_smc_rp5 <- fread('USA_100ksites_smc.csv',data.table = F)

full_us_100k_rp5 <- rbind(NorCal_100k_smc_rp5,SoCal_100k_smc_rp5,Florida_100k_smc_rp5,Exeter_100k_smc_rp5,Clovis_100k_smc_rp5,USA_100k_smc_rp5)

full_us_100k_rp5_plot <- full_rp5_africa %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=0,col='white',show.legend = F)+
  geom_line(data=full_us_100k_rp5,size=1) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log()) +
  scale_y_log10(breaks = breaks_log(n=6),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab(expression(italic("N"["e"]))) +
  xlab('Years ago') +
  theme_cowplot()+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

full_us_100k_rp5_plot <- full_us_100k_rp5 %>% ggplot(aes(x=x,y=y,col=label)) + 
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log(),limits=c(0.068,750000)) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=pal,name="population") +
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+ 
  theme_classic()

plot_grid(full_100k_plot,full_us_100k_rp5_plot,
          align = 'v',
          ncol = 1, 
          axis = 'b')

#-------------------------------------#
# msprime likelihood fits
#-------------------------------------#
options(digits=4)

#USA
USA_100k_sim_LL_rp4 <- fread('USA_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
USA_100k_obs_LL_rp4 <- fread('USA_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
USA_100k_delta_rp4 <- USA_100k_obs_LL_rp4$V1 - USA_100k_sim_LL_rp4$V1
USA_100k_sim_sfs_rp4 <- fread('USA_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

USA_100k_obs_sfs_rp4 <- fread('USA_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

USA_100k_sfs_rp4 <- rbind(USA_100k_sim_sfs_rp4,USA_100k_obs_sfs_rp4) 

USA_100k_resid_rp4 <- merge(USA_100k_sim_sfs_rp4,USA_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

USA_100k_sfs_rp4_plot <- USA_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("USA rp4,", Delta,"LL=",!!USA_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_resid_rp4_plot <- USA_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_rp4_plot <- cowplot::plot_grid(USA_100k_sfs_rp4_plot,USA_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

USA_100k_sim_LL_rp5 <- fread('USA_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
USA_100k_obs_LL_rp5 <- fread('USA_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
USA_100k_delta_rp5 <- USA_100k_obs_LL_rp5$V1 - USA_100k_sim_LL_rp5$V1
USA_100k_sim_sfs_rp5 <- fread('USA_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

USA_100k_obs_sfs_rp5 <- fread('USA_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

USA_100k_sfs_rp5 <- rbind(USA_100k_sim_sfs_rp5,USA_100k_obs_sfs_rp5) 

USA_100k_resid_rp5 <- merge(USA_100k_sim_sfs_rp5,USA_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

USA_100k_sfs_rp5_plot <- USA_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("USA rp5, ", Delta, " LL=",!!USA_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_resid_rp5_plot <- USA_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_rp5_plot <- cowplot::plot_grid(USA_100k_sfs_rp5_plot,USA_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

USA_100k_sim_LL_rp6 <- fread('USA_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
USA_100k_obs_LL_rp6 <- fread('USA_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
USA_100k_delta_rp6 <- USA_100k_obs_LL_rp6$V1 - USA_100k_sim_LL_rp6$V1
USA_100k_sim_sfs_rp6 <- fread('USA_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

USA_100k_obs_sfs_rp6 <- fread('USA_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

USA_100k_sfs_rp6 <- rbind(USA_100k_sim_sfs_rp6,USA_100k_obs_sfs_rp6) 

USA_100k_resid_rp6 <- merge(USA_100k_sim_sfs_rp6,USA_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

USA_100k_sfs_rp6_plot <- USA_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("USA rp6, ", Delta, " LL=",!!USA_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_resid_rp6_plot <- USA_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

USA_100k_rp6_plot <- cowplot::plot_grid(USA_100k_sfs_rp6_plot,USA_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(USA_100k_rp4_plot,USA_100k_rp5_plot,USA_100k_rp6_plot,ncol=3)

#Brazil
Brazil_100k_sim_LL_rp4 <- fread('Brazil_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Brazil_100k_obs_LL_rp4 <- fread('Brazil_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Brazil_100k_delta_rp4 <- Brazil_100k_obs_LL_rp4$V1 - Brazil_100k_sim_LL_rp4$V1
Brazil_100k_sim_sfs_rp4 <- fread('Brazil_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Brazil_100k_obs_sfs_rp4 <- fread('Brazil_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Brazil_100k_sfs_rp4 <- rbind(Brazil_100k_sim_sfs_rp4,Brazil_100k_obs_sfs_rp4) 

Brazil_100k_resid_rp4 <- merge(Brazil_100k_sim_sfs_rp4,Brazil_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Brazil_100k_sfs_rp4_plot <- Brazil_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Brazil rp4, ", Delta, " LL=",!!Brazil_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_resid_rp4_plot <- Brazil_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_rp4_plot <- cowplot::plot_grid(Brazil_100k_sfs_rp4_plot,Brazil_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Brazil_100k_sim_LL_rp5 <- fread('Brazil_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Brazil_100k_obs_LL_rp5 <- fread('Brazil_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Brazil_100k_delta_rp5 <- Brazil_100k_obs_LL_rp5$V1 - Brazil_100k_sim_LL_rp5$V1
Brazil_100k_sim_sfs_rp5 <- fread('Brazil_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Brazil_100k_obs_sfs_rp5 <- fread('Brazil_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Brazil_100k_sfs_rp5 <- rbind(Brazil_100k_sim_sfs_rp5,Brazil_100k_obs_sfs_rp5) 

Brazil_100k_resid_rp5 <- merge(Brazil_100k_sim_sfs_rp5,Brazil_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Brazil_100k_sfs_rp5_plot <- Brazil_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Brazil rp5, ", Delta, " LL=",!!Brazil_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_resid_rp5_plot <- Brazil_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_rp5_plot <- cowplot::plot_grid(Brazil_100k_sfs_rp5_plot,Brazil_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Brazil_100k_sim_LL_rp6 <- fread('Brazil_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Brazil_100k_obs_LL_rp6 <- fread('Brazil_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Brazil_100k_delta_rp6 <- Brazil_100k_obs_LL_rp6$V1 - Brazil_100k_sim_LL_rp6$V1
Brazil_100k_sim_sfs_rp6 <- fread('Brazil_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Brazil_100k_obs_sfs_rp6 <- fread('Brazil_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Brazil_100k_sfs_rp6 <- rbind(Brazil_100k_sim_sfs_rp6,Brazil_100k_obs_sfs_rp6) 

Brazil_100k_resid_rp6 <- merge(Brazil_100k_sim_sfs_rp6,Brazil_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Brazil_100k_sfs_rp6_plot <- Brazil_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Brazil rp6, ", Delta, " LL=",!!Brazil_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_resid_rp6_plot <- Brazil_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Brazil_100k_rp6_plot <- cowplot::plot_grid(Brazil_100k_sfs_rp6_plot,Brazil_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Brazil_100k_rp4_plot,Brazil_100k_rp5_plot,Brazil_100k_rp6_plot,ncol=3)

#Cali
Cali_100k_sim_LL_rp4 <- fread('Cali_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Cali_100k_obs_LL_rp4 <- fread('Cali_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Cali_100k_delta_rp4 <- Cali_100k_obs_LL_rp4$V1 - Cali_100k_sim_LL_rp4$V1
Cali_100k_sim_sfs_rp4 <- fread('Cali_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Cali_100k_obs_sfs_rp4 <- fread('Cali_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Cali_100k_sfs_rp4 <- rbind(Cali_100k_sim_sfs_rp4,Cali_100k_obs_sfs_rp4) 

Cali_100k_resid_rp4 <- merge(Cali_100k_sim_sfs_rp4,Cali_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Cali_100k_sfs_rp4_plot <- Cali_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Cali rp4, ", Delta, " LL=",!!Cali_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_resid_rp4_plot <- Cali_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_rp4_plot <- cowplot::plot_grid(Cali_100k_sfs_rp4_plot,Cali_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Cali_100k_sim_LL_rp5 <- fread('Cali_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Cali_100k_obs_LL_rp5 <- fread('Cali_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Cali_100k_delta_rp5 <- Cali_100k_obs_LL_rp5$V1 - Cali_100k_sim_LL_rp5$V1
Cali_100k_sim_sfs_rp5 <- fread('Cali_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Cali_100k_obs_sfs_rp5 <- fread('Cali_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Cali_100k_sfs_rp5 <- rbind(Cali_100k_sim_sfs_rp5,Cali_100k_obs_sfs_rp5) 

Cali_100k_resid_rp5 <- merge(Cali_100k_sim_sfs_rp5,Cali_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Cali_100k_sfs_rp5_plot <- Cali_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Cali rp5, ", Delta, " LL=",!!Cali_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_resid_rp5_plot <- Cali_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_rp5_plot <- cowplot::plot_grid(Cali_100k_sfs_rp5_plot,Cali_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Cali_100k_sim_LL_rp6 <- fread('Cali_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Cali_100k_obs_LL_rp6 <- fread('Cali_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Cali_100k_delta_rp6 <- Cali_100k_obs_LL_rp6$V1 - Cali_100k_sim_LL_rp6$V1
Cali_100k_sim_sfs_rp6 <- fread('Cali_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Cali_100k_obs_sfs_rp6 <- fread('Cali_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Cali_100k_sfs_rp6 <- rbind(Cali_100k_sim_sfs_rp6,Cali_100k_obs_sfs_rp6) 

Cali_100k_resid_rp6 <- merge(Cali_100k_sim_sfs_rp6,Cali_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Cali_100k_sfs_rp6_plot <- Cali_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Cali rp6, ", Delta, " LL=",!!Cali_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_resid_rp6_plot <- Cali_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Cali_100k_rp6_plot <- cowplot::plot_grid(Cali_100k_sfs_rp6_plot,Cali_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Cali_100k_rp4_plot,Cali_100k_rp5_plot,Cali_100k_rp6_plot,ncol=3)

#Rio_Claro
Rio_Claro_100k_sim_LL_rp4 <- fread('Rio_Claro_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Rio_Claro_100k_obs_LL_rp4 <- fread('Rio_Claro_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Rio_Claro_100k_delta_rp4 <- Rio_Claro_100k_obs_LL_rp4$V1 - Rio_Claro_100k_sim_LL_rp4$V1
Rio_Claro_100k_sim_sfs_rp4 <- fread('Rio_Claro_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_obs_sfs_rp4 <- fread('Rio_Claro_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_sfs_rp4 <- rbind(Rio_Claro_100k_sim_sfs_rp4,Rio_Claro_100k_obs_sfs_rp4) 

Rio_Claro_100k_resid_rp4 <- merge(Rio_Claro_100k_sim_sfs_rp4,Rio_Claro_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Rio_Claro_100k_sfs_rp4_plot <- Rio_Claro_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Río Claro rp4, ", Delta, " LL=",!!Rio_Claro_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_resid_rp4_plot <- Rio_Claro_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_rp4_plot <- cowplot::plot_grid(Rio_Claro_100k_sfs_rp4_plot,Rio_Claro_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Rio_Claro_100k_sim_LL_rp5 <- fread('Rio_Claro_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Rio_Claro_100k_obs_LL_rp5 <- fread('Rio_Claro_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Rio_Claro_100k_delta_rp5 <- Rio_Claro_100k_obs_LL_rp5$V1 - Rio_Claro_100k_sim_LL_rp5$V1
Rio_Claro_100k_sim_sfs_rp5 <- fread('Rio_Claro_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_obs_sfs_rp5 <- fread('Rio_Claro_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_sfs_rp5 <- rbind(Rio_Claro_100k_sim_sfs_rp5,Rio_Claro_100k_obs_sfs_rp5) 

Rio_Claro_100k_resid_rp5 <- merge(Rio_Claro_100k_sim_sfs_rp5,Rio_Claro_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Rio_Claro_100k_sfs_rp5_plot <- Rio_Claro_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Río Claro rp5, ", Delta, " LL=",!!Rio_Claro_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_resid_rp5_plot <- Rio_Claro_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_rp5_plot <- cowplot::plot_grid(Rio_Claro_100k_sfs_rp5_plot,Rio_Claro_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Rio_Claro_100k_sim_LL_rp6 <- fread('Rio_Claro_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Rio_Claro_100k_obs_LL_rp6 <- fread('Rio_Claro_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Rio_Claro_100k_delta_rp6 <- Rio_Claro_100k_obs_LL_rp6$V1 - Rio_Claro_100k_sim_LL_rp6$V1
Rio_Claro_100k_sim_sfs_rp6 <- fread('Rio_Claro_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_obs_sfs_rp6 <- fread('Rio_Claro_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Rio_Claro_100k_sfs_rp6 <- rbind(Rio_Claro_100k_sim_sfs_rp6,Rio_Claro_100k_obs_sfs_rp6) 

Rio_Claro_100k_resid_rp6 <- merge(Rio_Claro_100k_sim_sfs_rp6,Rio_Claro_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Rio_Claro_100k_sfs_rp6_plot <- Rio_Claro_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Río Claro rp6, ", Delta, " LL=",!!Rio_Claro_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_resid_rp6_plot <- Rio_Claro_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Rio_Claro_100k_rp6_plot <- cowplot::plot_grid(Rio_Claro_100k_sfs_rp6_plot,Rio_Claro_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Rio_Claro_100k_rp4_plot,Rio_Claro_100k_rp5_plot,Rio_Claro_100k_rp6_plot,ncol=3)

#Senegal
Senegal_100k_sim_LL_rp4 <- fread('Senegal_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Senegal_100k_obs_LL_rp4 <- fread('Senegal_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Senegal_100k_delta_rp4 <- Senegal_100k_obs_LL_rp4$V1 - Senegal_100k_sim_LL_rp4$V1
Senegal_100k_sim_sfs_rp4 <- fread('Senegal_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Senegal_100k_obs_sfs_rp4 <- fread('Senegal_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Senegal_100k_sfs_rp4 <- rbind(Senegal_100k_sim_sfs_rp4,Senegal_100k_obs_sfs_rp4) 

Senegal_100k_resid_rp4 <- merge(Senegal_100k_sim_sfs_rp4,Senegal_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Senegal_100k_sfs_rp4_plot <- Senegal_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Senegal rp4, ", Delta, " LL=",!!Senegal_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_resid_rp4_plot <- Senegal_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_rp4_plot <- cowplot::plot_grid(Senegal_100k_sfs_rp4_plot,Senegal_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Senegal_100k_sim_LL_rp5 <- fread('Senegal_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Senegal_100k_obs_LL_rp5 <- fread('Senegal_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Senegal_100k_delta_rp5 <- Senegal_100k_obs_LL_rp5$V1 - Senegal_100k_sim_LL_rp5$V1
Senegal_100k_sim_sfs_rp5 <- fread('Senegal_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Senegal_100k_obs_sfs_rp5 <- fread('Senegal_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Senegal_100k_sfs_rp5 <- rbind(Senegal_100k_sim_sfs_rp5,Senegal_100k_obs_sfs_rp5) 

Senegal_100k_resid_rp5 <- merge(Senegal_100k_sim_sfs_rp5,Senegal_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Senegal_100k_sfs_rp5_plot <- Senegal_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Senegal rp5, ", Delta, " LL=",!!Senegal_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_resid_rp5_plot <- Senegal_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_rp5_plot <- cowplot::plot_grid(Senegal_100k_sfs_rp5_plot,Senegal_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Senegal_100k_sim_LL_rp6 <- fread('Senegal_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Senegal_100k_obs_LL_rp6 <- fread('Senegal_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Senegal_100k_delta_rp6 <- Senegal_100k_obs_LL_rp6$V1 - Senegal_100k_sim_LL_rp6$V1
Senegal_100k_sim_sfs_rp6 <- fread('Senegal_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Senegal_100k_obs_sfs_rp6 <- fread('Senegal_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Senegal_100k_sfs_rp6 <- rbind(Senegal_100k_sim_sfs_rp6,Senegal_100k_obs_sfs_rp6) 

Senegal_100k_resid_rp6 <- merge(Senegal_100k_sim_sfs_rp6,Senegal_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Senegal_100k_sfs_rp6_plot <- Senegal_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Senegal rp6, ", Delta, " LL=",!!Senegal_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_resid_rp6_plot <- Senegal_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Senegal_100k_rp6_plot <- cowplot::plot_grid(Senegal_100k_sfs_rp6_plot,Senegal_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Senegal_100k_rp4_plot,Senegal_100k_rp5_plot,Senegal_100k_rp6_plot,ncol=3)

#Kenya
Kenya_100k_sim_LL_rp4 <- fread('Kenya_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Kenya_100k_obs_LL_rp4 <- fread('Kenya_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Kenya_100k_delta_rp4 <- Kenya_100k_obs_LL_rp4$V1 - Kenya_100k_sim_LL_rp4$V1
Kenya_100k_sim_sfs_rp4 <- fread('Kenya_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Kenya_100k_obs_sfs_rp4 <- fread('Kenya_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Kenya_100k_sfs_rp4 <- rbind(Kenya_100k_sim_sfs_rp4,Kenya_100k_obs_sfs_rp4) 

Kenya_100k_resid_rp4 <- merge(Kenya_100k_sim_sfs_rp4,Kenya_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Kenya_100k_sfs_rp4_plot <- Kenya_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Kenya rp4, ", Delta, " LL=",!!Kenya_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_resid_rp4_plot <- Kenya_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_rp4_plot <- cowplot::plot_grid(Kenya_100k_sfs_rp4_plot,Kenya_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Kenya_100k_sim_LL_rp5 <- fread('Kenya_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Kenya_100k_obs_LL_rp5 <- fread('Kenya_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Kenya_100k_delta_rp5 <- Kenya_100k_obs_LL_rp5$V1 - Kenya_100k_sim_LL_rp5$V1
Kenya_100k_sim_sfs_rp5 <- fread('Kenya_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Kenya_100k_obs_sfs_rp5 <- fread('Kenya_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Kenya_100k_sfs_rp5 <- rbind(Kenya_100k_sim_sfs_rp5,Kenya_100k_obs_sfs_rp5) 

Kenya_100k_resid_rp5 <- merge(Kenya_100k_sim_sfs_rp5,Kenya_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Kenya_100k_sfs_rp5_plot <- Kenya_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Kenya rp5, ", Delta, " LL=",!!Kenya_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_resid_rp5_plot <- Kenya_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_rp5_plot <- cowplot::plot_grid(Kenya_100k_sfs_rp5_plot,Kenya_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Kenya_100k_sim_LL_rp6 <- fread('Kenya_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Kenya_100k_obs_LL_rp6 <- fread('Kenya_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Kenya_100k_delta_rp6 <- Kenya_100k_obs_LL_rp6$V1 - Kenya_100k_sim_LL_rp6$V1
Kenya_100k_sim_sfs_rp6 <- fread('Kenya_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Kenya_100k_obs_sfs_rp6 <- fread('Kenya_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Kenya_100k_sfs_rp6 <- rbind(Kenya_100k_sim_sfs_rp6,Kenya_100k_obs_sfs_rp6) 

Kenya_100k_resid_rp6 <- merge(Kenya_100k_sim_sfs_rp6,Kenya_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Kenya_100k_sfs_rp6_plot <- Kenya_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Kenya rp6, ", Delta, " LL=",!!Kenya_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_resid_rp6_plot <- Kenya_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Kenya_100k_rp6_plot <- cowplot::plot_grid(Kenya_100k_sfs_rp6_plot,Kenya_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Kenya_100k_rp4_plot,Kenya_100k_rp5_plot,Kenya_100k_rp6_plot,ncol=3)

#Gabon
Gabon_100k_sim_LL_rp4 <- fread('Gabon_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1)
Gabon_100k_obs_LL_rp4 <- fread('Gabon_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Gabon_100k_delta_rp4 <- Gabon_100k_obs_LL_rp4$V1 - Gabon_100k_sim_LL_rp4$V1
Gabon_100k_sim_sfs_rp4 <- fread('Gabon_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Gabon_100k_obs_sfs_rp4 <- fread('Gabon_100k_sites_rp4_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Gabon_100k_sfs_rp4 <- rbind(Gabon_100k_sim_sfs_rp4,Gabon_100k_obs_sfs_rp4) 

Gabon_100k_resid_rp4 <- merge(Gabon_100k_sim_sfs_rp4,Gabon_100k_obs_sfs_rp4,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Gabon_100k_sfs_rp4_plot <- Gabon_100k_sfs_rp4 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Gabon rp4, ", Delta, " LL=",!!Gabon_100k_delta_rp4)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_resid_rp4_plot <- Gabon_100k_resid_rp4 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_rp4_plot <- cowplot::plot_grid(Gabon_100k_sfs_rp4_plot,Gabon_100k_resid_rp4_plot,nrow = 2,align='hv',axis='tblr')

Gabon_100k_sim_LL_rp5 <- fread('Gabon_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1)
Gabon_100k_obs_LL_rp5 <- fread('Gabon_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Gabon_100k_delta_rp5 <- Gabon_100k_obs_LL_rp5$V1 - Gabon_100k_sim_LL_rp5$V1
Gabon_100k_sim_sfs_rp5 <- fread('Gabon_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Gabon_100k_obs_sfs_rp5 <- fread('Gabon_100k_sites_rp5_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Gabon_100k_sfs_rp5 <- rbind(Gabon_100k_sim_sfs_rp5,Gabon_100k_obs_sfs_rp5) 

Gabon_100k_resid_rp5 <- merge(Gabon_100k_sim_sfs_rp5,Gabon_100k_obs_sfs_rp5,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Gabon_100k_sfs_rp5_plot <- Gabon_100k_sfs_rp5 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Gabon rp5, ", Delta, " LL=",!!Gabon_100k_delta_rp5)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_resid_rp5_plot <- Gabon_100k_resid_rp5 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_rp5_plot <- cowplot::plot_grid(Gabon_100k_sfs_rp5_plot,Gabon_100k_resid_rp5_plot,nrow = 2,align='hv',axis='tblr')

Gabon_100k_sim_LL_rp6 <- fread('Gabon_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1)
Gabon_100k_obs_LL_rp6 <- fread('Gabon_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=1)
Gabon_100k_delta_rp6 <- Gabon_100k_obs_LL_rp6$V1 - Gabon_100k_sim_LL_rp6$V1
Gabon_100k_sim_sfs_rp6 <- fread('Gabon_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=2) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="sim",rel_af=(rowid/max(rowid))/2)

Gabon_100k_obs_sfs_rp6 <- fread('Gabon_100k_sites_rp6_eval.txt',data.table = F, header=F,nrows=1,skip=3) %>% 
  transpose() %>% 
  na.omit() %>% 
  mutate(freq=V1/sum(V1)) %>%
  select(freq) %>%
  tibble::rowid_to_column() %>%
  mutate(source="obs",rel_af=(rowid/max(rowid))/2)

Gabon_100k_sfs_rp6 <- rbind(Gabon_100k_sim_sfs_rp6,Gabon_100k_obs_sfs_rp6) 

Gabon_100k_resid_rp6 <- merge(Gabon_100k_sim_sfs_rp6,Gabon_100k_obs_sfs_rp6,by='rel_af') %>%
  mutate(resid=freq.x-freq.y) %>% select(rel_af,resid)

Gabon_100k_sfs_rp6_plot <- Gabon_100k_sfs_rp6 %>% ggplot(aes(x=rel_af,y=freq,col=source))+
  geom_line() +
  geom_point() +
  ylab('Relative count')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  ggtitle(expr(paste("Gabon rp6, ", Delta, " LL=",!!Gabon_100k_delta_rp6)))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        plot.title = element_text(size=10))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_resid_rp6_plot <- Gabon_100k_resid_rp6 %>% ggplot(aes(x=rel_af,y=resid))+
  geom_line(col="darkgreen") +
  geom_point(col="darkgreen") +
  ylab('Residual')+
  xlab('Relative (folded) frequency')+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))

Gabon_100k_rp6_plot <- cowplot::plot_grid(Gabon_100k_sfs_rp6_plot,Gabon_100k_resid_rp6_plot,nrow = 2,align='hv',axis='tblr')


cowplot::plot_grid(Gabon_100k_rp4_plot,Gabon_100k_rp5_plot,Gabon_100k_rp6_plot,ncol=3)
#-------------------------------------#
# SMC++ Split 
#-------------------------------------#
all_splits <- data.frame(N= numeric(0), 
                         split_model= numeric(0), 
                         split_gen = numeric(0), 
                         split_year = numeric(0), 
                         pair = character(0),
                         region = character(0))

Brazil_Gabon_json <- fromJSON("Brazil-Gabon.model.final.json")
Brazil_Gabon_split <- as.numeric(Brazil_Gabon_json$model$split)
Brazil_Gabon_N <- as.numeric(Brazil_Gabon_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Brazil_Gabon_N,
                                     split_model=Brazil_Gabon_split,
                                     split_gen=2*Brazil_Gabon_N*Brazil_Gabon_split,
                                     split_year=2*Brazil_Gabon_N*Brazil_Gabon_split*0.067,
                                     pair="Brazil-Gabon",
                                     region="Americas-Africa")

Brazil_Kenya_json <- fromJSON("Brazil-Kenya.model.final.json")
Brazil_Kenya_split <- as.numeric(Brazil_Kenya_json$model$split)
Brazil_Kenya_N <- as.numeric(Brazil_Kenya_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Brazil_Kenya_N,
                                     split_model=Brazil_Kenya_split,
                                     split_gen=2*Brazil_Kenya_N*Brazil_Kenya_split,
                                     split_year=2*Brazil_Kenya_N*Brazil_Kenya_split*0.067,
                                     pair="Brazil-Kenya",
                                     region="Americas-Africa")

Brazil_Senegal_json <- fromJSON("Brazil-Senegal.model.final.json")
Brazil_Senegal_split <- as.numeric(Brazil_Senegal_json$model$split)
Brazil_Senegal_N <- as.numeric(Brazil_Senegal_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Brazil_Senegal_N,
                                     split_model=Brazil_Senegal_split,
                                     split_gen=2*Brazil_Senegal_N*Brazil_Senegal_split,
                                     split_year=2*Brazil_Senegal_N*Brazil_Senegal_split*0.067,
                                     pair="Brazil-Senegal",
                                     region="Americas-Africa")

Brazil_USA_json <- fromJSON("Brazil-USA.model.final.json")
Brazil_USA_split <- as.numeric(Brazil_USA_json$model$split)
Brazil_USA_N <- as.numeric(Brazil_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Brazil_USA_N,
                                     split_model=Brazil_USA_split,
                                     split_gen=2*Brazil_USA_N*Brazil_USA_split,
                                     split_year=2*Brazil_USA_N*Brazil_USA_split*0.067,
                                     pair="Brazil-USA",
                                     region="Americas")

Cali_Brazil_json <- fromJSON("Cali-Brazil.model.final.json")
Cali_Brazil_split <- as.numeric(Cali_Brazil_json$model$split)
Cali_Brazil_N <- as.numeric(Cali_Brazil_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Cali_Brazil_N,
                                     split_model=Cali_Brazil_split,
                                     split_gen=2*Cali_Brazil_N*Cali_Brazil_split,
                                     split_year=2*Cali_Brazil_N*Cali_Brazil_split*0.067,
                                     pair="Cali-Brazil",
                                     region="Americas")

Cali_Gabon_json <- fromJSON("Cali-Gabon.model.final.json")
Cali_Gabon_split <- as.numeric(Cali_Gabon_json$model$split)
Cali_Gabon_N <- as.numeric(Cali_Gabon_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Cali_Gabon_N,
                                     split_model=Cali_Gabon_split,
                                     split_gen=2*Cali_Gabon_N*Cali_Gabon_split,
                                     split_year=2*Cali_Gabon_N*Cali_Gabon_split*0.067,
                                     pair="Cali-Gabon",
                                     region="Americas-Africa")

Cali_Kenya_json <- fromJSON("Cali-Kenya.model.final.json")
Cali_Kenya_split <- as.numeric(Cali_Kenya_json$model$split)
Cali_Kenya_N <- as.numeric(Cali_Kenya_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Cali_Kenya_N,
                                     split_model=Cali_Kenya_split,
                                     split_gen=2*Cali_Kenya_N*Cali_Kenya_split,
                                     split_year=2*Cali_Kenya_N*Cali_Kenya_split*0.067,
                                     pair="Cali-Kenya",
                                     region="Americas-Africa")

Cali_Senegal_json <- fromJSON("Cali-Senegal.model.final.json")
Cali_Senegal_split <- as.numeric(Cali_Senegal_json$model$split)
Cali_Senegal_N <- as.numeric(Cali_Senegal_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Cali_Senegal_N,
                                     split_model=Cali_Senegal_split,
                                     split_gen=2*Cali_Senegal_N*Cali_Senegal_split,
                                     split_year=2*Cali_Senegal_N*Cali_Senegal_split*0.067,
                                     pair="Cali-Senegal",
                                     region="Americas-Africa")

Cali_USA_json <- fromJSON("Cali-USA.model.final.json")
Cali_USA_split <- as.numeric(Cali_USA_json$model$split)
Cali_USA_N <- as.numeric(Cali_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Cali_USA_N,
                                     split_model=Cali_USA_split,
                                     split_gen=2*Cali_USA_N*Cali_USA_split,
                                     split_year=2*Cali_USA_N*Cali_USA_split*0.067,
                                     pair="Cali-USA",
                                     region="Americas")

Kenya_Gabon_json <- fromJSON("Kenya-Gabon.model.final.json")
Kenya_Gabon_split <- as.numeric(Kenya_Gabon_json$model$split)
Kenya_Gabon_N <- as.numeric(Kenya_Gabon_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Kenya_Gabon_N,
                                     split_model=Kenya_Gabon_split,
                                     split_gen=2*Kenya_Gabon_N*Kenya_Gabon_split,
                                     split_year=2*Kenya_Gabon_N*Kenya_Gabon_split*0.067,
                                     pair="Kenya-Gabon",
                                     region="Africa")

Kenya_Senegal_json <- fromJSON("Kenya-Senegal.model.final.json")
Kenya_Senegal_split <- as.numeric(Kenya_Senegal_json$model$split)
Kenya_Senegal_N <- as.numeric(Kenya_Senegal_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Kenya_Senegal_N,
                                     split_model=Kenya_Senegal_split,
                                     split_gen=2*Kenya_Senegal_N*Kenya_Senegal_split,
                                     split_year=2*Kenya_Senegal_N*Kenya_Senegal_split*0.067,
                                     pair="Kenya-Senegal",
                                     region="Africa")

Kenya_USA_json <- fromJSON("Kenya-USA.model.final.json")
Kenya_USA_split <- as.numeric(Kenya_USA_json$model$split)
Kenya_USA_N <- as.numeric(Kenya_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Kenya_USA_N,
                                     split_model=Kenya_USA_split,
                                     split_gen=2*Kenya_USA_N*Kenya_USA_split,
                                     split_year=2*Kenya_USA_N*Kenya_USA_split*0.067,
                                     pair="Kenya-USA",
                                     region="Americas-Africa")

Rio_Claro_Brazil_json <- fromJSON("Rio_Claro-Brazil.model.final.json")
Rio_Claro_Brazil_split <- as.numeric(Rio_Claro_Brazil_json$model$split)
Rio_Claro_Brazil_N <- as.numeric(Rio_Claro_Brazil_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_Brazil_N,
                                     split_model=Rio_Claro_Brazil_split,
                                     split_gen=2*Rio_Claro_Brazil_N*Rio_Claro_Brazil_split,
                                     split_year=2*Rio_Claro_Brazil_N*Rio_Claro_Brazil_split*0.067,
                                     pair="Rio_Claro-Brazil",
                                     region="Americas")

Rio_Claro_Cali_json <- fromJSON("Rio_Claro-Cali.model.final.json")
Rio_Claro_Cali_split <- as.numeric(Rio_Claro_Cali_json$model$split)
Rio_Claro_Cali_N <- as.numeric(Rio_Claro_Cali_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_Cali_N,
                                     split_model=Rio_Claro_Cali_split,
                                     split_gen=2*Rio_Claro_Cali_N*Rio_Claro_Cali_split,
                                     split_year=2*Rio_Claro_Cali_N*Rio_Claro_Cali_split*0.067,
                                     pair="Rio_Claro-Cali",
                                     region="Americas")

Rio_Claro_Gabon_json <- fromJSON("Rio_Claro-Gabon.model.final.json")
Rio_Claro_Gabon_split <- as.numeric(Rio_Claro_Gabon_json$model$split)
Rio_Claro_Gabon_N <- as.numeric(Rio_Claro_Gabon_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_Gabon_N,
                                     split_model=Rio_Claro_Gabon_split,
                                     split_gen=2*Rio_Claro_Gabon_N*Rio_Claro_Gabon_split,
                                     split_year=2*Rio_Claro_Gabon_N*Rio_Claro_Gabon_split*0.067,
                                     pair="Rio_Claro-Gabon",
                                     region="Americas-Africa")

Rio_Claro_Kenya_json <- fromJSON("Rio_Claro-Kenya.model.final.json")
Rio_Claro_Kenya_split <- as.numeric(Rio_Claro_Kenya_json$model$split)
Rio_Claro_Kenya_N <- as.numeric(Rio_Claro_Kenya_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_Kenya_N,
                                     split_model=Rio_Claro_Kenya_split,
                                     split_gen=2*Rio_Claro_Kenya_N*Rio_Claro_Kenya_split,
                                     split_year=2*Rio_Claro_Kenya_N*Rio_Claro_Kenya_split*0.067,
                                     pair="Rio_Claro-Kenya",
                                     region="Americas-Africa")

Rio_Claro_Senegal_json <- fromJSON("Rio_Claro-Senegal.model.final.json")
Rio_Claro_Senegal_split <- as.numeric(Rio_Claro_Senegal_json$model$split)
Rio_Claro_Senegal_N <- as.numeric(Rio_Claro_Senegal_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_Senegal_N,
                                     split_model=Rio_Claro_Senegal_split,
                                     split_gen=2*Rio_Claro_Senegal_N*Rio_Claro_Senegal_split,
                                     split_year=2*Rio_Claro_Senegal_N*Rio_Claro_Senegal_split*0.067,
                                     pair="Rio_Claro-Senegal",
                                     region="Americas-Africa")

Rio_Claro_USA_json <- fromJSON("Rio_Claro-USA.model.final.json")
Rio_Claro_USA_split <- as.numeric(Rio_Claro_USA_json$model$split)
Rio_Claro_USA_N <- as.numeric(Rio_Claro_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Rio_Claro_USA_N,
                                     split_model=Rio_Claro_USA_split,
                                     split_gen=2*Rio_Claro_USA_N*Rio_Claro_USA_split,
                                     split_year=2*Rio_Claro_USA_N*Rio_Claro_USA_split*0.067,
                                     pair="Rio_Claro-USA",
                                     region="Americas")

Senegal_Gabon_json <- fromJSON("Senegal-Gabon.model.final.json")
Senegal_Gabon_split <- as.numeric(Senegal_Gabon_json$model$split)
Senegal_Gabon_N <- as.numeric(Senegal_Gabon_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Senegal_Gabon_N,
                                     split_model=Senegal_Gabon_split,
                                     split_gen=2*Senegal_Gabon_N*Senegal_Gabon_split,
                                     split_year=2*Senegal_Gabon_N*Senegal_Gabon_split*0.067,
                                     pair="Senegal-Gabon",
                                     region="Africa")

Senegal_USA_json <- fromJSON("Senegal-USA.model.final.json")
Senegal_USA_split <- as.numeric(Senegal_USA_json$model$split)
Senegal_USA_N <- as.numeric(Senegal_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Senegal_USA_N,
                                     split_model=Senegal_USA_split,
                                     split_gen=2*Senegal_USA_N*Senegal_USA_split,
                                     split_year=2*Senegal_USA_N*Senegal_USA_split*0.067,
                                     pair="Senegal-USA",
                                     region="Americas-Africa")
Gabon_USA_json <- fromJSON("Gabon-USA.model.final.json")
Gabon_USA_split <- as.numeric(Gabon_USA_json$model$split)
Gabon_USA_N <- as.numeric(Gabon_USA_json$model$model1$N0)
all_splits <- all_splits %>% add_row(N=Gabon_USA_N,
                                     split_model=Gabon_USA_split,
                                     split_gen=2*Gabon_USA_N*Gabon_USA_split,
                                     split_year=2*Gabon_USA_N*Gabon_USA_split*0.067,
                                     pair="Gabon-USA",
                                     region="Americas-Africa")

all_splits$pair <- all_splits$pair %>% str_replace("Rio_Claro","Río Claro")
  

all_splits_plot <- all_splits %>% arrange(region, pair) %>%
  mutate(pair = fct_inorder(pair)) %>%
  ggplot(aes(x=split_year,y=pair,col=region)) +
  geom_point(size=2) +
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log(),limits=c(0.068,750000)) +
  theme_bw() +
  xlab("Split Year")+
  scale_color_manual(values=c('forestgreen','goldenrod2','tan3'))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

top_row_smc <- plot_grid(africas_smc,americas_smc,nrow=1,labels = c('A','B'),hjust = -1,align = 'vh')
top_row_smc_legend <- plot_grid(top_row_smc,smc_legend,rel_widths = c(2, .4))
bottom_row_smc <- plot_grid(all_splits_plot,nrow=1,labels=c("C"))
plot_grid(top_row_smc_legend,bottom_row_smc,
          align = 'v',
          ncol = 1, 
          axis = 'b')


#USA
us_splits <- data.frame(N= numeric(0), 
                         split_model= numeric(0), 
                         split_gen = numeric(0), 
                         split_year = numeric(0), 
                         pair = character(0),
                         region = character(0))

NorCal_Cali_json <- fromJSON("NorCal_Cali.model.final.json")
NorCal_Cali_split <- as.numeric(NorCal_Cali_json$model$split)
NorCal_Cali_N <- as.numeric(NorCal_Cali_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Cali_N,
                                     split_model=NorCal_Cali_split,
                                     split_gen=2*NorCal_Cali_N*NorCal_Cali_split,
                                     split_year=2*NorCal_Cali_N*NorCal_Cali_split*0.067,
                                     pair="NorCal_Cali",
                                     region="NA-SA")

NorCal_Brazil_json <- fromJSON("NorCal-Brazil.model.final.json")
NorCal_Brazil_split <- as.numeric(NorCal_Brazil_json$model$split)
NorCal_Brazil_N <- as.numeric(NorCal_Brazil_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Brazil_N,
                                   split_model=NorCal_Brazil_split,
                                   split_gen=2*NorCal_Brazil_N*NorCal_Brazil_split,
                                   split_year=2*NorCal_Brazil_N*NorCal_Brazil_split*0.067,
                                   pair="NorCal_Brazil",
                                   region="NA-SA")

NorCal_Senegal_json <- fromJSON("NorCal-Senegal.model.final.json")
NorCal_Senegal_split <- as.numeric(NorCal_Senegal_json$model$split)
NorCal_Senegal_N <- as.numeric(NorCal_Senegal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Senegal_N,
                                   split_model=NorCal_Senegal_split,
                                   split_gen=2*NorCal_Senegal_N*NorCal_Senegal_split,
                                   split_year=2*NorCal_Senegal_N*NorCal_Senegal_split*0.067,
                                   pair="NorCal_Senegal",
                                   region="NA-Africa")

NorCal_Kenya_json <- fromJSON("NorCal-Kenya.model.final.json")
NorCal_Kenya_split <- as.numeric(NorCal_Kenya_json$model$split)
NorCal_Kenya_N <- as.numeric(NorCal_Kenya_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Kenya_N,
                                   split_model=NorCal_Kenya_split,
                                   split_gen=2*NorCal_Kenya_N*NorCal_Kenya_split,
                                   split_year=2*NorCal_Kenya_N*NorCal_Kenya_split*0.067,
                                   pair="NorCal_Kenya",
                                   region="NA-Africa")

NorCal_Gabon_json <- fromJSON("NorCal-Gabon.model.final.json")
NorCal_Gabon_split <- as.numeric(NorCal_Gabon_json$model$split)
NorCal_Gabon_N <- as.numeric(NorCal_Gabon_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Gabon_N,
                                   split_model=NorCal_Gabon_split,
                                   split_gen=2*NorCal_Gabon_N*NorCal_Gabon_split,
                                   split_year=2*NorCal_Gabon_N*NorCal_Gabon_split*0.067,
                                   pair="NorCal_Gabon",
                                   region="NA-Africa")

SoCal_Cali_json <- fromJSON("SoCal-Cali.model.final.json")
SoCal_Cali_split <- as.numeric(SoCal_Cali_json$model$split)
SoCal_Cali_N <- as.numeric(SoCal_Cali_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Cali_N,
                                   split_model=SoCal_Cali_split,
                                   split_gen=2*SoCal_Cali_N*SoCal_Cali_split,
                                   split_year=2*SoCal_Cali_N*SoCal_Cali_split*0.067,
                                   pair="SoCal_Cali",
                                   region="NA-SA")

SoCal_Brazil_json <- fromJSON("SoCal-Brazil.model.final.json")
SoCal_Brazil_split <- as.numeric(SoCal_Brazil_json$model$split)
SoCal_Brazil_N <- as.numeric(SoCal_Brazil_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Brazil_N,
                                   split_model=SoCal_Brazil_split,
                                   split_gen=2*SoCal_Brazil_N*SoCal_Brazil_split,
                                   split_year=2*SoCal_Brazil_N*SoCal_Brazil_split*0.067,
                                   pair="SoCal_Brazil",
                                   region="NA-SA")

SoCal_Kenya_json <- fromJSON("SoCal-Kenya.model.final.json")
SoCal_Kenya_split <- as.numeric(SoCal_Kenya_json$model$split)
SoCal_Kenya_N <- as.numeric(SoCal_Kenya_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Kenya_N,
                                   split_model=SoCal_Kenya_split,
                                   split_gen=2*SoCal_Kenya_N*SoCal_Kenya_split,
                                   split_year=2*SoCal_Kenya_N*SoCal_Kenya_split*0.067,
                                   pair="SoCal_Kenya",
                                   region="NA-Africa")

Florida_Cali_json <- fromJSON("Florida-Cali.model.final.json")
Florida_Cali_split <- as.numeric(Florida_Cali_json$model$split)
Florida_Cali_N <- as.numeric(Florida_Cali_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Cali_N,
                                   split_model=Florida_Cali_split,
                                   split_gen=2*Florida_Cali_N*Florida_Cali_split,
                                   split_year=2*Florida_Cali_N*Florida_Cali_split*0.067,
                                   pair="Florida_Cali",
                                   region="NA-SA")

Florida_Brazil_json <- fromJSON("Florida-Brazil.model.final.json")
Florida_Brazil_split <- as.numeric(Florida_Brazil_json$model$split)
Florida_Brazil_N <- as.numeric(Florida_Brazil_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Brazil_N,
                                   split_model=Florida_Brazil_split,
                                   split_gen=2*Florida_Brazil_N*Florida_Brazil_split,
                                   split_year=2*Florida_Brazil_N*Florida_Brazil_split*0.067,
                                   pair="Florida_Brazil",
                                   region="NA-SA")

Florida_Senegal_json <- fromJSON("Florida-Senegal.model.final.json")
Florida_Senegal_split <- as.numeric(Florida_Senegal_json$model$split)
Florida_Senegal_N <- as.numeric(Florida_Senegal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Senegal_N,
                                   split_model=Florida_Senegal_split,
                                   split_gen=2*Florida_Senegal_N*Florida_Senegal_split,
                                   split_year=2*Florida_Senegal_N*Florida_Senegal_split*0.067,
                                   pair="Florida_Senegal",
                                   region="NA-Africa")

Florida_Kenya_json <- fromJSON("Florida-Kenya.model.final.json")
Florida_Kenya_split <- as.numeric(Florida_Kenya_json$model$split)
Florida_Kenya_N <- as.numeric(Florida_Kenya_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Kenya_N,
                                   split_model=Florida_Kenya_split,
                                   split_gen=2*Florida_Kenya_N*Florida_Kenya_split,
                                   split_year=2*Florida_Kenya_N*Florida_Kenya_split*0.067,
                                   pair="Florida_Kenya",
                                   region="NA-Africa")

Florida_Gabon_json <- fromJSON("Clovis-Cali.model.final.json")
Florida_Gabon_split <- as.numeric(Florida_Gabon_json$model$split)
Florida_Gabon_N <- as.numeric(Florida_Gabon_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Gabon_N,
                                   split_model=Florida_Gabon_split,
                                   split_gen=2*Florida_Gabon_N*Florida_Gabon_split,
                                   split_year=2*Florida_Gabon_N*Florida_Gabon_split*0.067,
                                   pair="Florida_Gabon",
                                   region="NA-Africa")

Clovis_Cali_json <- fromJSON("Clovis-Cali.model.final.json")
Clovis_Cali_split <- as.numeric(Clovis_Cali_json$model$split)
Clovis_Cali_N <- as.numeric(Clovis_Cali_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Cali_N,
                                   split_model=Clovis_Cali_split,
                                   split_gen=2*Clovis_Cali_N*Clovis_Cali_split,
                                   split_year=2*Clovis_Cali_N*Clovis_Cali_split*0.067,
                                   pair="Clovis_Cali",
                                   region="NA-SA")

Clovis_Brazil_json <- fromJSON("Clovis-Brazil.model.final.json")
Clovis_Brazil_split <- as.numeric(Clovis_Brazil_json$model$split)
Clovis_Brazil_N <- as.numeric(Clovis_Brazil_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Brazil_N,
                                   split_model=Clovis_Brazil_split,
                                   split_gen=2*Clovis_Brazil_N*Clovis_Brazil_split,
                                   split_year=2*Clovis_Brazil_N*Clovis_Brazil_split*0.067,
                                   pair="Clovis_Brazil",
                                   region="NA-SA")

Clovis_Senegal_json <- fromJSON("Clovis-Senegal.model.final.json")
Clovis_Senegal_split <- as.numeric(Clovis_Senegal_json$model$split)
Clovis_Senegal_N <- as.numeric(Clovis_Senegal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Senegal_N,
                                   split_model=Clovis_Senegal_split,
                                   split_gen=2*Clovis_Senegal_N*Clovis_Senegal_split,
                                   split_year=2*Clovis_Senegal_N*Clovis_Senegal_split*0.067,
                                   pair="Clovis_Senegal",
                                   region="NA-Africa")

Clovis_Kenya_json <- fromJSON("Clovis-Kenya.model.final.json")
Clovis_Kenya_split <- as.numeric(Clovis_Kenya_json$model$split)
Clovis_Kenya_N <- as.numeric(Clovis_Kenya_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Kenya_N,
                                   split_model=Clovis_Kenya_split,
                                   split_gen=2*Clovis_Kenya_N*Clovis_Kenya_split,
                                   split_year=2*Clovis_Kenya_N*Clovis_Kenya_split*0.067,
                                   pair="Clovis_Kenya",
                                   region="NA-Africa")

Clovis_Gabon_json <- fromJSON("Clovis-Gabon.model.final.json")
Clovis_Gabon_split <- as.numeric(Clovis_Gabon_json$model$split)
Clovis_Gabon_N <- as.numeric(Clovis_Gabon_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Gabon_N,
                                   split_model=Clovis_Gabon_split,
                                   split_gen=2*Clovis_Gabon_N*Clovis_Gabon_split,
                                   split_year=2*Clovis_Gabon_N*Clovis_Gabon_split*0.067,
                                   pair="Clovis_Gabon",
                                   region="NA-Africa")

Exeter_Cali_json <- fromJSON("Exeter-Cali.model.final.json")
Exeter_Cali_split <- as.numeric(Exeter_Cali_json$model$split)
Exeter_Cali_N <- as.numeric(Exeter_Cali_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Cali_N,
                                   split_model=Exeter_Cali_split,
                                   split_gen=2*Exeter_Cali_N*Exeter_Cali_split,
                                   split_year=2*Exeter_Cali_N*Exeter_Cali_split*0.067,
                                   pair="Exeter_Cali",
                                   region="NA-SA")

Exeter_Brazil_json <- fromJSON("Exeter-Brazil.model.final.json")
Exeter_Brazil_split <- as.numeric(Exeter_Brazil_json$model$split)
Exeter_Brazil_N <- as.numeric(Exeter_Brazil_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Brazil_N,
                                   split_model=Exeter_Brazil_split,
                                   split_gen=2*Exeter_Brazil_N*Exeter_Brazil_split,
                                   split_year=2*Exeter_Brazil_N*Exeter_Brazil_split*0.067,
                                   pair="Exeter_Brazil",
                                   region="NA-SA")

Exeter_Senegal_json <- fromJSON("Exeter-Senegal.model.final.json")
Exeter_Senegal_split <- as.numeric(Exeter_Senegal_json$model$split)
Exeter_Senegal_N <- as.numeric(Exeter_Senegal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Senegal_N,
                                   split_model=Exeter_Senegal_split,
                                   split_gen=2*Exeter_Senegal_N*Exeter_Senegal_split,
                                   split_year=2*Exeter_Senegal_N*Exeter_Senegal_split*0.067,
                                   pair="Exeter_Senegal",
                                   region="NA-Africa")

Exeter_Kenya_json <- fromJSON("Exeter-Kenya.model.final.json")
Exeter_Kenya_split <- as.numeric(Exeter_Kenya_json$model$split)
Exeter_Kenya_N <- as.numeric(Exeter_Kenya_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Kenya_N,
                                   split_model=Exeter_Kenya_split,
                                   split_gen=2*Exeter_Kenya_N*Exeter_Kenya_split,
                                   split_year=2*Exeter_Kenya_N*Exeter_Kenya_split*0.067,
                                   pair="Exeter_Kenya",
                                   region="NA-Africa")

Exeter_Gabon_json <- fromJSON("Exeter-Gabon.model.final.json")
Exeter_Gabon_split <- as.numeric(Exeter_Gabon_json$model$split)
Exeter_Gabon_N <- as.numeric(Exeter_Gabon_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Gabon_N,
                                   split_model=Exeter_Gabon_split,
                                   split_gen=2*Exeter_Gabon_N*Exeter_Gabon_split,
                                   split_year=2*Exeter_Gabon_N*Exeter_Gabon_split*0.067,
                                   pair="Exeter_Gabon",
                                   region="NA-Africa")

SoCal_Clovis_json <- fromJSON("SoCal-Clovis.model.final.json")
SoCal_Clovis_split <- as.numeric(SoCal_Clovis_json$model$split)
SoCal_Clovis_N <- as.numeric(SoCal_Clovis_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Clovis_N,
                                   split_model=SoCal_Clovis_split,
                                   split_gen=2*SoCal_Clovis_N*SoCal_Clovis_split,
                                   split_year=2*SoCal_Clovis_N*SoCal_Clovis_split*0.067,
                                   pair="SoCal_Clovis",
                                   region="NA")

SoCal_Florida_json <- fromJSON("SoCal-Florida.model.final.json")
SoCal_Florida_split <- as.numeric(SoCal_Florida_json$model$split)
SoCal_Florida_N <- as.numeric(SoCal_Florida_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Florida_N,
                                   split_model=SoCal_Florida_split,
                                   split_gen=2*SoCal_Florida_N*SoCal_Florida_split,
                                   split_year=2*SoCal_Florida_N*SoCal_Florida_split*0.067,
                                   pair="SoCal_Florida",
                                   region="NA")

SoCal_Exeter_json <- fromJSON("SoCal-Exeter.model.final.json")
SoCal_Exeter_split <- as.numeric(SoCal_Exeter_json$model$split)
SoCal_Exeter_N <- as.numeric(SoCal_Exeter_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Exeter_N,
                                   split_model=SoCal_Exeter_split,
                                   split_gen=2*SoCal_Exeter_N*SoCal_Exeter_split,
                                   split_year=2*SoCal_Exeter_N*SoCal_Exeter_split*0.067,
                                   pair="SoCal_Exeter",
                                   region="NA")

NorCal_Clovis_json <- fromJSON("NorCal-Clovis.model.final.json")
NorCal_Clovis_split <- as.numeric(NorCal_Clovis_json$model$split)
NorCal_Clovis_N <- as.numeric(NorCal_Clovis_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Clovis_N,
                                   split_model=NorCal_Clovis_split,
                                   split_gen=2*NorCal_Clovis_N*NorCal_Clovis_split,
                                   split_year=2*NorCal_Clovis_N*NorCal_Clovis_split*0.067,
                                   pair="NorCal_Clovis",
                                   region="NA")

NorCal_Florida_json <- fromJSON("NorCal-Florida.model.final.json")
NorCal_Florida_split <- as.numeric(NorCal_Florida_json$model$split)
NorCal_Florida_N <- as.numeric(NorCal_Florida_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Florida_N,
                                   split_model=NorCal_Florida_split,
                                   split_gen=2*NorCal_Florida_N*NorCal_Florida_split,
                                   split_year=2*NorCal_Florida_N*NorCal_Florida_split*0.067,
                                   pair="NorCal_Florida",
                                   region="NA")

NorCal_Exeter_json <- fromJSON("NorCal-Exeter.model.final.json")
NorCal_Exeter_split <- as.numeric(NorCal_Exeter_json$model$split)
NorCal_Exeter_N <- as.numeric(NorCal_Exeter_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Exeter_N,
                                   split_model=NorCal_Exeter_split,
                                   split_gen=2*NorCal_Exeter_N*NorCal_Exeter_split,
                                   split_year=2*NorCal_Exeter_N*NorCal_Exeter_split*0.067,
                                   pair="NorCal_Exeter",
                                   region="NA")

NorCal_SoCal_json <- fromJSON("NorCal-SoCal.model.final.json")
NorCal_SoCal_split <- as.numeric(NorCal_SoCal_json$model$split)
NorCal_SoCal_N <- as.numeric(NorCal_SoCal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_SoCal_N,
                                   split_model=NorCal_SoCal_split,
                                   split_gen=2*NorCal_SoCal_N*NorCal_SoCal_split,
                                   split_year=2*NorCal_SoCal_N*NorCal_SoCal_split*0.067,
                                   pair="NorCal_SoCal",
                                   region="NA")

Florida_Exeter_json <- fromJSON("Florida-Exeter.model.final.json")
Florida_Exeter_split <- as.numeric(Florida_Exeter_json$model$split)
Florida_Exeter_N <- as.numeric(Florida_Exeter_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Exeter_N,
                                   split_model=Florida_Exeter_split,
                                   split_gen=2*Florida_Exeter_N*Florida_Exeter_split,
                                   split_year=2*Florida_Exeter_N*Florida_Exeter_split*0.067,
                                   pair="Florida_Exeter",
                                   region="NA")

Florida_Clovis_json <- fromJSON("Florida-Clovis.model.final.json")
Florida_Clovis_split <- as.numeric(Florida_Clovis_json$model$split)
Florida_Clovis_N <- as.numeric(Florida_Clovis_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Clovis_N,
                                   split_model=Florida_Clovis_split,
                                   split_gen=2*Florida_Clovis_N*Florida_Clovis_split,
                                   split_year=2*Florida_Clovis_N*Florida_Clovis_split*0.067,
                                   pair="Florida_Clovis",
                                   region="NA")

Florida_Exeter_json <- fromJSON("Florida-Exeter.model.final.json")
Florida_Exeter_split <- as.numeric(Florida_Exeter_json$model$split)
Florida_Exeter_N <- as.numeric(Florida_Exeter_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Exeter_N,
                                   split_model=Florida_Exeter_split,
                                   split_gen=2*Florida_Exeter_N*Florida_Exeter_split,
                                   split_year=2*Florida_Exeter_N*Florida_Exeter_split*0.067,
                                   pair="Florida_Exeter",
                                   region="NA")

Clovis_Exeter_json <- fromJSON("Clovis-Exeter.model.final.json")
Clovis_Exeter_split <- as.numeric(Clovis_Exeter_json$model$split)
Clovis_Exeter_N <- as.numeric(Clovis_Exeter_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Exeter_N,
                                   split_model=Clovis_Exeter_split,
                                   split_gen=2*Clovis_Exeter_N*Clovis_Exeter_split,
                                   split_year=2*Clovis_Exeter_N*Clovis_Exeter_split*0.067,
                                   pair="Clovis_Exeter",
                                   region="NA")

SoCal_Rio_Claro_json <- fromJSON("SoCal-Rio_Claro.model.final.json")
SoCal_Rio_Claro_split <- as.numeric(SoCal_Rio_Claro_json$model$split)
SoCal_Rio_Claro_N <- as.numeric(SoCal_Rio_Claro_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Rio_Claro_N,
                                   split_model=SoCal_Rio_Claro_split,
                                   split_gen=2*SoCal_Rio_Claro_N*SoCal_Rio_Claro_split,
                                   split_year=2*SoCal_Rio_Claro_N*SoCal_Rio_Claro_split*0.067,
                                   pair="SoCal_Rio_Claro",
                                   region="NA-SA")
NorCal_Rio_Claro_json <- fromJSON("NorCal-Rio_Claro.model.final.json")
NorCal_Rio_Claro_split <- as.numeric(NorCal_Rio_Claro_json$model$split)
NorCal_Rio_Claro_N <- as.numeric(NorCal_Rio_Claro_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=NorCal_Rio_Claro_N,
                                   split_model=NorCal_Rio_Claro_split,
                                   split_gen=2*NorCal_Rio_Claro_N*NorCal_Rio_Claro_split,
                                   split_year=2*NorCal_Rio_Claro_N*NorCal_Rio_Claro_split*0.067,
                                   pair="NorCal_Rio_Claro",
                                   region="NA-SA")
Florida_Rio_Claro_json <- fromJSON("Florida-Rio_Claro.model.final.json")
Florida_Rio_Claro_split <- as.numeric(Florida_Rio_Claro_json$model$split)
Florida_Rio_Claro_N <- as.numeric(Florida_Rio_Claro_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Florida_Rio_Claro_N,
                                   split_model=Florida_Rio_Claro_split,
                                   split_gen=2*Florida_Rio_Claro_N*Florida_Rio_Claro_split,
                                   split_year=2*Florida_Rio_Claro_N*Florida_Rio_Claro_split*0.067,
                                   pair="Florida_Rio_Claro",
                                   region="NA-SA")
Exeter_Rio_Claro_json <- fromJSON("Exeter-Rio_Claro.model.final.json")
Exeter_Rio_Claro_split <- as.numeric(Exeter_Rio_Claro_json$model$split)
Exeter_Rio_Claro_N <- as.numeric(Exeter_Rio_Claro_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Exeter_Rio_Claro_N,
                                   split_model=Exeter_Rio_Claro_split,
                                   split_gen=2*Exeter_Rio_Claro_N*Exeter_Rio_Claro_split,
                                   split_year=2*Exeter_Rio_Claro_N*Exeter_Rio_Claro_split*0.067,
                                   pair="Exeter_Rio_Claro",
                                   region="NA-SA")
Clovis_Rio_Claro_json <- fromJSON("Clovis-Rio_Claro.model.final.json")
Clovis_Rio_Claro_split <- as.numeric(Clovis_Rio_Claro_json$model$split)
Clovis_Rio_Claro_N <- as.numeric(Clovis_Rio_Claro_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=Clovis_Rio_Claro_N,
                                   split_model=Clovis_Rio_Claro_split,
                                   split_gen=2*Clovis_Rio_Claro_N*Clovis_Rio_Claro_split,
                                   split_year=2*Clovis_Rio_Claro_N*Clovis_Rio_Claro_split*0.067,
                                   pair="Clovis_Rio_Claro",
                                   region="NA-SA")
SoCal_Senegal_json <- fromJSON("SoCal-Senegal.model.final.json")
SoCal_Senegal_split <- as.numeric(SoCal_Senegal_json$model$split)
SoCal_Senegal_N <- as.numeric(SoCal_Senegal_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Senegal_N,
                                   split_model=SoCal_Senegal_split,
                                   split_gen=2*SoCal_Senegal_N*SoCal_Senegal_split,
                                   split_year=2*SoCal_Senegal_N*SoCal_Senegal_split*0.067,
                                   pair="SoCal_Senegal",
                                   region="NA-Africa")
SoCal_Gabon_json <- fromJSON("SoCal-Gabon.model.final.json")
SoCal_Gabon_split <- as.numeric(SoCal_Gabon_json$model$split)
SoCal_Gabon_N <- as.numeric(SoCal_Gabon_json$model$model1$N0)
us_splits <- us_splits %>% add_row(N=SoCal_Gabon_N,
                                   split_model=SoCal_Gabon_split,
                                   split_gen=2*SoCal_Gabon_N*SoCal_Gabon_split,
                                   split_year=2*SoCal_Gabon_N*SoCal_Gabon_split*0.067,
                                   pair="SoCal_Gabon",
                                   region="NA-Africa")

us_splits$pair <- us_splits$pair %>% str_replace("Rio_Claro","Río Claro") %>% str_replace("_","-")

us_splits_plot <- us_splits %>% arrange(region, pair) %>%
  mutate(pair = fct_inorder(pair)) %>%
  ggplot(aes(x=split_year,y=pair,col=region)) +
  geom_point(size=2) +
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log(),limits=c(0.068,750000)) +
  theme_bw() +
  xlab("Split Year")+
  scale_color_manual(values=c('forestgreen','goldenrod2','tan3'))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=9),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))

plot_grid(full_us_100k_rp5_plot,us_splits_plot,
          align = 'v',
          ncol = 1, 
          axis = 'b',
          labels = c("A","B"),
          rel_heights = c(0.8,1))




#------------------------------------#
# MSMC2 within
#------------------------------------#
setwd('~/aedes/Results/msmc2/')

gen=0.067
mu=4.85e-9
Brazil_msmc <- fread("Brazil_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Brazil")

Rio_Claro_msmc <- fread("Rio_Claro_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Rio_Claro")

Cali_msmc <- fread("Cali_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Cali")

USA_msmc <- fread("USA_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="USA")

Kenya_msmc <- fread("Kenya_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Kenya")

Senegal_msmc <- fread("Senegal_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Senegal")

Gabon_msmc <- fread("Gabon_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Gabon")

Clovis_msmc <- fread("Clovis_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Clovis")

Florida_msmc <- fread("Florida_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Florida")

NorCal_msmc <- fread("NorCal_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="NorCal")

SoCal_msmc <- fread("SoCal_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="SoCal")

Exeter_msmc <- fread("Exeter_within_msmc.final.txt",data.table=F) %>%
  mutate(time=(gen*left_time_boundary/mu)+3, ne=(1/lambda)/(2*mu),pop="Exeter")

all_msmc <- rbind(Brazil_msmc,Rio_Claro_msmc,Cali_msmc,Kenya_msmc,Senegal_msmc,Gabon_msmc,Clovis_msmc,Florida_msmc,NorCal_msmc,SoCal_msmc,Exeter_msmc)

pal1=tanagr_palette("chlorochrysa_nitidissima",n = 7,discrete = F)[1:6]
pal2=pal=met.brewer("OKeeffe1", 6, type="continuous")[1:5]
pal=c("#D4940E","#976407","#6B200C","#DA6C42","#FBC2A9","#36260A","#062C3B","#BAD6F9","#00504F","#286F29","#447FDD")
all_msmc %>% ggplot(aes(x=time,y=ne,col=pop)) +
  geom_line(size=1.3) + 
  scale_x_log10(breaks = breaks_log(n=8),labels = label_log(),limits=c(0.068,1500000)) +
  scale_y_log10(breaks = breaks_log(n=4),labels = label_log()) +
  scale_color_manual(values=pal)+
  ylab('Ne (x 1k individuals)') +
  xlab('Years ago (x 1k)') +
  theme_classic()
