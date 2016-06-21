library(reshape2)
library(ggplot2)
library(neotoma)
library(rioja)
library(scales)
library(ggmap)
library(maps)
library(maptools)
library(grid)
library(gridExtra)
library(raster)
library(GGally)

# source('utils/ggcorplot.R')

# coordinates for the three sites
# meta = data.frame(site=character(0), lat=numeric(0), long=numeric(0))
# meta = rbind(meta, data.frame('JRI', -64.2017, -57.685), #              
#              -64° 51' 36", -64° 12' 0")

coords = read.csv('data/coords.csv', header=FALSE)
colnames(coords) = c('site', 'source', 'lat', 'long', 'variable', 'type')


#----------- Functions -----------#
make_props <- function(x) {
  if (sum(x) > 0){
    return(x/sum(x))
  } else {
    return(x)
  }
}


#------------------- Mapping -------------#
#get baselayer for main map
myMap <- get_map(location=c(-66, -66, -55, -60), #long/lat, long/lat
                 source="stamen", maptype="watercolor", crop=FALSE)
#plot sample points
myMap = ggmap(myMap) + geom_point(aes(x = -long, y = -lat, shape=type), data = coords, alpha = .5, color="darkred", size = 4) + 
        geom_text(data=coords, aes(x=-long, y=-lat,label=type), hjust=-0.1, vjust=-0.5, size=5) + xlab('Longitude') + ylab('Latitude')

#get world map for inset
world <- map_data("world")

#plot inset map
insetMap = ggplot() + 
           geom_polygon(data = world, aes(x=long, y = lat, group = group))+
           geom_rect(aes(ymin=-66,ymax=-55,xmin=-60,xmax=-55), colour="red", fill="red")+
           coord_cartesian(x=c(-110, -30), y=c(-80, 10))+
           scale_x_discrete(expand=c(0, 0))+
           theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 panel.grid.major=element_blank(),
                 plot.background=element_blank())

#render map with inset map
pdf("figures/proxy_locations.pdf", width=8, height=8)
grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.25, height = 0.3, x = 0.2, y = .845) #plot area for the inset map
print(myMap,vp=v1) 
print(insetMap,vp=v2)
dev.off()

##################################################################################################################################################
## Parse Linked Earth data
##################################################################################################################################################

# linked earth data 
pc1 = read.csv('ODP1098B.data.PaleoData.csv', header=FALSE)
pc2 = pc1[,c(3,5)]
colnames(pc1) = c('age', 'sst')

# noaa data
pc2 = read.table('jr_temp.csv', header=TRUE)
pc2 = pc2[,c(1,4)]
colnames(pc2) = c('year', 'value')

pol = get_download(217)[[1]]
# pol = compile_taxa(pol_dat, list.name='WhitmoreFull')[[1]]
counts = pol$counts

counts[counts <= 0.5]	<- 0 	
props = t(apply(counts, 1, make_props))

props2 = props[,t(apply(props, 2, function(x) any(x > 0.3)))]
props2 = t(apply(props2, 1, make_props))

age = pol$sample.meta$age

pc1 = cbind(pc1, 'sst')
colnames(pc1) = c('year', 'value', 'record')
pc1$value = rescale(pc1$value)

pc2 = cbind(pc2, 'pc2')
colnames(pc2) = c('year', 'value', 'record')
pc2$value = rescale(pc2$value)

pol2 = cbind(rep(age, time=ncol(props2)), melt(props2)[,c(3,2)])
colnames(pol2) = c('year', 'value', 'record')

pol = data.frame(year=age, props2)

dat_all = data.frame(year=numeric(0), value=numeric(0), record=character(0))
dat_all = rbind(dat_all, pc1)
dat_all = rbind(dat_all, pol2)
dat_all = rbind(dat_all, pc2)

ggplot(data=subset(dat_all,year<5000)) + geom_line(aes(x=year, y=value)) + coord_flip() + scale_x_reverse() + facet_grid(~record, scales='free_y')

####################################################################################################################################################
## interpolate using a smoothing spline
####################################################################################################################################################

t_all = seq(0,5000, by=100)

# datsets = c(pc1

interp <- function(dataset, df, t_all){
  ss1 = smooth.spline(dataset$year, dataset$value, df=df)
  ss_pred = predict(ss1, x=t_all)
  
  return(ss_pred)
}

pc1_smooth = interp(pc1, 20, t_all)$y
pc2_smooth = interp(pc2, 20, t_all)$y

dat_smooth = data.frame(year=t_all)
ntypes = unique(pol2$record)
for (type in ntypes) {
  dat_sub = pol2[pol2$record %in% type,]
  dat_sub = dat_sub[order(dat_sub$year),]
  
  dat_sub_smooth = data.frame(interp(dat_sub, 20, t_all)$y)
  colnames(dat_sub_smooth) = type
  
  dat_smooth = data.frame(dat_smooth, dat_sub_smooth)
}

dat_smooth = data.frame(dat_smooth, SST = pc1_smooth, TS = pc2_smooth)

# func = splinefun(x=pal$year, y=pal$value, method="fmm",  ties = mean)
# plot(pal$year, pal$value, type='l', xlim=c(0,5000))
# lines(smooth.spline(pal$year, pal$value, df=50))

colnames(dat_smooth) = c('year', 'PSU', 'UNK1', 'UNK5', 'UNK6', 'PSB', 'SST', 'TS')

saveRDS(dat_smooth, file='data/paleo_dat.RDS')

pdf(file='figures/scatter.pdf')
p1 <- ggpairs(dat_smooth[,2:ncol(dat_smooth)]) + theme_bw()  
print(p1, left=.5, bottom=.5)
# ggsave(p1, file='figures/scatter.pdf')
dev.off()

pairs(dat_smooth[,2:ncol(dat_smooth)])

cor_mat = cor(dat_smooth[,2:ncol(dat_smooth)])
melted_cormat = melt(cor_mat)

#plot correlation matrix
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# 
# ss <- smooth.spline(pal$year, pal$value, df=20)
# ss_pred <- predict(ss, x=t_all)
# 
# lines(ss_pred, col='blue')
# 
# lines(seq(0,5000,by=1000), func(seq(0,5000,by=1000)), col='blue')
# 
# 
# # library(ggplot2)
# # library(GGally)
# # library(reshape)
# 
# library(dplR)
# 
# plot(pal$year, pal$value, type='l', xlim=c(0,5000))
# lines(pal$year, ffcsaps(y=pal$value, nyrs = 40, f = 0.5))
