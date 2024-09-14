library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(pheatmap)
library(ggh4x) #rbase like plots with truncated x-y axis
library(ggthemes) #rbase like boxed-background panel plots
library(ggdendro) #for clustering
library(scales)
library(RColorBrewer)

#rbase like-plots generating ggplot codes

ggplot(mpg, aes(model, hwy, fill=factor(cyl))) + 
    theme_classic() +
    geom_boxplot() + 
    guides(x = "axis_truncated", y= "axis_truncated") + #library(ggh4x)
    

    scale_y_continuous(name = "COOL", breaks=seq(10,50,10), limits=c(0,50), expand=c(0,0)) +
    labs(x="X-axis", y="COOL", fill="NEW LEGEND") + 

    theme(
        axis.ticks.length = unit(4, "pt"), #length of ticks

        axis.ticks = element_line(size = 3, color=c("red","blue","green")), #width of ticks
        #axis.ticks.x = element_line(size = 0.25, color="red"),
        #axis.ticks.y = element_line(size = 0.3),
        #axis.ticks.x = element_blank(), #rbase histogram style
        
        axis.line = element_line(color="black", size=0.25), #width of lines
        #axis.line.x = element_line(color="black", size=0.25),
        #axis.line.y = element_line(color="black", size=0.25),
        #axis.line.x = element_blank(),

        axis.text = element_text(angle=45, vjust = 0.5, hjust=1), #position of texts
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),

        axis.title = element_text(size=12, face="bold", color="black"),
        #axis.title.x=element_text(size=10),
        #axis.title.y=element_text(size=20, color="black"),

        legend.position="top", #bottom #right
        legend.title = element_text(color = "dodgerblue", size = 20),
        legend.text = element_text(color = "darkred"),
        legend.key.size = unit(1.2, "line")
        #legend.key = element_blank(),
        )

#################################################
#testing of dplyr mutate and summarize functions
#with groupby objects    

mpg1 <- mpg %>%
    group_by(cyl) %>%
    mutate(z_score = scale(hwy)) 

mpg1 %>% data.frame %>% head

custom <- function(x) (x-mean(x))/sd(x)
custom1 <- function(x) mean(x)
custom2 <- function(x,y) (mean(x) + sd(y))

mpg2 <- mpg %>%
    group_by(cyl) %>%
    mutate(z_score = custom(hwy)) 

    #for multiple arguments
    #mutate(z_score = custom(hwy, displ, cty))
mpg2 %>% data.frame %>% head


mpgtest = mpg 
mpgtest$z <- ave(mpgtest$hwy, mpgtest$cyl, FUN=scale)
mpgtest %>% data.frame %>% head

mpgnew = mpg 
mpgnew$z <- ave(mpgnew$hwy, mpgnew$cyl, FUN=function(x) (x-mean(x))/sd(x))
mpgnew %>% data.frame %>% head

all(mpgtest$z == mpg1$z_score)
all(mpgnew$z == mpg1$z_score)
all(mpgnew$z == mpg2$z_score)

# unload tidyverse and ggplot2 as an example first
detach("package:plyr")    #to solve dplyr group_by issues()
library(dplyr)
detach("package:tidyverse", unload = TRUE)
detach("package:ggplot2", unload = TRUE)
options(tidyverse.quiet = TRUE)
#library(tidyverse)

mpgfinal <- mpg %>% group_by(cyl) %>% 
 summarise_at(vars(cty, hwy, displ),
      list(mean = mean, sd = sd, min = min),
      na.rm = TRUE)

mpgfinal1 <- mpg  %>% group_by(cyl) %>%
    summarise(mean=mean(cty), sd= sd(hwy), min=min(displ), .groups = 'drop') %>% 
    ungroup() %>% 
    data.frame

#custom functions based summarisation 
custom1 <- function(x) mean(x)
custom2 <- function(x,y) (mean(x) + sd(y))
custom3 <- function(x,y,z) (mean(x) + sd(y))/min(z)*100 # equiv to ((mean(x) + sd(y))/min(z))*100

mpgfinal2 <- mpg  %>% group_by(cyl) %>%
    summarise(
        mean=mean(cty, na.rm = TRUE), sd= sd(hwy), min=min(displ), 
        zscore= custom1(cty), zscore2=custom2(cty,hwy),
        zscore3percent = custom3(cty, hwy, displ),
        number_counts = n(),
        .groups = 'drop') %>% 
    data.frame


# library(ggh4x)
# library(ggdendro)
# library(scales)
# library(RColorBrewer)
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord
[1]  2  3  1  4  8  5  6 10  7  9

# # The only thing you have to do then is transforming your Time-column to a factor where the factor levels are ordered by ord:

# pd <- as.data.frame( data )
# pd$Time <- sub("_.*", "", rownames(pd))
pd.m <- melt( pd, id.vars = "Time", variable.name = "Gene" )

# pd.m$Gene <- factor( pd.m$Gene, levels = colnames(data), labels = seq_along( colnames(data) ) )
pd.m$Time <- factor( pd.m$Time, levels = rownames(data)[ord],  labels = c("0h", "0.25h", "0.5h","1h","2h","3h","6h","12h","24h","48h") )
# # The rest is done by ggplot automatically:

# ggplot( pd.m, aes(Time, Gene) ) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient2(low=muted("blue"), high=muted("red"))


highlight <- c("New York", "California", "Alabama", "Hawaii")
clust <- hclust(dist(USArrests))
melt_df <- melt(USArrests, id.vars="state")
names(melt_df) <- c("state", "crime", "value")

melt_df1 <- melt_df %>% group_by(crime)  %>% mutate(zscore= scale(value))
ggplot(melt_df, aes(crime, state, fill = value)) +
    theme_minimal() +
  geom_raster() +
  scale_fill_viridis_c() +

  #scale_y_dendrogram(hclust = clust, labels = NULL) +
  scale_y_dendrogram(hclust = clust) +
  #guides(y.sec = guide_axis_manual(breaks = highlight, labels = highlight))
  guides(y.sec = guide_axis_manual(breaks = melt_df$state, labels = melt_df$state))

###############
#plot on zscore
zscore_df <- melt_df1 %>% dcast(state~crime, value.var="zscore")
rownames(zscore_df) <- zscore_df$state 
zscore_finaldf <- zscore_df %>% select(-state)

#clustering of axis
zscoreclust_xaxis <- hclust(dist(zscore_finaldf))
zscoreclust_yaxis <- hclust(dist(t(zscore_finaldf)))

#extract the order by $order
zscore_finaldf[zscoreclust$order,]

ggplot(melt_df1, aes(crime, state, fill = zscore)) +
    theme_classic() +
  geom_raster() +
  scale_fill_viridis_c() +

  #scale_y_dendrogram(hclust = clust, labels = NULL) +
  scale_y_dendrogram(hclust = zscoreclust) +
  #guides(y.sec = guide_axis_manual(breaks = highlight, labels = highlight))
  guides(y.sec = guide_axis_manual(breaks = melt_df$state, labels = NULL)) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
    )


#######################
#Final tested version
#######################

#different color themes
# plot using a palette from RColorBrewer:

# ggplot(data =  dat, aes(x = Row, y = Col)) + 
#       geom_tile(aes(fill = Y1), colour = "white") +
#       scale_fill_brewer(palette = "PRGn")

library(RColorBrewer)
color_theme <- colorRampPalette(rev(brewer.pal(n = 7, name =
       "RdYlBu")))(100)

# yclus <- hclust(dist(USArrests), "ave")
# xclus <- hclust(dist(t(USArrests)), "ave")

melt_df <- melt(USArrests, id.vars="state")
names(melt_df) <- c("state", "crime", "value")
melt_df1 <- melt_df %>% group_by(crime)  %>% mutate(zscore= scale(value))

zscore_df <- melt_df1 %>% dcast(state~crime, value.var="zscore")
rownames(zscore_df) <- zscore_df$state 
zscore_finaldf <- zscore_df %>% select(-state)

#clustering of axis
zscoreclust_xaxis <- hclust(dist(zscore_finaldf))
zscoreclust_yaxis <- hclust(dist(t(zscore_finaldf)))

ggplot(melt_df1, aes(crime, state, fill = zscore)) +
    
    theme_classic() +

    geom_raster() +

    #scale_fill_viridis_c() +
    #scale_fill_gradient(low="yellow", high="red") +
    scale_fill_gradientn(colors=color_theme) +

    scale_y_dendrogram(hclust = zscoreclust_xaxis) +
    #scale_y_dendrogram(hclust = zscoreclust_xaxis, labels=NULL) +

    scale_x_dendrogram(guide = guide_dendro(position = "top"), hclust = zscoreclust_yaxis) +
    #scale_x_dendrogram(guide = guide_dendro(position = "top"), hclust = zscoreclust_yaxis, labels=NULL) +

    ###guides(y.sec = guide_axis_manual(breaks = highlight, labels = highlight))
    guides(
        x.sec = guide_axis_manual(breaks = test, labels = test, position="bottom"),
        y.sec = guide_axis_manual(breaks = melt_df$state, labels = melt_df$state)
    ) +

    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        #axis.text = element_text(angle=0, vjust = 0.5, hjust=1, size=6)
        axis.text.x = element_text(angle=0,size=8),
        axis.text.y = element_text(angle=0,size=6)
    ) +

    theme(
        legend.position="right", #bottom #right
        legend.title = element_text(color = "dodgerblue", size = 6, face="bold"),
        legend.text = element_text(color = "darkred"),
        legend.key.size = unit(1.2, "line")
    ) +
    
    labs(x="X-axis", y="COOL", fill="NEW LEGEND") +

    geom_text(aes(label = round(zscore,1)), color = "black", 
            fontface = "bold", size = 1)
  