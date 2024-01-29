library("tidyverse")
library("viridis")

the_data <- read_delim("histogram_data_ipd_20240129_102312008565_0_whole_pop",delim=";")

# make histogram of xp
ggplot(data = the_data %>% filter(trait == "xp")
       ,mapping = aes(x=time,y=midpoints)) +
    geom_tile(mapping=aes(fill=count)) +
    scale_fill_viridis(option="turbo")
    