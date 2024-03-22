#!/usr/bin/env Rscript 
library("tidyverse")

library("data.table")

filename <- commandArgs(trailingOnly = T)[1]

nbin <- 100

# obtain the data using the ultrafast fread() function
the_data <- fread(input=filename)

# first column is the column specifying time
# all the other columns are traits. 
all_names <- names(the_data)
time_name <- all_names[1]
trait_names <- all_names[2:length(all_names)]

columns_to_omit <- "type"

trait_names <- trait_names[!trait_names %in% c(columns_to_omit)]


# now find the ranges of all the traits
# you can also supply those manually in case you aren't too happy with
# the ranges that are supplied
ranges <- data.frame(trait=trait_names,min=NA,max=NA)

for (trait_i in trait_names)
{
    current_range <- range(the_data[, get(trait_i)])
   
    # check whether name only occurs once 
    stopifnot(nrow(ranges[which(ranges$trait == trait_i),]) == 1)
    
    ranges[which(ranges$trait == trait_i),"min"] <- current_range[1]
    ranges[which(ranges$trait == trait_i),"max"] <- current_range[2]
}

## end finding ranges
# now start making histograms for each timestep and for each trait
unique_t <- sort(unique(the_data[, get(time_name)]))

# the full histogram
full_hist <- NULL

for (trait_i in trait_names)
{
    print(paste0("working on trait ", trait_i))
    
    # calculate the breaks for the histogram
    max_val <- ranges[which(ranges$trait == trait_i),"max"]
    min_val <- ranges[which(ranges$trait == trait_i),"min"]
    print(class(min_val))
    if (class(min_val) == class("character data"))
    {
        next
    }
    
    break_values <- seq(min_val,max_val,length.out=100)
    
    for (time_step in unique_t)
    {
        # get a vector with the data
        subset <- the_data[get(time_name) == time_step, get(trait_i)]
        
        # make a histogram
        the_hist <- hist(subset, breaks = break_values, plot=F)
       
        # get the counts 
        counts <- the_hist$counts
        
        # calculate the midpoints
        midpoints <- the_hist$mids
        
        # make a column indicating the trait
        trait <- rep(trait_i,length.out = length(counts))
        
        # make a column indicating the time point
        time <- rep(time_step, length.out = length(counts))
        
        # get it all together in a data.frame
        sub_hist <- data.frame(trait=trait,time=time,midpoints=midpoints,count=counts)
        
        # append data.frame to previous data.frames
        full_hist <- rbind(full_hist,sub_hist)
    }
}

# generate a file name, if the file is nested in folder just get its basename()
# and use that.
output_file_name <- paste0("histogram_",basename(filename))

# write a table to a file
write.table(x=full_hist,file=output_file_name,row.names=F,sep=";")