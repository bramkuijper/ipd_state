#!/usr/bin/env Rscript

library("tidyverse")
library("patchwork")

#file <- commandArgs(trailingOnly = T)[1]

file <- "test_run.csv"



# find params function
find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    par_list <- list(max_lines=length(f))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            par_list <- c(par_list, list(end_data_line=line_i))

            return(par_list)
        }
    }

    return(NA)
}

param.line <- find.params(filename=file)

print(param.line)

parameters <- NULL

if (with(param.line, max_lines > end_data_line))
{
    # read in the parameters
    parameters.t <- t(read_delim(file=file, delim=";", skip = param.line$end_data_line, col_names = F))

    parameters <- as.numeric(parameters.t[2,])
    names(parameters) <- parameters.t[1,]
}

# read in the data
data <- read_delim(file=file, skip_empty_rows=F, delim=";", n_max = ifelse(is.null(parameters), Inf, param.line$end_data_line - 1))

# selection of column names involving i
i_headers <- names(data)[grep(pattern="^i\\d+", x=names(data))]
c_headers <- names(data)[grep(pattern="^c\\d+", x=names(data))]

# headers to include in selection
i_headers_and_time <- c("time",i_headers)
c_headers_and_time <- c("time",c_headers)

i_data <- data %>% select(all_of(i_headers_and_time))
c_data <- data %>% select(all_of(c_headers_and_time))

i_data_l <- pivot_longer(i_data, 
                         cols=all_of(i_headers), 
                         names_to="investment_val", 
                         values_to = "frequency")

i_data_l <- mutate(i_data_l
                   ,ival=as.numeric(gsub(pattern="^i",x=i_data_l$investment_val, replacement=""))
                   ,phenval=ival/max(ival)
                   )

c_data_l <- pivot_longer(c_data, 
                         cols=all_of(c_headers), 
                         names_to="investment_val", 
                         values_to = "frequency")

c_data_l <- mutate(c_data_l
                   ,cval=as.numeric(gsub(pattern="^c",x=c_data_l$investment_val, replacement=""))
                   ,phenval=cval/max(cval)
                   )

p_freq1 <- ggplot(data=i_data_l
       ,mapping=aes(x=time,y=phenval)) +
    geom_tile(mapping=aes(fill=frequency)) 

p_freq2 <- ggplot(data=c_data_l
       ,mapping=aes(x=time,y=phenval)) +
    geom_tile(mapping=aes(fill=frequency)) 

mean_data_l <- data %>% 
    select(all_of(c("mean_i","mean_c","time"))) %>% 
    pivot_longer(cols=all_of(c("mean_i","mean_c"))
                 ,names_to="mean"
                 ,values_to="mean_values")

mean_dat <- ggplot(data=mean_data_l
                   ,mapping=aes(x=time, y=mean_values)) +
    geom_line(mapping=aes(colour=mean))

p_freq1 / p_freq2 / mean_dat

ggsave(filename=paste0("plot_",file,".pdf"))
