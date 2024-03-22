#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}


# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) 
{
    if (!is.null(x))
    {
       sym(x)
    }
}


find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^mu_x;",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

# finds the start of the statistics part of the file
find.first.line <- function(filename) {
    f <- readLines(filename)

    for (line_i in 1:length(f))
    {
        if (length(grep("^time;",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}


xvar <- "time"

# structure of the graph in JSON
plots_json <- paste0('[
    {"xvar" : "',xvar,'",
    "yvar" : ["x","y"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_x","var_y"]
    },
    {"xvar" : "',xvar,'",
    "yvar" : ["xp","yp"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_xp","var_yp"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["npaired","nsingle","nmort"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["mean_resources_single","mean_resources_paired","mean_resources"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_resources_single","var_resources_paired","var_resources"]
    }
    ]'
)

# if the file name not provided raise error
if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

firstline <- find.first.line(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,skip=firstline - 1
        ,col_names=T)

if (nrow(data.tibble) > 50000)
{
    data.tibble <- data.tibble[data.tibble$time %% 10 == 0,]
}

params <- NULL

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line - 1
        ,col_names=c("name","value")
        ,n_max=firstline - param.line 
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

plot.structure <- fromJSON(plots_json, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l))

plot.list.idx <- 1

# first plot stuff for the tibble as a whole
for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)
    
    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- data.tibble %>% pivot_longer(
                cols=any_of(yvar)
                ,names_to=yvar_name
                ,values_to=yvar_values)
        

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values),
                FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args))
                            
    } else {

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar
                        ),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line()
	}
    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

plot.title <- ""

names.print <- c("startup_cost","mort","dismiss_error","game_type")

for (name in names.print)
{
    if (!(name %in% names(params) ))
    {
        params[name] <- NA
    }
}

if (exists("params")  && "startup_cost" %in% names(params))
{
    plot.title <- paste0(
                "startup cost: ",params["startup_cost"],
                ", mort: ",params["mortality_prob"],
                ", dismiss error: ",params["dismiss_error"],
                ", game: ",params["game_type"]
                )
}


wrap_plots(plot.list,ncol=1) + plot_annotation(
        title=plot.title)

file.name <- paste0("graph_",basename(file.name),".png")

ggsave(file.name,height= 3 * plot.structure.l)

problems()
