#!/usr/bin/env Rscript

library(tidyverse)
library(tools)
args <- commandArgs(trailingOnly=TRUE)

filename <- args[[1]]

find.params <- function(filename)
{
    f = readLines(filename)

    print(length(f))
}

param.line <- find.params(filename)


#data <- read.table("sim_mutual_direct_17_3_2016_103342_585174049",sep=";",header=T, nrow=5000)
data <- read.table(filename,sep=";",header=T, nrow=5000)

#qline <- expression(d * r * (1 + 2 * bm * var_q ) * su / (bm * (4 + 8 * bm * var_q - d^2 * var_q * su)))

#params <- list(d=0.25, r=1.5, bm=0.001,var_q=mean(data[(nrow(data)-10):nrow(data), "var_q"]), su=1)

data.t <- pivot_longer(data,cols=c(mean_p,mean_t,mean_q),names_to="trait",values_to="trait_value")

ggplot(mapping=aes(x=generation,y=trait_value),
        data=data.t) +
    geom_line(mapping=aes(colour=trait))

output_file <- paste0("graph_",file_path_sans_ext(filename),".pdf")

ggsave(file=output_file)
        
