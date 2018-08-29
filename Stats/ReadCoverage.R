# Yuka Takemon
# 08/29/18
# Calculate average read coverage
library(tidyverse)
library(stringr)
library(data.table)
wd <- "/projects/korstanje-lab/ytakemon/Col4a5xDO/"
setwd(wd)

# Identify samples
input_sample <- commandArgs(trailingOnly = TRUE)[1] # input_sample <- "1415-0557_cat_R1.fastq.gz"
name <- str_sub(input_sample,,-17)
# find samtools depth output
sample_dir <- list.files( path = "./civet_run", pattern = name, full.name = TRUE)
input_file <- list.files( path = sample_dir, pattern = "samtools.depth", full.name = TRUE)

# Load data
df <- fread(input = input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, verbose = TRUE)
print(paste0("Done reading in data from: ", name))

# rename and restructure
names(df) <- c("ID","pos","depth")
# get gene and founder information to split
df1 <- df[1:(nrow(df)/2),]
df1 <- df1 %>% mutate(
  Founder = str_sub(ID,-1,),
  Name = str_sub(ID,,-3))
df2 <- df[((nrow(df)/2)+1):nrow(df),]
df2 <- df2 %>% mutate(
  Founder = str_sub(ID,-1,),
  Name = str_sub(ID,,-3))
# join together
df <- rbind(df1, df2)

# split by founder to rename depth column to join later
df_A <-  filter(df, Founder == "A") %>%
  select(Name, pos, depth) %>%
  rename(depth_A = depth)
df_B <-  filter(df, Founder == "B") %>%
  select(Name, pos, depth) %>%
  rename(depth_B = depth)
df_C <-  filter(df, Founder == "C") %>%
  select(Name, pos, depth) %>%
  rename(depth_C = depth)
df_D <-  filter(df, Founder == "D") %>%
  select(Name, pos, depth) %>%
  rename(depth_D = depth)
df_E <-  filter(df, Founder == "E") %>%
  select(Name, pos, depth) %>%
  rename(depth_E = depth)
df_F <-  filter(df, Founder == "F") %>%
  select(Name, pos, depth) %>%
  rename(depth_F = depth)
df_G <-  filter(df, Founder == "G") %>%
  select(Name, pos, depth) %>%
  rename(depth_G = depth)
df_H <-  filter(df, Founder == "H") %>%
  select(Name, pos, depth) %>%
  rename(depth_H = depth)

# Join all based on name and position
df_all_founders <- full_join(df_A, df_B, by = c("Name","pos")) %>%
  full_join(., df_C, by = c("Name","pos")) %>%
  full_join(., df_D, by = c("Name","pos")) %>%
  full_join(., df_E, by = c("Name","pos")) %>%
  full_join(., df_F, by = c("Name","pos")) %>%
  full_join(., df_G, by = c("Name","pos")) %>%
  full_join(., df_H, by = c("Name","pos"))

# Assign 0 reads to NA reads
df_all_founders[is.na(df_all_founders$depth_A),]$depth_A <- 0
df_all_founders[is.na(df_all_founders$depth_B),]$depth_B <- 0
df_all_founders[is.na(df_all_founders$depth_C),]$depth_C <- 0
df_all_founders[is.na(df_all_founders$depth_D),]$depth_D <- 0
df_all_founders[is.na(df_all_founders$depth_E),]$depth_E <- 0
df_all_founders[is.na(df_all_founders$depth_F),]$depth_F <- 0
df_all_founders[is.na(df_all_founders$depth_G),]$depth_G <- 0
df_all_founders[is.na(df_all_founders$depth_H),]$depth_H <- 0

# add reads per position
df_all_founders <- df_all_founders %>% mutate(
  depth_all = depth_A + depth_B + depth_C + depth_D + depth_E + depth_F + depth_G + depth_H
)

# calcualte coverage: (sum of total founder depth / # of positions) / # of founders
cover <- (sum(df_all_founders$depth_all)/nrow(df_all_founders))/8
output <- matrix(c(name, cover), nrow = 1, ncol = 2)

# append output to list
write.table(output, file = "./civet_run/AvgReadCoverage.txt", sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
