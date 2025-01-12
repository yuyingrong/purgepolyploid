# Analyzing under- & over-purging from BUSCO full_table.tsv

Yuying Rong
20240923

Due to signs of both under- and over-purging, I check the exact number of each Busco gene in the unpurged, purged, and single haplotype assemblies. 

## Create Busco count table per assembly

```r
## read in the BUSCO output full_table.tsv

# only count frequency for Complete and Duplicated Buscos

# unpurged
up <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_unpurged/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(up)[1] <- "Busco.id"
nrow(up)# 14764
unique(up$Status)# "Duplicated" "Complete" "Missing" "Fragmented"
nrow(up[up$Status == "Complete" | up$Status == "Duplicated",])# 14414
# count freq for all complete and duplicated buscos
up_freq <- as.data.frame(table(up[up$Status == "Complete" | up$Status == "Duplicated", 1]))
colnames(up_freq)[1] <- "Busco"

# purged
pg <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_to100000/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(pg)[1] <- "Busco.id"
nrow(pg)# 6320
nrow(pg[pg$Status == "Complete" | pg$Status == "Duplicated",])# 4142
pg_freq <- as.data.frame(table(pg[pg$Status == "Complete" | pg$Status == "Duplicated", 1]))
colnames(pg_freq)[1] <- "Busco"

# h0
h0 <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_h0/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(h0)[1] <- "Busco.id"
nrow(h0)# 7196
nrow(h0[h0$Status == "Complete" | h0$Status == "Duplicated",])# 6197
h0_freq <- as.data.frame(table(h0[h0$Status == "Complete" | h0$Status == "Duplicated", 1]))
colnames(h0_freq)[1] <- "Busco"

# h1
h1 <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_h1/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(h1)[1] <- "Busco.id"
nrow(h1)# 6417
nrow(h1[h1$Status == "Complete" | h1$Status == "Duplicated",])# 5054
h1_freq <- as.data.frame(table(h1[h1$Status == "Complete" | h1$Status == "Duplicated", 1]))
colnames(h1_freq)[1] <- "Busco"

# h2
h2 <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_h2/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(h2)[1] <- "Busco.id"
nrow(h2)# 6401
nrow(h2[h2$Status == "Complete" | h2$Status == "Duplicated",])# 4944
h2_freq <- as.data.frame(table(h2[h2$Status == "Complete" | h2$Status == "Duplicated", 1]))
colnames(h2_freq)[1] <- "Busco"

# h3
h3 <- read.delim("/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/potato/busco_h3/run_solanales_odb10/full_table.tsv", header=TRUE, sep="\t", skip=2)
colnames(h3)[1] <- "Busco.id"
nrow(h3)# 6207
nrow(h3[h3$Status == "Complete" | h3$Status == "Duplicated",])# 4809
h3_freq <- as.data.frame(table(h3[h3$Status == "Complete" | h3$Status == "Duplicated", 1]))
colnames(h3_freq)[1] <- "Busco"


## combine columns of frequency counts into one table

# the challenge now is that the number of Busco id in each occurrence col is not identical

# total number of Busco_ids:
length(unique(h1$Busco.id))# 5950
# number of Complete and Duplicate Busco_ids in h1
nrow(h1_freq)# 4587, the Fragmented and Missing ones are not there, so the number is below 5950

# to over the the challenge, use the total number of Busco_ids as rownames to create table
t <- data.frame(Busco = unique(up$Busco.id))

# merge all colmns into one table
t2 <- merge(t, up_freq, by = "Busco", all.x = TRUE, sort = TRUE)
# the first item is considered x, and the second item is considered y
# if all.x, then all keys will match the keys in x, the missing ones get NA

# this does not work: t2 <- merge(t, c(up_freq, pg_freq, h0_freq, h1_freq, h2_freq, h3_freq), by = "Busco", all.x = TRUE)

# update colnames by name of the assembly
colnames(t2)[2] <- 'unpurged'
t2 <- merge(t2, pg_freq, by = "Busco", all.x = TRUE, sort = TRUE)
colnames(t2)[3] <- 'purged'
t2 <- merge(t2, h0_freq, by = "Busco", all.x = TRUE, sort = TRUE)
colnames(t2)[4] <- 'h0'
t2 <- merge(t2, h1_freq, by = "Busco", all.x = TRUE, sort = TRUE)
colnames(t2)[5] <- 'h1'
t2 <- merge(t2, h2_freq, by = "Busco", all.x = TRUE, sort = TRUE)
colnames(t2)[6] <- 'h2'
t2 <- merge(t2, h3_freq, by = "Busco", all.x = TRUE, sort = TRUE)
colnames(t2)[7] <- 'h3'

View(t2)
```

Busco | unpurged | purged | h0 | h1 | h2 | h3
--- | --- | --- | --- | --- | --- | --- 
17319at4069 | 18 | 1 | NA | 1 | 1 | 1
1739at4069 | 14 | 14 | 14 | 1 | 7 | 9
33919at4069 | 10 | 12 | 2 | 1 | NA | 10
28343at4069 | 10 | 4 | 3 | 4 | 6 | 3
17228at4069 | 6 | 5 | 4 | 7 | 2 | 4
... | ... | ... | ... | ... | ... | ...
10050at4069 | 3 | 1 | 1 | 1 | 1 | 1
10062at4069 | 3 | 1 | 1 | 1 | 1 | 0
10086at4069 | 1 | 1 | 0 | 0 | 1 | 0
... | ... | ... | ... | ... | ... | ...
9415at4069 | 3 | NA | 1 | 1 | 1 | 1
947at4069 | 3 | NA | 1 | 1 | NA | 1
991at4069 | 2 | NA | NA | 1 | NA | 1
965at4069 | 1 | NA | 3 | 1 | 1 | 1
9483at4069 | NA | NA | NA | NA | 1 | 1

From the Busco gene count table: 
- First, it seems that purging did a good job removing some duplications (such as in the middle rows) but sometimes underpurged the regions containing multi-copy Busco's (upper rows) or overpurged (bottom rows).
- Second, I can see that BUSCO sometimes detected Busco's in the single haplotype genomes, but not in the unpurged genome, or the Busco counts don’t add up (17319at4069, 10050at4069, 9450at4069). This could be a statistical phenomenon, where detection power fluctuated with factors such as different assembly sizes.

## Examining Busco containing contigs

TBA

## Examining Busco containing alignments

TBA
