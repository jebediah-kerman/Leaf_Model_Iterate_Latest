library("dplyr", quietly = T)

# Parameters
weights <- c(1,2)
n_out <- 6

# Processing
dat <- read.table("./Output_Summary/good_output.txt", sep = "\t", header = T)
dat.mean <- dat %>%
  group_by(i, j, k) %>%
  summarise(Crit_AFInit_Mean = mean(Crit_AFInit),
            Crit_HWRatio_Mean= mean(Crit_HWRatio))
ranks <- rank(dat.mean$Crit_AFInit_Mean) * weights[1] + rank(dat.mean$Crit_HWRatio_Mean) * weights[2]
out <- dat.mean[order(ranks)[1:n_out],]

# Output
print(as.data.frame(out))

