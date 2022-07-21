# Script to address more than one aroclor reported per sample
# The total PCB is obtained from the average of the Aroclor mixtures
# reported

# Read data in pg/L
# WaterDataCongenerAroclor062322.csv
# When getting individual PCB congeners, the total PCB
# was calculated using the sum of the aroclor reported.
wdc.0 <- read.csv("WaterDataCongenerAroclor062322.csv")

# Select samples with Aroclor method
aroc <- wdc.0[wdc.0$AroclorCongener == "Aroclor", ]

# Select samples with more than one Aroclor reported
MultiAroclor <- subset(aroc, rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1016 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1221 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1232 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1242 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1248 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1254 |
                         rowSums(aroc[, c(115:122)], na.rm = T) > aroc$A1260)

# Remove metadata
MultiAroclor.2 <- subset(MultiAroclor, select = -c(ID:AroclorCongener))
# Remove Aroclor data
MultiAroclor.3 <- subset(MultiAroclor.2, select = -c(A1016:A1260))
# Sum Aroclors
sum.aroc <- rowSums(MultiAroclor.2[, c(105:111)], na.rm = TRUE)
# Average Aroclors
ave.aroc <- rowMeans(MultiAroclor.2[, c(105:111)], na.rm = TRUE)
# Proportion of individual PCBs
prop.aroc <- MultiAroclor.3/sum.aroc
new.cong <- prop.aroc*ave.aroc

