# Data Visualization of PIRCHE-II
# Used and modified script "hlaR_ethnic_pop_counts.R" by Grace Wager

if(!require("dtplyr"))
  install.packages("dtplyr")
if(!require("dplyr"))
  install.packages("dplyr")
if (!require("tidyverse"))
  install.packages("tidyverse")
if (!require("dplyr"))
  install.packages("dplyr")
if (!require("plyr"))
  install.packages("plyr")
if (!require("tidyr"))
  install.packages("tidyr")
if (!require("stringr"))
  install.packages("stringr")
if (!require("scales"))
  install.packages("scales")
if (!require("data.table"))
  installed.packages("data.table")
if (!require("rlang"))
  installed.packages("rlang")
if (!require("pastecs"))
  installed.packages("pastecs")
if (!require("tidyselect"))
  install.packages("tidyselect")
if (!require("Hmisc"))
  install.packages("Hmisc")
if (!require("mlbench"))
  install.packages("mlbench")
if (!require("caret"))
  install.packages("caret")
if (!require("leaps"))
  install.packages("leaps")
if (!require("tools"))
  install.packages("tools")
if (!require("purrr"))
  install.packages("purrr")
if (!require("fs"))
  install.packages("fs")
if (!require("matrixStats"))
  install.packages("matrixStats")
if (!require("modeest"))
  install.packages("modeest")
if (!require("reshape2"))
  install.packages("reshape2")
if (!require("questionr"))
  install.packages("questionr")
if (!require("nlme"))
  install.packages("nlme")
if (!require("survival"))
  install.packages("survival")
if (!require("ranger"))
  install.packages("ranger")
library(ggfortify)

# Set path
path <- paste0(getwd(),"/")


# Set risk categories to analyze
ethnicity <- c("AFA","CAU","HIS","UKN","HPI","MLT","NAM","ASN")

# Get PIRCHE file
molecule <- read.table("./clean_pirche100_GLOBAL_pops.csv", sep = ',',header = T)
head(molecule)
risk_category_list <- molecule %>% dplyr::select("PIRCHE_II", "PIRCHE_II_promiscuity_corrected")


table(molecule['CAN_RACE'])

# Bar plot of distribution of eplet mismatch counts for DR and for DQ.
Drc <- molecule %>% dplyr::select("PX_ID", "CAN_RACE", "DRB1_Originated_Epitopes_count")
# Drc$DR = rowSums(cbind(Drc$DRB1_1_Originated_Epitopes_count, Drc$DRB1_2_Originated_Epitopes_count), na.rm = TRUE) 
# Drc = Drc[,!(names(Drc) %in% c("DRB1_1_Originated_Epitopes_count","DRB1_2_Originated_Epitopes_count"))]
Drc$loci <- "DR"
names(Drc)[3] <- "mm_cnt"


DQc <- molecule %>% dplyr::select("PX_ID", "CAN_RACE", "DQB1_Originated_Epitopes_count")
# DQc$DQ = rowSums(cbind(DQc$DQB1_1_Originated_Epitopes_count, DQc$DQB1_2_Originated_Epitopes_count), na.rm = TRUE)
# DQc = subset(DQc, select = -c(DQB1_1_Originated_Epitopes_count, DQB1_2_Originated_Epitopes_count))
DQc$loci <- "DQ"
names(DQc)[3] <- "mm_cnt"


DrDq <- rbind(Drc, DQc)


#barbell plot
DrDq %>% 
  ggplot(aes(x=CAN_RACE,
             y=mm_cnt,
             color=loci))+
  geom_boxplot(lwd=1)+
  #geom_jitter(width=0.1,alpha=0.2)+
  theme(legend.position = "right")+
  theme_bw()+
  xlab("Transplant Recipient Ethnicity")+
  ylab("DR|DQ Originated Epitopes Count")+
  labs(title="DR|DQ Originated Epitopes Counts 
       by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="DRDQ_OE_Counts_PIRCHE_BW.pdf", width = 10, height = 10, units="in", path = "./")

