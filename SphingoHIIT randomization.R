#===============================================================================
# The SphingoHIIT study: randomization
# Author: Justin Carrard, & Denis Infanger
#===============================================================================

#import_library
library(ggplot2)
library(randomizeR)
#install.packages("writexl")
library("writexl")

#females and males
N <- 36 # Total participants
K <- 2  # Number of groups
# groups <- c("f_HIIT", "f_noHIIT", "m_HIIT", "m_noHIIT") # Name of the groups
groups <- c("HIIT", "noHIIT") # Name of the groups

# Initialize the randomization procedure
params_men <- rpbrPar(N/2, rb = c(2, 4, 6), K = K, groups = groups)
params_women <- rpbrPar(N/2, rb = c(2, 4, 6), K = K, groups = groups)

# Generate one sequence of randomization (change seed for different sequence)
mySeq_men <- genSeq(params_men, seed = 654983511)
mySeq_women <- genSeq(params_women, seed = 9198751)

# Display sequence
list_men <- getRandList(mySeq_men)
list_women <- getRandList(mySeq_women)

table(list_men)
table(list_women)

final_list <- data.frame(
  # id = 1:N
  group = c(list_men, list_women)
  , sex = rep(c("m", "f"), each = N/2)
)

write_xlsx(final_list,"")
