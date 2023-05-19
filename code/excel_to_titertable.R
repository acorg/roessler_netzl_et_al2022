# format titer table from excel sheet
rm(list = ls())
library(tidyverse)
library(stringr)

# vaccine status for sheet I
vacc_stat <- c("n" = "",
               "1" = "ChAdOx1-S/ChAdOx1-S", 
               "2" = "ChAdOx1-S/BNT162b2",
               "3" = "BNT162b2/BNT162b2",
               "4" = "mRNA-1273/mRNA-1273",
               "d" = "divers")

# function for sr_info

sr_info_function <- function(x, kimpel1) {
  stat <- kimpel1[[x, "Immune status"]]
  vacc <- kimpel1[[x,4]]
  var <- kimpel1[[x, "Variant"]]
    
  if(stat == "convalescent") {
    res <- paste0(var, " conv.")
    
  } else if(str_detect(stat, "con/vacc")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(var, "+", substring(vacc, 4, nchar(vacc)-1))
    } else {
      res <- paste0(var, "+", vacc_stat[vacc])
    }
    
  } else if(str_detect(stat, "vacc/con")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(substring(vacc, 4, nchar(vacc)-1), "+", var)
    } else {
      res <- paste0(vacc_stat[vacc],"+", var)
    }
    
  } else {
    res <- vacc_stat[vacc]
  }
  
  return(res)
}

shorten_vacc_names <- function(sr_group) {
  sr_group <- gsub("ChAdOx1-S", "AZ", sr_group)
  sr_group <- gsub("ChAdOx1", "AZ", sr_group)
  sr_group <- gsub("BNT162b2", "BNT", sr_group)
  sr_group <- gsub("BNT162b", "BNT", sr_group)
  sr_group <- gsub("Ad26.COV2.S", "J\\&J", sr_group)
  sr_group <- gsub("mRNA\\-1273", "mRNA1273", sr_group)
  
  sr_group <- gsub("B.1.617.2", "delta", sr_group)
  sr_group <- gsub("B.1.351", "beta", sr_group)
  sr_group <- gsub("B.1.1.7\\/B.1.1.7\\+E484K", "alpha\\/alpha\\+E484K", sr_group)
  
  return(sr_group)
}

set_threshold <- function(tab, thresh) {
  tab[as.numeric(tab) < as.numeric(thresh)] <- paste0("<", thresh)
  tab[is.na(tab)] <- "*"
  
  return(tab)
}


# =============================== Omicron I sheet, published in NEJM
kimpel1 <- readxl::read_excel("./data/titer_data/220309_Omicron I_raw data.xlsx", sheet = 1)
colnames(kimpel1) <- kimpel1[2,]
kimpel1 <- kimpel1[c(3:68),]
kimpel1$Variant <- gsub(" ", "", kimpel1$Variant)

# add sample info
# save sample ID, then sr group, then sr info including time since dose
kimpel1 <- kimpel1 %>%
  mutate(
    sr_info = unlist(lapply(1:nrow(kimpel1), function(x) sr_info_function(x, kimpel1)
  )),
  sr_group = shorten_vacc_names(sr_info)
  )

# set serum_info_column that contains all info for map
kimpel1 <- kimpel1 %>%
  mutate(sr_info_complete = paste(`sample ID`,sr_group, sr_info, `Blood collection\r\nDays after 2nd dose`, sep = "_"))

# pivot wider to titer table
kimpel1 %>%
  column_to_rownames("sr_info_complete") %>%
  select(c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", "BA.5")) -> kimpel1_table

#======================================  Omicron sheet II 

kimpel <- readxl::read_excel("./data/titer_data/220309_Omicron II_raw data.xlsx", sheet = 1)
colnames(kimpel) <- kimpel[2,]
kimpel <- kimpel[c(3:59),2:ncol(kimpel)]

kimpel$Cohort[1] <- "Vaccinated"
# fill the cohort rows
no_na_cohort <- c(which(!is.na(kimpel$Cohort)), nrow(kimpel)+1)

for(p in c(1:(length(no_na_cohort)-1))) {
  
  kimpel$Cohort[c(no_na_cohort[p]:(no_na_cohort[(p+1)]-1))] <- kimpel$Cohort[no_na_cohort[p]]
}

# add sample info
# save sample ID, then sr group, then sr info including time since dose
kimpel <- kimpel %>%
  mutate(n_vacc = unlist(lapply(kimpel$Vaccination, function(x) {
    if(x == "-") {
      0
    } else {
      length(str_split(x, "/")[[1]])
    }
  })),
  sr_info = unlist(lapply(1:nrow(kimpel), function(x){
    if(kimpel[[x, "Cohort"]] == "Vaccinated") {
      paste0(kimpel[[x, "Vaccination"]], "+BA.1")
    } else if(kimpel[[x, "Cohort"]] == "Unvaccinated") {
      "BA.1 conv."
    } else if(kimpel[[x, "Cohort"]] == "Vaccinated with prior infection") {
      if(kimpel[[x,"Variant of prior infection"]] == "B.1.617.2") {
        paste0(kimpel[[x, "Vaccination"]], "+", kimpel[[x,"Variant of prior infection"]], "+BA.1")
      } else {
        paste0(kimpel[[x,"Variant of prior infection"]], "+", kimpel[[x, "Vaccination"]], "+BA.1")
      }
    } else {
      paste0(kimpel[[x,"Variant of prior infection"]], "+BA.1")
    }
  })), 
  sr_group = unlist(lapply(1:nrow(kimpel), function(x){
    if(kimpel[[x, "Cohort"]] == "Vaccinated") {
      paste0("Vacc+BA.1")
    } else if(kimpel[[x, "Cohort"]] == "Unvaccinated") {
      "BA.1 conv."
    } else if(kimpel[[x, "Cohort"]] == "Vaccinated with prior infection") {
      "Vacc+BA.1 reinf."
    } else {
      "BA.1 reinf."
    }
  })), 
  # account for infection sequence
  sr_group = shorten_vacc_names(sr_group),
  sr_info_complete = paste(`sample ID`,sr_group, sr_info, `Blood collection\r\nDays after positive qPCR for BA.1`, `Study Participant`, sep = "_")
  )

# select titer columns
# pivot wider to titer table
kimpel %>%
  column_to_rownames("sr_info_complete") %>%
  select(c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", "BA.5")) -> kimpel2_table

#======================================  3rd Omicron table

kimpel <- readxl::read_excel("./data/titer_data/220308_further data for Omicron III_raw data.xlsx", sheet = 1)
colnames(kimpel) <- kimpel[2,]
kimpel <- kimpel[c(3:71),]

# fill the cohort rows
no_na_cohort <- c(which(!is.na(kimpel$Cohort)), nrow(kimpel)+1)

for(p in c(1:(length(no_na_cohort)-1))) {
  
  kimpel$Cohort[c(no_na_cohort[p]:(no_na_cohort[(p+1)]-1))] <- kimpel$Cohort[no_na_cohort[p]]
}

kimpel <- kimpel %>%
  mutate(n_vacc = unlist(lapply(kimpel$Vaccination, function(x) {
    if(is.na(x)) {
      3
    } else if(x == "-") {
      0
    } else if(x == "?"){
      NA
    } else {
      length(str_split(x, "/")[[1]])
    }
  })),
  sr_info = unlist(lapply(1:nrow(kimpel), function(x){
    if(kimpel[[x, "Cohort"]] == "BA.2 convalescent sera") {
      temp <- paste0(kimpel[[x, "Variant of prior infection"]],"+BA.2")
      temp <- gsub("\\-\\+", "", temp)
      temp <- gsub("\\?\\_", "", temp)
      
      if(kimpel[[x, "Immune status"]] == "vaccinated") {
        temp <- paste0(kimpel[[x, "Vaccination"]], "+BA.2")
        temp <- gsub("\\?", "Vacc", temp)
      } else if(kimpel[[x, "Immune status"]] == "3x vaccinated") {
        temp <- "Vacc+BA.2"
      }
    } else if(kimpel[[x, "Cohort"]] == "First wave (WT) convalescent"){
      temp <- "WT conv."
      
    } else if(kimpel[[x, "Cohort"]] == "BA.1 convalescent") {
      temp <- "Vacc/Vacc/Vacc+BA.1"
      
    } else if(kimpel[[x, "Cohort"]] == "naive + 3x vaccinated"){
      temp <- kimpel[[x, "Vaccination"]]
    } else if(kimpel[[x, "Cohort"]] == "BA.5 convalescent"){
      if(kimpel[[x, "Immune status"]] == "vacc + BA.5") {
        temp <- "Vacc+BA.5"
      } else if(kimpel[[x, "Immune status"]] == "unvaccinated"){
        temp <- "BA.5 conv."
      } else {
        temp <- "BA.2+BA.5 reinf."
      }
    } else { # this here for Vaccinated idvls with non-omicron breakthrough
      temp <- paste0(kimpel[[x, "Vaccination"]], "+", kimpel[[x, "Variant of prior infection"]])
      temp <- gsub(" \\?", "", temp)
      
    } 
    temp <- gsub(" \\_","\\_", temp )
  })), 
  sr_group = unlist(lapply(1:nrow(kimpel), function(x){
    if(kimpel[[x, "Cohort"]] == "BA.2 convalescent sera") {
      if(kimpel[[x, "Immune status"]] == "vaccinated") {
        "Vacc+BA.2"
      } else if(kimpel[[x, "Immune status"]] == "unvaccinated + reinfected") {
        "BA.2 reinf."
      } else if(kimpel[[x, "Immune status"]] == "3x vaccinated") {
        "Vacc+BA.2" # potentially combine with vacc + BA.2. look at landscapes to see how similar they are
      } else {
        "BA.2 conv."
      }
    }  else if(kimpel[[x, "Cohort"]] == "First wave (WT) convalescent"){
      "WT conv."
    } else if(kimpel[[x, "Cohort"]] == "BA.1 convalescent") {
      temp <- "Vacc+BA.1"
    } else if(kimpel[[x, "Cohort"]] == "naive + 3x vaccinated"){
      temp <- kimpel[[x, "Vaccination"]]
    } else if(kimpel[[x, "Cohort"]] == "BA.5 convalescent"){
      if(kimpel[[x, "Immune status"]] == "vacc + BA.5") {
        temp <- "Vacc+BA.5"
      } else if(kimpel[[x, "Immune status"]] == "unvaccinated"){
        temp <- "BA.5 conv."
      } else {
        temp <- "BA.2+BA.5 reinf."
      }
    } else {
      temp <- paste0(kimpel[[x, "Vaccination"]], "+", kimpel[[x, "Variant of prior infection"]])
      gsub(" \\?", "", temp)
    }
  })), 
  sr_group = shorten_vacc_names(sr_group),
  sr_info_complete = paste(`sample ID`,sr_group, sr_info, `Blood collection\r\nDays after positive qPCR for BA.1`, `Study Participant`, sep = "_")
  )

# select titer columns
# pivot wider to titer table
kimpel %>%
  column_to_rownames("sr_info_complete") %>%
  select(c("D614G", "B.1.1.7","B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", "BA.5")) -> kimpel3_table

# combine the three tables
kimpel_table <- rbind(kimpel1_table, kimpel2_table, kimpel3_table)
kimpel_table_wide <- t(kimpel_table)
# set not titrated to "*"
kimpel_table_wide[kimpel_table_wide == "-"] <- "*"

# add lower threshold at the end
kimpel_table_wide <- set_threshold(kimpel_table_wide, 16)

write.csv(kimpel_table_wide, "./data/titer_data/titer_table.csv")

