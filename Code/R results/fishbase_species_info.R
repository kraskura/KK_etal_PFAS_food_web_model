# Krista Kraskura
# fish base species info
# june 2 2024

# install.packages("rfishbase")
library(rfishbase)

fish <- c("Lepomis gibbosus", "Lepomis macrochirus",
          "Clinostomus funduloides", "Rhinichthys cataractae",
          "Rhinichthys atratulus", "Notropis procne",
          "Etheostoma fusiforme", "Etheostoma olmstedi",
          "Umbra pygmaea", "Noturus insignis",
          "Erimyzon oblongus", "Fundulus diaphanus",
          "Micropterus salmoides", "Semotilus corporalis")

# filter out only ecology and size metrics for each species
fish_table<-fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% fish) %>% 
  select(sci_name, FBname, Length, Weight, LongevityWild, DemersPelag)

# add acronyms
fish_table$SppAlias<-c("Dac", "Chu", "Dar", "Dar", "Kil", 
                       "Pum", "Bgl", "Bas", "Swa", "Mad",
                       "Dac", "Dac", "Fal", "Min")

pry_vals<-data.frame(sci_name = "Prey",
                     FBname = "Prey",
                     Length = NA, 
                     Weight = NA, 
                     LongevityWild = NA, 
                     DemersPelag = "benthopelagic", 
                     SppAlias = "Pry")
fish_table<-rbind(fish_table, pry_vals)
