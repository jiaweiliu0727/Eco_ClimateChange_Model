require(tidyverse)

#Read the annual defoliation summary in USA
read_csv("area_acres_stats_usa.csv") -> acreage
acres_per_ha <- 2.471053814672
acreage %>% 
  select(sum, SURVEY_YEAR) %>% 
  rename(acres_defol = sum,
         year = SURVEY_YEAR) %>% 
  bind_rows(data.frame(acres_defol = c(1.7e6, 1.3e6),
                       year = c(2022, 2023))) %>% 
  mutate(ha_defol = acres_defol / acres_per_ha,
         country = "USA") -> acreage

#Read the annual defoliation summary in Ontario, Canada
read_csv("area_sqm_stats_ontario.csv") -> sqm_ontario
area_defoliated <- sqm_ontario %>% 
  rename(sqm_defol = sum,
         year = EVENT_YEAR) %>% 
  mutate(ha_defol = sqm_defol / 1e4) %>% 
  select(-sqm_defol) %>% 
  mutate(acres_defol = ha_defol * acres_per_ha,
         country = "CAN") %>% 
  select(names(acreage)) %>% 
  bind_rows(acreage)

#Plotting the annual defoliation curve (USA+Ontario)
curve<-area_defoliated %>% 
  select(year, country, ha_defol) %>% 
  spread(key = "country", value = ha_defol) %>% 
  gather(key = "country", value = "ha_defol", -year) %>% 
  mutate(ha_defol = ifelse(is.na(ha_defol), 0, ha_defol)) %>% 
  group_by(year) %>% 
  summarise(sum_defol = sum(ha_defol)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = sum_defol)) +
  geom_line(col = "#65024b", lwd = 1.5) +
  geom_point(col = "#65024b", size = 4) +
  scale_y_continuous(expand = c(0,0), limits = c(0,3e6),labels = scales::comma) +
  scale_x_continuous(expand = c(0,0), limits = c(1995,2024))+
  theme_bw() +
  #ggtitle("Area defoliated (US + Ontario)")+
  labs(x= "Year", y = "Area defoliated (ha, US + Ontario)")+
  theme(plot.title=element_text(size=16,face="bold",hjust=0.5), axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.title.y=element_text(vjust=0.5),
        legend.position = "bottom",
        legend.title=element_text(size=16,face="bold"),
        legend.text=element_text(size=16),legend.key.size=unit(1, "cm"),
        legend.background = element_blank(),legend.key = element_blank(),
        plot.background = element_blank(),panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))

ggsave("Figure6_defo_cuvre.pdf",curve,width=8, height=6)
