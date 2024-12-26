library(tidyverse)
mrsa_2013 <- read.csv("data/MRSA_bloodstream_CA_2013.csv")
mrsa_2014 <- read.csv("data/MRSA_bloodstream_CA_2014.csv")
mrsa_2015 <- read.csv("data/MRSA_bloodstream_CA_2015.csv")
mrsa_2016 <- read.csv("data/MRSA_bloodstream_CA_2016.csv", fileEncoding = "latin1")
mrsa_2017 <- read.csv("data/MRSA_bloodstream_CA_2017.csv", fileEncoding = "latin1")
mrsa_2018 <- read.csv("data/MRSA_bloodstream_CA_2018.csv", fileEncoding = "latin1")
mrsa_2019 <- read.csv("data/MRSA_bloodstream_CA_2019.csv", fileEncoding = "latin1")
mrsa_2020_1 <- read.csv("data/MRSA_bloodstream_CA_2020.csv", fileEncoding = "latin1")
mrsa_2020_2 <- read.csv("data/MRSA_bloodstream_CA_2020-1.csv", fileEncoding = "latin1")
mrsa_2021 <- read.csv("data/MRSA_bloodstream_CA_2021.csv")
mrsa_2022 <- read.csv("data/MRSA_bloodstream_CA_2022.csv", fileEncoding = "latin1")
mrsa_2023 <- read.csv("data/MRSA_bloodstream_CA_2023.csv", fileEncoding = "latin1")

sum_2013 <- sum(mrsa_2013$Infection_Count, na.rm = TRUE)
sum_2013 <- as.numeric(sum_2013)

sum_2014 <- sum(mrsa_2014$Cases, na.rm = TRUE)
sum_2014 <- as.numeric(sum_2014)

sum_2015 <- sum(mrsa_2015$Hospital_Onset_Cases, na.rm= TRUE)
sum_2015 <- as.numeric(sum_2015)

sum_2016 <- mrsa_2016$Infections_Reported[1] + mrsa_2016$Infections_Reported[2] +  mrsa_2016$Infections_Reported[3] +  mrsa_2016$Infections_Reported[4]
sum_2016 <- as.numeric(sum_2016)

sum_2017 <- mrsa_2017$Infections_Reported[1] + mrsa_2017$Infections_Reported[2] +  mrsa_2017$Infections_Reported[3] +  mrsa_2017$Infections_Reported[4]
sum_2017 <- as.numeric(sum_2017)

sum_2018 <- mrsa_2018$Infections_Reported[1] + mrsa_2018$Infections_Reported[2] +  mrsa_2018$Infections_Reported[3] +  mrsa_2018$Infections_Reported[4]
sum_2018 <- as.numeric(sum_2018)

sum_2019 <- mrsa_2019$Infections_Reported[1] + mrsa_2019$Infections_Reported[2] +  mrsa_2019$Infections_Reported[3] +  mrsa_2019$Infections_Reported[4]
sum_2019 <- as.numeric(sum_2019)

sum_2020 <- mrsa_2020_1$Infections_Reported[1] +mrsa_2020_1$Infections_Reported[2] +  mrsa_2020_1$Infections_Reported[3] +  mrsa_2020_1$Infections_Reported[4] + mrsa_2020_2$Infections_Reported[1] + mrsa_2020_2$Infections_Reported[2] +  mrsa_2020_2$Infections_Reported[3] +  mrsa_2020_2$Infections_Reported[4]
sum_2020 <- as.numeric(sum_2020)


sum_2021 <- mrsa_2021$Infections_Reported[1] + mrsa_2021$Infections_Reported[2] +  mrsa_2021$Infections_Reported[3] +  mrsa_2021$Infections_Reported[4]
sum_2021 <- as.numeric(sum_2021)

sum_2022 <- mrsa_2022$Infections_Reported[1] + mrsa_2022$Infections_Reported[2] +  mrsa_2022$Infections_Reported[3] +  mrsa_2022$Infections_Reported[4]
sum_2022 <- as.numeric(sum_2022)

sum_2023 <- mrsa_2023$Infections_Reported[1] + mrsa_2023$Infections_Reported[2] +  mrsa_2023$Infections_Reported[3] +  mrsa_2023$Infections_Reported[4]
sum_2023 <- as.numeric(sum_2023)

mrsa_cases <- data.frame(
  year = c(2013:2023),
  cases = c(sum_2013,sum_2014, sum_2015, sum_2016, sum_2017, sum_2018, sum_2019, sum_2020, sum_2021, sum_2022, sum_2023)
)
library(ggplot2)

ggplot(data = mrsa_cases, aes(x = year, y = cases))+
  geom_point(size = 2, color = "black")+
  geom_smooth()+
  labs(
    title = "California MRSA bloodstream infections",
    subtitle = "Data from California Health and Human Services Open Data Portal (2013 to 2023)"
  )+
  xlab("Number of Cases")+
  ylab("Year")+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_continuous(
    limits = c(2013, 2023),
    breaks = seq(2013, 2023, by = 1)
  )+
  theme(

    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(face = "bold", color = "red", size = 12, hjust = 0.5)

  )
ggsave("figures/mrsa_timeline.png")
