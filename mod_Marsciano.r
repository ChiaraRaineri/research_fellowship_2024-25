# Transcription of mod Marsciano, 2024

ram_data <- read.csv(file="ram_data.csv", header=TRUE, sep=";",
                     colClasses = c("character", "numeric", "numeric", "numeric", "numeric","numeric"))

str(ram_data)


library(dplyr)
library(lifecycle)
library(zoo)
library(ggplot2)
library(scales)
#library(magrittr)
#library(tidyverse)


# Calculate the score for every criterion
ram_data <- ram_data %>%
  dplyr::mutate(
    # Score 1: rain e LW
    score1 = ifelse(rain == 0 & LW < 16, 0, 1),
    
    # Score 2: rain
    score2 = ifelse(rain <= 1.5, 1, 2),
    
    # Score 3: T < 5
    score3 = ifelse(T < 5, 0, 1),
    
    # Score 4: T is a range
    score4 = ifelse(T <= 15, 2, 3),
    
    # Final result
    final_score = ifelse(score1 == 0 | score3 == 0, 0,
                         ifelse(score2 == 2, score4, 1))
  )




# Add a column with the moving average of the preceding 7 values of final_score (7dd column)
ram_data <- ram_data %>%
  mutate(
    rolling_avg = zoo::rollmean(final_score, k = 7, fill = NA, align = "right")
  ) #%>%
   #filter(!is.na(rolling_avg))  # Remove the first 7 rows with NA values 



# Add the equations for the calculation of the infection rate
ram_data <- ram_data %>%
  mutate(
    Teq_ger = ifelse((T - 0) / (36 - 0) > 0, (T - 0) / (36 - 0), 0),
    
    GERrate = ((4.15 * Teq_ger^1.051 * (1 - Teq_ger))^3.283) * (1 - 0.848^LW),
    
    Teq_myc = ifelse((T - 2) / (40 - 2) > 0, (T - 2) / (40 - 2), 0),
    
    MYCrate = ((3.808 * Teq_myc^0.961 * (1 - Teq_myc))^5.97) / (1 + 19.925 * exp(-0.064 * LW)),
    
    INFrate = (GERrate + MYCrate) / 2
  )



# Add the BBCH values to the specific dates
ram_data <- ram_data %>%
  mutate(
    BBCH = case_when(
      date == "16/01/2024" ~ 15,
      date == "20/01/2024" ~ 21,
      date == "21/03/2024" ~ 30,
      date == "30/03/2024" ~ 32,
      date == "03/04/2024" ~ 35,
      date == "13/04/2024" ~ 41,
      date == "27/04/2024" ~ 51,
      date == "30/04/2024" ~ 55,
      date == "21/05/2024" ~ 73,
      date == "25/05/2024" ~ 85,
      date == "06/06/2024" ~ 90,
      TRUE ~ NA_real_
    )
  )

# Linear regression GS vs BBCH
bbch_wgsg_table <- ram_data %>%
  filter(!is.na(BBCH)) %>%
  select(BBCH, WGSg)

reg_model <- lm(BBCH ~ WGSg, data = bbch_wgsg_table)
a <- coef(reg_model)["WGSg"]
b <- coef(reg_model)["(Intercept)"]

# Calculate BBCHest
ram_data <- ram_data %>%
  mutate(
    BBCHest = a * WGSg + b
  )


# Calculate INFcum and Ponset based on the interval of dates from BBCH 15 to 32
# Definition of the dates of the beginning and the end of the interval in Date format
data_inizio <- as.Date("2024-01-16", format = "%Y-%m-%d")
data_fine <- as.Date("2024-03-30", format = "%Y-%m-%d")

ram_data <- ram_data %>%
  # Conversion of 'date' column into Date format
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  
  # Calculate INFcum e Ponset between the specified dates
  mutate(
    INFcum = ifelse(date >= data_inizio & date <= data_fine,
                    cumsum(coalesce(ifelse(date >= data_inizio, INFrate, 0), 0)),  # Cumulative sum of INFrate
                    NA_real_),  # NA al di fuori dell'intervallo delle date
    
    Ponset = ifelse(date >= data_inizio,
                    1 / (1 + 15075.146 * exp(-0.226 * BBCHest)),
                    NA_real_)  # NA al di fuori dell'intervallo delle date
  )



#### Graphs ####


# Graph 1: rolling_avg
ggplot(ram_data, aes(x = as.Date(date), y = rolling_avg)) +
  geom_col(na.rm = TRUE, fill = "steelblue") +  
  scale_y_continuous(limits = c(0, 6)) +        # Limits of y axis
  scale_x_date(date_breaks = "1 month", date_labels = "%d-%m-%y") +    # Prime date del mese
  labs(x = "Date",
       y = "Rolling Average (7dd)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Ruota le date


# Graph 2: INFrate e INFcum
ggplot(ram_data, aes(x = as.Date(date))) +
  geom_col(aes(y = INFrate), fill = "darkblue", na.rm = TRUE) +
  scale_y_continuous(name = "INFrate", limits = c(0, 1), 
                     sec.axis = sec_axis(~ . * 16, name = "INFcum")) +   # Asse secondario
  # Linea per INFcum sull'asse secondario
  geom_line(aes(y = INFcum / 16), color = "orange", size = 1, na.rm = TRUE) +  # Scala adattata per INFcum
  #geom_hline(yintercept = 0.9, color = "red", size = 1) +  # Linea orizzontale a 0.9 (asse primario)
  # annotate("text", x = as.Date("2024-03-01"), y = 0.92, label = "", color = "red") +  # Filtra la linea orizzontale solo nell'intervallo delle date
  scale_x_date(date_breaks = "1 month", date_labels = "%d-%m-%y") +
  labs(x = "Date",
       y = "INFrate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Graph 3: Ponset
ggplot(ram_data, aes(x = as.Date(date), y = Ponset)) +
  geom_line(color = "darkgreen", size = 1, na.rm = TRUE) +  
  scale_y_continuous(limits = c(0, 1)) +                 
  
  labs(x = "Date",
       y = "Ponset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Ruota le date





