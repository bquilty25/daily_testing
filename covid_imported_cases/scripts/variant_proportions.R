library(ggplot2)
#library(scales)
library(data.table)
#ibrary(wesanderson)

s_gene_dropout_data <- data.table(read.csv("data/sgtfvoc.csv")) %>% 
  .[, date := lubridate::ymd(date)] %>%
  .[, .(region = nhs_name, date, sgtfv)]

# s_gene_dropout_data %>% 
#   ggplot(aes(x = date, y = sgtfv, color = region)) + 
#   geom_line(show.legend = FALSE, alpha = 0.5) +
#   facet_wrap(~region) + 
#   theme_minimal() + 
#   scale_color_manual(values = wes_palette(7, name = "Darjeeling1", type = "continuous")) + 
#   scale_x_date(labels = date_format("%b")) 

s_gene_national_data <- s_gene_dropout_data[,  .(sgtfv = mean(sgtfv)), by = "date"] %>% 
  .[, "region" := "national"]
  
s_gene_data_final <- merge.data.table(s_gene_dropout_data, 
                                      s_gene_national_data,
                                      by = c("date", "region", "sgtfv"), 
                                      all = TRUE) %>%
  .[order(region, date)]

write.csv(s_gene_data_final, "data/sgtfv_clean.csv")
