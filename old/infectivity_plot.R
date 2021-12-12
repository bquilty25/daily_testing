source("packages.R")
source("utils.R")
source("parameters.R")
source("plot_functions.R")

df <- data.frame(x=-15:15)

ashcroft <- function(x,shape,rate,shift){
   dgamma(x=x+shift,shape = shape,rate=rate) 
}

ggplot(data=df,aes(x=x))+
  stat_function(fun=ashcroft, args=list(shape=infect_shape, rate=infect_rate,shift=infect_shift))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA))+
  labs(x="Days since symptom onset",y="Density")

ggsave("results/ashcroft.png",width=210,height=100,units="mm",dpi=320)
