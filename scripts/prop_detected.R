a <- traj_ %>% 
  filter.(test_label==TRUE) %>% 
  slice_min.(test_t,.by=c(variant,sim,idx)) %>% 
  ggplot()+
  geom_bar(aes(x=test_t,y=..prop..,fill=variant),
           position = position_dodge(width = 0.5))+
  labs(x="Days until detectable by LFT",fill="")

b <- traj_ %>% 
  filter.(test_label==TRUE) %>% 
  count.(variant,sim,idx) %>% 
  ggplot()+
  geom_bar(aes(x=N,y=..prop..,fill=variant),
           position = position_dodge(width = 0.5))+
  labs(x="Days of test positivity",fill="")

a/b&scale_x_continuous(breaks = breaks_width(1),limits = c(0,10))&scale_fill_brewer(palette = "Set2")
ggsave("duration_pre_post_detectable.png",width=120,height=100,dpi=600,units="mm")

traj_ %>% 
  summarise.(min_test_t=min(test_t==TRUE),
          max_test_t=max(test_t==TRUE),.by=c(variant,sim,idx)) %>% 
  mutate_rowwise.(
          test_1=rdunif(1:max_test_t,n = 1),
          test_2=test_1+3) %>% View()
          chance_detect=(max_test_t-min_test_t)/max_test_t,.by=c(variant,sim,idx)) %>% 
  group_by(variant) %>% 
  summarise_at(vars(chance_detect), p_funs)

