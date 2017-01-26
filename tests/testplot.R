
library(PKAEvis)
library(tidyverse)

load("data/ae.sim.profile.rda")
load("data/pk.sim.profile.rda")
load("data/pk.sim.profile.sparse.rda")
load("data/subject.sim.rda")

if(0){
  # Settings for test of pkaeplot

  pk <- pk.sim.profile.sparse
  ae <- ae.sim.profile
  subj<-subject.sim
  ae.data.first.day <- 1
  scale.y.log10=T
  x.range=NULL
  y.range=NULL
  ae.col.var= "AEGR"
  ae.col.name=NULL
  pk.col.var="DOSE"
  #pk.col.var=NULL
  pk.col.name=NULL

  ae.palette=c("#56B4E9","#0072B2","#D55E00")
  ggtheme=NULL
}

theme.to.use <-
  theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_rect(fill="#fff5f0"),
        axis.title.y=element_text(size=15, margin=margin(0, 10, 0, 0)),
        axis.title.x=element_text(size=15, margin=margin(10, 0, 0, 0)),
        legend.position="bottom"
  )



glist <-
  pkaeplot(pk = pk.sim.profile,
         ae = ae.sim.profile,
         subj=subject.sim,
         ae.data.first.day = 1,
         y.range=c(3,300),
         ae.col.var="AEGR",
         ae.col.name="AE grade",
         pk.col.var="DOSE",
         ggtheme=theme.to.use)

combine_pkaeplot(glist)

