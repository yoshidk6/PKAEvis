library(tidyverse)
library(PKPDsim)


#### Prepare PK dataset using PKPDsim package ####
n.sub.active <- 6
n.sub.placebo<- 12
dose    <- c(60,200,400)
n.dose  <- length(dose)

model <- new_ode_model(model = "pk_1cmt_oral")

param <-
  list(KA = 0.3, CL = 0.2, V = 5)

# Simulation
set.seed(3)
pk.sim.profile <- data.frame() %>% tbl_df()

for(kdose in 1:length(dose)){

  regim <-
    new_regimen(amt = dose[kdose],
                n = 3,
                interval = 0,
                type = "oral")

  ## Run ODE
  pk.sim.profile.each <- sim(
    ode = model, parameters = param, regimen = regim,
    t_obs = seq(0,70,by=0.5),
    n_ind = n.sub.active,
    only_obs=T,
    omega = c(0.1,
              0.05, 0.1,
              0.01, 0.05, 0.1)) %>%
    tbl_df() %>%
    select(ID=id,DAY=t,DV=y) %>%
    mutate(ID=ID + n.sub.active*(kdose-1))

  pk.sim.profile <- bind_rows(pk.sim.profile,pk.sim.profile.each)
}


# Sparse dataset to try calculating middle points
pk.sim.profile.sparse <-
  pk.sim.profile %>%
  filter(DAY%%2 == 0)

#ggplot(dat,aes(DAY,DV,group=ID)) + geom_line()


#### Prepare subject dataset ####
subject.sim <-
  data.frame(ID=1:(n.sub.active*n.dose+n.sub.placebo),
             DOSE=c(rep(dose,each=n.sub.active),
                    rep(0,n.sub.placebo)),
             ARM=c(rep("Drug",n.sub.active*n.dose),
                   rep("Placebo",n.sub.placebo)))


#### Prepare AE dataset ####
ae.sim.profile <-
  read_csv("data-raw/ae.manual.data.csv",
           col_types = cols(
             ID = col_integer(),
             AE = col_character(),
             AEGR = col_factor(levels=c(1,2,3)),
             AESTDY = col_integer(),
             AEEDDY = col_integer()
           ))


#### Export dataset ####
devtools::use_data(pk.sim.profile,
                   ae.sim.profile,
                   pk.sim.profile.sparse,
                   subject.sim,
                   overwrite = TRUE)

