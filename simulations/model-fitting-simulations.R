library(brms)
library(tidyverse)
library(future)

muller.model = readRDS("muller_sim.rds")
popnsAB.model = muller.model %>%
	group_by(day) %>%
	filter(time==max(time)&rep<=60) %>%
	ungroup() %>%
	mutate(state.simple = state) %>%
	mutate(antibiotic = treatment) %>%
	mutate(fday = as.factor(day)) %>%
	mutate(pmutS.text = recode_factor(p, `0`="none", `0.05`="low", `0.1`="intermediate", `0.3`="high")) %>%
	mutate(`double resistance` = as.numeric(as.character(state.simple)=="double resistance")) %>%
	mutate(state.double = factor(`double resistance`, levels=c(0,1), labels = c("no double resistance", "double resistance")))

popnsAB.model$pmutS.ordered = as.ordered(popnsAB.model$pmutS.text)

popnsday6.model = popnsAB.model %>%
	filter(day==6)

# Model SI
set.seed(19)

# This number of cores is set for a 12-core system, leaving some cores free for other tasks
# 3 cores (3 workers) leads to 3 sessions of 3 cores each (in brm calls), meaning 9 cores total!
# Adjust these values for your system (i.e. don't let mc.cores*3 > your # of cores).
mc.cores = 3
plan(multisession, workers = 3)

priors = c(set_prior ("student_t(7,-5,2.5)", class="Intercept"),
		set_prior ("student_t(7,0,2.5)", class="b"))

controls = list(adapt_delta = 0.99, max_treedepth=15)

today = format(Sys.time(), format="%Y-%m-%d-%H:%M:%S")

#brmname = "bernoulli1way"
#brmfile = paste0(today, "_", brmname)
#ber_1way %<-% {brm (state.double~(pmutS.text + antibiotic),
#	family=bernoulli("logit"), chains=3, iter=5000, warmup=1000,
#	prior=priors, control=controls,
#	data=popnsday6.model, cores=3, file=brmfile)}

#brmname = "bernoulli2way"
#brmfile = paste0(today, "_", brmname)
#ber_2way %<-% {brm (state.double~(pmutS.text + antibiotic)^2 ,
#	family=bernoulli("logit"), chains=3, iter=5000, warmup=1000,
#	prior=priors, control=controls,
#	data=popnsday6.model, cores=3, file=brmfile)}

brmname = "categorical2waymainsim"
brmfile = paste0(today, "_", brmname)
#cat_2way_sim %<-% {brm (state.simple~(pmutS.text + antibiotic),
cat_2way_sim = brm (state.simple~(pmutS.text + antibiotic),
	family=categorical("logit"), chains=3, iter=5000, warmup=1000,
	prior=priors, control=controls,
	data=popnsday6.model, cores=3, file=brmfile)

brmname = "categorical2wayinteraction"
brmfile = paste0(today, "_", brmname)
#cat_2way_interaction_sim %<-% {brm (state.simple~(pmutS.text + antibiotic)^2,
cat_2way_interaction_sim = brm (state.simple~(pmutS.text + antibiotic)^2,
	family=categorical("logit"), chains=3, iter=5000, warmup=1000,
	prior=priors, control=controls,
	data=popnsday6.model, cores=3, file=brmfile)

#brmname = "categorical3way_day"
#brmfile = paste0(today, "_", brmname)
#cat_3wayday %<-% {brm (state.simple~(pmutS.text + antibiotic + fday)^3,
#	family=categorical("logit"), chains=3, iter=5000, warmup=1000,
#	prior=priors, control=controls,
#	data=popnsAB.model, cores=3, file=brmfile)}

#brmname = "categorical2way_day"
#brmfile = paste0(today, "_", brmname)
#cat_2wayday %<-% {brm (state.simple~(pmutS.text + antibiotic + fday)^2,
#	family=categorical("logit"), chains=3, iter=5000, warmup=1000,
#	prior=priors, control=controls,
#	data=popnsAB.model, cores=3, file=brmfile)}


cat_2way_sim.waic <- waic(cat_2way_sim)
cat_2way_interaction_sim.waic <- waic(cat_2way_interaction_sim)
loo(cat_2way_sim,cat_2way_interaction_sim, compare=T, reloo = T)

