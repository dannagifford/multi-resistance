#install.packages(brms)
#install.packages(tidyverse)
#install.packages(future)

library(brms)
library(tidyverse)
library(future)

set.seed(19)

fixtable = function(matrix){
	matrix %>%
	as.data.frame() %>%
	rownames_to_column() %>%
	rename(parameter=rowname) %>%
	as_tibble()
}

# This number of cores is set for a 12-core system, leaving some cores free for other tasks
# 3 cores (3 workers) leads to 3 sessions of 3 cores each (in brm calls), meaning 9 cores total!
# Adjust these values for your system (i.e. don't let mc.cores*3 > your # of cores).
mc.cores = 4
plan(multisession, workers = 2)


# brms default behaviour is to load a file if it exists, rather than overwrite, and I find this confusing.
rerun.all.models = F
if(rerun.all.models == T) system("mv *brm_model.rds old/.")

version = "final"
# You could also append the time
#version = format(Sys.time(), format = "%Y-%m-%d-%H:%M:%S")

make.brmfile = function(name) paste(version, name, "brm_model", sep = "_")

# source("munging.R")
popns = readRDS("popns.Rds")

popnsAB = popns %>%
	mutate(state.simple = gsub(" ", "", state.simple)) %>% # Stan doesn't like spaces or underscores!!
	mutate(state.simple = factor(state.simple, c("noresistance","rifampicinresistance","nalidixicacidresistance","mixedresistance","doubleresistance"))) %>%
	filter(volume == 1) %>%
	mutate(row = gsub("[0-9]+", "", wellID), col = gsub("[A-Z]","", wellID)) %>%
	mutate(`no growth` = sensitive == 0) %>%
	mutate(plate = paste(pmutS.rank, antibiotic)) %>%
	mutate(pmutS_scaled = scale(pmutS)) %>%
	mutate(fday = as.factor(day)) %>%
	mutate(state.double = factor(`double resistance`, levels = c(0,1), labels = c("nodoubleresistance", "doubleresistance")))

popnsAB$pmutS.ordered = as.ordered(popnsAB$pmutS.text)

popnsday6 = popnsAB %>%
	filter(day == 6)

# Model M1: Detection of resistance arising during experimental evolution

priors = c(set_prior ("student_t(7, -5, 2.5)", class = "Intercept"),
		set_prior ("student_t(7, 0, 2.5)", class = "b"))

controls = list(adapt_delta = 0.99, max_treedepth = 15)

# Bernoulli models on double resistance were used to get a feel for the analysis
#brmname = "bernoulli1way"
#brmfile = make.brmfile(brmname)
#ber.1way %<-% {brm (state.double ~ (pmutS.text + antibiotic) + (1|row) + (1|col),
#	family = bernoulli("logit"), chains = 4, iter = 2000, warmup = 1000,
#	prior = priors, control = controls, sample_prior = "yes",
#	data = popnsday6, cores = 4, file = brmfile)}

#brmname = "bernoulli2way"
#brmfile = make.brmfile(brmname)
#ber.2way %<-% {brm (state.double ~ (pmutS.text + antibiotic)^2 + (1|row) + (1|col),
#	family = bernoulli("logit"), chains = 4, iter = 2000, warmup = 1000,
#	prior = priors, control = controls, sample_prior = "yes",
#	data = popnsday6, cores = 4, file = brmfile)}

brmname = "categorical2way"
brmfile = make.brmfile(brmname)
cat.2way %<-% {brm (state.simple ~ (pmutS.text + antibiotic)^2 + (1|row) + (1|col),
	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
	prior = priors, control = controls, sample_prior = "yes",
	data = popnsday6, cores = 4, file = brmfile)}

brmname = "categorical1way"
brmfile = make.brmfile(brmname)
cat.1way %<-% {brm (state.simple ~ (pmutS.text + antibiotic) + (1|row) + (1|col),
	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
	prior = priors, control = controls, sample_prior = "yes",
	data = popnsday6, cores = 4, file = brmfile)}

#WAIC(cat.1way, cat.2way) # WAIC on two objects is depricated in 2.8.0

cat.1way.waic <- waic(cat.1way)
cat.2way.waic <- waic(cat.2way)
loo_compare(cat.1way.waic, cat.2way.waic)

cat.1way.loo <- loo(cat.1way)
cat.2way.loo <- loo(cat.2way)
loo.compare1 = loo_compare(cat.1way.loo, cat.2way.loo)

# Note this is slow because it refits the model due to pareto_k>0.7
# It doesn't make a difference here
cat.1way.loo2 <- loo(cat.1way, reloo = T)
cat.2way.loo2 <- loo(cat.2way, reloo = T)
loo.compare2 = loo_compare(cat.1way.loo, cat.2way.loo)

# Beware: the models including day can take multiple hours/days to run
#brmname = "categorical3way_day"
#brmfile = make.brmfile(brmname)
#cat.day.3way %<-% {brm (state.simple ~ (pmutS.text + antibiotic + fday)^3 + (1|row) + (1|col),
#	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
#	prior = priors, control = controls, sample_prior = "yes",
#	data = popnsAB, cores = 4, file = brmfile)}

#brmname = "categorical2way_day"
#brmfile = make.brmfile(brmname)
#cat.day.2way %<-% {brm (state.simple ~ (pmutS.text + antibiotic + fday)^2 + (1|row) + (1|col),
#	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
#	prior = priors, control = controls, sample_prior = "yes",
#	data = popnsAB, cores = 4, file = brmfile)}


# Model M2: Fitness of fluctuation-derived strains
FTcurves.mean = readRDS("FTcurves_mean.Rds")

# brms doesn't like backticked names, so we rename these first
# brms 2.8.0 and later uses mvbind() instead of cbind()
brmname = "FT.model"
brmfile = make.brmfile(brmname)
controls = list(adapt_delta = 0.99, max_treedepth = 15)
priors = c(set_prior ("student_t(7, 10, 2.5)", class = "Intercept"),
		set_prior ("student_t(7, 0, 2.5)", class = "b"))

FT.model = FTcurves.mean %>%
	ungroup() %>%
	mutate(concentration_factor = factor(paste0("conc_", concentration))) %>%
	filter(!grepl("mutator", strain2)) %>%
	select(antibiotic, concentration_factor, id, strain, strain2, auc_e_mean) %>%
	spread(antibiotic, auc_e_mean) %>%
	mutate(nalidixicacid = `nalidixic acid`) %>%
	brm(mvbind(rifampicin, nalidixicacid, combination) ~ strain2 * concentration_factor,
	#brm(cbind(rifampicin, nalidixicacid, combination) ~ strain2 * concentration_factor,
		family = "student",		
		chains = 4, iter = 2000, warmup = 1000,
		control = controls, sample_prior = "yes", prior = priors,
		data = ., cores = 4, file = brmfile)
rifRvsWT = hypothesis(FT.model, hypothesis = "b_rifampicin_strain2rifampicinresistant = 0", class=NULL)
nalRvsWT = hypothesis(FT.model, hypothesis = "b_nalidixicacid_strain2nalidixicacidresistant = 0", class=NULL)
doublevsWT = hypothesis(FT.model, hypothesis = "b_combination_strain2doubleresistant = 0", class=NULL)

doubleRvsrifR = hypothesis(FT.model, hypothesis = "b_rifampicin_strain2doubleresistant - b_rifampicin_strain2rifampicinresistant = 0", class=NULL)
doubleRvsnalR = hypothesis(FT.model, hypothesis = "b_nalidixicacid_strain2doubleresistant - b_nalidixicacid_strain2nalidixicacidresistant = 0", class=NULL)


# Model M3: Fitness of double resistant strains in the presence and absence of the antibiotic combination
joinedcurves = readRDS("joinedcurves.Rds")

# brms doesn't like backticked names, so we rename these first
# brms 2.8.0 and later uses mvbind() instead of cbind()
brmname = "fitness.model"
brmfile = make.brmfile(brmname)
controls = list(adapt_delta = 0.99, max_treedepth = 15)
priors = c(set_prior ("student_t(7, 10, 2.5)", class = "Intercept"),
		set_prior ("student_t(7, 0, 2.5)", class = "b"))

fitness.model = joinedcurves %>%
	mutate(concentration = paste0("conc_", concentration)) %>%
	spread(concentration, auc_e) %>%
	rename(row = `Well Row`, col = `Well Col`) %>%
	brm(mvbind(conc_0, conc_2) ~ pmutS.text + (1|rep),
	#brm(cbind(conc_0, conc_2) ~ pmutS.text + (1|rep), # use if brms<2.8.0
		family = "student",		
		chains = 4, iter = 2000, warmup = 1000,
		control = controls, sample_prior = "yes", prior = priors,
		data = ., cores = 4, file = brmfile)

brmname = "fitness.model.intercept"
brmfile = make.brmfile(brmname)
controls = list(adapt_delta = 0.99, max_treedepth = 15)
priors = c(set_prior ("student_t(7, 10, 2.5)", class = "Intercept"))

fitness.model.intercept = joinedcurves %>%
	mutate(concentration = paste0("conc_", concentration)) %>%
	spread(concentration, auc_e) %>%
	rename(row = `Well Row`, col = `Well Col`) %>%
	brm(mvbind(conc_0, conc_2) ~ 1 + (1|rep),
	#brm(cbind(conc_0, conc_2) ~ 1 + (1|rep), # use if brms<2.8.0
		family = "student",		
		chains = 4, iter = 2000, warmup = 1000,
		control = controls, sample_prior = "yes", prior = priors,
		data = ., cores = 4, file = brmfile)

#WAIC(fitness.model, fitness.model.intercept) # WAIC on two objects is depricated in 2.8.0

waic1 <- waic(fitness.model)
waic2 <- waic(fitness.model.intercept)
loo_compare(waic1, waic2)

# Model M4: Compare data and simulation models
#source("simulations/plots.R") # Run this to generate muller.rds
muller.sim = readRDS("simulations/muller.rds")
popnsAB.sim = muller.sim %>%
	group_by(day) %>%
	filter(time==max(time)&rep<=60) %>%
	ungroup() %>%
	mutate(state.simple = state) %>%
	mutate(state.simple = gsub(" ", "", state.simple)) %>% # Stan doesn't like spaces or underscores!!
	mutate(state.simple = factor(state.simple, c("noresistance","rifampicinresistance","nalidixicacidresistance","mixedresistance","doubleresistance"))) %>%
	mutate(antibiotic = treatment) %>%
	mutate(fday = as.factor(day)) %>%
	mutate(pmutS.text = recode_factor(p, `0`="none", `0.05`="low", `0.1`="intermediate", `0.3`="high")) %>%
	mutate(`double resistance` = as.numeric(as.character(state.simple)=="double resistance")) %>%
	mutate(state.double = factor(`double resistance`, levels = c(0,1), labels = c("nodoubleresistance", "doubleresistance")))

popnsAB.sim$pmutS.ordered = as.ordered(popnsAB.sim$pmutS.text)

popnsday6.sim = popnsAB.sim %>%
	filter(day==6)
priors = c(set_prior ("student_t(7, -5, 2.5)", class = "Intercept"),
		set_prior ("student_t(7, 0, 2.5)", class = "b"))
controls = list(adapt_delta = 0.99, max_treedepth = 15)

brmname = "categorical2way_sim"
brmfile = make.brmfile(brmname)
cat.2way.sim %<-% {brm (state.simple ~ (pmutS.text + antibiotic)^2, 
	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
	prior = priors, control = controls, sample_prior = "yes",
	data = popnsday6.sim, cores = 4, file = brmfile)}

brmname = "categorical1way_sim"
brmfile = make.brmfile(brmname)
cat.1way.sim %<-% {brm (state.simple ~ pmutS.text + antibiotic, 
	family = categorical("logit"), chains = 4, iter = 2000, warmup = 1000,
	prior = priors, control = controls, sample_prior = "yes",
	data = popnsday6.sim, cores = 4, file = brmfile)}

compare_data_sim = left_join(fixtable(fixef(cat.1way)), fixtable(fixef(cat.1way.sim)),
	by = "parameter", suffix = c(" from experiment"," from simulation")) %>%
	separate(parameter, into = c("Response", "Predictor"), sep="_") %>%
	mutate(Response = gsub(x = Response, pattern = "mu", replacement = "")) %>%
	mutate(is.mutator = grepl("pmutS", Predictor)*1) %>%
	mutate(is.antibiotic = grepl("antibiotic", Predictor)*2) %>%
	mutate(is.interaction = grepl(":", Predictor)) %>%
	mutate(`Predictor type` = recode_factor(is.antibiotic + is.mutator,
		`0` = "intercept", `1` = "mutator frequency", `2` = "antibiotic treatment", `3` = "interaction")) %>%
	mutate(Response = gsub(x = Response, pattern = "resistance", replacement = " resistance")) %>%
	mutate(Response = gsub(x = Response, pattern = "nalidixicacid", replacement = "nalidixic acid")) %>%
	mutate(Response = recode_factor(Response, `rifampicin resistance`="rifampicin resistance",
				`nalidixic acid resistance`="nalidixic acid resistance",
				`mixed resistance`="mixed resistance",
				`double resistance`="double resistance")) %>%
	mutate(Predictor = gsub(x = Predictor, pattern = "pmutS.text", replacement = "")) %>%
	mutate(Predictor = gsub(x = Predictor, pattern = "antibiotic", replacement = "")) %>%
	mutate(Predictor = gsub(x = Predictor, pattern = "nalidixicacid", replacement = "nalidixic acid")) %>%
	mutate(Predictor = factor(Predictor, levels = c(
		"Intercept",
		"low", "intermediate", "high",
		"rifampicin", "nalidixic acid", "combination",
		"low:rifampicin", "low:nalidixic acid", "low:combination",
		"intermediate:rifampicin", "intermediate:nalidixic acid", "intermediate:combination",
		"high:rifampicin", "high:nalidixic acid", "high:combination")))

saveRDS(compare_data_sim, "compare_data_sim.Rds")
