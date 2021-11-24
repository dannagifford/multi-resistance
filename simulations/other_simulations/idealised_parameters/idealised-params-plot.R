require(tidyverse)
require(ggpubr)
require(tidyverse)
require(grid)
require(scales)


nreps = 1000
dt = 1
nsteps = 22*6*(1/dt) + 1 # 133
spd = 1/dt*22 #steps per day

theme_set(theme_bw())
theme_update(text = element_text(size = 10), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank()
)

quantiles = function(x, probs = c(0.25, 0.5, 0.75)) {
	setNames(quantile(x, probs = probs), nm = c("ymin", "y", "ymax"))}

# Munging
results = data_frame(files = list.files(pattern = ".txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
		mutate(time = rep(1:nsteps, times = nreps), rep = rep(1:nreps, each = nsteps)) %>%
		mutate(day = (time-1) %/% spd + 1) %>%
		mutate(concentration = 10*(1/16)*2^(day-1)))) %>%
	unnest() %>%
	separate(files, into = c("treatment", "p"), sep = "/") %>%
	mutate(p = gsub("([A-Za-z=_]|.txt)", "", p)) %>%
	type_convert() %>%
	mutate(St = S + mut_S, At = A + mut_A, Bt = B + mut_B, Dt = D + mut_D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (u = 0)", 
		`0.05` = "low (u = 0.05)", 
		`0.1` = "intermediate (u = 0.1)", 
		`0.3` = "high (u = 0.3)", 
		`0.5` = "extrahigh (u = 0.5)")) %>%
	filter(p <= 0.3, day <= 6) %>%
	mutate(treatment = recode_factor(treatment, S = "none", A = "drug A", B = "drug B", D = "combination"))

states = list(St = "no resistance", 
		At = "A resistance", 
		Bt = "B resistance", 
		Mt = "mixed resistance", 
		Dt = "D resistance")

results.gathered = results %>%
	select(-c(St, At, Bt, Dt)) %>%
	gather(key = "strain", value = "n", -c(treatment, rep, p, pmutS.text, day, concentration, time)) %>%
	mutate(mutator = as.factor(ifelse(!grepl("mut_", strain), "wild-type", "mutator"))) %>%
	mutate(mutator = relevel(mutator, "wild-type")) %>%
#	mutate(strain = recode_factor(strain, S = "S sensitive", A = "A resistant", B = "B resistant", D = "D double resistant", mut_S = "S sensitive", mut_A = "A resistant", mut_B = "B resistant", mut_D = "D double resistant"))
	mutate(strain = recode_factor(strain, S = "S", A = "A", B = "B", D = "D", mut_S = "S", mut_A = "A", mut_B = "B", mut_D = "D"))



muller = results %>%
	select(treatment, day, concentration, p, pmutS.text, rep, time, St, At, Bt, Dt) %>%
	mutate(total = St + At + Bt + Dt) %>%
	gather(state, value, -c(treatment, p, pmutS.text, time, day, rep, concentration, total)) %>%
	mutate(detected = map2_lgl(total, value, ~rbinom(n = 1, size = .x, prob = .y/(.x*200))>0)) %>%
	select(-value, -total) %>%
	spread(state, detected) %>%
	mutate(Mt = At&Bt, St = ifelse(At|Bt|Mt|Dt, F, 1), At = ifelse(Mt|Dt, F, At), Bt = ifelse(Mt|Dt, F, Bt), Mt = ifelse(Dt, F, Mt)) %>%
	gather(state, value, -c(treatment, day, concentration, p, pmutS.text, rep, time)) %>%
	filter(value == T) %>%
	select(-value) %>%
	mutate(state = factor(unlist(states[state]), levels = states)) %>%
	arrange(treatment, day, time, rep) %>%
	droplevels()

saveRDS(muller, "muller.rds")



# Plots
#bestcolours = c("white", "grey80", "grey50", "grey30", "black")
#bestcolours2 = c("black", bestcolours[-1])

#bestcolours.arbitrary = c("white", "grey80", "grey50", "grey40", "black")
#bestcolours2.arbitrary = c("black", bestcolours.arbitrary[-1])

bestcolours.arbitrary = c("white", "#6d8755", "#48839f", "#a2639b", "#463255")
bestcolours2.arbitrary = c("black", "#6d8755", "#48839f", "#a2639b", "#463255")



n.plot = results.gathered %>%
	filter() %>%
	ggplot(aes(x = time, y = n + 1, fill = strain, colour = strain, linetype = mutator)) + 
		stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = 0.4) + 
		facet_grid(pmutS.text~treatment) + 
		scale_y_log10(name = "Number of individuals", sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initial mutator frequency")) + 
		scale_colour_manual(name = "Strain", values = c(bestcolours2.arbitrary[-4], bestcolours2.arbitrary[-4])) + 
		scale_fill_manual(name = "Strain", values = c(bestcolours2.arbitrary[-4], bestcolours2.arbitrary[-4])) + 
		guides(linetype = guide_legend(title = "Mutation rate", override.aes = list(colour = "black", fill = NA))) + 
		scale_x_continuous(name = "Time (# transfers)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0), 

			sec.axis = sec_axis(~ ., breaks = (unique(results$day)-0.5)*spd, 

			labels = unique(muller$concentration), name = "Simulated antibiotic concentration (mg/l)")) + 
		geom_vline(xintercept = spd*(1:6) + 2, color = "grey80") + 
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0))

#n.plot = ggplotGrob(n.plot)
#n.plot[["grobs"]][[18]][["children"]][[2]] = nullGrob()
#n.plot = as_ggplot(n.plot)

results.gathered.mean = results.gathered %>% group_by(treatment, p, time, day, pmutS.text, strain, mutator) %>% summarise(mean = mean(n + 1), sd = sd(n + 1))

#median absolute deviation

plotinator = function(data, suppress.strain = F){
output = ggplot(data = data, aes(x = time, y = n + 1, fill = mutator, colour = mutator)) + 
#		stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = 0.4) + 
		stat_summary(fun.data = quantiles, geom = "ribbon", alpha = 0.4) + 
		facet_grid(strain~treatment) + 
		scale_y_log10(name = "Number of simulated bacteria\n(interquartile range)", labels = trans_format("log10", math_format(10^.x))) + 
		scale_fill_manual(name = "Genetic background", values = c("grey60","black")) +
		scale_colour_manual(name = "Genetic background", values = c("grey60","black")) +
#		scale_colour_manual(name = "Strain", values = c(bestcolours2.arbitrary[-4], bestcolours2.arbitrary[-4])) + 
#		scale_fill_manual(name = "Strain", values = c(bestcolours2.arbitrary[-4], bestcolours2.arbitrary[-4])) + 
#		guides(fill = F, color = F, linetype = F) + 
		scale_x_continuous(name = "Time (# transfers)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0)) + 
		scale_linetype_manual(values = c("44", "solid")) + 
#		geom_vline(xintercept = spd*(1:6) + 2, color = "grey80") + 
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0))

if(suppress.strain) output = output + theme(strip.text.y = element_blank())


#output = ggplotGrob(output)
#output[["grobs"]][[18]][["children"]][[2]] = nullGrob()
#output = as_ggplot(output)
return(output)
}

n.plot.separated.iqr = results.gathered %>%
	mutate(treatment = recode_factor(treatment, 
		none = "none", `drug A` = "drug A" , `drug B` = "drug B", combination = "combo")) %>%
	group_by(p) %>%
	nest() %>%
	mutate(plots = map(data, plotinator))
	
n.plot.all.iqr = ggarrange(plotlist = n.plot.separated.iqr$plots, labels=LETTERS[3:6], common.legend = T, legend = "right")




# Muller-like plot equivalent to experiment
muller.sim.plot = muller %>%
	mutate(treatment = recode_factor(treatment, 
		none = "none", `drug A` = "drug A" , `drug B` = "drug B", combination = "combo")) %>%
	ggplot(aes(x = time, fill = state)) +
		geom_bar(width = 1.2, size = 0) +
		facet_grid(pmutS.text~treatment) +
		scale_fill_manual(values = bestcolours.arbitrary) +
		scale_x_continuous(name = "Time (# transfers)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0)) +
		scale_y_continuous(name = "Number of simulated populations", 
			sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initial mutator frequency")) +
		geom_vline(xintercept = (1:6)*spd+1, size = 0.2) +
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0)) +
		guides(fill = guide_legend(title = "Detection of", override.aes = list(color = "black", size = 0.2)))

muller.sim.plot = ggplotGrob(muller.sim.plot)
muller.sim.plot[["grobs"]][[18]][["children"]][[2]] = nullGrob()
muller.sim.plot = as_ggplot(muller.sim.plot)

ggsave("arbitrary_muller_sim.pdf", muller.sim.plot, width = 8, height = 6)

params = data_frame(files = list.files(pattern = "[SABD].ini$", recursive = T), 
	inis = map(files, ~read_table(file = ., col_names = "parameter"))) %>%
	separate(files, into = "treatment") %>%
	unnest() %>%
	separate(parameter, into = c("parameter", "value"), sep = " = ") %>%
	type_convert() %>%
	filter(!is.na(value)) 

params$day = rep(c(rep(NA, 9), rep(1:6, each = 16)), times = 4)

params = params %>%
	filter(!is.na(day), day>0) %>%
	spread(key = parameter, value = value) %>%
	mutate(treatment = recode_factor(treatment, S = "none", A = "drug A", B = "drug B", D = "combination"))

params_plot = params %>%
	gather(key, value, -c(treatment, day)) %>%
	mutate(param_class = gsub("[SABD]","",key), strain = gsub("[kr12]","",key)) %>%
	mutate(strain = recode_factor(strain, S = "S sensitive", A = "A resistant", B = "B resistant", D = "D resistant")) %>%
	filter(grepl("1", param_class)) %>%
	mutate(param_class = recode(param_class, k1="k1=k2", r1="r1=r2")) %>%
	select(-key) %>%
	unique() %>%
	ggplot(aes(x=day, y=value, fill = strain, shape = strain)) +
	geom_line(aes(colour = strain), position = position_dodge(width=0.3)) +
	geom_point(size = 2, position = position_dodge(width=0.3)) +
	facet_grid(param_class~treatment, scales = "free_y") + 
	scale_colour_manual(values = bestcolours2.arbitrary[c(1,2,3,5)], name = "Simulated strain with\n arbitrary parameters") +
	scale_fill_manual(values = bestcolours.arbitrary[c(1,2,3,5)], name = "Simulated strain with\n arbitrary parameters") +
	scale_shape_manual(values = c(21, 22, 23, 25, 24), name = "Simulated strain with\n arbitrary parameters") +
	labs(x = "Time (# transfers)", y = "Parameter value")

#ggsave("arbitrary_params.pdf", params_plot, width = 6, height = 3)	


#n.plot.withparam = ggarrange(params_plot, n.plot.all.iqr, labels = c("A",""), nrow = 2, heights = c(.5,1))
#ggsave("arbitrary_nplot.pdf", n.plot.withparam, width = 7, height = 8)

arbitrary.sim.plot = ggarrange(
ggarrange(params_plot, muller.sim.plot, nrow = 2, labels = c("A","B","C"), heights = c(.60,.9)),
n.plot.all.iqr, nrow = 2)

ggsave("arbitrary_simplot.pdf", arbitrary.sim.plot, width = 7, height = 10)

write_csv(params, "arbitrary_simulation_parameters.csv")




