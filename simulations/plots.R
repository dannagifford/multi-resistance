#install.packages(tidyverse)
#install.packages(ggpubr)
#install.packages(tidyverse)
#install.packages(grid)
#install.packages(scales)

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


scientific_10x <- function(values, digits = 1) {
 if(!is.numeric(values)){
 stop("values must be numbers")
 }
 if(grepl("^\\d{2}$", digits)){
 stop("digits must a one or two digit whole number")
 }

 x_zeros <- which(values == 0)

 x <- sprintf(paste0("%.", digits, "e"), values)

 x <- gsub("^(.*)e", "'\\1'e", x)

 longestExponent <- max(sapply(gregexpr("\\d{1, }$", x), attr, 'match.length'))
 zeroTrimmed <- ifelse(longestExponent > 2, 
 paste0("\\1", paste(rep("~", times = longestExponent-1), collapse = "")), 
 "\\1")
 x <- gsub("(e[ + |-])[0]", zeroTrimmed, x)

 x <- gsub("e", "~x~10^", x)

 if(any(grepl("\\^\\-", x))){
 x <- gsub("\\^\\ + ", "\\^~~", x)
 } else {
 x <- gsub("\\^\\ + ", "\\^", x)
 }

 x[x_zeros] <- "0"
 # return this as an expression
 parse(text = x)
} 


mean_sd = function(x, mult = 1) { 
 x <- na.omit(x)
 sd <- sd(x)
 mean <- mean(x)
 data.frame(y = mean, ymin = max(mean - sd, 0), ymax = mean + sd)
}

quantiles = function(x, probs = c(0.25, 0.5, 0.75)) {
	setNames(quantile(x, probs = probs), nm = c("ymin", "y", "ymax"))}

mads = function(x) {
	tibble(y = median(x), ymin = y - mad(x), ymax = y + mad(x))}

# Munging
results = tibble(files = list.files(pattern = ".txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
		mutate(time = rep(1:nsteps, times = nreps), rep = rep(1:nreps, each = nsteps)) %>%
		mutate(day = (time-1) %/% spd + 1) %>%
		mutate(concentration = 10*(1/16)*2^(day-1)))) %>%
	unnest() %>%
	separate(files, into = c("treatment", "p"), sep = "/") %>%
	mutate(p = gsub("([a-z]+_p=|.txt)", "", p)) %>%
	type_convert() %>%
	mutate(St = S + mut_S, At = A + mut_A, Bt = B + mut_B, Dt = D + mut_D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (p = 0)", 
		`0.05` = "low (p = 0.05)", 
		`0.1` = "intermediate (p = 0.1)", 
		`0.3` = "high (p = 0.3)", 
		`0.5` = "extrahigh (p = 0.5)")) %>%
	filter(p <= 0.3, day <= 6) %>%
	mutate(treatment = recode_factor(treatment, none = "none", rif = "rifampicin", nal = "nalidixic acid", both = "combination"))

results.gathered = results %>%
	select(-c(St, At, Bt, Dt)) %>%
	gather(key = "strain", value = "n", -c(treatment, rep, p, pmutS.text, day, concentration, time)) %>%
	mutate(mutator = as.factor(ifelse(!grepl("mut_", strain), "wild-type", "mutator"))) %>%
	mutate(mutator = relevel(mutator, "wild-type")) %>%
	mutate(strain = recode_factor(strain, S = "sensitive", A = "rifampicin resistant", B = "nalidixic acid resistant", D = "double resistant", mut_S = "sensitive", mut_A = "rifampicin resistant", mut_B = "nalidixic acid resistant", mut_D = "double resistant"))


states = list(St = "no resistance", 
		At = "rifampicin resistance", 
		Bt = "nalidixic acid resistance", 
		Mt = "mixed resistance", 
		Dt = "double resistance")

muller = results %>%
	select(treatment, day, concentration, p, pmutS.text, rep, time, St, At, Bt, Dt) %>%
	mutate(total = St + At + Bt + Dt) %>%
	gather(state, value, -c(treatment, p, pmutS.text, time, day, rep, concentration, total)) %>%
	mutate(detected = map2_lgl(total, value, ~rbinom(n = 1, size = .x, prob = .y/(.x*200))>0)) %>%
	select(-value, -total) %>%
	spread(state, detected) %>%
	mutate(Mt = At&Bt, St = ifelse(At|Bt|Mt|Dt, F, St), At = ifelse(Mt|Dt, F, At), Bt = ifelse(Mt|Dt, F, Bt), Mt = ifelse(Dt, F, Mt)) %>%
	gather(state, value, -c(treatment, day, concentration, p, pmutS.text, rep, time)) %>%
	filter(value == T) %>%
	select(-value) %>%
	mutate(state = factor(unlist(states[state]), levels = states)) %>%
	arrange(treatment, day, time, rep) %>%
	droplevels()



# Plots
bestcolours = c("grey90", "#469067", "#235161", "#6f2c3d", "#ef0358")
bestcolours2 = c("black", bestcolours[-1])


n.plot = results.gathered %>%
	filter() %>%
	ggplot(aes(x = time, y = n + 1, fill = strain, colour = strain, linetype = mutator)) + 
		stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = 0.4) + 
		facet_grid(pmutS.text~treatment) + 
		scale_y_log10(name = "Number of individuals", sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initial mutatorr frequency")) + 
		scale_colour_manual(name = "Strain", values = c(bestcolours2[-4], bestcolours2[-4])) + 
		scale_fill_manual(name = "Strain", values = c(bestcolours2[-4], bestcolours2[-4])) + 
		guides(linetype = guide_legend(title = "Mutation rate", override.aes = list(colour = "black", fill = NA))) + 
		scale_x_continuous(name = "Time (days)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0), 

			sec.axis = sec_axis(~ ., breaks = (unique(results$day)-0.5)*spd, 

			labels = unique(results.gathered$concentration), name = "Simulated antibiotic concentration (mg/l)")) + 
		geom_vline(xintercept = spd*(1:6) + 2, color = "grey80") + 
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0))

n.plot = ggplotGrob(n.plot)
n.plot[["grobs"]][[18]][["children"]][[2]] = nullGrob()
n.plot = as_ggplot(n.plot)

results.gathered.mean = results.gathered %>% group_by(treatment, p, time, day, pmutS.text, strain, mutator) %>% summarise(mean = mean(n + 1), sd = sd(n + 1))

#median absolute deviation

plotinator = function(data, suppress.strain = T){
output = ggplot(data = data, aes(x = time, y = n + 1, fill = strain, color = strain, linetype = mutator, alpha = mutator)) + 
#		stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = 0.4) + 
		stat_summary(fun.data = quantiles, geom = "ribbon", alpha = 0.5) + 
		facet_grid(strain~treatment) + 
		scale_y_log10(name = "Number of simulated bacteria\n(interquartile range)", labels = trans_format("log10", math_format(10^.x))) + 
		scale_colour_manual(name = "Strain", values = c(bestcolours2[-4], bestcolours2[-4])) + 
		scale_fill_manual(name = "Strain", values = c(bestcolours2[-4], bestcolours2[-4])) + 
		guides(fill = F, color = F, linetype = F) + 
		scale_x_continuous(name = "Time (days)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0), 

			sec.axis = sec_axis(~ ., breaks = (unique(results$day)-0.5)*spd, 

			labels = unique(muller$concentration), name = "Simulated antibiotic concentration (mg/l)")) + 
		scale_linetype_manual(values = c("44", "solid")) + 
#		geom_vline(xintercept = spd*(1:6) + 2, color = "grey80") + 
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0))

if(suppress.strain) output = output + theme(strip.text.y = element_blank())


output = ggplotGrob(output)
output[["grobs"]][[18]][["children"]][[2]] = nullGrob()
output = as_ggplot(output)
return(output)
}

#n.plot.combo = results.gathered %>%
#	filter(treatment == "combination") %>%
#	plotinator(.)

n.plot.separated.iqr = results.gathered %>%
	group_by(p) %>%
	nest() %>%
	mutate(plots = map(data, plotinator))

n.plot.maintext = results.gathered %>%
	filter(p==0.1) %>%
	plotinator(., suppress.strain = T)

n.plot.all.iqr = ggarrange(plotlist = n.plot.separated.iqr$plots, labels = "AUTO")



# Muller-like plot equivalent to experiment

# Muller-like plot equivalent to experiment
muller.sim.plot = muller %>%
	ggplot(aes(x = time, fill = state)) +
		geom_bar(width = 1.2, size = 0) +
		facet_grid(pmutS.text~treatment) +
		scale_fill_manual(values = bestcolours) +
		scale_x_continuous(name = "Time (days)", breaks = (unique(results$day)-0.5)*spd, labels = unique(results$day), limits = c(0, nsteps-1), expand = expand_scale(mult = 0, add = 0), 
			sec.axis = sec_axis(~ ., breaks = (unique(results$day)-0.5)*spd, 
			labels = unique(muller$concentration), name = "Simulated antibiotic concentration (mg/l)")) +
		scale_y_continuous(name = "Number of simulated populations", 
			sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initial mutator frequency")) +
		geom_vline(xintercept = (1:6)*spd+1, size = 0.2) +
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0)) +
		guides(fill = guide_legend(title = "Detection of", override.aes = list(color = "black", size = 0.2)))

muller.sim.plot = ggplotGrob(muller.sim.plot)
muller.sim.plot[["grobs"]][[18]][["children"]][[2]] = nullGrob()
muller.sim.plot = as_ggplot(muller.sim.plot)


day6.sim.plot = muller %>%
	filter(time == max(muller$time)) %>%
	group_by(treatment, p, state) %>%
	summarise(count = n()) %>%
	complete(treatment, p, state, fill = list(count = 0)) %>%
	ggplot(aes(x = p, y = count, colour = state, fill = state)) + 
		geom_point(aes(shape = state), size = 2, colour = "black") + 
		geom_line() + 
		facet_grid(~treatment) + 
		scale_colour_manual(name = "Detection of", values = bestcolours2) + 
		scale_fill_manual(name = "Detection of", values = bestcolours2) + 
		scale_shape_manual(values = c(21, 22, 23, 25, 24), name = "Detection of") + 
		labs(x = 'Initial mutator frequency', y = "Number of populations") + 
		geom_hline(yintercept = 1000, linetype = "dashed", color = "grey40")

expevol = ggarrange(muller.sim.plot, day6.sim.plot, nrow = 2, heights = c(3, 1.25), labels = "AUTO", font.label = list(size = 10, face = "plain"))

ggsave("muller_sim.pdf", muller.sim.plot, width = 8, height = 6)
ggsave("nplot_sim_main.pdf", n.plot.maintext, width = 5, height = 5)
ggsave("nplot_sim_supplement.pdf", n.plot.all.iqr, width = 8, height = 8)
#ggsave("nplot_combo_sim.pdf", n.plot.combo, width = 8, height = 5)
ggsave("day6_sim.pdf", day6.sim.plot, width = 8, height = 2.5)
ggsave("expevol.pdf", expevol, width = 8, height = 7.5)



params = tibble(files = list.files(pattern = "[a-z].ini$", recursive = T), 
	inis = map(files, ~read_table(file = ., col_names = "parameter"))) %>%
	separate(files, into = "treatment") %>%
	unnest() %>%
	separate(parameter, into = c("parameter", "value"), sep = " = ") %>%
	type_convert() %>%
	filter(!grepl("\\[", parameter)) 

params$day = rep(c(rep(NA, 9), rep(1:6, each = 20)), times = 4)

params = params %>%
	filter(!is.na(day), day>0) %>%
	spread(key = parameter, value = value) %>%
	mutate(treatment = recode_factor(treatment, none = "none", rif = "rifampicin", nal = "nalidixic acid", both = "combination"))

write_csv(params, "simulation_parameters.csv")
