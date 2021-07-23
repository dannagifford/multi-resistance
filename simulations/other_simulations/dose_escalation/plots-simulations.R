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

# Munging
results = data_frame(files = list.files(pattern = "both_[1-5]_p.*.txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
		mutate(time = rep(1:nsteps, times = nreps), rep = rep(1:nreps, each = nsteps)) %>%
		mutate(day = (time-1) %/% spd + 1))) %>%
	unnest() %>%
	separate(files, into = c("treatment","ramp", "p"), sep = "_") %>%
	mutate(p = gsub("([A-Za-z=_]|.txt)", "", p)) %>%
	type_convert() %>%
	mutate(St = S + mut_S, At = A + mut_A, Bt = B + mut_B, Dt = D + mut_D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (q = 0)", 
		`0.05` = "low (q = 0.05)", 
		`0.1` = "intermediate (q = 0.1)", 
		`0.3` = "high (q = 0.3)", 
		`0.5` = "extrahigh (q = 0.5)")) %>%
	filter(p <= 0.3, day <= 6) %>%
	mutate(treatment = recode_factor(treatment,  both = "combination")) %>%
	mutate(MICday = 6-ramp) %>%
	select(-ramp)

nsteps2 = 22*12*(1/dt) + 1 # 133

results_double = data_frame(files = list.files(pattern = "both_10.+\\.txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
	mutate(time = rep(1:nsteps2, times = nreps), rep = rep(1:nreps, each = nsteps2)) %>%
	mutate(day = (time-1) %/% spd + 1))) %>%
	unnest() %>%
	separate(files, into = c("treatment","ramp", "p"), sep = "_") %>%
	mutate(p = gsub("([A-Za-z=_]|.txt)", "", p)) %>%
	type_convert() %>%
	mutate(St = S + mut_S, At = A + mut_A, Bt = B + mut_B, Dt = D + mut_D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (q = 0)", 
		`0.05` = "low (q = 0.05)", 
		`0.1` = "intermediate (q = 0.1)", 
		`0.3` = "high (q = 0.3)", 
		`0.5` = "extrahigh (q = 0.5)")) %>%
	filter(p <= 0.3, day <= 12) %>%
	mutate(treatment = recode_factor(treatment,  both = "combination")) %>%
	mutate(MICday = 10) %>%
	select(-ramp)


nsteps4 = 22*24*(1/dt) + 1 

results_quadruple = data_frame(files = list.files(pattern = "both_20.+\\.txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
	mutate(time = rep(1:nsteps4, times = nreps), rep = rep(1:nreps, each = nsteps4)) %>%
	mutate(day = (time-1) %/% spd + 1))) %>%
	unnest() %>%
	separate(files, into = c("treatment","ramp", "p"), sep = "_") %>%
	mutate(p = gsub("([A-Za-z=_]|.txt)", "", p)) %>%
	type_convert() %>%
	mutate(St = S + mut_S, At = A + mut_A, Bt = B + mut_B, Dt = D + mut_D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (q = 0)", 
		`0.05` = "low (q = 0.05)", 
		`0.1` = "intermediate (q = 0.1)", 
		`0.3` = "high (q = 0.3)", 
		`0.5` = "extrahigh (q = 0.5)")) %>%
	filter(p <= 0.3, day <= 24) %>%
	mutate(treatment = recode_factor(treatment,  both = "combination")) %>%
	mutate(MICday = 20) %>%
	select(-ramp)


	
	
states = list(St = "no resistance", 
		At = "rifampicin resistance", 
		Bt = "nalidixic acid resistance", 
		Mt = "mixed resistance", 
		Dt = "double resistance")
results.all = bind_rows(results,results_double, results_quadruple)

muller = results.all %>%
	group_by(day) %>%
	filter(day == MICday, time == max(time)) %>%
	ungroup() %>%
	select(treatment, MICday, day, p, pmutS.text, rep, time, St, At, Bt, Dt) %>%
	mutate(total = St + At + Bt + Dt) %>%
	gather(state, value, -c(treatment, MICday, p, pmutS.text, time, day, rep, total)) %>%
	mutate(detected = map2_lgl(total, value, ~ifelse(.y>0&.x>0, rbinom(n = 1, size = .x, prob = .y/(.x*200))>0,F))) %>%
	select(-value, -total) %>%
	spread(state, detected) %>%
	mutate(Mt = At&Bt, St = ifelse(At|Bt|Mt|Dt, F, 1), At = ifelse(Mt|Dt, F, At), Bt = ifelse(Mt|Dt, F, Bt), Mt = ifelse(Dt, F, Mt)) %>%
	gather(state, value, -c(treatment, MICday, day, p, pmutS.text, rep, time)) %>%
	filter(value == T) %>%
	select(-value) %>%
	mutate(state = factor(unlist(states[state]), levels = states)) %>%
	arrange(treatment, MICday, day, time, rep) %>%
	droplevels()

muller.sum = muller %>%
	count(treatment, MICday, day, pmutS.text, state, .drop=F)



# Plots
bestcolours = c("grey90", "#469067", "#235161", "#6f2c3d", "#ef0358")

dose.plot = muller.sum %>%
	group_by(treatment, MICday, day) %>%
	filter(state == "double resistance")  %>%
	ggplot(aes(x=MICday, y=n, fill = pmutS.text, shape = pmutS.text)) +
	geom_vline(xintercept = 5, linetype = "dotted") +
#	geom_line(aes(colour = pmutS.text)) +
	geom_point(size = 2) +
#	facet_grid(~state) +
	labs(x = "Transfer # when MIC is achieved", y = "Number of double resistant populations\nwhen MIC is achieved") +
	scale_shape_manual(values = c(21,22,23,24), name = "Initial proportion of mutators") +
	scale_fill_brewer(palette="RdBu", name = "Initial proportion of mutators", direction = -1) +
	scale_colour_brewer(palette="RdBu", name = "Initial proportion of mutators", direction = -1)

MICday.sim.plot = muller %>%
	filter(day == MICday) %>%
	group_by(day) %>%
	filter(time == max(time)) %>%
	ggplot(aes(x=as.factor(MICday), fill=state)) +
	geom_bar() +
	facet_grid(pmutS.text~.) +
	labs(x="Time MIC achieved (# transfers)",
		y = "Number of populations\nwhen MIC achieved") +
	scale_fill_manual(name = "Detection of", values = bestcolours)

ggsave("dose_plot.pdf", dose.plot, width = 5, height = 2.75)

