require(tidyverse)
require(ggpubr)
require(tidyverse)
require(grid)
require(scales)


nreps = 100
dt = 1
nsteps = 22*6*(1/dt) + 1 # 133
spd = 1/dt*22 #steps per day

theme_set(theme_bw())
theme_update(text = element_text(size = 10), 

	panel.grid.major = element_blank(), 

	panel.grid.minor = element_blank()
)


labelzero = function(x){format(x, zero.print = "0")}

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
R
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
	data_frame(y = median(x), ymin = y - mad(x), ymax = y + mad(x))}

# Munging
results = data_frame(files = list.files(pattern = ".txt", recursive = T)) %>%
	mutate(data = map(files, ~read_tsv(., col_names = c("S", "A", "B", "D", "mut_S", "mut_A", "mut_B", "mut_D"), col_types = "dddddddd") %>%
		mutate(time = rep(1:nsteps, times = nreps), rep = rep(1:nreps, each = nsteps)) %>%
		mutate(day = (time-1) %/% spd + 1) %>%
		mutate(concentration = 10*(1/16)*2^(day-1)))) %>%
	unnest() %>%
	separate(files, into = c("treatment","treatment2","N","p"), sep="[/_]") %>%
	select(-treatment2) %>%
	mutate(p = gsub("(p=|.txt)", "", p)) %>%
	type_convert() %>%
	select(-starts_with("mut_")) %>% # there are no mutators!
	mutate(St = S , At = A, Bt = B, Dt = D) %>%
	mutate(pmutS.text = recode_factor(p, 
		`0` = "none (p = 0)")) %>%
	mutate(treatment = recode_factor(treatment, none = "none", rif = "rifampicin", nal = "nalidixic acid", both = "combination"))

results.gathered = results %>%
	select(-c(St, At, Bt, Dt)) %>%
	gather(key = "strain", value = "n", -c(treatment, rep, N, p, pmutS.text, day, concentration, time)) %>%
	mutate(mutator = as.factor(ifelse(!grepl("mut_", strain), "wild-type", "mutator"))) %>%
	mutate(mutator = relevel(mutator, "wild-type")) %>%
	mutate(strain = recode_factor(strain, S = "sensitive", A = "rifampicin resistant", B = "nalidixic acid resistant", D = "double resistant", mut_S = "sensitive", mut_A = "rifampicin resistant", mut_B = "nalidixic acid resistant", mut_D = "double resistant"))


states = list(St = "no resistance", 
		At = "rifampicin resistance", 
		Bt = "nalidixic acid resistance", 
		Mt = "mixed resistance", 
		Dt = "double resistance")

proportions = results %>%
	select(treatment, day, concentration, N, pmutS.text, rep, time, St, At, Bt, Dt) %>%
	mutate(total = St + At + Bt + Dt) %>%
	gather(state, value, -c(treatment, N, pmutS.text, time, day, rep, concentration, total)) %>%
	mutate(proportion = value / total)

proportions.mean = proportions %>%
	group_by(day, concentration, pmutS.text, time, N, treatment, state) %>%
	summarise(proportion.mean = mean(proportion), proportion.sd = sd(proportion, na.rm=T)) %>%
	mutate(strain = recode_factor(state, St = "sensitive", At = "rifampicin resistant", Bt = "nalidixic acid resistant", Dt = "double resistant"))

# Plots
bestcolours = c("grey90", "#469067", "#235161", "#6f2c3d", "#ef0358")
bestcolours2 = c("black", bestcolours[-1])

Nproportions.plot = proportions.mean %>%
	ungroup() %>%
	filter(time == max(time)) %>%
	ggplot(aes(x=5.71*10^N, y=proportion.mean, fill = strain, shape = strain)) +
#		geom_bar(stat="identity", position="dodge") +
		geom_errorbar(aes(ymin = proportion.mean - proportion.sd,
			ymax = proportion.mean + proportion.sd,
			colour = strain),
#			position = position_dodge(0.9),
			width = 0) +
		geom_line(aes(colour = strain)) +
		geom_point(size = 2) +
		scale_shape_manual(values = c(21, 22, 23, 24), name = "Strain") +
		facet_grid(~treatment) +
		scale_colour_manual(values = bestcolours2[c(1:3,5)], name = "Strain") +
		scale_fill_manual(values = bestcolours[c(1:3,5)], name = "Strain") +
		labs(x="Simulated maximum population size", y="Mean proportion\nwithin populations") +
		scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
		scale_y_continuous(labels = labelzero, limits = c(0,1))



ggsave("Nproportions.pdf", Nproportions.plot, height = 2, width = 8)

