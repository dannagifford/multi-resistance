### PLOTTING
#install.packages("tidyverse")
#install.packages("ggpubr")
#install.packages("grid")
#install.packages("egg")
#install.packages("scales")
#install.packages("rflan")

require(tidyverse)
require(ggpubr)
require(grid)
#require(egg)
require(scales)
require(flan)


#source("munging.R") # This code will recreate all Rds files, which is slow, so we load them via readRDS instead after they are created.

tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {

  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}

theme_set(theme_bw())
theme_update(text = element_text(size = 10),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.justification = "top",
#	strip.background = element_blank()
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

  longestExponent <- max(sapply(gregexpr("\\d{1,}$", x), attr, 'match.length'))
  zeroTrimmed <- ifelse(longestExponent > 2,
                        paste0("\\1", paste(rep("~", times = longestExponent-1), collapse = "")),
                        "\\1")
  x <- gsub("(e[+|-])[0]", zeroTrimmed, x)
  x <- gsub("e", "~x~10^", x)
  if(any(grepl("\\^\\-", x))){
    x <- gsub("\\^\\+", "\\^~~", x)
  } else {
    x <- gsub("\\^\\+", "\\^", x)
  }

  x[x_zeros] <- "0"
  # return this as an expression
  parse(text=x)
} 

bestcolours = c("grey90","#469067", "#235161","#6f2c3d", "#ef0358")
bestcolours2 = c("grey40", bestcolours[-1])

popns = readRDS("popns.Rds")
muller = popns %>%
	filter(volume == 1) %>%
	ggplot(aes(x = day, fill = state.simple)) +
		geom_bar(width = 1, stat = "count", color = "black", size = 0.2) +
		facet_grid(pmutS.text~antibiotic) +
		scale_fill_manual(values = bestcolours, name = "Detection of") +
		labs(x = "Time (days)", y = "Number of populations") +
		scale_x_continuous("Time (days)", breaks = unique(popns$day), labels = unique(popns$day), sec.axis = sec_axis(~ ., breaks = unique(popns$day), labels = unique(popns$concentration)*10, name = "Antibiotic concentration (mg/l)")) +
		scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initial mutator frequency")) +
		theme(legend.key = element_rect(color = "white")) +
		theme(axis.text.x.top = element_text(angle = 45, hjust = 0))

mullert = ggplotGrob(muller)
mullert[["grobs"]][[18]][["children"]][[2]] = nullGrob()
muller = as_ggplot(mullert)

#day6 = popns %>%
#	group_by(volume, day, antibiotic, pmutS.text, pmutS, concentration, state.simple) %>%
#	summarise(count = n()) %>%
#	filter(day == 6, volume == 1) %>%
#	complete(state.simple, fill = list(count = 0))

#day6Tot = day6 %>%
#	mutate(state.sum = gsub("double","any", state.simple)) %>%
#	mutate(state.sum = map2_chr(antibiotic, state.sum, ~gsub(.x, "any", .y))) %>%
#	mutate(state.sum = map2_chr(antibiotic, state.sum, ~ifelse(.x != "combination",
#		gsub("mixed", "any", .y), .y))) %>%
#	filter(state.sum == "any resistance") %>%
#	group_by(volume, day, antibiotic, pmutS.text, pmutS, concentration, state.sum) %>%
#	summarise(count = sum(count)) %>%
#	filter(day == 6, volume == 1, antibiotic != "no antibiotic") %>%
#	complete(state.sum, fill = list(count = 0))

pointcolours = bestcolours
pointcolours[1] = "grey60"

#day6p = day6 %>% ggplot(aes(x = pmutS, y = count)) +
#	geom_line(aes(group = state.simple, color = state.simple)) +
#	geom_point(aes(fill = state.simple, shape = state.simple), size = 2) +
#	facet_grid(~antibiotic, scale = "free_x") +
#	geom_hline(yintercept = 60, linetype = "dashed", color = "grey40") +
#	labs(x = 'Initial mutator frequency', y = "Number of populations") +
#	scale_color_manual(values = muller.cols.day6, guide = F) +
#	scale_fill_manual(values = muller.cols.day6, name = "Detection of") +
#	scale_shape_manual(values = c(21,22,23,25,24), name = "Detection of") +
#	theme(plot.margin = margin(5,14,7,15))

ggsave("expevol.pdf", muller, device = "pdf", height = 6, width = 8)

### Growth curve plots
#### Fluctuation test generated strains
# Now shown in the supplemental material:
FTsamples.mean = readRDS("FTsamples_mean.Rds")
FT.gc.plot = FTsamples.mean %>%
	ungroup() %>%
	filter(!grepl("mutator", strain2), time > 1, time <= 30) %>%
	mutate(concentration = concentration *10) %>%
	ggplot(aes(x = time, y = OD_mean, color = strain2, linetype = as.factor(id))) +
	geom_line(aes(group = paste(id, strain2))) +
	facet_grid(antibiotic~concentration) +
	scale_color_manual(values = pointcolours[-4]) +
	geom_vline(xintercept = 22, linetype = "11", color = "grey40") +
	guides(color = guide_legend(title = "Strain"), linetype = F) +
	scale_x_continuous(name = "Time (h)", sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Antibiotic concentration (mg/l)")) +
	scale_y_continuous(name = "Blanked OD") +
	theme(legend.position = "bottom")

ggsave("ftgc.pdf", FT.gc.plot, device = "pdf", height = 4, width = 8)

FTcurves.mean = readRDS("FTcurves_mean.Rds")
FT.auc.plot = FTcurves.mean %>%
	filter(!grepl("mutator", strain2)) %>%
	ggplot(aes(x = as.factor(concentration*10), y = auc_e_mean, fill = strain2, shape = strain2)) +
	geom_line(aes(group = paste(strain2,id), color = strain2)) +
	geom_point(size = 2) +
#	geom_errorbar(aes(ymin = auc_e_mean-auc_e_sd, ymax = auc_e_mean+auc_e_sd), width = 0) +
	facet_grid(~antibiotic) +
	scale_fill_manual(values = pointcolours[-4], name = "Strain") +
	scale_color_manual(values = pointcolours[-4], guide = F) +
	scale_shape_manual(values = c(21,22,23,24), name = "Strain") +
	scale_x_discrete(labels = unique(FTcurves.mean$concentration*10), breaks = unique(FTcurves.mean$concentration*10), name = "Antibiotic concentration (mg/l)") +
	scale_y_continuous(name = "Growth (AUC of OD600)", limits = c(0,12)) +
	geom_vline(aes(xintercept = ifelse(antibiotic == "combination",5,6)), linetype = "22", color = "grey60") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ftauc.pdf", FT.auc.plot, device = "pdf", height = 2.5, width = 7)

##### Selection strains
#Not shown in the text, but plots raw growth curves (i.e. OD vs time):
#EV.gc.plot = EVdata.mean %>%
#	filter(volume == 1) %>%
#	ggplot(aes(x = time, y = OD, group = paste(pmutS, concentration), color = lineage)) +
#	geom_line(aes(group = paste(pmutS.text, lineage))) +
#	facet_grid(concentration~pmutS.text) +
#	geom_vline(xintercept = 22, linetype = "11", color = "grey40") +
#	guides(color = F, linetype = F) +
#	scale_y_continuous( sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Concentration of combination treatment \n(mg/l of each antibiotic)"))

#Not shown in text, plots mutator proportion against AUC, but we treat mutators as a categorical in models so not strictly appropriate...
#EV.auc.plot = EVcurves %>%
#		filter(volume == 1) %>%
#		ggplot(aes(x = pmutS, y = auc_e, fill = as.factor(concentration*10), shape = as.factor(concentration*10), linetype = as.factor(concentration*10))) +
#			geom_smooth(method = "lm", alpha = 0, aes(, color = as.factor(concentration*10))) +
#			geom_point(position = position_dodge(width = 0.2)) +
#		scale_x_continuous(name = "Initial mutator frequency", trans = "log2")+
#		scale_y_continuous(name = "Fitness (AUC of OD600)", limits = c(0,16)) +
#		guides(color = guide_legend(title = "Combination treatment concentration
#(mg/l of each antibiotic)"), shape = guide_legend(title = "Combination treatment concentration
#(mg/l of each antibiotic)"), fill = guide_legend(title = "Combination treatment concentration
#(mg/l of each antibiotic)"), linetype = guide_legend(title = "Combination treatment concentration
#(mg/l of each antibiotic)")) +
#	scale_shape_manual(values = c(21,24)) +
#	scale_color_manual(values = c("#ef0358","#ef0358")) +
#	scale_fill_manual(values = c("#ef0358","#ef0358")) +
#		theme(legend.key = element_rect(color = "white"))

# Instead we plot AUC in the presence and absence of antibiotics, coloured by mutator proportion.
joinedcurves.mean = readRDS("joinedcurves_mean.Rds")
joined.auc.plot = joinedcurves.mean %>%
	group_by(antibiotic, concentration, pmutS.text) %>%
	nest(auc_e_mean, auc_e_sd, .key = auc_e) %>%
	spread(concentration, auc_e) %>%
	unnest(`0`, `2`, .sep = "_") %>%
	ggplot(aes(`0_auc_e_mean`,`2_auc_e_mean`, fill = pmutS.text, shape = pmutS.text)) +
		geom_errorbarh(aes(xmin = `0_auc_e_mean`-`0_auc_e_sd`, xmax = `0_auc_e_mean`+`0_auc_e_sd`, color = pmutS.text)) +
		geom_errorbar(aes(ymin = `2_auc_e_mean`-`2_auc_e_sd`, ymax = `2_auc_e_mean`+`2_auc_e_sd`, color = pmutS.text)) +
		geom_point(size = 2) +
		geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60") +
		scale_y_continuous(name = "Growth in 20 mg/l", limits = c(5,10)) +
		scale_x_continuous(name = "Growth in 0 mg/l", limits = c(5,10)) +
		scale_colour_manual(values = c("grey60","#f6acc7","#ef0358","#82002f"))+
		scale_fill_manual(values = c("white","#f6acc7","#ef0358","#82002f"))+
		scale_shape_manual(values = c(1:4+20)) +
		guides(color = "none", fill = guide_legend(title = "Initial mutator frequency"), shape = guide_legend(title = "Initial mutator frequency"))

ggsave("joinedauc.pdf", joined.auc.plot, device = "pdf", height = 3, width = 5)


compare_data_sim = readRDS("compare_data_sim.Rds")
compare.plot = compare_data_sim %>%
	ggplot(aes(`Estimate from experiment`, `Estimate from simulation`, fill = Predictor, shape = `Predictor type`)) +
		facet_grid(~Response) +
		xlab("Estimate from experiment") + # bugfix
		ylab("Estimate from simulation") +
		geom_hline(yintercept = 0, color = "grey80", linetype = "11") +
		geom_vline(xintercept = 0, color = "grey80", linetype = "11") +
		geom_abline(intercept = 0, slope = 1, color = "grey60", linetype = "dashed") +
		geom_errorbarh(aes(xmin = `Q2.5 from experiment`, xmax = `Q97.5 from experiment`)) +
		geom_errorbar(aes(ymin = `Q2.5 from simulation`, ymax = `Q97.5 from simulation`)) +
		geom_point(size = 2.5) +
		scale_shape_manual(values = c(22:24)) +
		scale_fill_manual(values = c("darkgreen","grey20","grey60","white","red", "blue", "purple",
			"#990000", "#002b80", "#602040", "#ff3333", "#3399ff", "#ac39ac", "#ffb3b3", "#99ccff", "#df9fbf")) +
		guides(shape = FALSE, fill = guide_legend(override.aes = list(shape = c(22,23,23,23,24,24,24))))

ggsave("sim_vs_experiment.pdf", compare.plot, device = "pdf", height = 2, width = 8)
dailyOD = readRDS("dailyOD.Rds")
odp = dailyOD %>%
	filter(volume == 1) %>%
	ggplot(aes(day, OD, fill = state.simple, shape = state.simple)) +
	geom_jitter(alpha = 1, size = 2, width = 0.3) +
	facet_grid(pmutS.text~antibiotic) +
	scale_fill_manual(values = bestcolours2, name = "") +
	scale_shape_manual(values = c(21,22,23,25,24), name = "") +
	labs(y = "Optical density (600 nm)") +
	scale_x_continuous("Time (days)", breaks = unique(dailyOD$day), labels = unique(dailyOD$day), sec.axis = sec_axis(~ ., breaks = unique(dailyOD$day), labels = unique(dailyOD$concentration)*10, name = "Concentration (mg/l)")) +
	scale_y_continuous(limits = c(0.1,2), sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL, name = "Initiial mutator frequency")) +
	theme(axis.text.x.top = element_text(angle = 45, hjust = 0))


# Lose the concentration axis from the "no" antibiotic treatment
odpt = ggplotGrob(odp)
odpt[["grobs"]][[18]][["children"]][[2]] = nullGrob()
pdf("OD.pdf", height = 7 , width = 8)
grid.newpage()
grid.draw(odpt)
dev.off()

# Old NT plot.
###### OD versus N
#ODvsNT = readRDS("ODvsNT.Rds")
#odnp = ODvsNT %>%
#	ggplot(aes(x=correctedOD, y=NT)) +
#	geom_point() +
#	scale_x_continuous(name = " ") +
#	scale_y_continuous(name = " ", labels = scientific_10x) +
#	theme(plot.background = element_rect(fill = "transparent", colour = NA)) +
#	geom_vline(xintercept = 1, color = "grey40", linetype = "22") +
#	geom_hline(yintercept = 5.71e8, color = "grey40", linetype = "22")

#odnp_log = odnp +
#	scale_x_log10(name = "Optical density") +
#	scale_y_log10(name = "Number of bacteria", labels = trans_format("log10", math_format(10^.x))) +
#	annotate("text", x = 0.001, y = 1.4*5.71e8, label = scientific_10x(5.71e8, digits = 2), size = 3)

#odnp_inset <-
#  cowplot::ggdraw() +
#  cowplot::draw_plot(odnp_log) +
#  cowplot::draw_plot(odnp, x = 0.46, y = 0.1, width = 0.5, height = 0.4)

#ggsave("ODvsNT.pdf", odnp_inset, width = 3, height = 3)	


#### New plots to address reviewers' comments
MHcols = "grey60"

rifcols = rev(c("#1C3929", "#316448", "#469067", "#7DB194", "#B5D2C2", "#eaefec"))
nalcols = rev(c("#183843", "#235161", "#658590", "#A7B9BF", "#BDCACF","#FFFFFF"))
doubcols = rev(c("#5F0123", "#A7023E", "#EF0359", "#F34E8A", "#F89ABC", "#FFEFEF"))

odvnt.key = read_csv("odkill/odvnt-key.csv")
odvnt.count = read_csv("odkill/odvnt.csv")

odvnt.od = tibble(files = list.files(path="odkill/", pattern = ".+OD-NT.+CSV"),
	data = map(files, ~read_csv(file = paste0("odkill/",.), skip = 2, n_max = 1))) %>%
	separate(files, into = letters[1:3], sep = "_") %>%
	mutate(series = gsub(".CSV","",c)) %>%
	select(series, data) %>%
	type_convert() %>%
	arrange(series) %>%
	unnest() %>%
	gather(well, OD, `Sample X1`:`Sample X96`) 

odvnt = left_join(left_join(odvnt.key, odvnt.od), odvnt.count)

odvnt.blanks = odvnt %>%
	filter(is.blank) %>%
	group_by(series) %>%
	summarise(blank = mean(OD))

odvnt = odvnt %>%
	left_join(., odvnt.blanks) %>%
	mutate(OD.blanked = OD - blank,
		experiment = "odvnt")


kill.key = read_csv("odkill/kill-key.csv")
kill.count = read_csv("odkill/killcurves.csv")

kill.od = tibble(files = list.files(pattern = ".+Kill.+CSV", path = "odkill/"),
	data = map(files, ~read_csv(file = paste0("odkill/",.), skip = 2, n_max = 1))) %>%
	separate(files, into = letters[1:4], sep = "_") %>%
	mutate(series = gsub("[A-Za-z.]","",d)) %>%
	select(series, data) %>%
	type_convert() %>%
	unnest() %>%
	gather(well, OD, `Sample X1`:`Sample X96`) %>%
	arrange(series)

kill = left_join(left_join(kill.key, kill.od), kill.count)

kill.blanks = kill %>%
	filter(is.blank) %>%
	group_by(series) %>%
	summarise(blank = mean(OD))

kill = kill %>%
	left_join(., kill.blanks) %>%
	mutate(OD.blanked = OD - blank,
		experiment = "kill")


odvnt_original = read_csv("2018-06-08_ODvsNT.csv") %>%
	mutate(concentration=0, antibiotic = "MH", experiment = "odvnt_original", strain = "sensitive", series = 0) %>%
	select(experiment, row=`Well Row`, col=`Well Col`, series, strain, concentration, antibiotic, OD.blanked=correctedOD, NT)


alldata = bind_rows(odvnt, kill) %>%
	mutate(NT = 4*count*10^(-plating)) %>%
	bind_rows(., odvnt_original) %>%
	mutate(abconc = paste0(antibiotic, ", ", as.factor(10*concentration))) %>%
	filter(!is.na(NT)) %>%
	 select(experiment, row, col, series, strain, abconc, antibiotic, concentration, OD.blanked, NT)


abconc = sort(unique(alldata$abconc))
abconc = c(abconc[7],
	abconc[c(14,15,17,19,16, 18)],
	abconc[c(14,15,17,19,16, 18)-6],
	abconc[c(14,15,17,19,16, 18)-13])

alldata$abconc = factor(alldata$abconc, levels = abconc)

alldata = alldata %>%
	mutate(strain = recode_factor(strain,
	sensitive="sensitive",
	mutS = "mutator",
	rifampicin = "rifampicin resistant",
	`nalidixic acid`="nalidixic acid resistant",
	combination = "double resistant")) %>%
	mutate(antibiotic = recode_factor(antibiotic,
	MH="MH",
	rifampicin = "rifampicin",
	`nalidixic acid`="nalidixic acid",
	combination = "combination"))

original = alldata %>%
	filter(NT>0, experiment=="odvnt_original") %>%
	ggplot(aes(x = OD.blanked,
	y = NT,
	shape = strain,
	fill = abconc)) +
	geom_point(size = 2) +
	scale_shape_manual(values = c(21,25,22:24), name = "Strain") +
	scale_x_log10() +
	scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
	guides(fill = guide_legend(override.aes=list(shape=21))) +
	scale_fill_manual(values = c(MHcols, rifcols, nalcols, doubcols), name = "Antibiotic concentration") + 
	labs(x="Blanked OD", y=" (200μl)") + 
	facet_wrap(~antibiotic) + 
	geom_hline(yintercept = 5.71e8, linetype="22") +
	geom_vline(xintercept = 1, linetype="22")

odp = alldata %>%
	filter(NT>0, series!=0.01, experiment!="odvnt_original") %>%
	ggplot(aes(x = OD.blanked,
	y = NT,
	shape = strain,
	fill = abconc)) +
	geom_point(size = 2) +
	scale_shape_manual(values = c(21,25,22:24), name = "Strain") +
	scale_x_log10() +
	scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
	guides(fill = guide_legend(override.aes=list(shape=22, size=5))) +
	scale_fill_manual(values = c(MHcols, rifcols, nalcols, doubcols), name = "Antibiotic concentration (mg/l)") + 
	labs(x="Blanked OD", y="Number of bacteria (in 200μl)") + 
	facet_wrap(~antibiotic, nrow=4) + 
	theme(legend.key.size = unit(0.8, 'lines')) +
	theme(plot.margin = unit(c(0,0.5,0,0), "lines"))

killp = alldata %>%
	filter(concentration>-1, experiment == "kill") %>%
	ggplot(aes(x = as.factor(series),
	y = NT,
	fill = abconc)) +
	geom_jitter(size = 2, shape = 21, width=0.2) +
#	scale_shape_manual(values = c(21,25,22:24), name = "Strain") +
	guides(fill = guide_legend(override.aes=list(shape=23, size=6), ncol=3)) +
	scale_fill_manual(values = c(MHcols, rifcols[1:6], nalcols[1:6], doubcols[1:6]), name = "Antibiotic concentration (mg/l)") + 
	labs(x="Time (h)", y="Number of bacteria (in 200μl)") +
	scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
	facet_wrap(~antibiotic, nrow=4) +
	theme(plot.margin = unit(c(0,0.5,0,0), "lines"))

lmdata = alldata %>%
	filter(OD.blanked>0.05, NT>0) %>%
	mutate(antibiotic2 = gsub("combination", "double", antibiotic)) %>%
	mutate(strain2 = factor(ifelse(grepl("resistant",strain),"resistant","sensitive"),
		levels = c("sensitive","resistant"))) %>%
	mutate(grows_in = map2_lgl(strain, antibiotic2,
		~ grepl(.y, .x)  | .y=="MH")) %>%
	select(-antibiotic2) %>%
	filter((experiment == "kill" & series == 24) | experiment == "odvnt" | experiment == "odvnt_original", (grows_in | concentration <=10)) %>%
	droplevels()

m1 = lm(log10(NT) ~ log10(OD.blanked) * (antibiotic/concentration)*strain, data = lmdata)
m2 = step(m1)
summary(m2)

odp_rev = alldata %>%
	mutate(fitted_vals = 10^predict(m2, newdata = alldata)) %>%
	filter(NT > 0) %>%
	ggplot(aes(x = OD.blanked,
	y = NT,
	shape = strain,
	fill = abconc)) +
	scale_shape_manual(values = c(21,25,22:24), name = "Strain") +
	scale_x_log10(limits = c(0.01,1)) +
	scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
	guides(fill = guide_legend(override.aes=list(shape=22, size=4), colour=NA)) +
#	scale_fill_manual(values = c(MHcols, rifcols, nalcols[-3], doubcols[-3]), name = "Antibiotic concentration (mg/l)") + 
	scale_fill_manual(values = c(MHcols, rifcols, nalcols, doubcols), name = "Antibiotic concentration (mg/l)") + 
	labs(x="Blanked OD", y="Number of bacteria (in 200μl)") + 
	facet_wrap(~antibiotic, nrow = 4) + 
	theme(legend.key.size = unit(0.8, 'lines')) +
	theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
	geom_point(size = 2) +
	geom_line(aes(x = OD.blanked, y = fitted_vals, colour = abconc), inherit.aes=F) +
	scale_colour_manual(values = c(MHcols, rifcols, nalcols, doubcols), name = "Antibiotic concentration (mg/l)")


joinedODp = ggarrange(odp_rev, killp, ncol = 2, legend.grob = get_legend(odp), legend = "right", labels = "AUTO", font.label = list(face = "plain"), align = "v")
cairo_pdf("ODvsNT.pdf", width = 6, height = 6)
joinedODp
dev.off()

# Reviewer suggestion on p^2
## Replicates the first 4 days of the experiment (i.e. wild-types could still grow)
## 10^6 times and counts populations where at least one mutant was observed.
## Mutation rate mu is product of the two mutation rates to simulate the probability
## of two independent resistance mutations occurring simultaneously.

set.seed(19)
flansim = tibble(
	rep = 1:1e6,
	pops.mutators = map_dbl(rep,
		~sum(rflan(n=4*60, mutprob = 3.17e-14, mfn = 5.71e8)$mc > 0)),
	pops.wt = map_dbl(rep,
		~sum(rflan(n=4*60, mutprob = 5.0e-18, mfn = 5.71e8)$mc > 0)),
	prop.mutators = pops.mutators / 60,
	prop.wt = pops.wt / 60)

mean(flansim$prop.mutators)
sd(flansim$prop.mutators)

mean(flansim$prop.wt)
sd(flansim$prop.wt)
