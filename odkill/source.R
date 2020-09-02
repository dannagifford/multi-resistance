require(tidyverse)
require(ggpubr)
require(grid)
require(scales)

theme_set(theme_bw())
theme_update(text = element_text(size = 10),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	strip.background = element_blank()
	legend.justification = "top",
	panel.background = element_blank()
)
MHcols = "grey60"

rifcols = rev(c("#1C3929", "#316448", "#469067", "#7DB194", "#B5D2C2", "#eaefec"))

nalcols = rev(c("#183843", "#235161", "#658590", "#A7B9BF", "#BDCACF","#FFFFFF"))

doubcols = rev(c("#5F0123", "#A7023E", "#EF0359", "#F34E8A", "#F89ABC", "#FFEFEF"))




odvnt.key = read_csv("odvnt-key.csv")
odvnt.count = read_csv("odvnt.csv")

odvnt.od = tibble(files = list.files(pattern = ".+OD-NT.+CSV"),
	data = map(files, ~read_csv(file = ., skip = 2, n_max = 1))) %>%
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
		experiment = "odvnt",
		concentration = concentration * 10)


kill.key = read_csv("kill-key.csv")
kill.count = read_csv("killcurves.csv")

kill.od = tibble(files = list.files(pattern = ".+Kill.+CSV"),
	data = map(files, ~read_csv(file = ., skip = 2, n_max = 1))) %>%
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
		experiment = "kill",
		concentration = concentration * 10)

odvnt_original = read_csv("../2018-06-08_ODvsNT.csv") %>%
	mutate(concentration = 0, antibiotic = "MH", experiment = "odvnt_original", strain = "sensitive", series = 0) %>%
	select(experiment, row=`Well Row`, col=`Well Col`, series, strain, concentration, antibiotic, OD.blanked=correctedOD, NT)


alldata = bind_rows(odvnt, kill) %>%
	mutate(NT = 4*count*10^(-plating)) %>%
	bind_rows(., odvnt_original) %>%
	mutate(abconc = paste(antibiotic, concentration)) %>%
	filter(!is.na(NT)) %>%
	 select(experiment, row, col, series, strain, abconc, antibiotic, concentration, OD.blanked, NT)


#abconc = sort(unique(alldata$abconc))
#abconc = c(abconc[7], abconc[14:19], abconc[8:13], abconc[1:6])

#alldata$abconc = factor(alldata$abconc, levels = abconc)

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
	
abconc.levels = alldata %>%
	select(abconc, antibiotic, concentration) %>%
	arrange(antibiotic, concentration) %>%
	unique() %>%
	.$abconc

alldata = alldata %>%
	mutate(abconc = factor(abconc, levels = abconc.levels))


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
	theme(plot.margin = unit(c(0,0.5,0,0), "lines"))+
	geom_line(aes(x = OD.blanked, y = fitted_vals, colour = abconc), inherit.aes=F) +
	geom_point(size = 2) +
	scale_colour_manual(values = c(MHcols, rifcols, nalcols, doubcols), name = "Antibiotic concentration (mg/l)")


joinedODp = ggarrange(odp_rev, killp, ncol = 2, legend.grob = get_legend(odp), legend = "right", labels = "AUTO", font.label = list(face = "plain"), align = "v")
cairo_pdf("ODvsNT.pdf", width = 6, height = 6)
joinedODp
dev.off()

