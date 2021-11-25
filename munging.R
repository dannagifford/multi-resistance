#install.packages("tidyverse")
#install.packages("lubridate")
#install.packages("growthcurver")
#install.packages("minpack.lm")
#install.packages("reshape2")

library(tidyverse)
require(lubridate)
require(growthcurver)
require(minpack.lm)

### Custom functions #############################################
read_csv_drop = function(...) read_csv(...) %>%
	select(-ncol(.))

# Some functions for importing BMG-style CSV files in a tidy way

BMGtime = function(otime){
# A function to convert BMG's default time format to fractional hours
	missing.s = grep("[0-9] s", otime, value = F, invert = T)
	missing.m = grep("[0-9] min", otime, value = F, invert = T)
	missing.h = grep("[0-9] h", otime, value = F, invert = T)
	otime[missing.s] = sub("$", " 0 s", otime[missing.s])
	otime[missing.h] = sub("^", "0 h ", otime[missing.h])
	otime[missing.m] = sub("h", "h 0 min ", otime[missing.m])
	otime = period_to_seconds(hms(gsub("[a-z ]+",":", otime)))/3600
return(otime)
}

readBMG = function(conn){
		data = read_csv_drop(conn, skip = 2, col_names = c("Well Row", "Well Col", "Content")) %>%
			select(-ends_with("_1")) 
		n_rows = nrow(data)	

		time_vec = read_csv(conn, n_max = 1) %>% as.vector() %>% t() %>% BMGtime()
		time_vec = time_vec[4:(length(time_vec)-1)]
		time_vec = rep(time_vec, times = n_rows)
		
		data = data %>% pivot_longer(-c(`Well Row`, `Well Col`, `Content`),
				names_to = "time", values_to = "OD") %>%
			mutate(time = time_vec)
		
		return(data)
		}


##################################################################

original.labels = as.vector(mapply(paste0, LETTERS[2:11], MoreArgs = list(2:7), SIMPLIFY = T,USE.NAMES = F))
original.cols = paste(original.labels, collapse = ";")
standard.labels = as.vector(mapply(paste0, LETTERS[2:7], MoreArgs = list(2:11), SIMPLIFY = T,USE.NAMES = F))
standard.cols = paste(standard.labels, collapse = ";")

filesaggregate = list.files(pattern = "*aggregate.*.csv")

### Data from the selection experiments

#### The aggregate files include the measured proportions of mutS strains
states.aggregated = tibble(file = filesaggregate, contents = map(file, ~read_csv(.))) %>%
	select(-file) %>%
	unnest(contents) %>%
		gather(pmutS, count, -c(volume, day, concentration, antibiotic, resistance)) %>%
	type_convert() %>%
	mutate(resistance = recode_factor(resistance, MH = "surviving",
					Rif = "rifampicin resistance",
					Nal = "nalidixic acid resistance",
					`Nal + Rif` = "double resistance"), 
	antibiotic = recode_factor(antibiotic, none = "no antibiotic",rifampicin = "rifampicin",
				`nalidixic acid` = "nalidixic acid",
				`rifampicin+nalidixic acid` = "combination")) %>%
	filter(!is.na(count)) %>%
	group_by(volume, antibiotic, pmutS) %>%
	arrange(volume, antibiotic) %>%
	nest() %>%
	group_by(volume, antibiotic) %>%
	mutate(pmutS.rank = as.integer(c(0,10,25,50)[rank(pmutS)])) %>%
		mutate(pmutS.text = recode_factor(pmutS.rank, `0` = "none",
		`10` = "low",
		`25` = "intermediate",
		`50` = "high")) %>%
	unnest()

#### The individual files include the fate of each well
filesindividual = list.files(pattern = "*individual.*.csv")

popns = tibble(filesindividual) %>%
	mutate(data = map(filesindividual, ~read_delim(., delim = ","))) %>%
	unnest() %>%
	mutate(antibiotic = ifelse(grepl("none",filesindividual),"none", ifelse(grepl("double", filesindividual), "combination", antibiotic))) %>%
	select(-filesindividual, -expt) %>%
	mutate(antibiotic = recode_factor(antibiotic, none = "no antibiotic", rifampicin = "rifampicin", `nalidixic acid` = "nalidixic acid", combination = "combination")) %>%
	mutate(medium = recode_factor(medium, MH = "sensitive", rif = "rifampicin resistance", nal = "nalidixic acid resistance", `nal+rif` = "double resistance")) %>%
	group_by(antibiotic, volume) %>%
	rename(pmutS.rank = pmutS) %>%
	type_convert() %>%
	left_join(., states.aggregated %>% select(., volume, antibiotic, pmutS, pmutS.rank, pmutS.text) %>% distinct()) %>%
	mutate(wells = ifelse(wells == "ALL", original.cols, wells)) %>%
	mutate(concentration = 2^(day-5)) %>%
	mutate(complement = grepl("-", wells)) %>%
	mutate(wells = gsub("0", NA, wells)) %>%
	mutate(wells = gsub("-", "", wells)) %>%
	mutate(wells = map(wells, ~unlist(strsplit(., split = ";")))) %>%
	mutate(wellsc = map(wells, ~original.labels[!(original.labels %in% .)])) %>%
	mutate(well = ifelse(complement, wellsc, wells)) %>%
	select(-wells, -wellsc, -complement) %>%
	unnest(well) %>%
	mutate(value = 1) %>%
	spread(well, value, fill = 0) %>%
	gather(wellID, value, -c(volume, antibiotic, day, concentration, pmutS, pmutS.text, pmutS.rank, medium)) %>%
	spread(medium, value, fill = 0) %>%
	filter(wellID != "<NA>") %>%
	mutate(wellID = standard.labels[match(wellID, original.labels)]) %>%
	unite(state, sensitive, `rifampicin resistance`, `nalidixic acid resistance`, `double resistance`, sep = "", remove = F) %>%
	mutate(state = recode(state, `0000` = "no resistance",
		`1000` = "no resistance",

		`1100` = "rifampicin resistance",
		`0100` = "rifampicin resistance (low density)",

		`1010` = "nalidixic acid resistance",
		`0010` = "nalidixic acid resistance (low density)",

		`1110` = "mixed resistance",
		`0110` = "mixed resistance (low density)",

		`1111` = "double resistance",
		`1101` = "double resistance (low density)",
		`1011` = "double resistance (low density)",
		`1001` = "double resistance (low density)",
		`0111` = "double resistance (low density)",
		`0101` = "double resistance (low density)",
		`0011` = "double resistance (low density)",
		`0001` = "double resistance (low density)")) %>%
	mutate(state.simple = gsub(" (low density)","",state, fixed = T)) %>%
	mutate(state = factor(state, levels = c("no resistance","rifampicin resistance", "rifampicin resistance (low density)", "nalidixic acid resistance", "nalidixic acid resistance (low density)", "mixed resistance", "mixed resistance (low density)", "double resistance", "double resistance (low density)"))) %>%
	mutate(state.simple = factor(state.simple, levels = c("no resistance","rifampicin resistance", "nalidixic acid resistance", "mixed resistance", "double resistance"))) %>% 
	type_convert()

saveRDS(popns, file = "popns.Rds")

##################################################################

### Fluctuation strains growth curve data ########################
concentrations = c(0,2^(-4:1),2^(-4:1),2^(-4:1))
antibiotics = factor(c("no antibiotic",rep("combination", times = 6), rep("rifampicin", times = 6), rep("nalidixic acid", times = 6)), levels = c("no antibiotic","rifampicin","nalidixic acid", "combination"))


FTfiles = list.files(path = "gc_fluctuation/", pattern = ".csv", recursive = T, full.names = F)
FTdata = tibble(files = FTfiles) %>%
	mutate(data = map(files, ~readBMG(paste0("gc_fluctuation/",.))), files = as.numeric(gsub(".csv","",files))) %>%
	arrange(files)
FTdata$antibiotic = antibiotics[FTdata$files]
FTdata$concentration = concentrations[FTdata$files]

FTdata = FTdata %>%
		filter(antibiotic == "no antibiotic") %>%
		bind_rows(.,.,.,.) %>%
		mutate(antibiotic = factor(levels(FTdata$antibiotic), levels = levels(FTdata$antibiotic))) %>%
		bind_rows(., FTdata) %>%
		filter(antibiotic != "no antibiotic") %>%
		arrange(antibiotic) %>%
		unnest()

FTsamples = FTdata %>%
	filter(!`Well Row` %in% c("A","H"), !`Well Col` %in% c(1,12), !(`Well Row`%in%c("E","F","G")&`Well Col`%in%c(6:11))) %>%
	mutate(id = as.integer(`Well Col`-1-5*((`Well Col`-1)%/%6)), rep = as.integer((`Well Col`-1)%/%6+1), strain = recode_factor(`Well Row`, E = "S",B = "A", C = "B",D = "D",F = "Sp1",G = "Sp2")) %>%
	mutate(rep = ifelse(grepl("S", strain), id, rep), id = ifelse(grepl("S", strain), 1, id)) %>%
	mutate(strain2 = recode_factor(strain, S = "sensitive", Sp1 = "mutator (1)", Sp2 = "mutator (2)",A = "rifampicin resistant", B = "nalidixic acid resistant", D = "double resistant")) %>%
	mutate(rep = paste0("rep", rep)) %>%
	type_convert()


blanks = FTdata %>%
	filter(grepl("Blank",Content)) %>% 
	type_convert() %>%
	group_by(antibiotic, concentration, time) %>%
	summarise(blank = mean(OD), blank.sd = sd(OD))

	
#### Remove outliers
outliers = FTsamples %>%
	group_by(antibiotic, concentration, strain, id, time) %>%
	mutate(OD_mean = mean(OD), OD_sd = sd(OD)) %>%
	mutate(dOD = OD-OD_mean) %>%
	ungroup() %>%
	filter(dOD>0.1) %>%
	select(antibiotic, concentration, id, rep, strain) %>%
	distinct() %>%
	unite(strains, antibiotic, concentration, id, rep, strain) %>%
	.$strains

FTsamples.outliers = FTsamples %>%
	unite(strains, antibiotic, concentration, id, rep, strain, remove = F) %>%
	mutate(outlier = strains%in%outliers) %>%
#	filter(strains%in%outliers) %>%
	select(-strains)

FTsamples = FTsamples %>%
	unite(strains, antibiotic, concentration, id, rep, strain, remove = F) %>%
	filter(!strains%in%outliers) %>%
	select(-strains)

### Blank correct
FTsamples.blanked = left_join(FTsamples, blanks) %>%
	mutate(blankOD = OD-blank)

#### Recalculate means
FTsamples.mean = FTsamples.blanked %>%
	group_by(antibiotic, concentration, strain, strain2, id, time) %>%
	summarise(OD_mean = mean(OD), OD_sd = sd(OD), blankOD_mean = mean(OD), blankOD_sd = sd(OD))

saveRDS(FTsamples.mean, file = "FTsamples_mean.Rds")

### OD vs NT
ODvsNT = read_csv("2018-06-08_ODvsNT.csv")
NTmodel = ODvsNT %>%
	mutate(blankOD = correctedOD) %>%
	filter(`Well Row` == "A" ) %>%
	lm(NT~blankOD, data = .)
saveRDS(ODvsNT, file = "ODvsNT.Rds")

#### Fit growth curves
FTcurves = FTsamples.blanked %>%
	filter(time <= 22) %>%
		group_by(antibiotic,`Well Row`, `Well Col`, concentration, id, strain, strain2, rep) %>%
		nest() %>%
		mutate(gc_fit = map(data, ~SummarizeGrowth(.$time, .$OD))) %>%
		mutate(vals = map(gc_fit, ~unlist(.$vals) %>%
			t() %>%
			as_tibble() %>%
			type_convert())) %>%
		unnest(vals) %>%
	arrange(id) %>%
	mutate(id = as.character(id))


FTcurves.mean = FTcurves %>%
	group_by(antibiotic, concentration, strain, strain2, id) %>%
	summarise_at(.vars = vars(k:auc_e), .funs = c("mean","sd"))

saveRDS(FTcurves.mean, file = "FTcurves_mean.Rds")

blankscurves = blanks %>%
#		select(-files) %>%
		filter(time <= 22) %>%
		group_by(antibiotic, concentration) %>%
#		group_by(antibiotic, concentration, `Well Row`, `Well Col`, Content) %>%
		nest() %>%
		mutate(gc_fit = map(data, ~SummarizeGrowth(.$time, .$blank))) %>%
		mutate(vals = map(gc_fit, ~unlist(.$vals) %>%
			t() %>%
			as_tibble() %>%
			type_convert())) %>%
		unnest(vals)

##################################################################


### Evolved strains growth curve data ###########################

#### These lines define the plate layout
pm = states.aggregated %>%
	filter(concentration == 1, antibiotic == "combination",	resistance == "surviving", pmutS>0) %>%
		select(volume, pmutS, pmutS.text) %>%
	group_by(volume) %>%
	mutate(rank = rank(pmutS))

rows = tibble(`Well Row` = rep(LETTERS[1:8],times = 2),
	volume = rep(c(1,20), each = 8),
	rank = c(NA,3,3,3,2,2,1,NA,NA,3,3,2,2,1,NA,NA)) %>%
	left_join(.,pm)

col1 = reshape2::melt(list(B = 2:11,C = 2:11,D = 2:10,E = 2:11,F = 2:3,G = 2:10)) %>%
	mutate(volume = 1)

col20 = reshape2::melt(list(B = 2:11,C = 2:5,D = 2:11,E = 2:9,F = 2:7)) %>%
	mutate(volume = 20)

cols = rbind(col1, col20) %>%
	rename(`Well Col` = value, `Well Row` = L1) %>%
	as_tibble()

EV.layout = left_join(rows, cols)
####

EVfiles = list.files(pattern = "*.csv", path = "gc_evolved", recursive = T, full.names = T)

EVdata = tibble(files = EVfiles, parameters = gsub("gc_evolved/","",files)) %>%
	mutate(data = map(files, ~readBMG(.))) %>%
	separate(col = parameters, sep = "[/_.]", into = c("concentration","volume","rep")) %>%
	mutate(rep = paste0("plate", rep)) %>%
	select(-files) %>%
	unnest(data) %>%
	type_convert() %>%
	left_join(., EV.layout) %>%
	filter(!`Well Row` %in% c("A","H"), !`Well Col` %in% c(1,12), !is.na(antibiotic)) %>%
	mutate(concentration.text = recode_factor(concentration, `0` = "0 mg/L rifampicin+nalidixic acid", `2` = "20 mg/L rifampicin+nalidixic acid")) %>%
	rename(lineage = Content)

EVdata.mean = EVdata %>%
	group_by(concentration, volume, antibiotic, time, pmutS, pmutS.text, lineage) %>%
	summarise(OD = mean(OD))



#### Fit growth curves

EVcurves = EVdata %>%
	filter(time <= 22) %>%
	group_by(rep, concentration, volume, antibiotic, `Well Row`, `Well Col`, pmutS, pmutS.text, lineage) %>%
		nest() %>%
		mutate(gc_fit = map(data, ~SummarizeGrowth(.$time, .$OD))) %>%
		mutate(vals = map(gc_fit, ~unlist(.$vals) %>%
			t() %>%
			as_tibble() %>%
			type_convert())) %>%
		unnest(vals) %>%
	arrange(pmutS)

EVcurves.mean = EVcurves %>%
	group_by(concentration, volume, antibiotic, pmutS, pmutS.text, lineage) %>%
	filter(auc_e>0.1) %>% #things smaller didn't grow...
	summarise_at(.vars = vars(k:auc_e), .funs = c("mean","sd")) %>%
	mutate(id = lineage)
	
	

joinedcurves = FTcurves %>%
	filter(antibiotic == "combination",
		concentration %in% c(0,2),
		strain == "D") %>%
	bind_rows(., EVcurves) %>%
	replace_na(list(pmutS.text="none")) %>%
	mutate(pmutS.text = recode_factor(pmutS.text, none="none (from fluctuation test)")) %>%
	mutate(strain = "D", strain2 = "double resistant") %>%
	mutate(id_EV = paste(`Well Row`, `Well Col`)) %>%
	mutate(id = ifelse(is.na(id), id_EV, id)) %>%
	mutate(rep = ifelse(grepl("plate[789]", rep), 
		paste0("plate",as.numeric(gsub("[^0-9.-]", "", rep))-6 ),rep)) %>%
	select(-id_EV)
saveRDS(joinedcurves, file = "joinedcurves.Rds")
	

joinedcurves.mean = FTcurves.mean %>%
	filter(antibiotic == "combination",
		concentration %in% c(0,2),
		strain == "D") %>%
	bind_rows(., EVcurves.mean) %>%
	replace_na(list(pmutS.text="none")) %>%
	mutate(pmutS.text = recode_factor(pmutS.text, none="none (from fluctuation test)")) %>%
	mutate(strain = "D", strain2 = "double resistant")

saveRDS(joinedcurves.mean, file = "joinedcurves_mean.Rds")


#### Daily OD measurements

dailyODfiles = list.files(path = "./daily_ods/", pattern = "(none|rifampicin|acid|double).csv")
dailyOD = tibble(files = dailyODfiles, data = map(files, ~read_csv_drop(paste0("./daily_ods/",.)))) %>%
	mutate(files = gsub(".csv","",files)) %>%
	separate(files, sep = "_", into = c("date", "day", "pmutS", "volume", "antibiotic")) %>%
	unnest() %>% 
	filter(!`Well Row`%in%c("A","H"), !`Well Col`%in%c(1,12)) %>%
	rename(OD = `Raw Data (600)`) %>%
	type_convert() %>%
	mutate(wellID = paste0(`Well Row`, `Well Col`),
		pmutS = 100-pmutS) %>%
	mutate(antibiotic = recode_factor(antibiotic,
		none = "no antibiotic",
		rifampicin = "rifampicin",
		`nalidixic acid` = "nalidixic acid",
		double = "combination")) %>%
	mutate(expt = ifelse(antibiotic == "combination", "combination", "single")) %>%
	mutate(pmutS.text = recode_factor(pmutS, `0` = "none",
		`10` = "low",
		`25` = "intermediate",
		`50` = "high")) %>%
	select(-c(pmutS, `Well Row`, `Well Col`,Content, expt)) %>%
	left_join(popns %>% ungroup() %>% select(volume, day, antibiotic, pmutS, pmutS.text, wellID, state, state.simple), .) %>%
	mutate(N = OD/1.38*10^9) %>%
	mutate(concentration = concentrations[day]) %>%
	mutate(state2 = as.factor(recode(state.simple,
		sensitive = "S",
		`rifampicin resistant` = "R",
		`nalidixic acid resistant` = "N",
		`mixed resistant` = "M",
		`double resistant` = "D",
		`no growth` = "E"))) %>%
	mutate(state2 = recode_factor(state2, S = "S",
			R = "A",
			N = "B",
			M = "A+B",
			D = "D",
			E = "E"))
saveRDS(dailyOD, file = "dailyOD.Rds")

