library(ggmuller)
library(tidyverse)

#bestcolours3 = c("#dddddddd", "#46906788", "#23516188", "#ef035888", "#008088", "#999999ff", "#469067ff", "#235161ff", "#ef0358ff", "#0080ff") # AD and BD different

bestcolours3 = c("#dddddddd", "#46906788", "#23516188", "#ef035888", "#ef035888", "#999999ff", "#469067ff", "#235161ff", "#ef0358ff", "#ef0358ff") #AD and BD same


population.data.original = results %>%
	select(treatment, p, rep, times = time, S:mut_BD) %>%
	pivot_longer(cols = S:mut_BD, names_to = "names", values_to = "abundances") %>%
	group_by(treatment, p, rep) %>%
	mutate(strains = recode_factor(names,
		S="sensitive", A="rifampicin resistant",B="nalidixic acid resistant",AD="double resistant (from A)",BD="double resistant (from B)",
		mut_S="mutator (sensitive)", mut_A="rifampicin resistant (mutator)",mut_B="nalidixic acid resistant (mutator)",mut_AD="double resistant (from A, mutator)",mut_BD="double resistant (from B, mutator)"))

population.data = population.data.original %>%
	nest()
	
attributes = data.frame(names = unique(population.data.original$names), parents = c(NA, "S", "S", "A", "B", "S", "mut_S", "mut_S", "mut_A", "mut_B"), colors = bestcolours3)

edges = attributes[2:10,2:1]
names(edges) = c("Parents", "Identity")	
edges = edges %>%
	mutate(Identity = recode_factor(Identity,
		S="sensitive", A="rifampicin resistant",B="nalidixic acid resistant",AD="double resistant (from A)",BD="double resistant (from B)",
		mut_S="mutator (sensitive)", mut_A="rifampicin resistant (mutator)",mut_B="nalidixic acid resistant (mutator)",mut_AD="double resistant (from A, mutator)",mut_BD="double resistant (from B, mutator)")) %>%
	mutate(Parents = recode_factor(Parents,
		S="sensitive", A="rifampicin resistant",B="nalidixic acid resistant",
		mut_S="mutator (sensitive)", mut_A="rifampicin resistant (mutator)",mut_B="nalidixic acid resistant (mutator)")) 




population.data.means = population.data %>%
	ungroup() %>%
	group_by(treatment, p, names, times) %>%
	summarise(abundances = mean(abundances))
	
	
has.doubles = population.data.original %>%
	filter(names %in% c("AD", "BD", "mut_AD", "mut_BD") & times == max(times) & abundances > 0) %>%
	select(treatment, p, rep, -times) %>%
	unique()
	
population.data.has.doubles = population.data %>%
	left_join(has.doubles, .)
	
get_mdf = function(edges, data){
	data = data %>% select(Generation = times, Identity = strains, Population = abundances)
	get_Muller_df(edges, data)  %>%
	mutate(Identity = recode_factor(Identity,
		S="sensitive", A="rifampicin resistant",B="nalidixic acid resistant",AD="double resistant (from A)",BD="double resistant (from B)",
		mut_S="mutator (sensitive)", mut_A="rifampicin resistant (mutator)",mut_B="nalidixic acid resistant (mutator)",mut_AD="double resistant (from A, mutator)",mut_BD="double resistant (from B, mutator)"))
	}


# Examples of outcomes from the p=0.1 combination treatment
j = c(46, 2, 23, 39)
k = 0.1
pdf(paste0("example_popn_","p=",k,".pdf"), width = 8, height = 2.5)
for(i in "combination"){
for(j in j){
data_temp = population.data %>% filter(treatment == i, p == k)

mdf = get_mdf(edges, data_temp$data[[j]])

mdp = Muller_pop_plot(mdf, add_legend=T, conceal_edges = F) +
theme(axis.line=element_blank(),
#     axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
#     axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
#     panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank()) +
      scale_fill_manual(values = bestcolours3, na.value = "#FFFFFF00") + 
      scale_x_continuous(name = "Time (# transfers)", breaks = ((1:6)-0.5)*spd+2, labels = 1:6, limits = c(1, nsteps-1), sec.axis = sec_axis(~ ., breaks = ((1:6)-0.5)*spd+2, labels = unique(results$concentration), name = "Simulated antibiotic concentration (mg/l)"), expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) + 
      geom_vline(xintercept = (1:6)*spd+2, colour = "grey60", size = 0.5, linetype = "dashed")


arise_times = data_temp$data[[j]] %>%
	group_by(names, strains) %>%
	mutate(lagged_ab = lag(abundances), lead_ab = lead(abundances)) %>%
	filter(abundances<=0&lead_ab>0, !names%in%c("S","mut_S"), grepl("mut", names)) %>%
	mutate(Identity = strains, Generation = times) %>%
	ungroup() %>%
	left_join(mdp$data) %>%
	select(Identity, Generation, Group_id) %>%
	unique()

extinct_times = data_temp$data[[j]] %>%
	group_by(names, strains) %>%
	mutate(lagged_ab = lag(abundances), lead_ab = lead(abundances)) %>%
	filter(abundances>0&lead_ab<=0, !names%in%c("S","mut_S"), grepl("mut", names)) %>%
	mutate(Identity = strains, Generation = times) %>%
	ungroup() %>%
	left_join(mdp$data) %>%
	select(Identity, Generation, Group_id) %>%
	unique()


#arise_times = data_temp$data[[j]] %>%
#	filter(!names%in%c("S","mut_S"), abundances >0) %>%
#	group_by(Identity = strains) %>%
#	summarise(Generation = min(times)) %>%
#	left_join(mdp$data) %>%
#	select(Identity, Generation, Group_id) %>%
#	unique()

if(nrow(arise_times)>0){		
mdp = mdp +
	geom_star(data = arise_times, 
		aes(x = Generation, y = max(mdp$data$Population), fill = Identity, colour = Identity), size = 2, starshape = 9) +
		scale_colour_manual(values = bestcolours3, na.value = "#FFFFFF00", drop = F)
}

if(nrow(extinct_times)>0){	
mdp = mdp +
	geom_point(data = extinct_times,
		aes(x = Generation, y = max(mdp$data$Population), colour = Identity), size = 4, shape = 4) +
		scale_colour_manual(values = bestcolours3, na.value = "#FFFFFF00", drop = F)
}
plot(mdp)
}
}
dev.off()

## All plots for grid
#for(i in levels(population.data$treatment)){
#for(k in unique(population.data$p)){

#pdf(paste0("for_grid_popn_",i,"_","p=",k,".pdf"))
#data_temp = population.data %>% filter(treatment == i, p == k)

#for(j in 1:100){
#mdf = get_mdf(edges, data_temp$data[[j]])
#mdp = Muller_pop_plot(mdf, add_legend=T, conceal_edges = F) +
#theme(axis.line=element_blank(),
#      axis.text.x=element_blank(),
#      axis.text.y=element_blank(),
#      axis.ticks.y=element_blank(),
#      axis.title.x=element_blank(),
#      axis.title.y=element_blank(),
#      legend.position="none",
#      panel.background=element_blank(),
##     panel.border=element_blank(),
#      panel.grid.major=element_blank(),
#      panel.grid.minor=element_blank(),
#      plot.background=element_blank()) +
#      scale_fill_manual(values = bestcolours3, na.value = "#FFFFFF00") + 
#      scale_colour_manual(values = bestcolours3) +
##     labs(title = paste0(data_temp$treatment[j], ", p = ", data_temp$p[j], ", replicate = ", data_temp$rep[j])) +
#      scale_x_continuous(name = "Time (# transfers)", breaks = ((1:6)-0.5)*spd+2, labels = 1:6, limits = c(1, nsteps-1), sec.axis = sec_axis(~ ., breaks = ((1:6)-0.5)*spd+2, labels = unique(results$concentration), name = "Simulated antibiotic concentration (mg/l)"), expand = c(0,0)) +
#      scale_y_continuous(expand = c(0,0)) + 
#      geom_vline(xintercept = (1:6)*spd+2, colour = "grey60", size = 0.5, linetype = "dashed")

#arise_times = data_temp$data[[j]] %>%
#	filter(!names%in%c("S","mut_S"), abundances >0) %>%
#	group_by(Identity = strains) %>%
#	summarise(Generation = min(times)) %>%
#	left_join(mdp$data)

#mdp = mdp +
#	geom_star(data = arise_times, 
#	aes(x = Generation, y = max(mdp$data$Population), fill = Identity, colour = Identity, starshape = grepl("mutator", Identity)), size = 2) +
#	scale_starshape_manual(values = c(9,9)) +
#	scale_colour_manual(values = bestcolours3, na.value = "#FFFFFF00", drop = F)
#plot(mdp)
#}
#dev.off()

#}
#}

