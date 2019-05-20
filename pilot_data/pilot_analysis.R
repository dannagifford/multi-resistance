library(tidyverse)
library(grid)
library(LaplacesDemon)
library(Cairo)
require(ggpubr)

theme_set(theme_bw())
theme_update(text = element_text(size=10),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
)

# Pilot experiment with either all BW25113 or mutators to set priors

pilot_survival = read_csv("2017-12-01_population-survival-datda.csv") %>%
	mutate(antibiotic = factor(antibiotic,
		levels=c("none", "rifampicin", "nalidixic acid", "combination"))) %>%
	mutate(strain = recode_factor(strain, BW25113 = "wild-type", `ΔmutS`="mutator"))


pilot_plot = pilot_survival %>%
	group_by(strain, day, antibiotic) %>%
	summarise(palive = sum(OD>=0.15)/n()) %>%
	ggplot(aes(x = day-1, y = palive, linetype = strain)) +
		geom_line() +
		facet_grid(~antibiotic) +
		labs(x="Time (days)", y="Proportion alive", tag="A") +
	scale_x_continuous(breaks=unique(pilot_survival$day-1),
			labels = unique(pilot_survival$day-1),
			sec.axis = sec_axis(~ ., breaks = unique(pilot_survival$day-1),
			labels = unique(pilot_survival$concentration)*10, name="Antibiotic concentration (mg/l)")) +
		theme(axis.text.x.top= element_text(angle = 45, hjust = 0)) +
		scale_linetype(name="Strain (pure culture)")

pilot_plot = ggplotGrob(pilot_plot)
pilot_plot[["grobs"]][[18]][["children"]][[2]] = nullGrob()

# Model S1 Prior illustration
prior_plotS1 = ggplot(data.frame(x=c(-20, 20)), aes(x)) +	
	stat_function(fun=dst, args=list(nu=7,mu=0,sigma=2.5), aes(linetype="Main effects\nand interactions")) +
	stat_function(fun=dst, args=list(nu=7,mu=-5,sigma=2.5), aes(linetype="Intercepts")) +
	geom_vline(xintercept = c(-10,logit(1/54), logit(39/54), 10), linetype="33", color="grey60")  +
	geom_text(data = tibble(x = c(-10,logit(1/54),logit(39/54), 10), 
		text = paste("p =",c("0.00005", "1/54", "39/54", 0.99995))), mapping = aes(x, y=0.14, label = text), hjust="left", angle = -90, nudge_x=-0.75, size=3) + 
	labs(x="Parameter", y="Density", tag="B") +
scale_linetype_manual("Priors", values = c("11", "solid"))

prior_explainS1 = ggarrange(pilot_plot, prior_plotS1, nrow = 2, heights=c(1.3, 1))

ggsave("prior_justificationS1.pdf", prior_explainS1, width=6, height=4, device=cairo_pdf)


# Model S2 and S3 Prior illustration
prior_plotS2 = ggplot(data.frame(x=c(-10, 30)), aes(x)) +	
	stat_function(fun=dst, args=list(nu=7,mu=0,sigma=2.5), aes(linetype="Main effects")) +
	stat_function(fun=dst, args=list(nu=7,mu=10,sigma=2.5), aes(linetype="Intercepts")) +
	geom_vline(xintercept = c(0,24), linetype="33", color="grey60")  +
	geom_text(data = tibble(x = c(0,24), 
		text = paste("AUC =",c("0", "24"))), mapping = aes(x, y=0.11, label = text), hjust="left", angle = -90, nudge_x=-0.75, size=3) + 
	labs(x="Parameter", y="Density") +
scale_linetype_manual("Priors", values = c("11", "solid"))

ggsave("prior_justificationS2.pdf", prior_plotS2, width=6, height=2, device=cairo_pdf)


