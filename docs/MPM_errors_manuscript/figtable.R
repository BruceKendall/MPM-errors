# Test script for making figures/tables from case studies
setwd("~/Documents/Github/MPM-errors/docs/MPM_errors_manuscript")
# Read in data
library(dplyr)
library(magrittr)
lf <- read.csv("lionfish9models.csv")
lf <- lf %>% mutate(MatMod = rep(c("FAS", "SAS", "AAS"), 3),
                    FertSurv = c(rep(TRUE, 6), rep(FALSE, 3)),
                    JuvRep = c(rep(TRUE, 3), rep(FALSE, 6)),
                    Population = "Lionfish")
lf$MatMod <- factor(lf$MatMod, levels=c("AAS", "SAS", "FAS"))

aN <- read.csv("alligatorNorth.csv")
aN <- aN %>% mutate(MatMod = rep(c("FAS", "SAS", "AAS"), 3),
                    FertSurv = c(rep(TRUE, 6), rep(FALSE, 3)),
                    JuvRep = c(rep(TRUE, 3), rep(FALSE, 6)),
                    Population = "Alligator North")
aN$MatMod <- factor(aN$MatMod, levels=c("AAS", "SAS", "FAS"))

aS <- read.csv("alligatorSouth.csv")
aS <- aS %>% mutate(MatMod = rep(c("FAS", "SAS", "AAS"), 3),
                    FertSurv = c(rep(TRUE, 6), rep(FALSE, 3)),
                    JuvRep = c(rep(TRUE, 3), rep(FALSE, 6)),
                    Population = "Alligator South")
aS$MatMod <- factor(aS$MatMod, levels=c("AAS", "SAS", "FAS"))
alldata <- rbind(lf[,-(2:17)], aN[,-(2:21)], aS[,-(2:21)])
alldata$Population <- factor(alldata$Population, 
                             levels=c("Lionfish", "Alligator South",
                                      "Alligator North"))
# Figure of lionfish lambdas
library(ggplot2)
ggplot(alldata, aes(x = MatMod, y = Lambda)) +
  geom_point(aes(shape = JuvRep, fill = FertSurv), size = 4) + 
  scale_shape_manual(values = 21:22) +
  scale_fill_manual(values = c("black", "white")) +
  theme_classic() + theme(legend.position = "none") +
  ylab(expression(paste("Aymptotic growth rate (", lambda[1], ")"))) + 
  xlab("Maturation model") + #ylim(1, NA) +
  facet_wrap("Population", scales = "free")

# Figure of lionfish elasticities
# reshape the data
library(reshape2)
lf2 <- melt(lf, c("MatMod", "FertSurv", "JuvRep"), 
            c("Elasticity_1", "Elasticity_2", "Elasticity_3"),
            "Stage", value.name = "Elasticity")
levels(lf2$Stage) <- c("Larvae", "Juveniles", "Adults")

# Contrast maturation models
p1 <- ggplot(subset(lf2, FertSurv & JuvRep),
       aes(x = Stage, y = Elasticity)) +
  geom_point(aes(shape = MatMod), size = 4) +
  geom_line(aes(group = MatMod)) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.3)) +
  labs(shape = "Maturation Model") +
  ylab("") + xlab("")

# Contrast fertility coefficients
fert_names <- c("correct",
                "missing survival",
                "missing survival and maturing juveniles")
lf2 <- mutate(lf2,
              fert_index = 1+ (!FertSurv) * (!JuvRep) + !JuvRep,
              Fertility = fert_names[fert_index])

p2 <- ggplot(subset(lf2, MatMod == "AAS"),
       aes(x = Stage, y = Elasticity)) +
  geom_point(aes(shape = Fertility), size = 4) +
  geom_line(aes(group = Fertility)) +
  theme_classic()+
  theme(legend.position = c(0.7, 0.3)) +
  labs(shape = "Fertility Coefficient") +
  ylab("")

library(ggpubr)

pp <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = "auto", label.x = 0.95)
# THere seems to be no way to use plotmath with annotate_figure
annotate_figure(pp, left = "Elasticity of asymptotic growth rate to stage-specific survival")

# "anova" table of lambdas
lambda_stats <- function(st) {
  lambda <- filter(st, MatMod == "AAS" & FertSurv & JuvRep)$Lambda
  lambda_SAS <- filter(st, MatMod == "SAS" & FertSurv & JuvRep)$Lambda
  lambda_FAS <- filter(st, MatMod == "FAS" & FertSurv & JuvRep)$Lambda
  lambda_good_fert <- mean(c(lambda, lambda_FAS, lambda_SAS))
  lambda_no_juv <- mean(filter(st, (FertSurv) & (!JuvRep))$Lambda)
  lambda_no_surv <- mean(filter(st, (!FertSurv) & (!JuvRep))$Lambda)
  c(lambda, lambda_SAS - lambda, lambda_FAS - lambda, 
    lambda_no_juv - lambda_good_fert, lambda_no_surv - lambda_no_juv)
}

lambda_a <- tibble(lionfish = lambda_stats(lf),
                       "alligator (north)" = lambda_stats(aN),
                       `alligator (south)` = lambda_stats(aS))
rownames(lambda_a) <- c("$\\lambda_1$", "SAS", "FAS", "delayed maturity", "no survival")
knitr::kable(round(lambda_a, 3))
