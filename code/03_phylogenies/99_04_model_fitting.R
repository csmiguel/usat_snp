###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: fit models to bootstrap data
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)

bs_treemetrics <- readRDS("data/intermediate/bs_treemetrics.rds")
#i tried to fit glm with binomial distirbution and others but they do not fit well the data.
m1 <- glm(bs_tm ~ dnp + dnr, data = bs_treemetrics$dart_pelo, family = "binomial")

summary(m1)
plot(bs_tm~dnp, data = bs_treemetrics$dart_pelo)

#also gam do not fit well:

#http://environmentalcomputing.net/intro-to-gams/
library(mgcv)
gam_y <- mgcv::gam(bs_tm ~ s(dnp), data = h, method = "REML")
y_pred <- predict(gam_y, data.frame(dnp = seq(0, max(h$dnr), length.out = 1000)))

ggplot(data.frame(h), aes(dnp, bs_tm)) +
  geom_point() +
  geom_smooth(method = gam, formula = y ~ s(x))

gam.check(gam_y)
summary(gam_y)
