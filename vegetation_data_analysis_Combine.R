## Vegetation analysis to support '02-vegetationanalysis.Rmd'.

## Data manipulation --> analysis (this doc) --> Rmarkdown

## This document replaces prior documents including the
## 'OfficeDevelopmentVegetationX.Rmd' document series and the
## 'veg_data_analysisVX.R' series of scripts.



## -- SETUP -------------------------------------------------------------

date <- format(Sys.Date(), "%m%d%Y")

## Load required libraries
    library(labdsv)
    library(vegan)
    library(cluster)
    library(indicspecies)
    library(tidyr)
    library(data.table)
    library(MVN)
    library(stringr)
    library(dplyr)
    library(lmerTest)


## Check if these are necessary
    #library(dendextend)
    #library(Formula)
    #library(corrgram)


## Source data + helper R script

    source("vegetation_data_processing_Combined.R")
    source('../../../../RCode/R_Scripts/AICc_PERMANOVA.R')
    source('../../../../RCode/R_Scripts/repeat_multipatt.R')
    source("../../../../RCode/R_Scripts/group_PERMANOVA.R")
    source("../../../../RCode/R_Scripts/AICc_table_generation.R")

sample.covariates <- sample.covariates[ order(sample.covariates$SiteName), ]

matrify.tree.sp.only <- matrify.tree.sp.only[ order(rownames(matrify.tree.sp.only)), ]
matrify.tree.sp.size <- matrify.tree.sp.size[ order(rownames(matrify.tree.sp.size)), ]
dens.matrify.tree.sp.only <- dens.matrify.tree.sp.only[ order(rownames(dens.matrify.tree.sp.only)), ]

matrify.gc <- matrify.gc[ order(rownames(matrify.gc)), ]
dens.matrify.gc <- dens.matrify.gc[ order(rownames(dens.matrify.gc)), ]

matrify.shrub.site <- matrify.shrub.site[ order(rownames(matrify.shrub.site)), ]
dens.matrify.shrub.site <- dens.matrify.shrub.site[ order(rownames(dens.matrify.shrub.site)), ]

## -- DESCRIPTIVE STATISTICS ---------------------------------------------

matrify.tree.sp.only.native <-
    matrify.tree.sp.only[, colnames(matrify.tree.sp.only) %in%
                             tree.native$tree.species[tree.native$tree.origin == "native"]]
matrify.shrub.site.native <-
    matrify.shrub.site[, colnames(matrify.shrub.site) %in%
                           shrub.native$shrub.scientific.update[shrub.native$shrub.origin == "native"]]


    ## Create a table of relevant descriptive statistics for trees and shrubs
    ## first collapsed on SITES


descriptive.stats.site <- tibble(
    site = rownames(matrify.tree.sp.only),
    site.area = sample.covariates$acres,
    tree.site.abundance = rowSums(matrify.tree.sp.only),
    tree.site.density = rowSums(dens.matrify.tree.sp.only),
    tree.sp.richness = apply(matrify.tree.sp.only, MARGIN = 1, specnumber),
    tree.native.richness = apply(matrify.tree.sp.only.native, MARGIN = 1, specnumber),
    tree.Shannon = apply(matrify.tree.sp.only, MARGIN = 1, diversity),
    tree.effective.sp = exp(tree.Shannon),
    tree.native.Shannon = apply(matrify.tree.sp.only.native, MARGIN = 1, diversity),
    tree.native.effective.sp = exp(tree.native.Shannon),
    tree.prop.sp.rich.of.tot = tree.sp.richness / ncol(matrify.tree.sp.only),

    tree.conifer.native.abund = matrify.tree.sp.only$douglas.fir + matrify.tree.sp.only$western.hemlock +
        matrify.tree.sp.only$western.red.cedar,
    tree.conifer.nat.dens = tree.conifer.native.abund / sample.covariates$acres,
    tree.conifer.L.native = matrify.tree.sp.size$douglas.fir.L + matrify.tree.sp.size$western.red.cedar.L +
        matrify.tree.sp.size$western.hemlock.L,
    tree.conifer.L.nat.dens = tree.conifer.L.native / sample.covariates$acres,
    tree.conifer.ML.native = tree.conifer.L.native + matrify.tree.sp.size$douglas.fir.M +
        matrify.tree.sp.size$western.red.cedar.M + matrify.tree.sp.size$western.hemlock.M,
    tree.conifer.ML.nat.dens = tree.conifer.ML.native / sample.covariates$acres,
    tree.all.L.native = tree.conifer.L.native + matrify.tree.sp.size$big.leaf.maple.L +
        matrify.tree.sp.size$pacific.madrone.L +
        matrify.tree.sp.size$red.alder.L + matrify.tree.sp.size$black.cottonwood.L,
    tree.native.abundance = apply(matrify.tree.sp.only.native, MARGIN = 1, sum),
    tree.native.density = tree.native.abundance/sample.covariates$acres,
    
    
    shrub.site.abundance = rowSums(matrify.shrub.site),
    shrub.site.density = rowSums(dens.matrify.shrub.site),
    shrub.sp.richness = apply(matrify.shrub.site, MARGIN = 1, specnumber),
    shrub.native.richness = apply(matrify.shrub.site.native, MARGIN = 1, specnumber),
    shrub.Shannon = apply(matrify.shrub.site, MARGIN = 1, diversity),
    shrub.effective.sp = exp(shrub.Shannon),
    shrub.native.Shannon = apply(matrify.shrub.site.native, MARGIN = 1, diversity),
    shrub.native.effective.sp = exp(shrub.native.Shannon),
    shrub.prop.sp.rich.of.tot = shrub.sp.richness / ncol(matrify.shrub.site),
    shrub.native.abundance = apply(matrify.shrub.site.native, MARGIN = 1, sum),
    shrub.native.density = shrub.native.abundance/sample.covariates$acres,

    gc.pct.impervious = round(dens.matrify.gc$impervious.sqft * 100, 2),
    gc.pct.pervious = round(dens.matrify.gc$pervious.sqft * 100, 2),
    gc.pct.dense.veg = round(dens.matrify.gc$dense.veg * 100, 2),
    gc.pct.dirt.litter = round(dens.matrify.gc$dirt.litter * 100, 2),
    gc.pct.grass = round(dens.matrify.gc$grass * 100, 2),
    gc.pct.ivy = round(dens.matrify.gc$ivy * 100, 2),
    gc.pct.mulch = round(dens.matrify.gc$mulch * 100, 2),
    gc.pct.other = round((
        dens.matrify.gc$gravel + dens.matrify.gc$water
    ) * 100, 2)
)


management.landscaping$tree.conifer.nat.dens <- descriptive.stats.site$tree.conifer.nat.dens



        temp <- summarise_all(descriptive.stats.site[ , -1],
                              list(min, max, mean, sd, median)) %>%
            round(1)
        pretty.descriptive.stats.site <- tibble(
                    Metric = colnames(descriptive.stats.site)[-1],
                    Minimum = unlist(temp[1:38]),
                    Maximum = unlist(temp[39:76]),
                    Mean = unlist(temp[77:114]),
                    `Standard Deviation` = unlist(temp[115:152]),
                    Median = unlist(temp[153:190])

                )
        remove(temp)

        pretty.descriptive.stats.site$Metric <- c(
            "Site Area (acres)",
            "Tree Abundance",
            "Tree Density",
            "Tree Species Richness",
            "Native Tree Species Richness",
            "Tree Shannon Diversity",
            "Tree Effective Species Richness",
            "Native Tree Shannon Diversity",
            "Native Tree Effective Species Richness",
            "Proportion of Total Tree Species Richness",
            
            "Native Conifer Abundance",
            "Native Conifer Density",
            "Native Conifer Abundance Large Trees Only",
            "Native Conifer Density Large Trees Only",
            "Native Confier Abundance Medium & Large Trees",
            "Native Conifer Density Medium & Large Trees",
            "Native Tree Abundance Large Trees Only",
            "Native Tree Abundance",
            "Native Tree Density",
            
            "Shrub Abundance",
            "Shrub Density",
            "Shrub Species Richness",
            "Native Shrub Species Richness",
            "Shrub Shannon Diversity",
            "Shrub Effective Species Richness",
            "Native Shrub Shannon Diversity",
            "Native Shrub Effective Species Richness",
            "Proportion of Total Shrub Species Richness",
            "Native Shrub Abundance",
            "Native Shrub Density",
            
            "Percent Impervious",
            "Percent Pervious",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Dense Vegetation",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Dirt & Litter",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Grass",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Ivy",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Mulch",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Percent Other"
        )
        colnames(pretty.descriptive.stats.site)[5] <- "S.D."






    ## Create a table of relevant descriptive statistics for trees and shrubs
    ## first collapsed on TREE/SHRUB SPECIES

        descriptive.stats.tree.sp <- tibble(
            species.name = colnames(matrify.tree.sp.only),
            species.min = apply(matrify.tree.sp.only, 2, min),
            species.max = apply(matrify.tree.sp.only, 2, max),
            species.total.abundance = colSums(matrify.tree.sp.only),
            species.mean.abundance = colMeans(matrify.tree.sp.only),
            species.median.abundance = apply(matrify.tree.sp.only, 2, median),
            species.abundance.sd = apply(matrify.tree.sp.only, 2, sd),
            species.incidence = apply(matrify.tree.sp.only, 2, specnumber)

        )


        descriptive.stats.shrub.sp <- tibble(
            species.name = colnames(matrify.shrub.site),
            species.min = apply(matrify.shrub.site, 2, min),
            species.max = apply(matrify.shrub.site, 2, max),
            species.total.abundance = colSums(matrify.shrub.site),
            species.mean.abundance = colMeans(matrify.shrub.site),
            species.median.abundance = apply(matrify.shrub.site, 2, median),
            species.abundance.sd = apply(matrify.shrub.site, 2, sd),
            species.incidence = apply(matrify.shrub.site, 2, specnumber)

        )




    ## Compare observed vs. expected vegetation classes

        # Set up anova and contrasts testing:
        aov.vars <-
            data.frame(
                veg.class = sample.covariates$vegetation_class,
                conifer.dens = descriptive.stats.site$tree.conifer.ML.nat.dens,
                shrub.rich = descriptive.stats.site$shrub.sp.richness,
                shrub.effsr = descriptive.stats.site$shrub.effective.sp
            )

        # Test for significant differences in conifer density using one-way anova:
            anova(lm(conifer.dens ~ veg.class, data = aov.vars))
                # or equivalently:
                summary(aov(conifer.dens ~ veg.class, data = aov.vars))

        # And between group differences using a priori contrasts LL, MM, MD vs. MC, HH
            contrasts(aov.vars$veg.class) <- cbind(c(-3/2, 1, -3/2, 1, 1),
                                                   c(1, 0, -1, 0, 0),
                                                   c(0, 1, 0, -1, 0),
                                                   c(0, 0, 0, 1, -1))
            tree.aov <- summary.lm(aov(conifer.dens ~ veg.class, data = aov.vars))
            summary.aov(aov(conifer.dens ~ veg.class, data = aov.vars),
                        split = list(veg.class = list(
                            "HH/MC vs LL/MM/MD" = 1,
                            "HH vs MC" = 2,
                            "LL vs MD" = 3,
                            "MD vs MM" = 4

                        )))

        # Test for significant differences in shrub species richness
            summary(aov(shrub.effsr ~ veg.class, data = aov.vars))


            contrasts(aov.vars$veg.class) <- cbind(c(-3/2, 1, 1, -3/2, 1),
                                                   c(1, 0, 0, -1, 0),
                                                   c(0, 1, -1, 0, 0),
                                                   c(0, 0, 1, 0, -1))
            shrub.aov <- summary.lm(aov(shrub.effsr ~ veg.class, data = aov.vars))
            summary.aov(aov(shrub.effsr ~ veg.class, data = aov.vars),
                        split = list(veg.class = list(
                            "HH/MD vs LL/MM/MC" = 1,
                            "HH vs MD" = 2,
                            "LL vs MC" = 3,
                            "MC vs MM" = 4

                        )))






## -- INDEPENDENT VARIABLE SUMMARY --------------------------------------------------

        # Create a summary table to match the ones done for trees/shrubs across all
        # sites. Need to mention both sample and population statistics where
        # possible...

        # Ok, sample.covariates and population.covariates already have metrics
        # computed, so we can skip the first step (metric calculation) that was
        # necessary for trees and shrubs (descriptive.stats.site).


        # Instead, subset sample. and population.covariates to only include
        # relevant metrics, then run through summarise_all.

        # Sample:
            temp <- dplyr::select(sample.covariates,
                                  sample.acres = acres,
                                  sample.yrbuilt = YrBuiltAvg,
                                  sample.age2017 = age.2017,
                                  
                                  sample.parc500m.imp = dissolved_parcel_500m_buffer_impervious_500m_mean,
                                  sample.medincome = MedianHouseholdIncome_B19013e1,
                                  sample.propimmig = Proportion_ForeignBorn_B99051e5) %>%
                summarise_all(list(min, max, mean, sd, median)) %>%
                round(2)

            descriptive.sample.cov <- tibble(
                Metric = c("Sample: Area (acres)", "Sample: Year Built",
                           "Sample: Percent Impervious w/in 500m", "Sample: Median Income ($)",
                           "Sample: Percent Foreign-Born"),
                Minimum = unlist(temp[1:5]),
                Maximum = unlist(temp[6:10]),
                Mean = unlist(temp[11:15]),
                `Standard Deviation` = unlist(temp[16:20]),
                Median = unlist(temp[21:25])
            )
            remove(temp)

            # Population:
            temp <- dplyr::select(
                population.covariates,
                popn.acres = acres,
                popn.yrbuilt = YrBuiltAvg,
                popn.parc500m.imp = dissolved_parcel_500m_buffer_impervious_500m_mean,
                popn.medincome = MedianHouseholdIncome_B19013e1,
                popn.propimmig = Proportion_ForeignBorn_B99051e5
            ) %>%
                summarise_all(list(
                    ~min(., na.rm = TRUE),
                    ~max(., na.rm = TRUE),
                    ~mean(., na.rm = TRUE),
                    ~sd(., na.rm = TRUE),
                    ~median(., na.rm = TRUE)
                )) %>%
                round(2)
            
            descriptive.popn.cov <- tibble(
                Metric = c(
                    "Population: Area (acres)",
                    "Population: Year Built",
                    "Population: Percent Impervious w/in 500m",
                    "Population: Median Income ($)",
                    "Population: Percent Foreign-Born"
                ),
                Minimum = unlist(temp[1:5]),
                Maximum = unlist(temp[6:10]),
                Mean = unlist(temp[11:15]),
                `Standard Deviation` = unlist(temp[16:20]),
                Median = unlist(temp[21:25])
            )
            remove(temp)

            pretty.descriptive.covariates <- rbind(descriptive.sample.cov, descriptive.popn.cov)
            pretty.descriptive.covariates <- pretty.descriptive.covariates[ c(1,6,3,8,2,7,4,9,5,10), ]

            colnames(pretty.descriptive.covariates)[5] <- "S.D."





    # Need to handle grouped variables slightly differently...
            popn.categorical.covariates <- dplyr::select(population.covariates,
                                                         DistrictName, BldgQuality_lookup,
                                                         vegetation_class, acres, in.sample, SiteName) %>%
                filter(str_detect(vegetation_class, "HH|MM|MD|MC|LL")) %>%
                mutate(Frequency = 1)
            pop.site.count <- sum(popn.categorical.covariates$Frequency)

            samp.categorical.covariates <- filter(popn.categorical.covariates, str_detect(in.sample, "Y")) %>%
                filter(!str_detect(SiteName, "Forest.Park.X"))


            VegClass.p <- aggregate(list(Count = popn.categorical.covariates$Frequency),
                                    by = list(VegClass = popn.categorical.covariates$vegetation_class),
                                    FUN = sum) %>%
                mutate(bin = "Population") %>%
                mutate(Proportion = Count/pop.site.count)
            VegClass.s <- aggregate(list(Count = samp.categorical.covariates$Frequency),
                                    by = list(VegClass = samp.categorical.covariates$vegetation_class),
                                    FUN = sum) %>%
                mutate(bin = "Sample") %>%
                mutate(Proportion = Count/20)
            VegClass <- rbind(VegClass.p, VegClass.s); remove(VegClass.p, VegClass.s)
            VegClass$VegClass <- factor(VegClass$VegClass, levels = c("HH", "MC", "MD", "MM", "LL"),
                                        ordered = TRUE)


            Town.p <- aggregate(list(Count = popn.categorical.covariates$Frequency),
                                by = list(Town = popn.categorical.covariates$DistrictName), FUN = sum) %>%
                mutate(bin = "Population") %>%
                mutate(Proportion = Count/pop.site.count)
            Town.s <- aggregate(list(Count = samp.categorical.covariates$Frequency),
                                by = list(Town = samp.categorical.covariates$DistrictName), FUN = sum) %>%
                mutate(bin = "Sample") %>%
                mutate(Proportion = Count/20)
            Town <- rbind(Town.p, Town.s); remove(Town.p, Town.s)


            BldQual.p <- aggregate(list(Count = popn.categorical.covariates$Frequency),
                                   by = list(BldQual = popn.categorical.covariates$BldgQuality_lookup),
                                   FUN = sum) %>%
                mutate(bin = "Population") %>%
                mutate(Proportion = Count/sum(Count))
            BldQual.s <- aggregate(list(Count = samp.categorical.covariates$Frequency),
                                   by = list(BldQual = samp.categorical.covariates$BldgQuality_lookup),
                                   FUN = sum) %>%
                mutate(bin = "Sample") %>%
                mutate(Proportion = Count/20)
            BldQual <- rbind(BldQual.p, BldQual.s); remove(BldQual.p, BldQual.s)
            BldQual$BldQual <- factor(BldQual$BldQual,
                                      levels = c("GOOD/EXCELLENT", "GOOD",
                                                 "AVERAGE/GOOD", "AVERAGE",
                                                 "LOW/AVERAGE", "LOW COST"),
                                      ordered = TRUE)



















## -- FLEXIBLE BETA CLUSTER ANALYSIS ----------------------------------------------------

            tree.mahalanobis <- mahalanobis(dens.matrify.tree.sp.only,
                                            center = colMeans(dens.matrify.tree.sp.only),
                                            cov = cov(dens.matrify.tree.sp.only),
                                            inverted = TRUE)


            shrub.mahalanobis <- mahalanobis(dens.matrify.shrub.site,
                                             center = colMeans(dens.matrify.shrub.site),
                                             cov = cov(dens.matrify.shrub.site),
                                             inverted = TRUE)

            # Both tree and shrub communities contain candidate outliers based on mahalanobis distances.


    ## Choose value of beta based on indicator species analysis.

        beta.val <- c(-.1, -.2, -.25, -.3, -.4, -.5, -.6, -.7)


    ## Tree density
        ## Create candidate clusters
            i = 1

            for (i in 1:length(beta.val)) {
                alpha <- (beta.val[i] - 1)/-2

                beta.var <- agnes(vegdist(dens.matrify.tree.sp.only, method = "bray"),
                                  diss = TRUE, method = "flexible",
                                  par.method = c(alpha,alpha,beta.val[i]))
                plot(beta.var,  main = paste("SeqVar Beta =", beta.val[i]), which.plots = 2)

                print(paste(beta.val[i], "with Agglomerative Coefficient of ", round(beta.var$ac, 3)))
                if (beta.val[i] == -.25) {
                    beta.25.tree <- beta.var
                }
                if (beta.val[i] == -.4) {
                    beta.4.tree <- beta.var
                }
                if (beta.val[i] == -.5) {
                    beta.5.tree <- beta.var
                }
                if (beta.val[i] == -.6) {
                    beta.6.tree <- beta.var
                }

            }

        # All have pretty good agglomerative coefficients! (Above .7 is strong structure)
        # beta = -.1, -.2 are pretty chained, equvalent
        # -.25 and -.3 equivalent, still chained
        # I feel -.4 starts having less chained structure, but only top level changes
        # but -.5 thru -.7 equivalent and have better structure. 


        # based on the Agglomerative coefficient (which can be compared to
        # silhouette coefficient), beta = -.25 and above have strong structure
        # (see http://www.stat.berkeley.edu/~s133/Cluster2a.html)


        # tanglegram(as.dendrogram(beta.25.tree), as.dendrogram(beta.4.tree))


        ## Determine number of groups for each beta level.

            beta.25.tree.groups <-
                as.data.frame.matrix(cutree(as.hclust(beta.25.tree),
                                            k = 2:5))
            
            beta.5.tree.groups <-
                as.data.frame.matrix(cutree(as.hclust(beta.5.tree),
                                            k = 2:5))
            


        for (k in 2:4) {

            summary(multipatt(x = dens.matrify.tree.sp.only, max.order = 2,
                              cluster = beta.5.tree.groups[[k - 1]], func = "IndVal.g"),
                    indvalcomp = TRUE, alpha = .05, minstat = .8)

        }

        ## -.25 has consistently confused indicator species structure... few species, mostly in group 1/ 1+2 (same, depends on k)
        ## -.5 is better, has group specific species for both 1 and 2 (k = 2). k = 3 and 4 have combined groups (1+2, 2+3).

        ## From this it seems that -.5 with k = 2 is a good grouping solution.


        ## Indicator species analysis

        tree.multipatt <- repeat.multipatt(matrix.name = dens.matrify.tree.sp.only,
                                           p.cutoff = .05,
                                           freq.cutoff = .5,
                                           cluster.name = beta.5.tree.groups[[1]],
                                           func.name = "IndVal.g",
                                           xlab.input = "Tree cluster indicator species",
                                           repeats = 100,
                                           quiet = FALSE)



    ## Shrub density
        ## Create candidate clusters
        i = 1

        for (i in 1:length(beta.val)) {
            alpha <- (beta.val[i] - 1)/-2

            beta.var <- agnes(vegdist(dens.matrify.shrub.site, method = "bray"),
                              diss = TRUE, method = "flexible",
                              par.method = c(alpha,alpha,beta.val[i]))
            plot(beta.var,  main = paste("SeqVar Beta =", beta.val[i]), which.plots = 2)

            print(paste(beta.val[i], "with Agglomerative Coefficient of ", round(beta.var$ac, 3)))
            if (beta.val[i] == -.25) {
                beta.25.shrub <- beta.var
            }
            if (beta.val[i] == -.4) {
                beta.4.shrub <- beta.var
            }
            if (beta.val[i] == -.5) {
                beta.5.shrub <- beta.var
            }
            if (beta.val[i] == -.6) {
                beta.6.shrub <- beta.var
            }

        }

        # First, note that the agglomerative coefficients are not as good as
        # trees.

        # beta = -.2, -.25, -.3 are largely equvalent
        # -.4 and above are equivalent.
        # Not sure which one is better; compare -.25 and -.5 to match trees.


        # based on the Agglomerative coefficient (which can be compared to
        # silhouette coefficient), beta = -.4 and to -.7 have strong structure
        # (see http://www.stat.berkeley.edu/~s133/Cluster2a.html)


        # tanglegram(as.dendrogram(beta.4.shrub), as.dendrogram(beta.6.shrub))


        ## Determine number of groups for each beta level.

        # beta.25.shrub.groups <- as.data.frame.matrix(cutree(as.hclust(beta.25.shrub), k = 2:5))
        beta.5.shrub.groups <- as.data.frame.matrix(cutree(as.hclust(beta.5.shrub), k = 2:5))



        for (k in 2:4) {

            summary(multipatt(x = dens.matrify.shrub.site, max.order = 2,
                              cluster = beta.5.shrub.groups[[k - 1]], func = "IndVal.g"),
                    indvalcomp = TRUE, alpha = .05, minstat = .8)

        }

        ## -.25 suggests 2 groups, higher k leads to more combined groups.
        ## -.5 is better in that it has more indicator species. As with -.25, k = 2 avoids cross groups.

        ## From this it seems that -.5 with k = 2 is a good grouping solution.
        ## This is convient since it matches the tree grouping.


        ## Indicator species analysis

        shrub.multipatt <- repeat.multipatt(matrix.name = dens.matrify.shrub.site,
                                            p.cutoff = .05,
                                            freq.cutoff = .5,
                                           cluster.name = beta.5.shrub.groups[[1]],
                                           func.name = "IndVal.g",
                                           xlab.input = "Shrub cluster indicator species",
                                           repeats = 100,
                                           quiet = FALSE)






    # Create some summary tables for each of the clusters.

        vegetation.clusters <- data.frame(SiteName = rownames(beta.5.shrub.groups),
                                          tree.cluster = beta.5.tree.groups[[1]],
                                          shrub.cluster = beta.5.shrub.groups[[1]])
            vegetation.clusters$tree.cluster.name <- as.factor(ifelse(vegetation.clusters$tree.cluster == 1,
                                                            "Native", "Ornamental Deciduous"))
            vegetation.clusters$shrub.cluster.name <- as.factor(ifelse(vegetation.clusters$shrub.cluster == 1,
                                                             "Cultivated","Native"))


        tree.cluster.native <- matrify.tree.sp.only[
                                      which(rownames(matrify.tree.sp.only) %in%
                                                vegetation.clusters$SiteName[which(
                                                    vegetation.clusters$tree.cluster == 1)]) , ] %>%
                                colSums()
        tree.cluster.native <- sort(tree.cluster.native, decreasing = TRUE)

        tree.cluster.orndec <- matrify.tree.sp.only[
                                        which(rownames(matrify.tree.sp.only) %in%
                                                  vegetation.clusters$SiteName[which(
                                                      vegetation.clusters$tree.cluster == 2)]) , ] %>%
                                colSums()
        tree.cluster.orndec <- sort(tree.cluster.orndec, decreasing = TRUE)

        shrub.cluster.cultiv <- matrify.shrub.site[
                                        which(rownames(matrify.shrub.site) %in%
                                                  vegetation.clusters$SiteName[which(
                                                      vegetation.clusters$shrub.cluster == 1)]) , ] %>%
                                colSums()
        shrub.cluster.cultiv <- sort(shrub.cluster.cultiv, decreasing = TRUE)

        shrub.cluster.native <- matrify.shrub.site[
                                        which(rownames(matrify.shrub.site) %in%
                                                  vegetation.clusters$SiteName[which(
                                                      vegetation.clusters$shrub.cluster == 2)]) , ] %>%
                                colSums()
        shrub.cluster.native <- sort(shrub.cluster.native, decreasing = TRUE)


        # zones <- str_sub(rownames(matrify.shrub.zones), 
        #          start = 1 + str_locate(rownames(matrify.shrub.zones), pattern = "\\.")[,1])
        # zone.multipatt <- multipatt(matrify.shrub.zones, zones)


        # zone.multipatt <- repeat.multipatt(matrix.name = matrify.shrub.zones,
        #                                     p.cutoff = .05,
        #                                     freq.cutoff = .5,
        #                                     cluster.name = as.factor(zones),
        #                                     func.name = "IndVal.g",
        #                                     xlab.input = "Shrub cluster indicator species",
        #                                     repeats = 100,
        #                                     quiet = FALSE)
        # 



## -- GROUP DIFFERENCES --------------------
        # Want to see if impervious surface is equal between tree groups.
        temp <-
            data.frame(
                cbind(beta.5.tree.groups[[1]], dens.matrify.gc$impervious.sqft),
                rownames(dens.matrify.gc)
            )
        colnames(temp) <- c("group", "imp", "site")
        
        bartlett.test(temp$imp, temp$group)
        fligner.test(temp$imp, temp$group) #homoskedasticity
        
        
        anova(lm(imp ~ group, temp))
        anova(lm(group ~ imp, temp))
        t.test(temp[temp$group == 1, "imp"], temp[temp$group == 2, "imp"])
        
        means <-
            temp %>% group_by(group) %>% summarise_at(vars(imp), list(mean, sd))
        group1 <- temp[temp$group == 1,]
        group2 <- temp[temp$group == 2,]
        
        MSE <-
            (sum((group1$imp - means[[1, 2]]) ^ 2) + sum((group2$imp - means[[2, 2]]) ^
                                                             2)) / 18
        # this matches the anova(lm) so think its good
        
        adonis2(temp$imp ~ as.factor(temp$group))
        
        remove(temp, means, MSE)
        #Garbage <- rnorm(nrow(temp))
        # anova(lmer( imp ~ Garbage + (1|group), temp), ddf="Kenward-Roger") so lmer requires a mixed effect model, which I don't have...
        
 # ok testing done, now to see if tree and shrub groups differ in their means for certain variables
        
        
        # grass
        
        grass.tree <-
            adonis2(
                dens.matrify.gc$grass ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        grass.shrub <-
            adonis2(
                dens.matrify.gc$grass ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # impervious
        
        imp.tree <-
            adonis2(
                dens.matrify.gc$impervious.sqft ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        imp.shrub <-
            adonis2(
                dens.matrify.gc$impervious.sqft ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # dirt/litter
        
        dirt.tree <-
            adonis2(
                dens.matrify.gc$dirt.litter ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        dirt.shrub <-
            adonis2(
                dens.matrify.gc$dirt.litter ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # mulch
        
        
        mulch.tree <-
            adonis2(
                dens.matrify.gc$mulch ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        mulch.shrub <-
            adonis2(
                dens.matrify.gc$mulch ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # abundance
        count.tree <-
            adonis2(
                rowSums(matrify.tree.sp.only) ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        count.shrub <-
            adonis2(
                rowSums(matrify.shrub.site) ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # density
        dens.tree <-
            adonis2(
                rowSums(dens.matrify.tree.sp.only) ~ vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        dens.shrub <-
            adonis2(
                rowSums(dens.matrify.shrub.site) ~ vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # median conifer height
        
        mheight.tree <-
            adonis2(
                management.landscaping$height.m.median ~
                    vegetation.clusters$tree.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        mheight.shrub <-
            adonis2(
                management.landscaping$height.m.median ~
                    vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # area
        
        area.tree <- adonis2(
            sample.covariates$acres ~
                vegetation.clusters$tree.cluster.name,
            permutations = 99999,
            method = "eucl"
        )
        area.shrub <- adonis2(
            sample.covariates$acres ~
                vegetation.clusters$shrub.cluster.name,
            permutations = 99999,
            method = "eucl"
        )
        
        
        # dead wood
        
        dw.tree <- adonis2(
            management.landscaping$dead.wood ~
                vegetation.clusters$tree.cluster.name,
            permutations = 99999,
            method = "eucl"
        )
        dw.shrub <- adonis2(
            management.landscaping$dead.wood ~
                vegetation.clusters$shrub.cluster.name,
            permutations = 99999,
            method = "eucl"
        )
        
        
        # median conifer dbh
        
        mdbh.tree <- adonis2(
            management.landscaping$DBH.in.median ~
                vegetation.clusters$tree.cluster.name,
            permutations = 99999,
            method = "eucl"
        )
        mdbh.shrub <-
            adonis2(
                management.landscaping$DBH.in.median ~
                    vegetation.clusters$shrub.cluster.name,
                permutations = 99999,
                method = "eucl"
            )
        
        # irrigation--only 3 nos, so unreliable
        
        irr.tree <-
            chisq.test(management.landscaping$Irrigation[-1],
                       vegetation.clusters$tree.cluster.name[-1])
        irr.shrub <-
            chisq.test(management.landscaping$Irrigation[-1],
                       vegetation.clusters$shrub.cluster.name[-1])
        
        # mulch-- only 3 nos
        
        mulch.tree <- chisq.test(management.landscaping$Mulch,
                                 vegetation.clusters$tree.cluster.name)
        mulch.shrub <- chisq.test(management.landscaping$Mulch,
                                  vegetation.clusters$shrub.cluster.name)
        
        # herb-- only 3 nos
        
        herb.tree <- chisq.test(management.landscaping$Herbicide,
                                vegetation.clusters$tree.cluster.name)
        herb.shrub <- chisq.test(management.landscaping$Herbicide,
                                 vegetation.clusters$shrub.cluster.name)
        
        # fert-- only 3 nos
        
        fert.tree <- chisq.test(management.landscaping$Fertilizer,
                                vegetation.clusters$tree.cluster.name)
        fert.shrub <- chisq.test(management.landscaping$Fertilizer,
                                 vegetation.clusters$shrub.cluster.name)
        

        
            adonis2(
                rowSums(matrify.tree.sp.only) ~
                    sample.covariates$MedianHouseholdIncome_B19013e1,
                permutations = 99999,
                method = "eucl"
            )
            
            adonis2(
                dens.matrify.tree.sp.only ~
                    management.landscaping$height.m.median,
                permutations = 99999,
                method = "bray"
            )
        
        for (i in 1:20){
            x <- adonis2(
                matrify.tree.sp.only[-i,] ~
                    sample.covariates$MedianHouseholdIncome_B19013e1[-i],
                permutations = 99999,
                method = "bray"
            )    
            print(x)
        }
            
            
## -- NMDS --------------------------------------------------------------------------

    # Tree density NMDS:

        ## Two dimensions are sufficient (Stress = 0.12)
            NMDS.tree.2d <- metaMDS(dens.matrify.tree.sp.only, distance = "bray", k = 2, trymax = 9999,
                                    autotransform = FALSE)
            for (i in 1:100) {

                NMDS.tree.2d <- metaMDS(dens.matrify.tree.sp.only, distance = "bray", k = 2, trymax = 9999,
                                        autotransform = FALSE, previous.best = NMDS.tree.2d)

            }

            # Diagnostics:
                stressplot(NMDS.tree.2d)
                plot(goodness(NMDS.tree.2d))
                plot(NMDS.tree.2d, display = "sites")

                # Interesting note: k = 1 NMDS has high stress (0.25) but perfectly splits along ordination groups.
                plot(metaMDS(dens.matrify.tree.sp.only, distance = "bray", k = 1, trymax = 9999,
                        autotransform = FALSE))

    # Shrub density NMDS:
        NMDS.shrub.2d <- metaMDS(dens.matrify.shrub.site, distance = "bray", k = 2, trymax = 9999,
                                autotransform = FALSE)
        for (i in 1:100) {

            NMDS.shrub.2d <- metaMDS(dens.matrify.shrub.site, distance = "bray", k = 2, trymax = 9999,
                                    autotransform = FALSE, previous.best = NMDS.shrub.2d)

        }
#
#         NMDS.shrub.3d <- metaMDS(dens.matrify.shrub.site, distance = "bray", k = 3, trymax = 9999,
#                                  autotransform = FALSE)
#         for (i in 1:100){
#
#             NMDS.shrub.3d <- metaMDS(dens.matrify.shrub.site, distance = "bray", k = 3, trymax = 9999,
#                                      autotransform = FALSE, previous.best = NMDS.shrub.3d)
#
#         }

        # The two dimensional is probably sufficient...
            stressplot(NMDS.shrub.2d)
            plot(goodness(NMDS.shrub.2d))
            plot(NMDS.shrub.2d, display = "sites")




## -- PERMANOVA --------------------------------------------------------------------


            ## Group Membership


            vegetation.clusters$vegetation_class <-
                as.factor(sample.covariates$vegetation_class)
            write.csv(vegetation.clusters, "vegetationClusters.csv")

            # Tree Density
            PERMANOVA.tree.groups <-
                group.PERMANOVA(
                    var.names = c("tree.cluster.name", "vegetation_class"),
                    var.table =  vegetation.clusters,
                    var.table.c = "vegetation.clusters",
                    species.table = dens.matrify.tree.sp.only,
                    species.table.c = "dens.matrify.tree.sp.only",
                    num.control.vars = 0,
                    by.adonis2 = "terms"
                )



            # Shrub Density
            PERMANOVA.shrub.groups <-
                group.PERMANOVA(
                    var.names = c("shrub.cluster.name", "vegetation_class"),
                    var.table =  vegetation.clusters,
                    var.table.c = "vegetation.clusters",
                    species.table = dens.matrify.shrub.site,
                    species.table.c = "dens.matrify.shrub.site",
                    num.control.vars = 0,
                    by.adonis2 = "terms"
                )


            ## Based on above, vegetation class is significant for both trees
            ## and shrubs; betadispersion test not significant so this is based
            ## on location, not multivariate spread

            ## However, the clusters appear to be better predictors, and
            ## multivariate spread is not significant here, either.





    # Socio-economic variables...
        socio.economic.vars <- c("acres",
                                 "Proportion_ForeignBorn_B99051e5",
                                 "DistrictName",
                                 "MedianHouseholdIncome_B19013e1",
                                 "dissolved_parcel_500m_buffer_impervious_500m_mean",
                                 "age.2017",
                                 "BldgQuality_lookup",
                                 "landrentperacre",
                                 "SMV_mean",
                                 "TV_mean")

            # Tree density

        PERMANOVA.socio.economic.tree.dens <-
            group.PERMANOVA(
                var.names = socio.economic.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = dens.matrify.tree.sp.only,
                species.table.c = "dens.matrify.tree.sp.only",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )


            # Shrub density

            PERMANOVA.socio.economic.shrub.dens <-
                group.PERMANOVA(
                    var.names = socio.economic.vars,
                    var.table =  sample.covariates,
                    var.table.c = "sample.covariates",
                    species.table = dens.matrify.shrub.site,
                    species.table.c = "dens.matrify.shrub.site",
                    num.control.vars = 0,
                    by.adonis2 = "terms"
                )


    ## Landscaping variables
            management.landscaping$tree.cluster.name <-
                vegetation.clusters$tree.cluster.name

            landscaping.vars <- c(
                "tree.cluster.name",
                "dead.wood",
                "height.m.median",
                "stand.predate.development",
                "tree.conifer.nat.dens"
            )

        # Tree Density--this doesn't make a ton of sense, given that these are
        # mostly descripters of the community composition itself...

            # PERMANOVA.landscaping.tree.dens <-
            #     group.PERMANOVA(
            #         var.names = landscaping.vars[-1],
            #         var.table =  management.landscaping,
            #         var.table.c = "management.landscaping",
            #         species.table = dens.matrify.tree.sp.only,
            #         species.table.c = "dens.matrify.tree.sp.only",
            #         num.control.vars = 0,
            #         by.adonis2 = "terms"
            #     )

            # Shrub Density
            PERMANOVA.landscaping.shrub.dens <-
                group.PERMANOVA(
                    var.names = landscaping.vars,
                    var.table =  management.landscaping,
                    var.table.c = "management.landscaping",
                    species.table = dens.matrify.shrub.site,
                    species.table.c = "dens.matrify.shrub.site",
                    num.control.vars = 0,
                    by.adonis2 = "terms"
                )




            ## Median DF height, Density of native conifers, tree cluster name,
            ## stands predating develpment explain variation in shrub community
            ## composition.


    ## Management variables
        # create the matrix for the management variables:


            management.vars <-
                c(
                    "Herbicide",
                    "Fertilizer",
                    "Irrigation",
                    "Mulch",
                    "Mushroom"
                )



            # Tree Density--this doesn't make a ton of sense, given that these
            # are mostly descripters of the community composition itself...


        PERMANOVA.management.tree.dens <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  management.landscaping,
                var.table.c = "management.landscaping",
                species.table = dens.matrify.tree.sp.only,
                species.table.c = "dens.matrify.tree.sp.only",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )

        # Shrub Density
        PERMANOVA.management.shrub.dens <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  management.landscaping,
                var.table.c = "management.landscaping",
                species.table = dens.matrify.shrub.site,
                species.table.c = "dens.matrify.shrub.site",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )

            ## None are significant.








    # Compare models using AICc.

            # Trees: Tested non-management and tree cluster variables. Only
            # vegetation class and tree cluster was significant. Vegetation
            # class should probably be the "baseline" since it is how I sampled!


                AICc.tree.vc <- AICc.PERMANOVA2(
                    adonis2.model = adonis2(
                        dens.matrify.tree.sp.only ~
                            factor(sample.covariates$vegetation_class),
                        permutations = 99999,
                        method = "bray"
                    )
                )

                AICc.tree.gp <- AICc.PERMANOVA2(
                    adonis2.model = adonis2(
                        dens.matrify.tree.sp.only ~
                            factor(vegetation.clusters$tree.cluster.name),
                        permutations = 99999,
                        method = "bray"
                    )
                )

                # Tree cluster is a significantly better model than vegetation
                # class for explaining tree community composition.

                # and the least bad of the socio-economic variables:
                AICc.tree.mhhi <- AICc.PERMANOVA2(
                    adonis2.model = adonis2(
                        dens.matrify.tree.sp.only ~
                            sample.covariates$MedianHouseholdIncome_B19013e1,
                        permutations = 99999,
                        method = "bray"
                    )
                )




            # Shrubs: Tested non-management, shrub/tree cluster variables, and
            # landscaping variables (tree stuff). Significant variables included
            # shrub/tree cluster variables, median DF height, Stands predating
            # development, and native conifer density. The tree group/tree
            # variables are likely highly correlated but AICc should indicate
            # best model.

                AICc.shrub.vc <- AICc.PERMANOVA2(
                    adonis2.model = adonis2(
                        dens.matrify.shrub.site ~
                            factor(sample.covariates$vegetation_class),
                        permutations = 99999,
                        method = "bray"
                    )
                )

                AICc.shrub.gp <- AICc.PERMANOVA2(
                    adonis2.model = adonis2(
                        dens.matrify.shrub.site ~
                            factor(vegetation.clusters$shrub.cluster.name),
                        permutations = 99999,
                        method = "bray"
                    )
                )


        # Shrub cluster is also a substantial improvement over the vegetation
        # class.

        # Development and landscaping variables:
        #



                AICc.shrub.table <- AICc.table.all(
                    sig.vars = paste0("management.landscaping$", landscaping.vars[-2]),
                    matrix.char = "dens.matrify.shrub.site",
                    extra.var = TRUE,
                    extra.var.char = c(
                        "sample.covariates$MedianHouseholdIncome_B19013e1",
                        "vegetation.clusters$shrub.cluster.name",
                        "sample.covariates$vegetation_class"
                    ),
                    comb.incl = 1:2,
                    type = "AICc",
                    perm = 99999
                )


            ## need to fix Delta AICc
            AICc.shrub.table$`Delta AICc`[1:11] <- AICc.shrub.table$AICc.values[1:11] -
                min(AICc.shrub.table$AICc.values[1:11])
            AICc.shrub.table$`Delta AICc`[12:13] <- AICc.shrub.table$AICc.values[12:13] -
                min(AICc.shrub.table$AICc.values[12:13])

            shrub.model.summary <- AICc.shrub.table

            # None of the different AICc calculations really changes Delta AICc.
            
            AICc.all.shrub <- AICc.table.all(
                sig.vars = c(paste0("sample.covariates$", socio.economic.vars), 
                             paste0("management.landscaping$", landscaping.vars)),
                matrix.char = "dens.matrify.shrub.site",
                extra.var = FALSE,
                comb.incl = 1,
                type = "AICc",
                perm = 99999
            )

            AICc.all.tree <- AICc.table.all(
                sig.vars = c(paste0("sample.covariates$", socio.economic.vars), 
                             paste0("management.landscaping$", landscaping.vars)[2]),
                matrix.char = "dens.matrify.tree.sp.only",
                extra.var = FALSE,
                comb.incl = 1,
                type = "AICc",
                perm = 99999
            )

        # Note that sqrt() makes distribution for tree conifer density slightly more normal and increases significance slightly.


    # Create summary tables of PERMANOVA/AICc results

        # Trees

                tree.model.summary <- tibble(
                    Model = c(
                        "Vegetation Class",
                        "Tree Group Membership",
                        "Socio-economic Vars",
                        "Site Area (acres)"
                    ),
                    `Delta AICc` = c(
                        round(abs(AICc.tree.vc$AICc - AICc.tree.gp$AICc), 2),
                        0,
                        "--",
                        "--"
                        #paste0("> ", abs(AICc.tree.mhhi$AICc - AICc.tree.gp$AICc))
                    ),
                    `Pseudo-_F_` = c(
                        round(PERMANOVA.tree.groups$pseudo.F[2], 2),
                        round(PERMANOVA.tree.groups$pseudo.F[1], 2),
                        paste0("< ", round(PERMANOVA.socio.economic.tree.dens$pseudo.F[5], 2)),
                        round(PERMANOVA.socio.economic.tree.dens$pseudo.F[1], 2)
                    ),
                    `p-value` = c(
                        round(PERMANOVA.tree.groups$adj.F.pval[2], 5),
                        round(PERMANOVA.tree.groups$adj.F.pval[1], 5),
                        paste0("> ", round(PERMANOVA.socio.economic.tree.dens$adj.F.pval[5], 2)),
                        round(PERMANOVA.socio.economic.tree.dens$adj.F.pval[1], 2)
                    ),
                    `% Variance Explained per df` = c(
                        PERMANOVA.tree.groups$avg.var.explnd[2],
                        PERMANOVA.tree.groups$avg.var.explnd[1],
                        paste0("< ", PERMANOVA.socio.economic.tree.dens$avg.var.explnd[5]),
                        PERMANOVA.socio.economic.tree.dens$avg.var.explnd[1]
                    )
                )



        # Shrubs

        shrub.model.summary$Model <-  c("Tree Group Membership",
                                        "Median Douglas Fir Height",
                                        "Stands Predate Development",
                                        "Native Conifer Density",
                                        #"Dead Wood (cube root)",
                                        "Tree Group + Median DF Height",
                                        "Tree Group + Stands Predate Development",
                                        "Tree Group + Native Conifer Density",
                                        #"Tree Group + Dead Wood (cube root)",
                                        "Median DF Height + Stands Predate Development",
                                        "Median DF Height + Native Conifer Density",
                                        "Stands Predate Development + Native Conifer Density",
                                        #"Stands Predate Development + Dead Wood (cube root)",
                                        #"Stands Predate Development + Median DF Height (m^2)",
                                        #"Native Conifer Density + Dead Wood (cube root)",
                                        #"Native Conifer Density + Median DF Height (m^2)",
                                        #"Dead Wood (cube root) and Median DF Height (m^2)",
                                         "Median Household Income",
                                         "Shrub Cluster Name",
                                         "Vegetation Class"
                                         )














         

## -- EXTRAPOLATION -------------------------------------------------------------------------


    # Need to look at means for different veg classes for extrapolation

    extrapolate.site.stats <- merge(descriptive.stats.site[ , 1:18],
                                    vegetation.clusters[ , c(1,4:6)],
                                    by.x = "site", by.y = "SiteName", sort = F)

            extrapolate.site.stats$tree.cluster.name <- 
                as.character(extrapolate.site.stats$tree.cluster.name)
            
        cluster.site.stats <- group_by(extrapolate.site.stats[ , 2:19], 
                                       tree.cluster.name) %>%
                                summarise_all(list(min = min, max = max,
                                                   mean = mean, sd = sd,
                                                   median = median)) 
        
        cluster.site.stats[2:86] <- round(cluster.site.stats[2:86] , 2)

    pretty.cluster.site.stats <- tibble(
            tree.cluster.name = cluster.site.stats$tree.cluster.name,
            site.area.range = paste(
                cluster.site.stats$site.area_min,
                "-",
                cluster.site.stats$site.area_max
            ),
            site.area.mean = paste(
                cluster.site.stats$site.area_mean,
                "(+/-",
                cluster.site.stats$site.area_sd,
                ")"
            ),
            tree.site.density.range = paste(
                cluster.site.stats$tree.site.density_min,
                "-",
                cluster.site.stats$tree.site.density_max
            ),
            tree.site.density.mean = paste(
                cluster.site.stats$tree.site.density_mean,
                "(+/-",
                cluster.site.stats$tree.site.density_sd,
                ")"
            ),
            tree.con.nat.dens.range = paste(
                cluster.site.stats$tree.conifer.nat.dens_min,
                "-",
                cluster.site.stats$tree.conifer.nat.dens_max
            ),
            tree.con.nat.dens.mean = paste(
                cluster.site.stats$tree.conifer.nat.dens_mean,
                "(+/-",
                cluster.site.stats$tree.conifer.nat.dens_sd,
                ")"
            )
        )



    
    
    vegclass.site.stats <-
        group_by(extrapolate.site.stats[, c(2, 4:5, 13, 21)], vegetation_class) %>%
        summarise_all(list(min = min, max = max, mean = mean, sd = sd, median = median))
    
    pretty.vegclass.site.stats <- tibble(
        vegetation_class = vegclass.site.stats$vegetation_class,
        site.area.range = paste(
            vegclass.site.stats$site.area_min,
            "-",
            vegclass.site.stats$site.area_max
        ),
        site.area.mean = paste(
            vegclass.site.stats$site.area_mean,
            "(+/-",
            vegclass.site.stats$site.area_sd,
            ")"
        ),
        tree.site.density.range = paste(
            vegclass.site.stats$tree.site.density_min,
            "-",
            vegclass.site.stats$tree.site.density_max
        ),
        tree.site.density.mean = paste(
            vegclass.site.stats$tree.site.density_mean,
            "(+/-",
            vegclass.site.stats$tree.site.density_sd,
            ")"
        ),
        tree.con.nat.dens.range = paste(
            vegclass.site.stats$tree.conifer.nat.dens_min,
            "-",
            vegclass.site.stats$tree.conifer.nat.dens_max
        ),
        tree.con.nat.dens.mean = paste(
            vegclass.site.stats$tree.conifer.nat.dens_mean,
            "(+/-",
            vegclass.site.stats$tree.conifer.nat.dens_sd,
            ")"
        )
    )


        vegclass.site.stats <- merge(vegclass.site.stats, VegClass[1:5 , c(1,2,4)],
                                     by.x = "vegetation_class", by.y = "VegClass")

    # Need to double check how these are calculated

        weighted.average.vegclass <- tibble(
            ` ` = "Population Weighted Mean",
            `Site Area` = sum(
                vegclass.site.stats$site.area_mean * vegclass.site.stats$Proportion
            ),
            `Species Richness` = sum(
                vegclass.site.stats$tree.sp.richness_mean *
                    vegclass.site.stats$Proportion
            ),
            `Tree Density` = sum(
                vegclass.site.stats$tree.site.density_mean *
                    vegclass.site.stats$Proportion
            ),
            `Native Conifer Density` = sum(
                vegclass.site.stats$tree.conifer.nat.dens_mean *
                    vegclass.site.stats$Proportion
            )
        )




## -- Save everything... ----------------------------------------------------------

save.image(file = paste0("vegetation_data_analysis_combined_", date, ".RData"))



## QUICK TEST
## 
## 
        df <- tibble(a = rep(0, 20),
                     b = a,
                     c = a,
                     d = a,
                     e = a,
                     f = a,
                     g = a,
                     x = a,
                     y = a,
                     z = a,
                     h = sample(1:10, 20, replace = T),
                     i= sample(1:10, 20, replace = T),
                     j= sample(1:10, 20, replace = T),
                     k= sample(1:10, 20, replace = T),
                     l= sample(1:10, 20, replace = T),
                     m= sample(1:10, 20, replace = T),
                     n= sample(1:10, 20, replace = T),
                     o= sample(1:10, 20, replace = T),
                     p= sample(1:10, 20, replace = T)
                    )
adonis(df ~ sample(1:10, 20, replace = T))
adonis(c(df$a, df$h) ~ sample(1:10, 40, replace = T), method = "eucl")

