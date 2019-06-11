# Script to annmze vegetatn data to prepare for publication. 

source("../../../../RCode/R_Scripts/anon_data.R")

# population and site covariates

population.covariates <-
    read.csv(
        "../../DataAnalysis/SiteStats/sitecovariateswtax.csv",
        stringsAsFactors = FALSE,
        strip.white = TRUE
    )

    population.covariates <- population.covariates[, -1]
    population.covariates$BldgQuality_lookup <-
        trimws(population.covariates$BldgQuality_lookup, which = "right")
    population.covariates$acres <-
        population.covariates$DissolvedArea * .00002296
    
    
neighborhood.veg <- 
    read.csv(
        "../../DataRepository/SiteData/dissolved_parcel_neighborhood_veg.csv",
        stringsAsFactors = FALSE,
        strip.white = TRUE
    )


population.covariates <-
    merge(
        population.covariates,
        neighborhood.veg[, c(3, 8:10)],
        by.x = "Grp_PIN",
        by.y = "PIN",
        all.x = TRUE
    )

population.covariates <- population.covariates[ , -c(1,3)]

population.covariates <- anon.data(
        input.table = population.covariates, 
        lookup.table = "../../DataAnalysis/site_lookup.csv",
        site.column = "SiteName",
        match.column = "site.name",
        anon.column = "site.abbr",
        output.path = "popcovariateswtax.csv",
        write.to.disk = TRUE
    )

## Maintenance data

maintenance.data <- anon.data(
        input.table = "../../DataRepository/VegetationData/maintenance_data.csv",
        lookup.table = "../../DataAnalysis/site_lookup.csv",
        site.column = "Site",
        match.column = "site.name",
        anon.column = "site.abbr",
        output.path = "maintenance_data.csv",
        write.to.disk = TRUE
    )


## Douglas fir data

DF.measurements <- anon.data(
        input.table = "../../DataRepository/VegetationData/DF_measurements2.csv",
        lookup.table = "../../DataAnalysis/site_lookup.csv",
        site.column = "Site",
        match.column = "site.name",
        anon.column = "site.abbr",
        output.path = "DF_measurements2.csv",
        write.to.disk = TRUE
    )


## Site tree data

tree.data <- anon.data(
        input.table = "../../DataRepository/VegetationData/TreeData/treedata.csv",
        lookup.table = "../../DataAnalysis/site_lookup.csv",
        site.column = "site",
        match.column = "site.name",
        anon.column = "site.abbr",
        output.path = "treedata.csv",
        write.to.disk = TRUE
    )


## Site shrub data

shrub.data <- anon.data(
    input.table = "../../DataRepository/VegetationData/ShrubsCSV.csv",
    lookup.table = "../../DataAnalysis/site_lookup.csv",
    site.column = "Site",
    match.column = "site.name",
    anon.column = "site.abbr",
    write.to.disk = FALSE
)

shrub.data$Site.Standard.Group <- do.call(paste,
                                         c(shrub.data[c("Site", "Standard.Group")],
                                           sep = ".")) 

write.csv(shrub.data, file = "shrubdata.csv")

## Site ground cover data

groundcover <- anon.data(
    input.table = "../../DataRepository/VegetationData/GroundCover/groundcover_area2.csv",
    lookup.table = "../../DataAnalysis/site_lookup.csv",
    site.column = "Site",
    match.column = "site.name",
    anon.column = "site.abbr",
    write.to.disk = FALSE
)

groundcover$SiteSamTyp <- do.call(paste,
                                  c(groundcover[c("Site", "Sample.Typ")],
                                    sep = ".")) 

write.csv(groundcover, file = "groundcover_area2.csv")



