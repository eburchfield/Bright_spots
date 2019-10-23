run_model <- function(formula, data, name, save=T) {
  out <- inla(formula, data=data, family="gaussian", 
              control.predictor = list(compute=T), 
              control.compute = list(dic=T, cpo=T))
  if (save==T) {
    saveRDS(out, paste0("./out/", name, ".RDS")) 
  } 
  return(out)
}

model_diagnostics <- function(model_run, data) {
  # add plots and additional checks in model.results()
  # Gelman R2 in cross.val function
  
  if (is.character(model_run)) {
    r <- readRDS(model_run)
  } else (r <- model_run)
  
  fdf <- cbind(data, r$summary.fitted.values$mean)  # posterior distribution mean 
  colnames(fdf)[dim(fdf)[2]] <- c("FittedVals")
  
  # Deviance information criterion, useful for cross-model comparison (p. 170)
  # Models with smaller DIC are better supported by the data
  DIC <- r$dic$dic
  
  # Cross-validation (5.6.1, p. 166)
  
  # Conditional predictive ordinate, LOOCV
  # Compare copetitive models in terms of predictive performance, with larger values denoting better fit
  CPO <- sum(log(r$cpo$cpo), na.rm=T)
  # Probability Integral Transform, uniform distribution means the predictive distribution is coherent with the data
  PIT <- ggplot() +
    geom_histogram(aes(x = r$cpo$pit)) +
    ggtitle("PIT")
  
  # Posterior predictive checks (5.6.1, p. 168)
  
  # Posterior predictive distribution
  ppv <- c()
  for (i in 1:nrow(fdf)) {
    ppv[i] <- inla.pmarginal(q = fdf$YIELD[i], 
                             marginal = r$marginals.fitted.values[[i]])
  }
  # Linear indicates close fit between observed and predicted
  PPC_PT <- ggplot(fdf) +
    geom_point(aes(x = YIELD, y = fdf$FittedVals), alpha=0.2) +
    xlab("Observed") + ylab("Mean posterior predictive distribution") +
    theme_classic()
  # Reasonable fit should have few high P-values
  PPC_HIST <- ggplot() +
    geom_histogram(aes(x = ppv)) +
    ggtitle("Posterior predictive p-values")
  
  # MSE <- 1/(length(fdf$YIELD)*sum((fdf$YIELD - fdf$FittedVals)^2, na.rm=T))
  MSE <- MSE(fdf$FittedVals, fdf$YIELD)
  pred_res2<-(fdf$FittedVals[!is.na(fdf$YIELD)] - mean(fdf$YIELD, na.rm=T))^2
  obs_res2<-(fdf$YIELD[!is.na(fdf$YIELD)] - mean(fdf$YIELD, na.rm=T))^2
  R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
  
  # R2 <- function(x,y) cor(x,y) ^ 2
  # 
  # rss <- sum((fdf$FittedVals - fdf$YIELD) ^ 2)
  # tss <- sum((fdf$YIELD - mean(fdf$YIELD))^2)
  # rsq <- 1 - (rss/tss)
  
  fit <- list(DIC, CPO, MSE, R2, PIT, PPC_PT, PPC_HIST)
  names(fit) <- c("DIC", "CPO", "MSE", "R2", "PIT", "PPC_PT", "PPC_HIST")
  return(fit)
}

spatial_effects <- function(model_run, data, level) {
  # assumes set-up code has been run loading cty shapefile
  
  r <- model_run
  n <- length(unique(data[,level]))
  # extract area-specific residuals
  asr <- r$summary.random[[level]][1:n, c(1, 2, 3)]  # ID, mean, sd
  if (level == "LRR") {
    asr <- r$summary.random[[level]][1:(n + 2), c(1, 2, 3)]  # ID, mean, sd
    
    colnames(asr) <- c("LRR", "mean", "SD")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, LRR)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID", all=T)
  } 
  
  if (level == "ECO") {
    colnames(asr) <- c("ECO", "mean", "SD")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, ECO)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID")
  }
  map_asr <- sf::st_as_sf(merge(county_sub, asr, by = level), all=T)
  # include indicator of significance (not just mean)
  
  require(viridis)
  
  se <- ggplot(map_asr) +
    geom_sf(color = "transparent", size = 0.05, aes(fill = mean)) +
    theme_minimal() +
    scale_fill_viridis(option = "magma")
  
  return(se)
  
}

spatial_residuals <- function(model_run, data, level) {
  # assumes set-up code has been run loading cty shapefile
  
  r <- model_run
  n <- length(unique(data[,level]))
  # extract area-specific residuals
  asr <- r$summary.random[[level]][(n+1):(n*2), c(1, 2, 3)]  # ID, mean, sd
  asr$ID <- 1:n
  if (level == "STATE") {
    colnames(asr) <- c("STATE", "mean", "SD")
  } 
  
  if (level == "ECO") {
    colnames(asr) <- c("ECO", "mean", "SD")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, ECO)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID", all=T)
  }
  map_asr <- sf::st_as_sf(merge(county_sub, asr, by = level))
  # include indicator of significance (not just mean)
  
  require(viridis)
  
  se <- ggplot(map_asr) +
    geom_sf(color = "transparent", size = 0.05, aes(fill = mean)) +
    theme_minimal() +
    scale_fill_viridis(option = "magma")
  
  return(se)
  
}

nonlinear_effect <- function(model_run, variable) {
  r <- model_run
  nle <- ggplot(data = r$summary.random[[variable]][,c(1,4:6)], 
                aes(x = ID, y = `0.5quant`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
    #geom_line(aes(x = ID, y = `0.025quant`) , col="#0d98ba") +
    #geom_line(aes(x=ID, y=`0.975quant`), col="#0d98ba") + 
    theme_minimal() +
    theme(legend.position="none") +
    xlab(variable) +
    ylab("Effect on yield (bu/ac)")
  return(nle)
}

generate_DIC <- function(fm, mod_id) {
  modr <- run_model(fm, data=null, name = mod_id, save=F)
  modd <- model_diagnostics(modr, null)
  return(list(modd$DIC, modd$R2))
}

find_bs <- function(model_run, data, level, th) {
  
  # assumes set-up code has been run loading cty shapefile
  
  r <- model_run
  b0 <- r$summary.fixed[1,1]
  
  # extract area-specific residuals
  n <- length(unique(data[,"ID"]))
  asr <- r$summary.random[["ID"]][1:n, c(1, 2)]  # E_ij = u_0jk + v_0jk
  asr$mean <- asr$mean # mean spatial effects, county
  colnames(asr) <- c("ID", "COUNTY_MEAN")
  
  # extract area-specific residuals
  n2 <- length(unique(data[,level]))
  asr2 <- r$summary.random[[level]][ , c(1, 2)]  # STATE_FP, mean
  asr2$mean <- asr2$mean # v_00k
  
  if (level == "LRR") {
    colnames(asr2) <- c("LRR", "REGION_MEAN")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, LRR)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID", all=T)
    
    m1 <- merge(county_sub, asr, by = "ID", all=T)
    m2 <- merge(m1, asr2, by = level, all=T)
    m2$DIFF <- m2$REGION_MEAN - m2$COUNTY_MEAN
    m2$BS_CAT <- NA
    
    lrr_ids <- unique(m2@data[,level])
    
    for (i in 1:length(lrr_ids)) {
      id <- lrr_ids[[i]]
      df <- m2@data %>% filter(LRR == id)
      diff <- df$REGION_MEAN - df$COUNTY_MEAN
      sdh <- mean(diff) + (th*sd(diff))
      sdl <- mean(diff) - (th*sd(diff)) 
      
      bs <- ifelse(diff > sdh, "Dark Spot", NA)
      bs <- ifelse(diff < sdl, "Bright Spot", bs)
      bs[is.na(bs)] <- "Average"
      m2$BS_CAT[m2$LRR == id] <- bs
      
    }
    
  } 
  
  if (level == "ECO") {
    colnames(asr2) <- c("ECO", "REGION_MEAN")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, ECO)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID")
    
    m1 <- merge(county_sub, asr, by = "ID", all=T)
    m2 <- merge(m1, asr2, by = level, all=T)
    m2$DIFF <- m2$REGION_MEAN - m2$COUNTY_MEAN
    m2$BS_CAT <- NA
    
    lrr_ids <- unique(m2@data[,level])
    
    for (i in 1:length(lrr_ids)) {
      id <- lrr_ids[[i]]
      df <- m2@data %>% filter(ECO == id)
      diff <- df$REGION_MEAN - df$COUNTY_MEAN
      sdh <- mean(diff) + (th*sd(diff))
      sdl <- mean(diff) - (th*sd(diff)) 
      
      bs <- ifelse(diff > sdh, "Dark Spot", NA)
      bs <- ifelse(diff < sdl, "Bright Spot", bs)
      bs[is.na(bs)] <- "Average"
      m2$BS_CAT[m2$ECO == id] <- bs
      
    }
  }
  
  if (level == "STATE") {
    colnames(asr2) <- c("STATE", "REGION_MEAN")
    ne <- null %>% group_by(GEOID) %>%
      filter(row_number() == 1) %>%
      select(GEOID, STATE)
    county_sub@data <- merge(county_sub@data, ne, by = "GEOID")
    county_sub$STATE <- county_sub$STATE.x
    
    m1 <- merge(county_sub, asr, by = "ID", all=T)
    m2 <- merge(m1, asr2, by = level, all=T)
    m2$DIFF <- m2$REGION_MEAN - m2$COUNTY_MEAN
    m2$BS_CAT <- NA
    
    lrr_ids <- unique(m2@data[,level])
    
    for (i in 1:length(lrr_ids)) {
      id <- lrr_ids[[i]]
      df <- m2@data %>% filter(STATE == id)
      diff <- df$REGION_MEAN - df$COUNTY_MEAN
      sdh <- mean(diff) + (th*sd(diff))
      sdl <- mean(diff) - (th*sd(diff)) 
      
      bs <- ifelse(diff > sdh, "Dark Spot", NA)
      bs <- ifelse(diff < sdl, "Bright Spot", bs)
      bs[is.na(bs)] <- "Average"
      m2$BS_CAT[m2$STATE == id] <- bs
      
    }
  }
  
  m2$BS_CAT <- as.factor(m2$BS_CAT)
  return(m2)
  
}

plot_bsds <- function(spdf) {
  spdf <- st_as_sf(spdf)
  # spdf$BS_CAT[spdf$BS_CAT == "Average"] <- NA
  se <- ggplot() +
    geom_sf(data = spdf, color = "transparent", size = 0.05, aes(fill = BS_CAT)) +
    geom_sf(data = lrr_shp, color = "black", fill = NA, size=.5) +
    theme_minimal() +
    #theme(legend.position = "none") +
    #scale_fill_discrete(na.value = "#f2f2f2") +
    #scale_fill_manual(values = c("#4C9A2A", "#000080"), na.value ="#e5e5e5", name = "")
    scale_fill_manual(values = c("#e5e5e5", "#4C9A2A", "#000080"), name = "")
  se
  
}




# plot_bs <- function(spdf) {
#   spdf$BS <- ifelse(spdf$BS_NUM > 0, 1, NA)
#   map_asr <- sf::st_as_sf(spdf)
#   se <- ggplot() +
#     geom_sf(data = spdf, color = "transparent", size = 0.05, aes(fill = BS)) +
#     geom_sf(data = lrr_shp, color = "black", fill = NA, size=.5) +
#     theme_minimal() +
#     theme(legend.position = "none") +
#     scale_fill_discrete(na.value = "#f2f2f2")
#   return(se)
# }

# plot_ds <- function(spdf, th) {
#   spdf$BS <- ifelse(spdf$ZSCORE < -th, spdf$ZSCORE, NA)
#   map_asr <- sf::st_as_sf(spdf)
#   require(viridis)
#   se <- ggplot(map_asr) +
#     geom_sf(color = "transparent", size = 0.05, aes(fill = BS)) +
#     theme_minimal() +
#     scale_fill_viridis(option = "magma") +
#     theme(legend.position = "none")
#   return(se)
# }

national_bs <- function(model_run, data, th) {
  r <- model_run
  b0 <- r$summary.fixed[1,1]
  n <- length(unique(data[,"ID"]))
  # extract area-specific residuals
  asr <- r$summary.random[["ID"]][1:n, c(1, 2)]  # ID, mean
  asr$mean <- asr$mean # mean spatial effects, county
  colnames(asr) <- c("ID", "COUNTY_MEAN")
  
  m1 <- merge(county_sub, asr, by = "ID", all=T)
  
  diff <- m1$COUNTY_MEAN - b0
  sdh <- mean(diff) + (th*sd(diff))
  sdl <- mean(diff) - (th*sd(diff)) 
  
  bs <- ifelse(diff > sdh, "Bright Spot", NA)
  bs <- ifelse(diff < sdl, "Dark Spot", bs)
  bs[is.na(bs)] <- "Average"
  m1$BS_CAT <- bs
  m1$DIFF <- diff
  
  return(m1)
}  

viz_difference <- function(voi, ylab) {
  ggplot(data = nulli) +
    geom_boxplot(aes(x = BS, y = voi, fill = BS)) +
    theme_minimal() +
    ylab(ylab) +
    xlab("")
}