setwd('C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt_utilities')
library(here)
library(tinyplot)
library(data.table)

rm(list=ls())

i_am('R/test_count_models.R')

# Source the new count-based models
source(here('R','mortality_models_counts.R'))
source(here('R','lt_functions.R'))

# =========================================================================
# TEST 1: BASIC-TESTING WITH HMD DATA
# =========================================================================

# Abridged data for test
exposures_abridged_data <- fread(here('data','exposures_abridged_data.txt'))
exposures_abridged_data[,c('female','male'):=NULL]
setnames(exposures_abridged_data,c('total'),c('Ex'))

deaths_abridged_data <- fread(here('data','deaths_abridged_data.txt'))
deaths_abridged_data[,c('female','male'):=NULL]
setnames(deaths_abridged_data,c('total'),c('Dx'))

lt_abridged_data <- merge(exposures_abridged_data,deaths_abridged_data,by=c('year','age'))

# Some cleaning
# 1. get start age for each interval
lt_abridged_data[, x_start := as.integer(stringr::str_extract(age,'\\d+'))]
# 2. get mid age for each interval (assume its 0.14 for 0-1 and 111.5 for 110+)
lt_abridged_data[, x_mid := x_start + 2.5]
lt_abridged_data[x_start == 0, x_mid := 0.14]
lt_abridged_data[x_start == 1, x_mid := 2.5]
lt_abridged_data[x_start == 110, x_mid := 111.5]

# 3. keep age, deaths, and exposures for the new model
# (Assuming your updated text file contains Dx and Ex)
lt_abridged_data <- lt_abridged_data[, .(x = x_mid, Dx, Ex)]

# 4. We still need the empirical mx for plotting and classical life tables
lt_abridged_data[, mx := Dx / Ex]
setorder(lt_abridged_data,x)

# Single-year standard (The standard curve still expects rates, so this remains unchanged)
lt_single_standard <- fread(here('data','lt_single_standard.txt'))

# Some cleaning
lt_single_standard[, x_start := as.integer(stringr::str_extract(age,'\\d+'))]
lt_single_standard[, x_mid := x_start + 0.5]
lt_single_standard[x_start == 0, x_mid := 0.14]
lt_single_standard[x_start == 110, x_mid := 111.5]
lt_single_standard <- lt_single_standard[, .(x = x_mid, mx_std = mx)]

# Singularize now takes Dx and Ex, and fits the Negative Binomial model.
# It still outputs mx_s (predicted continuous rates).
lt_singularized <- singularize(lt_abridged_data, lt_single_standard)

# Combine data for plot
lt_singularized[, type := 'Singularized']
lt_abridged_data[, type := 'Abridged']
lt_single_standard[, type := 'Standard']

lt_data <- rbindlist(
  list(
    lt_abridged_data[, .(type, x, mx)],
    lt_singularized[, .(type, x, mx = mx_s)],
    lt_single_standard[, .(type, x, mx = mx_std)]
  )
)

setorder(lt_data,type,x)

# Compare age-specific mortality rates
tinyplot(log(mx) ~ x | type, 
         data = lt_data, 
         type = "l", 
         lwd = 2,
         xlab = "Age (x)", 
         ylab = "Log(mx)",
         legend = list(title = "Data Type", bty = "n"),
         grid = TRUE,
         frame = FALSE)

# Compute abridged life table with the original empirical rates
my_lt_abridged <- lt_abridged(nMx = lt_abridged_data$mx)
setDT(my_lt_abridged)
my_lt_abridged[, type := 'Abridged']

# Compute single-year life table from singularized data
my_lt_singularized <- lt_single(Mx = lt_singularized$mx_s)
setDT(my_lt_singularized)
my_lt_singularized[, type := 'Singularized']

# Load the corresponding single-year data and compute single-year life table
lt_single_data <- fread(here('data','lt_single_data.txt'))
my_lt_single <- lt_single(Mx = lt_single_data$mx)
setDT(my_lt_single)
my_lt_single[, type := 'Single']

# Combine data for plot
comparison_data <- rbindlist(
  list(
    my_lt_abridged[, .(type, x, lx)],
    my_lt_single[, .(type, x, lx)],
    my_lt_singularized[, .(type, x, lx)]
  )
)

# Compare survival function
tinyplot(lx ~ x | type, 
         data = comparison_data, 
         type = "l", 
         lwd = 2,
         xlab = "Age (x)", 
         ylab = "Survival function (lx)",
         legend = list(title = "Data Type", bty = "n"),
         grid = TRUE,      
         frame = FALSE)    

# =========================================================================
# TEST 2: STRESS-TESTING WITH PROBLEMATIC DATA
# =========================================================================

# 1. Load the problematic abridged data from the CSV
dt_prob <- fread(here('outputs','data','problematic_abridged_mx.csv'))

# Set mx to zero for age 10 to test if it the code breaks
dt_prob[age==10,mx:=0]

# 2. Clean to match the expected format for 'singularize'
dt_prob[, x_mid := age + 2.5]
dt_prob[age == 0, x_mid := 0.14]
dt_prob[age == 1, x_mid := 2.5]
dt_prob[age == 90, x_mid := 95.5] 

# Add exposures
fake_exposures <- copy(lt_abridged_data)[,x:=fifelse(x>=90,95.5,x)][,.(Ex=sum(Ex)),by=.(x)]
dt_prob <- merge(dt_prob,fake_exposures[,.(x_mid=x,Ex)],by=c('x_mid'),all.x=T)
dt_prob[,Dx:=floor(mx*Ex)]

# Assume the problematic dataset now has Dx and Ex columns
prob_abridged <- dt_prob[, .(x = x_mid, Dx, Ex)]

# Derive empirical mx for the life table functions
prob_abridged[, mx := Dx / Ex]

# 3. Run the model on the problematic data
prob_singularized <- singularize(
  prob_abridged,
  lt_single_standard,
  dof = 15,
  skip_open = TRUE,
  skip_ages = list(c(0, 5),c(80,100)), # Now takes more than one interval if needed
  smoooth_old_ages = TRUE # This adds a further Kannisto smoothing, with might
                          # be a good idea if the results at the older age don't
                          # look plausible.
)

# 4. Combine data for Mortality Rates plot
prob_singularized[, type := 'Singularized']
prob_abridged[, type := 'Abridged Problematic']

prob_data <- rbindlist(
  list(
    prob_abridged[, .(type, x, mx)],
    prob_singularized[, .(type, x, mx = mx_s)],
    lt_single_standard[, .(type, x, mx = mx_std)]
  )
)

# Plot 1: Compare log(mx) for problematic data
tinyplot(log(mx) ~ x | type, 
         data = prob_data, 
         type = "l", 
         lwd = 2,
         xlab = "Age (x)", 
         ylab = "Log(mx)",
         main = "Mortality Rates: Problematic vs. Singularized",
         legend = list(title = "Data Type", bty = "n"),
         grid = TRUE,
         frame = FALSE)

# 5. Compute Life Tables
my_prob_lt_abridged <- lt_abridged(nMx = prob_abridged$mx)
setDT(my_prob_lt_abridged)
my_prob_lt_abridged[, type := 'Abridged Problematic']

my_prob_lt_singularized <- lt_single(Mx = prob_singularized$mx_s)
setDT(my_prob_lt_singularized)
my_prob_lt_singularized[, type := 'Singularized']

# 6. Combine data for Survival plot 
prob_comparison_data <- rbindlist(
  list(
    my_prob_lt_abridged[, .(type, x, lx)],
    my_prob_lt_singularized[, .(type, x, lx)]
  )
)

# Plot 2: Compare Survival Function (lx)
tinyplot(lx ~ x | type, 
         data = prob_comparison_data, 
         type = "l", 
         lwd = 2,
         xlab = "Age (x)", 
         ylab = "Survival function (lx)",
         main = "Survival: Problematic vs. Singularized",
         legend = list(title = "Data Type", bty = "n"),
         grid = TRUE,      
         frame = FALSE)
