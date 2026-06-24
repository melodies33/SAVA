# Main script for running all algorithms under real data analysis in Section 4.

# The script generates CSV files that record the performance of the competing methods. Once the CSV files are generated, execute the corresponding code blocks in realdata_plot.R to reproduce the figures reported in the paper.

library(truncnorm)
library(tidyverse)
library(Rcpp)
library(onlineFDR)
library(SAVA)
library(foreach)
library(doParallel)
library(VGAM)
source('functions_sava.R')
library(lubridate)
library(data.table)

###### Functions for algorithms #####

# Function: amazon_sava
#
# Implement the SAVA algorithm in Section 4.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   bound            : size of support. 
#   q                : target FSR level.
#   k                : tuning parameter k.
#   w0               : initial alpha wealth.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   FSP     : false selection proportion sequence at each decision time.
#   realsig : cumulated numbers of tasks at each decision time.
#   posi    : number of true selection for items whose ratings are greater than 3.
#   nega    : number of true selection for items whose ratings are less than 3.
# 
# Used in:
#   SAVA algorithm in Section 4
#
# See Section 4 for details.

amazon_sava <- function(item_metrics, rating_sequences, thetas, bound = 4, q = 0.1, k = 25, w0 = q, itemnumbner = itemnumbner, decision_times = decision_times) {
  
  
  deltaesti <- rep(0, itemnumbner)  
  stoptime <- rep(Inf, itemnumbner) 
  thresholdvec <- rep(0, itemnumbner)
  thresholdvec[1:min(k, itemnumbner)] <- w0/k
  
  selectlarger1 <- 0  
  num_decisions <- length(decision_times)
  selectposi <- integer(0) 
  selectnega <- integer(0) 
  falsetotal <- numeric(num_decisions)
  decidetotal <- numeric(num_decisions)
  realtotal <- numeric(num_decisions)
  positotal <- numeric(num_decisions)
  negatotal <- numeric(num_decisions)
  pb <- txtProgressBar(min = 0, max = num_decisions, style = 3)
  for (ti in 1:num_decisions) {
    
    current_time <- decision_times[ti]
    active_tasks <- which(item_metrics$arrival_time <= current_time & stoptime > current_time)
    if (length(active_tasks) > 0) {
      for (task_idx in active_tasks) {
        if (deltaesti[task_idx] != 0) next 
        current_item <- item_metrics$item_id[task_idx]
        current_ratings <- rating_sequences %>%
          filter(item_id == current_item, review_time <= current_time) %>%
          pull(rating)
        
        if (length(current_ratings) == 0) next 
        threshold <- thresholdvec[task_idx]
        
        cs.below <- falphaplus(threshold, current_ratings, bound, length(current_ratings), q)
        cs.above <- falphaminus(threshold, current_ratings, bound, length(current_ratings), q)
        
        deltaplus <- cs.below >= 0  
        deltaminus <- cs.above < 0 
        deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
        
        if (deltatotal != 0) {
          deltaesti[task_idx] <- deltatotal
          stoptime[task_idx] <- current_time
          selectlarger1 <- selectlarger1 + 1
          if (deltatotal > 0) {
            selectposi <- c(selectposi, task_idx)
          } else {
            selectnega <- c(selectnega, task_idx)
          }
          endindex <- min(task_idx + k, itemnumbner)
          if (endindex > task_idx) {
            if (selectlarger1 == 1) {
              thresholdvec[(task_idx + 1):endindex] <- thresholdvec[(task_idx + 1):endindex] + (q - w0) / k
            } else {
              thresholdvec[(task_idx + 1):endindex] <- thresholdvec[(task_idx + 1):endindex] + q / k
            }
          }
        }
      }
    }
    
    selected_tasks <- which(stoptime <= current_time & deltaesti != 0)
    decidetotal[ti] <- length(selected_tasks)
    real_tasks <- which(thetas != 0 & item_metrics$arrival_time <= current_time)
    realtotal[ti] <- length(real_tasks)
    
    positotal[ti] <- length(selectposi)
    negatotal[ti] <- length(selectnega)
    if (decidetotal[ti] > 0) {
      false_selections <- sum(deltaesti[selected_tasks] != thetas[selected_tasks])
      falsetotal[ti] <- false_selections
    } else {
      if (ti > 1) {
        falsetotal[ti] <- falsetotal[ti-1]
        realtotal[ti] <- realtotal[ti-1]
      }
    }
    setTxtProgressBar(pb, ti)
  }
  close(pb)
  fsp <- falsetotal / pmax(decidetotal, 1)
  
  return(list(
    FSP = fsp,
    realsig = realtotal,
    posi = positotal,
    nega = negatotal
  ))
}

# Function: amazon_lordpp
#
# Implement the LORD++ algorithm in Section 4.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   FSP     : false selection proportion sequence at each decision time.
#   realsig : cumulated numbers of tasks at each decision time.
#   posi    : number of true selection for items whose ratings are greater than 3.
#   nega    : number of true selection for items whose ratings are less than 3.
# 
# Used in:
#   LORD++ algorithm in Section 4
#
# See Section 4 for details.

amazon_lordpp <- function(item_metrics, rating_sequences, thetas, q = 0.1, itemnumbner = itemnumbner, decision_times = decision_times) {
  selectposi <- integer(0) 
  selectnega <- integer(0)  
  deltaesti <- rep(0, itemnumbner)  
  falsetotal <- numeric(itemnumbner)
  decidetotal <- numeric(itemnumbner)
  positotal <- numeric(itemnumbner)
  negatotal <- numeric(itemnumbner)
  realtotal <- numeric(itemnumbner)
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
    
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
    if (length(current_ratings) >= 3) { 
      p.posi <- wilcox.test(current_ratings, mu = 0, alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, alternative = 'less', exact = FALSE)$p.value
    } else {
      p.posi <- 1  
      p.nega <- 1
    }
    thresholdposi <- threldpp(q, q/10, i, selectposi)
    thresholdnega <- threldpp(q, q/10, i, selectnega)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
      } else {
        selectnega <- c(selectnega, i)
      }
    }
    selected_tasks <- which(deltaesti[1:i] != 0)
    decidetotal[i] <- length(selected_tasks)
    real_tasks <- which(thetas[1:i] != 0)
    realtotal[i] <- length(real_tasks)
    positotal[i] = length(selectposi)
    negatotal[i] = length(selectnega)
    if (decidetotal[i] > 0) {
      false_selections <- sum(deltaesti[selected_tasks] != thetas[selected_tasks])
      falsetotal[i] <- false_selections
    } else {
      if (i > 1) {
        falsetotal[i] <- falsetotal[i-1]
        realtotal[i] <- realtotal[i-1]
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  fsp <- falsetotal / pmax(decidetotal, 1)
  return(list(
    FSP = fsp,
    realsig = realtotal,
    posi = positotal,
    nega = negatotal
  ))
}

# Function: amazon_saffron
#
# Implement the SAFFRON algorithm in Section 4.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.
#   Clamb            : parameter lambda in SAFFRON.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   FSP     : false selection proportion sequence at each decision time.
#   realsig : cumulated numbers of tasks at each decision time.
#   posi    : number of true selection for items whose ratings are greater than 3.
#   nega    : number of true selection for items whose ratings are less than 3.
# 
# Used in:
#   SAFFRON algorithm in Section 4
#
# See Section 4 for details.

amazon_saffron <- function(item_metrics, rating_sequences, thetas, q = 0.1, Clamb = 0.5, itemnumbner = itemnumbner, decision_times = decision_times) {
  selectposi <- integer(0)  
  selectnega <- integer(0)  
  cplus <- rep(0, itemnumbner)  
  cminus <- rep(0, itemnumbner) 
  deltaesti <- rep(0, itemnumbner) 
  falsetotal <- numeric(itemnumbner)
  decidetotal <- numeric(itemnumbner)
  realtotal <- numeric(itemnumbner)
  positotal <- numeric(itemnumbner)
  negatotal <- numeric(itemnumbner)
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
   
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
    if (length(current_ratings) >= 3) {
      p.posi <- wilcox.test(current_ratings, mu = 0, alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, alternative = 'less', exact = FALSE)$p.value
    } else {
      p.posi <- 1 
      p.nega <- 1
    }
    if (p.posi <= Clamb) {  
      cplus[i] <- 1
    }
    if (p.nega <= Clamb) {
      cminus[i] <- 1
    }
    thresholdposi <- thresaff(q, q/2, i, selectposi, cplus, Clamb)
    thresholdnega <- thresaff(q, q/2, i, selectnega, cminus, Clamb)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
      } else {
        selectnega <- c(selectnega, i)
      }
    }
    selected_tasks <- which(deltaesti[1:i] != 0)
    decidetotal[i] <- length(selected_tasks)
    real_tasks <- which(thetas[1:i] != 0)
    realtotal[i] <- length(real_tasks)
    positotal[i] = length(selectposi)
    negatotal[i] = length(selectnega)
    if (decidetotal[i] > 0) {
      false_selections <- sum(deltaesti[selected_tasks] != thetas[selected_tasks])
      falsetotal[i] <- false_selections
    } else {
      if (i > 1) {
        falsetotal[i] <- falsetotal[i-1]
        realtotal[i] <- realtotal[i-1]
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  fsp <- falsetotal / pmax(decidetotal, 1)
  return(list(
    FSP = fsp,
    realsig = realtotal,
    posi = positotal,
    nega = negatotal
  ))
}

# Function: amazon_addis
#
# Implement the ADDIS algorithm in Section 4.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.           
#   lamb             : parameter lambda in ADDIS.
#   tau              : parameter tau in ADDIS. 
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   FSP     : false selection proportion sequence at each decision time.
#   realsig : cumulated numbers of tasks at each decision time.
#   posi    : number of true selection for items whose ratings are greater than 3.
#   nega    : number of true selection for items whose ratings are less than 3.
# 
# Used in:
#   ADDIS algorithm in Section 4
#
# See Section 4 for details.

amazon_addis <- function(item_metrics, rating_sequences, thetas, q = 0.1, lamb = 0.25, tau = 0.5, itemnumbner = itemnumbner, decision_times = decision_times) {
  selectposi <- integer(0)  
  selectnega <- integer(0) 
  kstarposi <- integer(0)   
  kstarnega <- integer(0)  
  splus <- rep(0, itemnumbner) 
  sminus <- rep(0, itemnumbner)
  cplus <- rep(0, itemnumbner)  
  cminus <- rep(0, itemnumbner)
  deltaesti <- rep(0, itemnumbner) 
  falsetotal <- numeric(itemnumbner)
  decidetotal <- numeric(itemnumbner)
  realtotal <- numeric(itemnumbner)
  positotal = numeric(itemnumbner)
  negatotal = numeric(itemnumbner)
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
    if (length(current_ratings) >= 3) {
      p.posi <- wilcox.test(current_ratings, mu = 0, alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, alternative = 'less', exact = FALSE)$p.value
    } else {
      p.posi <- 1  
      p.nega <- 1
    }
    if (p.posi <= tau) {
      splus[i] <- 1
    }
    if (p.nega <= tau) {
      sminus[i] <- 1
    }
    
    if (p.posi <= lamb) {
      cplus[i] <- 1
    }
    if (p.nega <= lamb) {
      cminus[i] <- 1
    }
    sindexposi <- sum(splus[1:i])
    sindexnega <- sum(sminus[1:i])
    thresholdposi <- threaddis(q, q/2, i, selectposi, cplus, sindexposi, kstarposi, lamb = lamb, tau = tau)
    thresholdnega <- threaddis(q, q/2, i, selectnega, cminus, sindexnega, kstarnega, lamb = lamb, tau = tau)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
        kstarposi <- c(kstarposi, sum(splus[1:i]))
      } else {
        selectnega <- c(selectnega, i)
        kstarnega <- c(kstarnega, sum(sminus[1:i]))
      }
    }
    selected_tasks <- which(deltaesti[1:i] != 0)
    decidetotal[i] <- length(selected_tasks)
    real_tasks <- which(thetas[1:i] != 0)
    realtotal[i] <- length(real_tasks)
    positotal[i] = length(selectposi)
    negatotal[i] = length(selectnega)
    if (decidetotal[i] > 0) {
      false_selections <- sum(deltaesti[selected_tasks] != thetas[selected_tasks])
      falsetotal[i] <- false_selections
    } else {
      if (i > 1) {
        falsetotal[i] <- falsetotal[i-1]
        realtotal[i] <- realtotal[i-1]
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  fsp <- falsetotal / pmax(decidetotal, 1)
  return(list(
    FSP = fsp,
    realsig = realtotal,
    posi = positotal,
    nega = negatotal
  ))
}



###### pre-process data from Amazon Fashion #####
amazon_data <- fread("real-amazon25/AMAZON_FASHION.csv",
                     col.names = c("item_id", "user_id", "rating", "timestamp"))

processed_data <- amazon_data %>%
  filter(rating >= 1 & rating <= 5) %>%
  arrange(item_id, timestamp)

item_metrics <- processed_data %>%
  group_by(item_id) %>%
  summarise(
    first_review_time = min(timestamp),
    total_reviews = n(),
    avg_rating = mean(rating)
  ) %>%
  arrange(first_review_time) %>%
  filter(total_reviews > 50)  

min_time <- min(item_metrics$first_review_time)
item_metrics <- item_metrics %>%
  mutate(arrival_day = as.integer((first_review_time - min_time) / 86400)) %>%
  arrange(arrival_day, item_id)

item_metrics <- item_metrics %>%
  group_by(arrival_day) %>%
  slice(1) %>% 
  ungroup() %>%
  mutate(arrival_time = arrival_day) %>%  
  select(-arrival_day, -first_review_time)

thetas <- numeric(nrow(item_metrics))
thetas[item_metrics$avg_rating >= 3] <- 1  
thetas[item_metrics$avg_rating < 3] <- -1  
arrival_map <- item_metrics %>% select(item_id, arrival_time)

rating_sequences <- processed_data %>%
  inner_join(arrival_map, by = "item_id") %>%
  arrange(item_id, timestamp) %>%
  group_by(item_id) %>%
  mutate(review_seq = row_number(),
         review_time = as.integer((timestamp - min_time) / 86400)) %>%
  ungroup()

rating_sequences$rating = rating_sequences$rating - 3

item_metrics <- item_metrics %>% arrange(arrival_time)
itemnumbner <- nrow(item_metrics)

time_max = max(rating_sequences$review_time)
decision_times <- c(item_metrics$arrival_time[-1],time_max)
num_decision = length(decision_times)

#### Compare methods under data from Amazon Fashion (for figure 5) ####
q = 0.2  # target FSR level
startpeek = floor(num_decision/10)
time.peek = seq(startpeek, num_decision, startpeek) # time points for recording FSP, TSP

# SAVA
re = amazon_sava(item_metrics, rating_sequences, thetas, q, k = 100, w0 = q, bound = 4, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.sava = re$FSP[time.peek]
posi.sava = re$posi[time.peek]
nega.sava = re$nega[time.peek]

# LORD++
re = amazon_lordpp(item_metrics, rating_sequences, thetas, q, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.lordpp = re$FSP[time.peek]
posi.lordpp = re$posi[time.peek]
nega.lordpp = re$nega[time.peek]

# SAFFRON
re = amazon_saffron(item_metrics, rating_sequences, thetas, q, Clamb = 0.5,   itemnumbner = itemnumbner, decision_times = decision_times)
FSR.saffron = re$FSP[time.peek]
posi.saffron = re$posi[time.peek]
nega.saffron = re$nega[time.peek]

# ADDIS
re = amazon_addis(item_metrics, rating_sequences, thetas, q, lamb = 0.25, tau = 0.5, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.addis = re$FSP[time.peek]
posi.addis = re$posi[time.peek]
nega.addis = re$nega[time.peek]


select.sava = posi.sava + nega.sava
select.lordpp = posi.lordpp + nega.lordpp
select.saffron = posi.saffron + nega.saffron
select.addis = posi.addis + nega.addis

# summary of performances 
re.dat = tibble(
  decision.time = rep(time.peek, 4),
  FSP = c(FSR.sava, FSR.lordpp, FSR.saffron, FSR.addis), 
                Positive = c(posi.sava, posi.lordpp, posi.saffron, posi.addis), 
                Negative = c(nega.sava, nega.lordpp, nega.saffron, nega.addis), 
                Totalselection = c(select.sava, select.lordpp, select.saffron, select.addis),
                Method = rep(c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS'), each = length(time.peek)))
write.csv(re.dat, 'amazon-fashion-comprison.csv')

#### pre-process data from All beauty data ####
amazon_data <- fread("real-amazon25/AMAZON_All_Beauty.csv",
                     col.names = c("item_id", "user_id", "rating", "timestamp"))

processed_data <- amazon_data %>%
filter(rating >= 1 & rating <= 5) %>%
arrange(item_id, timestamp)
item_metrics <- processed_data %>%
group_by(item_id) %>%
summarise(
  first_review_time = min(timestamp),
  total_reviews = n(),
  avg_rating = mean(rating)
) %>%
arrange(first_review_time) %>%
filter(total_reviews > 50) 
min_time <- min(item_metrics$first_review_time)
item_metrics <- item_metrics %>%
mutate(arrival_day = as.integer((first_review_time - min_time) / 86400)) %>%
arrange(arrival_day, item_id)
item_metrics <- item_metrics %>%
group_by(arrival_day) %>%
slice(1) %>% 
ungroup() %>%
mutate(arrival_time = arrival_day) %>%  
select(-arrival_day, -first_review_time)
thetas <- numeric(nrow(item_metrics))
thetas[item_metrics$avg_rating >= 3] <- 1  
thetas[item_metrics$avg_rating < 3] <- -1  
arrival_map <- item_metrics %>% select(item_id, arrival_time)

rating_sequences <- processed_data %>%
  inner_join(arrival_map, by = "item_id") %>%
  arrange(item_id, timestamp) %>%
  group_by(item_id) %>%
  mutate(review_seq = row_number(),
         review_time = as.integer((timestamp - min_time) / 86400)) %>%
  ungroup()

rating_sequences$rating = rating_sequences$rating - 3
item_metrics <- item_metrics %>% arrange(arrival_time)
itemnumbner <- nrow(item_metrics)
time_max = max(rating_sequences$review_time)
decision_times <- c(item_metrics$arrival_time[-1],time_max)
num_decision = length(decision_times)

#### Compare methods under All beauty data (for figure 5) ####
q = 0.2
startpeek = floor(num_decision/10)
time.peek = seq(startpeek, num_decision, startpeek) # time points for recording FSP, TSP
# SAVA
re = amazon_sava(item_metrics, rating_sequences, thetas, q, k = 100, w0 = q, bound = 4, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.sava = re$FSP[time.peek]
posi.sava = re$posi[time.peek]
nega.sava = re$nega[time.peek]

# LORD++
re = amazon_lordpp(item_metrics, rating_sequences, thetas, q, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.lordpp = re$FSP[time.peek]
posi.lordpp = re$posi[time.peek]
nega.lordpp = re$nega[time.peek]
# SAFFRON

re = amazon_saffron(item_metrics, rating_sequences, thetas, q, Clamb = 0.5,  itemnumbner = itemnumbner, decision_times = decision_times)
FSR.saffron = re$FSP[time.peek]
posi.saffron = re$posi[time.peek]
nega.saffron = re$nega[time.peek]
# ADDIS

re = amazon_addis(item_metrics, rating_sequences, thetas, q,  lamb = 0.25, tau = 0.5, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.addis = re$FSP[time.peek]
posi.addis = re$posi[time.peek]
nega.addis = re$nega[time.peek]

select.sava = posi.sava + nega.sava
select.lordpp = posi.lordpp + nega.lordpp
select.saffron = posi.saffron + nega.saffron
select.addis = posi.addis + nega.addis

# summary of performances 
re.dat = tibble(
  decision.time = rep(time.peek, 4),
  FSP = c(FSR.sava, FSR.lordpp, FSR.saffron, FSR.addis), 
                Positive = c(posi.sava, posi.lordpp, posi.saffron, posi.addis), 
                Negative = c(nega.sava, nega.lordpp, nega.saffron, nega.addis), 
                Totalselection = c(select.sava, select.lordpp, select.saffron, select.addis),
                Method = rep(c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS'), each = length(time.peek)))
write.csv(re.dat, 'amazon-allbeauty-comprison.csv')




#### pre-process data from Luxury beauty data ####
amazon_data <- fread("real-amazon25/AMAZON_Luxury_Beauty.csv",
                     col.names = c("item_id", "user_id", "rating", "timestamp"))
processed_data <- amazon_data %>%
  filter(rating >= 1 & rating <= 5) %>%
  arrange(item_id, timestamp)
item_metrics <- processed_data %>%
  group_by(item_id) %>%
  summarise(
    first_review_time = min(timestamp),
    total_reviews = n(),
    avg_rating = mean(rating)
  ) %>%
  arrange(first_review_time) %>%
  filter(total_reviews > 50) 
min_time <- min(item_metrics$first_review_time)
item_metrics <- item_metrics %>%
  mutate(arrival_day = as.integer((first_review_time - min_time) / 86400)) %>%
  arrange(arrival_day, item_id)
item_metrics <- item_metrics %>%
  group_by(arrival_day) %>%
  slice(1) %>% 
  ungroup() %>%
  mutate(arrival_time = arrival_day) %>%  
  select(-arrival_day, -first_review_time)

thetas <- numeric(nrow(item_metrics))
thetas[item_metrics$avg_rating >= 3] <- 1  
thetas[item_metrics$avg_rating < 3] <- -1  
arrival_map <- item_metrics %>% select(item_id, arrival_time)

rating_sequences <- processed_data %>%
  inner_join(arrival_map, by = "item_id") %>%
  arrange(item_id, timestamp) %>%
  group_by(item_id) %>%
  mutate(review_seq = row_number(),
         review_time = as.integer((timestamp - min_time) / 86400)) %>%
  ungroup()

rating_sequences$rating = rating_sequences$rating - 3

item_metrics <- item_metrics %>% arrange(arrival_time)
itemnumbner <- nrow(item_metrics)

time_max = max(rating_sequences$review_time)
decision_times <- c(item_metrics$arrival_time[-1],time_max)
num_decision = length(decision_times)

#### Compare methods under Luxury beauty data (for figure 5) ####
q = 0.2
startpeek = floor(num_decision/10)
time.peek = seq(startpeek, num_decision, startpeek) # time points for recording FSP, TSP
# SAVA
re = amazon_sava(item_metrics, rating_sequences, thetas, q, k = 100, w0 = q, bound = 4, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.sava = re$FSP[time.peek]
posi.sava = re$posi[time.peek]
nega.sava = re$nega[time.peek]

# LORD++
re = amazon_lordpp(item_metrics, rating_sequences, thetas, q, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.lordpp = re$FSP[time.peek]
posi.lordpp = re$posi[time.peek]
nega.lordpp = re$nega[time.peek]
# SAFFRON

re = amazon_saffron(item_metrics, rating_sequences, thetas, q, Clamb = 0.5,   itemnumbner = itemnumbner, decision_times = decision_times)
FSR.saffron = re$FSP[time.peek]
posi.saffron = re$posi[time.peek]
nega.saffron = re$nega[time.peek]
# ADDIS

re = amazon_addis(item_metrics, rating_sequences, thetas, q, lamb = 0.25, tau = 0.5, itemnumbner = itemnumbner, decision_times = decision_times)
FSR.addis = re$FSP[time.peek]
posi.addis = re$posi[time.peek]
nega.addis = re$nega[time.peek]

select.sava = posi.sava + nega.sava
select.lordpp = posi.lordpp + nega.lordpp
select.saffron = posi.saffron + nega.saffron
select.addis = posi.addis + nega.addis

# summary of performances 
re.dat = tibble(
  decision.time = rep(time.peek, 4),
  FSP = c(FSR.sava, FSR.lordpp, FSR.saffron, FSR.addis), 
                Positive = c(posi.sava, posi.lordpp, posi.saffron, posi.addis), 
                Negative = c(nega.sava, nega.lordpp, nega.saffron, nega.addis), 
                Totalselection = c(select.sava, select.lordpp, select.saffron, select.addis),
                Method = rep(c('SAVA', 'LORD++', 'SAFFRON', 'ADDIS'), each = length(time.peek)))
write.csv(re.dat, 'amazon-luxurybeauty-comprison.csv')


######## Compare test levels for different methods (for figure 6) #######


# Function: amazon_sava_testlevel
#
# Implement the SAVA algorithm in Section 4 and record the test levels at specific decision times.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   bound            : size of support. 
#   q                : target FSR level.
#   k                : tuning parameter k.
#   w0               : initial alpha wealth.
#   obs_index        : pre-specified tasks for recording test levels.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   test_level_obs     : test levels for the pre-specified tasks.
# 
# Used in:
#   Check test levels of SAVA algorithm in Section 4.
#
# See Section 4 for details.

amazon_sava_testlevel <- function(item_metrics, rating_sequences, thetas, bound = 4, q = 0.2, k = 100, w0 = q,  obs_index, itemnumbner = itemnumbner, decision_times = decision_times) {
  num_decisions <- length(decision_times)
  obs_index <- obs_index[obs_index <= num_decision]
  if (length(obs_index) == 0) {
    stop("No valid observation indices provided")
  }
  time_obs = seq(100, 900, 50)
  
  test_obs <- matrix(NA, length(obs_index), length(time_obs))
  deltaesti <- rep(0, itemnumbner)  
  stoptime <- rep(Inf, itemnumbner) 
  thresholdvec <- rep(0, itemnumbner)
  thresholdvec[1:min(k, itemnumbner)] <- w0/k
  selectlarger1 <- 0  
  pb <- txtProgressBar(min = 0, max = num_decisions, style = 3)
  for (ti in 1:num_decisions) {
    
    current_time <- decision_times[ti]
    active_tasks <- which(item_metrics$arrival_time <= current_time & stoptime > current_time)
    if (ti %in% time_obs){
      test_obs[, which(time_obs == ti)] = thresholdvec[obs_index]
    }
    if (length(active_tasks) > 0) {
      for (task_idx in active_tasks) {
        if (deltaesti[task_idx] != 0) next
        current_item <- item_metrics$item_id[task_idx]
        current_ratings <- rating_sequences %>%
          filter(item_id == current_item, review_time <= current_time) %>%
          pull('rating')
        
        if (length(current_ratings) == 0) next
        threshold <- thresholdvec[task_idx]
        cs.below <- falphaplus(threshold, current_ratings, bound, length(current_ratings), q)
        cs.above <- falphaminus(threshold, current_ratings, bound, length(current_ratings), q)
        deltaplus <- cs.below >= 0 
        deltaminus <- cs.above < 0
        deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
        
        if (deltatotal != 0) {
          deltaesti[task_idx] <- deltatotal
          stoptime[task_idx] <- current_time
          selectlarger1 <- selectlarger1 + 1
          endindex <- min(task_idx + k, itemnumbner)
          if (endindex > task_idx) {
            if (selectlarger1 == 1) {
              thresholdvec[(task_idx + 1):endindex] <- thresholdvec[(task_idx + 1):endindex] + (q - w0) / k
            } else {
              thresholdvec[(task_idx + 1):endindex] <- thresholdvec[(task_idx + 1):endindex] + q / k
            }
          }
        }
      }
    }
    
    setTxtProgressBar(pb, ti)
  }
  close(pb)
  return(list(
    test_level_obs = test_obs
  ))
}


# Function: amazon_lordpp_testlevel
#
# Implement the LORD++ algorithm in Section 4 and record the test levels at specific decision times.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.
#   obs_index        : pre-specified tasks for recording test levels.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   test_level_obs_posi : test levels for arm A from the pre-specified tasks.
#   test_level_obs_nega :test levels for arm B from the pre-specified tasks.
# 
# Used in:
#   Check test levels of LORD++ algorithm in Section 4.
#
# See Section 4 for details.

amazon_lordpp_testlevel <- function(item_metrics, rating_sequences, thetas, q = 0.1,  obs_index, itemnumbner = itemnumbner, decision_times = decision_times) {
  obs_index <- obs_index[obs_index <= itemnumbner]
  if (length(obs_index) == 0) {
    stop("No valid observation indices provided")
  }
  test_obs_posi <- rep(NA, length(obs_index))
  test_obs_nega <- rep(NA, length(obs_index))
  obs_counter <- 1
  selectposi <- integer(0) 
  selectnega <- integer(0) 
  deltaesti <- rep(0, itemnumbner) 
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
    
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
    if (length(current_ratings) >= 3) {
      p.posi <- wilcox.test(current_ratings, mu = 0, 
                            alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, 
                            alternative = 'less', exact = FALSE)$p.value
    } else {
      p.posi <- 1  
      p.nega <- 1
    }
    thresholdposi <- threldpp(q, q/10, i, selectposi)
    thresholdnega <- threldpp(q, q/10, i, selectnega)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
      } else {
        selectnega <- c(selectnega, i)
      }
    }
    if (i %in% obs_index) {
      test_obs_posi[obs_counter] <- thresholdposi
      test_obs_nega[obs_counter] <- thresholdnega
      obs_counter <- obs_counter + 1
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(
    test_level_obs_posi = test_obs_posi,
    test_level_obs_nega = test_obs_nega))
}

# Function: amazon_saffron_testlevel
#
# Implement the SAFFRON algorithm in Section 4 and record the test levels at specific decision times.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.
#   Clamb            : parameter lambda in SAFFRON.
#   obs_index        : pre-specified tasks for recording test levels.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   test_level_obs_posi : test levels for arm A from the pre-specified tasks.
#   test_level_obs_nega :test levels for arm B from the pre-specified tasks.
# 
# Used in:
#   Check test levels of SAFFRON algorithm in Section 4.
#
# See Section 4 for details.

amazon_saffron_testlevel <- function(item_metrics, rating_sequences, thetas, q = 0.1, Clamb = 0.5, obs_index, itemnumbner = itemnumbner, decision_times = decision_times) {
  obs_index <- obs_index[obs_index <= itemnumbner]
  if (length(obs_index) == 0) {
    stop("No valid observation indices provided")
  }
  test_obs_posi <- rep(NA, length(obs_index))
  test_obs_nega <- rep(NA, length(obs_index))
  obs_counter <- 1
  selectposi <- integer(0)  
  selectnega <- integer(0)  
  cplus <- rep(0, itemnumbner)  
  cminus <- rep(0, itemnumbner) 
  deltaesti <- rep(0, itemnumbner)  
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
      if (length(current_ratings) >= 3) {
      p.posi <- wilcox.test(current_ratings, mu = 0, alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, alternative = 'less', exact = FALSE)$p.value
      if (p.posi <= Clamb) {
        cplus[i] <- 1
      }
      if (p.nega <= Clamb) {
        cminus[i] <- 1
      }
    } else {
      p.posi <- 1 
      p.nega <- 1
      cplus[i] <- 0
      cminus[i] <- 0
    }
    thresholdposi <- thresaff(q, q/2, i, selectposi, cplus, Clamb)
    thresholdnega <- thresaff(q, q/2, i, selectnega, cminus, Clamb)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
      } else {
        selectnega <- c(selectnega, i)
      }
    }
    if (i %in% obs_index) {
      test_obs_posi[obs_counter] <- thresholdposi
      test_obs_nega[obs_counter] <- thresholdnega
      obs_counter <- obs_counter + 1
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(
    test_level_obs_posi = test_obs_posi,
    test_level_obs_nega = test_obs_nega
  ))
}

# Function: amazon_saffron_testlevel
#
# Implement the ADDIS algorithm in Section 4 and record the test levels at specific decision times.
#
# Inputs:
#   item_metrics     : the pre-processed dataframe of items.
#   rating_sequences : the sequence of ratings.
#   thetas           : true states of tasks.
#   q                : target FSR level.
#   lamb             : parameter lambda in ADDIS.
#   tau              : parameter tau in ADDIS.
#   obs_index        : pre-specified tasks for recording test levels.
#   itemnumbner      : number of items.
#   decision_times   : sequence of decision times.
#
# Returns:
#   test_level_obs_posi : test levels for arm A from the pre-specified tasks.
#   test_level_obs_nega :test levels for arm B from the pre-specified tasks.
# 
# Used in:
#   Check test levels of ADDIS algorithm in Section 4.
#
# See Section 4 for details.

amazon_addis_testlevel <- function(item_metrics, rating_sequences, thetas, q = 0.1, lamb = 0.25, tau = 0.5, obs_index, itemnumbner = itemnumbner, decision_times = decision_times) {
  
  obs_index <- obs_index[obs_index <= itemnumbner]
  if (length(obs_index) == 0) {
    stop("No valid observation indices provided")
  }
  test_obs_posi <- rep(NA, length(obs_index))
  test_obs_nega <- rep(NA, length(obs_index))
  obs_counter <- 1
  selectposi <- integer(0)  
  selectnega <- integer(0)  
  kstarposi <- integer(0)  
  kstarnega <- integer(0)   
  splus <- rep(0, itemnumbner) 
  sminus <- rep(0, itemnumbner) 
  cplus <- rep(0, itemnumbner)  
  cminus <- rep(0, itemnumbner) 
  deltaesti <- rep(0, itemnumbner)  
  pb <- txtProgressBar(min = 0, max = itemnumbner, style = 3)
  for (i in 1:itemnumbner) {
    current_time <- decision_times[i]
    current_item <- item_metrics$item_id[i]
    current_ratings <- rating_sequences %>%
      filter(item_id == current_item, review_time <= current_time) %>%
      pull(rating)
    if (length(current_ratings) >= 3) {
      p.posi <- wilcox.test(current_ratings, mu = 0, alternative = 'greater', exact = FALSE)$p.value
      p.nega <- wilcox.test(current_ratings, mu = 0, alternative = 'less', exact = FALSE)$p.value
      if (p.posi <= tau) {
        splus[i] <- 1
      }
      if (p.nega <= tau) {
        sminus[i] <- 1
      }
      
      if (p.posi <= lamb) {
        cplus[i] <- 1
      }
      if (p.nega <= lamb) {
        cminus[i] <- 1
      }
    } else {
      p.posi <- 1 
      p.nega <- 1
      splus[i] <- 0
      sminus[i] <- 0
      cplus[i] <- 0
      cminus[i] <- 0
    }
    
    sindexposi <- sum(splus[1:i])
    sindexnega <- sum(sminus[1:i])
    
    thresholdposi <- threaddis(q, q/2, i, selectposi, cplus, sindexposi, kstarposi, lamb = lamb, tau = tau)
    thresholdnega <- threaddis(q, q/2, i, selectnega, cminus, sindexnega, kstarnega, lamb = lamb, tau = tau)
    deltaplus <- p.posi < thresholdposi
    deltaminus <- p.nega < thresholdnega
    deltatotal <- as.numeric(deltaplus) - as.numeric(deltaminus)
    
    if (deltatotal != 0) {
      deltaesti[i] <- deltatotal
      if (deltatotal > 0) {
        selectposi <- c(selectposi, i)
        kstarposi <- c(kstarposi, sum(splus[1:i]))
      } else {
        selectnega <- c(selectnega, i)
        kstarnega <- c(kstarnega, sum(sminus[1:i]))
      }
    }
    if (i %in% obs_index) {
      test_obs_posi[obs_counter] <- thresholdposi
      test_obs_nega[obs_counter] <- thresholdnega
      obs_counter <- obs_counter + 1
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(
    test_level_obs_posi = test_obs_posi,
    test_level_obs_nega = test_obs_nega
  ))
}

#### Compare test levels of different methods in Section 4 (figure 6) #####
amazon_data <- fread("real-amazon25/AMAZON_FASHION.csv",
                     col.names = c("item_id", "user_id", "rating", "timestamp"))

processed_data <- amazon_data %>%
  filter(rating >= 1 & rating <= 5) %>%
  arrange(item_id, timestamp)

item_metrics <- processed_data %>%
  group_by(item_id) %>%
  summarise(
    first_review_time = min(timestamp),
    total_reviews = n(),
    avg_rating = mean(rating)
  ) %>%
  arrange(first_review_time) %>%
  filter(total_reviews > 50)  

min_time <- min(item_metrics$first_review_time)
item_metrics <- item_metrics %>%
  mutate(arrival_day = as.integer((first_review_time - min_time) / 86400)) %>%
  arrange(arrival_day, item_id)

item_metrics <- item_metrics %>%
  group_by(arrival_day) %>%
  slice(1) %>% 
  ungroup() %>%
  mutate(arrival_time = arrival_day) %>%  
  select(-arrival_day, -first_review_time)

thetas <- numeric(nrow(item_metrics))
thetas[item_metrics$avg_rating >= 3] <- 1  
thetas[item_metrics$avg_rating < 3] <- -1  
arrival_map <- item_metrics %>% select(item_id, arrival_time)

rating_sequences <- processed_data %>%
  inner_join(arrival_map, by = "item_id") %>%
  arrange(item_id, timestamp) %>%
  group_by(item_id) %>%
  mutate(review_seq = row_number(),
         review_time = as.integer((timestamp - min_time) / 86400)) %>%
  ungroup()

rating_sequences$rating = rating_sequences$rating - 3

item_metrics <- item_metrics %>% arrange(arrival_time)
itemnumbner <- nrow(item_metrics)

time_max = max(rating_sequences$review_time)
decision_times <- c(item_metrics$arrival_time[-1],time_max)
num_decision = length(decision_times)

obs_index = seq(100,500,100)
time_obs = seq(100, 900, 50)
num_obs = length(obs_index)
re = amazon_sava_testlevel(item_metrics, rating_sequences, thetas, bound = 4, q, k = 100, w0 = q, obs_index, itemnumbner = itemnumbner, decision_times = decision_times)
threposi.sava = re$test_level_obs
threnega.sava = re$test_level_obs

re = amazon_lordpp_testlevel(item_metrics, rating_sequences, thetas, q, obs_index, itemnumbner = itemnumbner, decision_times = decision_times)
threposi.lordpp = re$test_level_obs_posi
threnega.lordpp = re$test_level_obs_nega

re = amazon_saffron_testlevel(item_metrics, rating_sequences, thetas, q, Clamb = 0.5,  obs_index, itemnumbner = itemnumbner, decision_times = decision_times)
threposi.saffron = re$test_level_obs_posi
threnega.saffron = re$test_level_obs_nega

re = amazon_addis_testlevel(item_metrics, rating_sequences, thetas, q,  lamb = 0.25, tau = 0.5, obs_index, itemnumbner = itemnumbner, decision_times = decision_times)
threposi.addis = re$test_level_obs_posi
threnega.addis = re$test_level_obs_nega

threposi.dat = tibble(threshold = c(threposi.lordpp[1:5], threposi.saffron[1:5], threposi.addis[1:5]), Method = rep(c('LORD++', 'SAFFRON', 'ADDIS'), each = num_obs),
                      Index = rep(seq(100,500,100), 3))

threnega.dat = tibble(threshold = c( threnega.lordpp[1:5], threnega.saffron[1:5], threnega.addis[1:5]), Method = rep(c('LORD++', 'SAFFRON', 'ADDIS'), each = num_obs),
                      Index = rep(seq(100,500,100), 3))

write.csv(threposi.dat, 'threposi-online.csv')
write.csv(threnega.dat, 
          'threnega-online.csv')

thresava = threposi.sava
thresava.dat = tibble(threshold = as.vector(thresava), index = rep(seq(100,500,100), length(time_obs)), decision_time = rep(time_obs, each = 5))
write.csv(thresava.dat, 'threshold-sava.csv')

thresholdcompare = tibble(Test_level = c(rep(as.vector(thresava),2), threposi.dat$threshold, threnega.dat$threshold), 
                          Index = c(rep(thresava.dat$index, 2), threposi.dat$Index, threnega.dat$Index),
                          Decision_time = c(rep(thresava.dat$decision_time, 2), threposi.dat$Index, threnega.dat$Index),
                          Method = c(rep('SAVA', 2*length(thresava.dat$index)), threposi.dat$Method, threnega.dat$Method),
                          Arm = c(rep('A', length(thresava.dat$index)), rep('B', length(thresava.dat$index)), rep('A', length(threposi.dat$Index)), rep('B', length(threnega.dat$Index)))
                          )
write.csv(thresholdcompare, 
          'threshold-compare.csv', row.names = F)