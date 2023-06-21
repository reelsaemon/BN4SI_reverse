# Preliminaries -----------------------------------------------------------

gc()

# extended functionality
require(lubridate)
require(data.table)
require(tidyr)
require(plyr)
require(dplyr)
require(stringr)
require(rlist)
require(DescTools)

# model packages
require(bnlearn)

# evaluation metrics
require(caret)
require(Metrics)
require(StatRank)
require(comparator)

# parallelization
require(foreach)
require(doParallel)

source('src/BN4SI_functions.R')


# Model setup -------------------------------------------------------------

set.seed(26477)

# constants
totalCores <- detectCores()
usedCores <- totalCores-1
n_samples <- 10000
training_size <- 2/3
test_size <- 1-training_size

# Data import -------------------------------------------------------------

event_log_data <- as.data.table(read.csv("data/sim_event_log.csv", header = T, sep = ","))
id_colname <- 'order_id'
activity_colname <- 'station'
timestamp_colname <- 'timestamp_in'
event_log_data <- event_log_data[, c(id_colname, activity_colname, timestamp_colname), with=FALSE]

# Data preprocessing ------------------------------------------------------

event_log_data[, (activity_colname) := as.character(get(activity_colname))]
event_log_data[, FirstTimestamp := min(get(timestamp_colname)), by=get(id_colname)]
event_log_data <- event_log_data[order(FirstTimestamp), ]

# chronological train test split
trace_ids <- event_log_data[, unique(get(id_colname))]
train_ids <- trace_ids[1:floor((training_size)*length(trace_ids))]
test_ids <- trace_ids[(floor((training_size)*length(trace_ids))+1):length(trace_ids)]

# # random train test split
# train_ids <- sample(event_log_data[, unique(get(id_colname))], floor(training_size*event_log_data[, length(unique(get(id_colname)))]), replace=FALSE)
# test_ids <- setdiff(event_log_data[, unique(get(id_colname))], train_ids)

training_set <- event_log_data[get(id_colname) %in% train_ids, ]
test_set <- event_log_data[get(id_colname) %in% test_ids, ]

# widening datasets for parameter learning
training_set_wide <- make_wide_format_reverse(training_set, activity_colname, id_colname, handle_loops=TRUE)
test_set_wide <- make_wide_format_reverse(test_set, activity_colname, id_colname, handle_loops=TRUE, preserve_exec=TRUE)

# add non-occurring (test set) columns to the test set
no_cols <- setdiff(colnames(training_set_wide), colnames(test_set_wide))
occ_cols <- setdiff(colnames(test_set_wide), c(id_colname, "execution_order"))
train_cols <- setdiff(colnames(training_set_wide), c(id_colname, "execution_order"))
full_levels <- unique(unlist(lapply(train_cols, function(x) levels(unlist(training_set_wide[, x, with=FALSE])))))

if(length(no_cols)>1) {
  
  test_set_wide[, (no_cols) := factor(NA, levels=full_levels)]
}

# Network construction and fitting ----------------------------------------

network_structure <- make_sparse_binary_structure_reverse_transitivity_check_max_parents(training_set, activity_colname, id_colname, handle_loops=TRUE)

fitted_network <- bn.fit(network_structure, training_set_wide[, nodes(network_structure), with=FALSE], replace.unidentifiable=TRUE)

# Next Activity/Remaining Trace Prediction --------------------------------

test_data <- test_set_wide
test_data_traces <- test_data[, unique(get(id_colname))]

cluster <- makeCluster(usedCores, outfile="parallel_log.txt")
registerDoParallel(cluster)

remaining_traces_list <- foreach(i=test_data_traces, .packages=c("data.table", "bnlearn", "tidyr", "DescTools")) %dopar% {

  set.seed(26477)
  
  current_data <- test_data[get(id_colname)==i, ]
  
  current_data[, pred_next_activity := character(0)]
  current_data[, pred_remaining_trace := list(character(0))]
  
  print(which(test_data_traces==i)/length(test_data_traces))
  
  for(j in 1:nrow(current_data)) {
    
    non_na_activities <- unlist(current_data[j, lapply(.SD, function(x) !is.na(x)), .SDcols=setdiff(colnames(current_data), c(id_colname, "execution_order", "pred_remaining_trace", "pred_next_activity"))])
    
    known_activities <- names(non_na_activities[non_na_activities==TRUE])
    
    unknown_activities <- names(non_na_activities[non_na_activities==FALSE])
    
    current_step <- j
    
    prediction <- predict_case_activities_NO_only_above_50(current_data[j, ],
                                                           network_structure=network_structure,
                                                           fitted_network=fitted_network,
                                                           id_colname=id_colname,
                                                           activity_colname=activity_colname,
                                                           n_samples=n_samples)
    
    next_activity <- names(prediction$predicted_activities_without_NO[[1]][which(prediction$predicted_activities_without_NO[[1]]==current_step+1)])
    remaining_trace <- prediction$predicted_activities_without_NO[[1]][which(prediction$predicted_activities_without_NO[[1]]>=current_step+1)]
    
    
    if(any(table(remaining_trace) > 1)) {
      
      for(place in names(table(remaining_trace)[which(table(remaining_trace)>1)])) {
        
        place_activities <- names(remaining_trace[remaining_trace==place])
        picked_activity_for_place <- names(which.max(prediction$predicted_activities_without_NO_probs[[1]][place_activities]))
        
        remaining_trace <- c(remaining_trace[remaining_trace!=place], setNames(place, picked_activity_for_place))
      }
    }
    
    predicted_remaining_trace <- names(sort(remaining_trace, decreasing=FALSE))
    
    if(length(next_activity)>1) {
      next_activity <- next_activity[which.max(prediction$predicted_activities_without_NO_probs[[1]][next_activity])]
    }
    
    predicted_next_activity <- next_activity
    
    if(!is.null(predicted_remaining_trace)) {
      
      current_data[j, pred_remaining_trace:=list(predicted_remaining_trace)]
    }
    
    if(length(predicted_next_activity)!=0) {
      
      current_data[j, pred_next_activity:=predicted_next_activity]
    }
  }
  
  return(current_data)
}

stopCluster(cluster)

remaining_traces <- do.call(rbind, remaining_traces_list)[, c(id_colname, "execution_order", "pred_remaining_trace", "pred_next_activity"), with=FALSE]

# Evaluation of predicted next activities/remaining traces ----------------

overall_comp_frame <- copy(test_data)
act_colnames <- setdiff(colnames(test_data), c(id_colname, "execution_order"))

test_data_complete <- copy(overall_comp_frame)
test_data_complete[, (act_colnames) := lapply(.SD, function(x) unique(na.omit(x))), .SDcols=act_colnames, by=get(id_colname)]

test_data_remaining_activities <- suppressWarnings(test_data_complete[, lapply(.SD, function(x) as.numeric(as.character(x)) > as.numeric(as.character(execution_order))), .SDcols=act_colnames])
test_data_remaining_activities[, act_remaining_trace := lapply(apply(test_data_remaining_activities, 1, function(x) which(x==TRUE)), names)]

test_data_next_activities <- suppressWarnings(test_data_complete[, lapply(.SD, function(x) as.numeric(as.character(x)) == as.numeric(as.character(execution_order))+1), .SDcols=act_colnames])
test_data_next_activities[, act_next_activity := lapply(apply(test_data_next_activities, 1, function(x) which(x==TRUE)), names)]

test_data_complete[, act_next_activity:=test_data_next_activities[, act_next_activity]]
test_data_complete[, act_remaining_trace:=test_data_remaining_activities[, act_remaining_trace]]

overall_comp_frame <- left_join(overall_comp_frame, remaining_traces, by=c(id_colname, "execution_order"))
overall_comp_frame <- left_join(overall_comp_frame, test_data_complete[, c(id_colname, "execution_order", "act_remaining_trace", "act_next_activity"), with=FALSE], by=c(id_colname, "execution_order"))

# with all prefixes except last (1 <= k < |T|)
eval_frame <- overall_comp_frame[execution_order != 1 & lengths(act_remaining_trace) != 0,]
eval_frame[, act_next_activity:=unlist(act_next_activity)]

# next activities
eval_frame[, act_next_activity := act_next_activity]
eval_frame[, pred_next_activity := pred_next_activity]
print(paste0("NAP Accuracy: ", sum(na.omit(eval_frame[, act_next_activity == pred_next_activity]))/nrow(eval_frame)))

# remaining traces
eval_frame[, act_remaining_trace_char := as.character(act_remaining_trace)]
eval_frame[, pred_remaining_trace_char := as.character(pred_remaining_trace)]
print(paste0("RTP Accuracy: ", sum(eval_frame[, act_remaining_trace == pred_remaining_trace_char])/nrow(eval_frame)))

similarities <- c()

for(index in 1:nrow(eval_frame)) {
  
  similarities <- c(similarities, manual_ndls(eval_frame[index, act_remaining_trace], eval_frame[index, pred_remaining_trace]))
}

print(paste0("Suffix NDLS: ", mean(similarities)))