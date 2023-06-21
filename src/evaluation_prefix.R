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
test_set_wide[, (occ_cols) := lapply(.SD, function(x) unique(na.omit(x))), .SDcols=occ_cols, by=get(id_colname)]


test_data_complete <- copy(test_set_wide)
keep_info_table <- copy(test_set_wide)
keep_info_table[, execution_order := as.numeric(as.character(execution_order))]
keep_info_table[, (occ_cols) := lapply(.SD, function(x) replace(x, x=="NO", NA)), .SDcols=occ_cols]
keep_info_table[, (occ_cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols=occ_cols]
keep_info_table[, (occ_cols) := lapply(.SD, function(x) x>=execution_order), .SDcols=occ_cols]

keep_info_table[, (occ_cols) := lapply(.SD, function(x) replace(x, is.na(x), FALSE)), .SDcols=occ_cols]

keep_matrix <- as.matrix(keep_info_table[, occ_cols, with=FALSE])
test_set_wide_matrix <- as.matrix(test_set_wide[, occ_cols, with=FALSE])
test_set_wide_matrix[keep_matrix==FALSE] <- NA
test_set_wide_info_table <- as.data.table(test_set_wide_matrix)
test_set_wide_info_table[, (colnames(test_set_wide_info_table)) := lapply(.SD, function(x) factor(x, levels=full_levels)), .SDcols=colnames(test_set_wide_info_table)]
test_set_wide[, (occ_cols) := NULL]
test_set_wide <- cbind(test_set_wide, test_set_wide_info_table)

# add non-occurring (test set) columns to the test set
if(length(no_cols)>1) {
  
  test_set_wide[, (no_cols) := factor(NA, levels=full_levels)]
}

# Network construction and fitting ----------------------------------------

network_structure <- make_sparse_binary_structure_reverse_transitivity_check_max_parents(training_set, activity_colname, id_colname, handle_loops=TRUE)

fitted_network <- bn.fit(network_structure, training_set_wide[, nodes(network_structure), with=FALSE], replace.unidentifiable=TRUE)

# Previous Activity/Preceding Trace Prediction ----------------------------

test_data <- test_set_wide

test_data_traces <- test_data[, unique(get(id_colname))]

cluster <- makeCluster(usedCores, outfile="parallel_log.txt")
registerDoParallel(cluster)

# for(i in test_data_traces) {
preceding_traces_list <- foreach(i=test_data_traces, .packages=c("data.table", "bnlearn", "tidyr", "DescTools")) %dopar% {
  
  set.seed(26477)
  
  current_data <- test_data[get(id_colname)==i, ]
  
  current_data[, pred_previous_activity := character(0)]
  current_data[, pred_preceding_trace := list(character(0))]
  
  print(which(test_data_traces==i)/length(test_data_traces))
  
  for(j in 2:nrow(current_data)) {
    
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
    
    previous_activity <- names(prediction$predicted_activities_without_NO[[1]][which(prediction$predicted_activities_without_NO[[1]]==current_step-1)])
    preceding_trace <- prediction$predicted_activities_without_NO[[1]][which(prediction$predicted_activities_without_NO[[1]]<=current_step-1)]
    
    if(any(table(preceding_trace) > 1)) {
      
      for(place in names(table(preceding_trace)[which(table(preceding_trace)>1)])) {
        
        place_activities <- names(preceding_trace[preceding_trace==place])
        picked_activity_for_place <- names(which.max(prediction$predicted_activities_without_NO_probs[[1]][place_activities]))
        
        preceding_trace <- c(preceding_trace[preceding_trace!=place], setNames(place, picked_activity_for_place))
      }
    }
    
    predicted_preceding_trace <- names(sort(preceding_trace, decreasing=FALSE))
    
    if(length(previous_activity)>1) {
      previous_activity <- previous_activity[which.max(prediction$predicted_activities_without_NO_probs[[1]][previous_activity])]
    }
    
    predicted_previous_activity <- previous_activity
    
    if(!is.null(predicted_preceding_trace)) {
      
      current_data[j, pred_preceding_trace:=list(predicted_preceding_trace)]
    }
    
    if(length(predicted_previous_activity)!=0) {
      
      current_data[j, pred_previous_activity:=predicted_previous_activity]
    }
  }
  return(current_data)
}

stopCluster(cluster)

preceding_traces <- do.call(rbind, preceding_traces_list)[, c(id_colname, "execution_order", "pred_preceding_trace", "pred_previous_activity"), with=FALSE]


# Evaluation of predicted previous activities/preceding traces ------------

overall_comp_frame <- copy(test_data)
act_colnames <- setdiff(colnames(test_data), c(id_colname, "execution_order"))

test_data_complete <- copy(overall_comp_frame)
test_data_complete[, (act_colnames) := lapply(.SD, function(x) unique(na.omit(x))), .SDcols=act_colnames, by=get(id_colname)]

test_data_preceding_activities <- suppressWarnings(test_data_complete[, lapply(.SD, function(x) as.numeric(as.character(x)) < as.numeric(as.character(execution_order))), .SDcols=act_colnames])
test_data_preceding_activities[, act_preceding_trace := lapply(apply(test_data_preceding_activities, 1, function(x) which(x==TRUE)), names)]

test_data_previous_activities <- suppressWarnings(test_data_complete[, lapply(.SD, function(x) as.numeric(as.character(x)) == as.numeric(as.character(execution_order))-1), .SDcols=act_colnames])
test_data_previous_activities[, act_previous_activity := lapply(apply(test_data_previous_activities, 1, function(x) which(x==TRUE)), names)]

test_data_complete[, act_previous_activity:=test_data_previous_activities[, act_previous_activity]]
test_data_complete[, act_preceding_trace:=test_data_preceding_activities[, act_preceding_trace]]

overall_comp_frame <- left_join(overall_comp_frame, preceding_traces, by=c(id_colname, "execution_order"))
overall_comp_frame <- left_join(overall_comp_frame, test_data_complete[, c(id_colname, "execution_order", "act_preceding_trace", "act_previous_activity"), with=FALSE], by=c(id_colname, "execution_order"))

eval_frame <- copy(overall_comp_frame)
eval_frame[, first := (execution_order==1)]
eval_frame[, last := execution_order==last(execution_order), by=get(id_colname)]

# with all prefixes except last (1 <= k < |T|)
eval_frame <- eval_frame[first==FALSE, ]
eval_frame[, first := NULL]
eval_frame[, last := NULL]

# previous activities
eval_frame[, act_previous_activity := act_previous_activity]
eval_frame[, pred_previous_activity := pred_previous_activity]
print(paste0("PAP Accuracy: ", sum(na.omit(eval_frame[, act_previous_activity == pred_previous_activity]))/nrow(eval_frame)))

# preceding traces
eval_frame[, act_preceding_trace_char := as.character(act_preceding_trace)]
eval_frame[, pred_preceding_trace_char := as.character(pred_preceding_trace)]
print(paste0("PTP Accuracy: ", sum(eval_frame[, act_preceding_trace == pred_preceding_trace_char])/nrow(eval_frame)))

similarities <- c()

for(index in 1:nrow(eval_frame)) {
  
  similarities <- c(similarities, manual_ndls(eval_frame[index, act_preceding_trace], eval_frame[index, pred_preceding_trace]))
}

print(paste0("Prefix NDLS: ", mean(similarities)))