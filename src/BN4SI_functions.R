# function for follows graphs (eventually follows, directly follows)
get_follows_graph <- function(log, activity_colname, id_colname, plot=TRUE, add_text=TRUE, sort=TRUE, horizon=1, only_plot=FALSE, transition_threshold=1) {
  
  # params
  # horizon     if horizon=1 a directly follows graph is computed
  #             if 1 < horizon <= n a follows graph with n steps into the future trace is computed. If horizon>=n with n being the maximum trace length in the event log
  #             we get a eventually follows graph for all possible traces and trace lengths
  # 
  # sort=TRUE   sorts the activities with respect to the most followers, i.e. activities that have many followers are assigned a high spot in the matrix
  #             with eventually follows graphs those activities are most likely the ones earlier in the process
  #             with directly follows graphs those activities are most likely choice node activities, i.e. activities after which many different activities can be carried out
  
  activities <- sort(unlist(unique(log[[activity_colname]])))
  n_activities <- length(activities)
  
  ids <- sort(unlist(unique(log[[id_colname]])))
  
  follows_matrix <- matrix(0, nrow=n_activities, ncol=n_activities, dimnames=list(activities, activities))
  
  for(current_id in ids) {
    
    current_trace_indices <- which(log[[id_colname]] == current_id)
    current_trace <- log[current_trace_indices, ][[activity_colname]]
    trace_length <- length(current_trace)
    
    i <- 1
    
    while(i <= trace_length) {
      
      current_activity <- current_trace[i]
      current_remaining_trace <- current_trace[-(1:i)][1:min(horizon, length(current_trace) - i)]
      
      if(!any(is.na(current_remaining_trace))) {
        for(act in current_remaining_trace) {
          
          follows_matrix[as.character(current_activity), as.character(act)] <- follows_matrix[as.character(current_activity), as.character(act)] + 1
        }
      }
      
      i <- i+1
    }
  }
  
  # apply threshold
  follows_matrix[apply(follows_matrix, 2, function(x) x/log[, table(get(activity_colname))]) < (1-transition_threshold)] <- 0
  
  if(sort==TRUE) {
    
    non_zero_entries <- follows_matrix
    non_zero_entries[non_zero_entries > 0] <- 1
    non_zero_followers <- apply(non_zero_entries, 1, sum)
    
    sorted_order <- names(sort(non_zero_followers, decreasing=TRUE))
    
    # reorder matrix rows/columns
    follows_matrix <- follows_matrix[sorted_order, sorted_order]
  }
  
  if(plot==TRUE) {
    
    image(t(follows_matrix[nrow(follows_matrix):1,]), xlab="to", ylab="from", axes=FALSE, col=heat.colors(dim(follows_matrix)[1]**2), useRaster=TRUE, main=paste0("Follows Graph of ", deparse(substitute(log)), " with horizon ", horizon))
    abline(1, -1)
    axis(1, at=seq(0,1,1/(ncol(follows_matrix)-1)), labels=rownames(follows_matrix))
    axis(2, at=seq(0,1,1/(ncol(follows_matrix)-1)), labels=rev(colnames(follows_matrix)))
    
    # values as text
    if(add_text==TRUE) {
      
      text(expand.grid(seq(0,1,1/(ncol(follows_matrix)-1)), seq(0,1,1/(ncol(follows_matrix)-1))), labels=t(follows_matrix[nrow(follows_matrix):1,]), cex=0.66)
    }
    
  }
  
  if(plot & only_plot) {
    return(invisible(NULL))
  } else {
    
    return(follows_matrix)
  }
  
}

# function for footprint matrix from a given follows graph
footprint <- function(follows_graph) {

  fg_dimnames <- dimnames(follows_graph)
  
  # initialize matrix with same dimensions as a directly/eventually follows graph
  footprint_matrix <- matrix(NA, nrow=nrow(follows_graph), ncol=ncol(follows_graph), dimnames=fg_dimnames)
  
  # check the follows graph for existing connections/transitions and fill them into the footprint matrix
  for(row in 1:nrow(footprint_matrix)) {
    
    for(col in row:ncol(footprint_matrix)) {
      
      if(follows_graph[row, col] > 0 & follows_graph[col, row] == 0) {
        
        footprint_matrix[row, col] <- ">"
        footprint_matrix[col, row] <- "<"
      } else if(follows_graph[row, col] == 0 & follows_graph[col, row] > 0) {
        
        footprint_matrix[row, col] <- "<"
        footprint_matrix[col, row] <- ">"
      } else if(follows_graph[row, col] > 0 & follows_graph[col, row] > 0){
        
        footprint_matrix[row, col] <- "||"
        footprint_matrix[col, row] <- "||" 
      } else {
        
        footprint_matrix[row, col] <- "-"
        footprint_matrix[col, row] <- "-"
      }
    }
  }
  
  return(footprint_matrix)
}

# function for converting event log to wide format
make_wide_format_reverse <- function(log_data, activity_colname, id_colname, preserve_exec=FALSE, handle_loops=FALSE) {
  
  require(mltools)
  require(tidyr)
  
  log_data <- copy(log_data[, c(id_colname, activity_colname), with=FALSE])
  original_max_trace_length <- max(log_data[, .N, by=get(id_colname)][, N])
  max_trace_length <- max(log_data[, trace_length := .N, by=get(id_colname)][, trace_length])
  
  
  log_data[, place := factor(1:.N), by=get(id_colname)]
  
  if(handle_loops) {
    
    recoded_activities <- log_data[, .(activities = activity_count_helper(get(activity_colname))), by=get(id_colname)][, activities]
    log_data[, (activity_colname) := recoded_activities]
  }
  
  if(preserve_exec) {
    
    log_data[, execution_order := place]
    
    # widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "place", "execution_order"), with=FALSE], names_from=activity_colname, values_from="place"))
    
    wide_colnames <- as.character(log_data[, unique(get(activity_colname))])
    wide_log_data[, (wide_colnames) := lapply(.SD, carry_helper), .SDcols=wide_colnames, by=get(id_colname)]
    
    wide_log_data[, last_helper := c(rep(FALSE, .N-1), TRUE), by=get(id_colname)]
    
    last_indices <- which(wide_log_data[, last_helper])
    
    for(j in wide_colnames) {
      
      set(wide_log_data, i=intersect(which(is.na(wide_log_data[[j]])), last_indices), j=j, value="NO")
    }
    
    wide_log_data[, last_helper := NULL]
    
  } else {
    
    # widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "place"), with=FALSE], names_from=activity_colname, values_from="place"))
    
    for(j in as.character(log_data[, unique(get(activity_colname))])) {
      
      set(wide_log_data, i=which(is.na(wide_log_data[[j]])), j=j, value="NO")
    }
    
  }
  
  return(wide_log_data)
}

# function for converting event log to wide format for prefix prediction
make_wide_format_reverse_prefix <- function(log_data, activity_colname, id_colname, preserve_exec=FALSE, handle_loops=FALSE) {
  
  require(mltools)
  require(tidyr)
  
  log_data <- copy(log_data[, c(id_colname, activity_colname), with=FALSE])
  original_max_trace_length <- max(log_data[, .N, by=get(id_colname)][, N])
  max_trace_length <- max(log_data[, trace_length := .N, by=get(id_colname)][, trace_length])
  
  
  log_data[, place := factor(1:.N), by=get(id_colname)]
  
  if(handle_loops) {
    
    recoded_activities <- log_data[, .(activities = activity_count_helper(get(activity_colname))), by=get(id_colname)][, activities]
    log_data[, (activity_colname) := recoded_activities]
  }
  
  if(preserve_exec) {
    
    log_data[, execution_order := place]
    
    # widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "place", "execution_order"), with=FALSE], names_from=activity_colname, values_from="place"))
    
    wide_colnames <- as.character(log_data[, unique(get(activity_colname))])
    wide_log_data[, (wide_colnames) := lapply(.SD, carry_helper_prefix), .SDcols=wide_colnames, by=get(id_colname)]
    
    wide_log_data[, last_helper := c(rep(FALSE, .N-1), TRUE), by=get(id_colname)]
    
    last_indices <- which(wide_log_data[, last_helper])
    
    for(j in wide_colnames) {
      
      set(wide_log_data, i=intersect(which(is.na(wide_log_data[[j]])), last_indices), j=j, value="NO")
    }
    
    wide_log_data[, last_helper := NULL]
    
  } else {
    
    # widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "place"), with=FALSE], names_from=activity_colname, values_from="place"))
    
    for(j in as.character(log_data[, unique(get(activity_colname))])) {
      
      set(wide_log_data, i=which(is.na(wide_log_data[[j]])), j=j, value="NO")
    }
    
  }
  
  return(wide_log_data)
}

# function to switch types of a vector with the type given as argument
type_switcher <- function(vec, type, levels) {
  
  if(type=="character") {
    
    switched_vec <- as.character(vec)
    
  } else if(type=="factor" & !is.na(levels)) {
    
    switched_vec <- factor(vec, levels=unlist(levels))
    
  } else if(type=="numeric") {
    
    switched_vec <- as.numeric(vec)
    
  } else if(type=="integer") {
    
    switched_vec <- as.integer(vec)
    
  } else if(type=="logical") {
    
    switched_vec <- as.logical(vec)
    
  } else if(type=="complex") {
    
    switched_vec <- as.complex(vec)
  }
  
  return(switched_vec)
}

# helper function to carry the first non-NA element of a vector onto the remaining elements of the vector
carry_helper <- function(vector) {
  
  if(length(vector)==1) {
    
    return(vector)
  }
  
  current_element <- NA
  
  for(i in 1:(length(vector)-1)) {
    current_element <- vector[i]
    if(!is.na(current_element)) {
      vector[i+1] <- current_element
    }
  }
  
  return(vector)
}

# helper function to carry the first non-NA element of a vector onto the remaining elements of the vector
carry_helper_prefix <- function(vector) {
  
  if(length(vector)==1) {
    
    return(vector)
  }
  
  current_element <- NA
  
  for(i in (length(vector)-1):1) {
    current_element <- vector[i]
    if(!is.na(current_element)) {
      vector[i-1] <- current_element
    }
  }
  
  return(vector)
}

# function for automatic generation of a sparse binary approach BN structure with the reverse approach while keeping transitivity intact in case of loops
make_sparse_binary_structure_reverse_transitivity_check_max_parents <- function(log, activity_colname, id_colname, keep_non_occurring=TRUE, handle_loops=FALSE, transition_threshold=1) {
  
  log <- copy(log[, c(id_colname, activity_colname), with=FALSE])
  max_trace_length <- max(log[, trace_length := .N, by=get(id_colname)][, trace_length])
  
  max_parent_count <- floor((2^31)^(1/(max_trace_length+1)))
  
  if(handle_loops) {
    
    recoded_activities <- log[, .(activities = activity_count_helper(get(activity_colname))), by=get(id_colname)][, activities]
    log[, (activity_colname) := recoded_activities]
  }
  
  follows_graph <- get_follows_graph(log, activity_colname=activity_colname, id_colname=id_colname, plot=FALSE, add_text=FALSE, sort=FALSE, horizon=max_trace_length, transition_threshold=transition_threshold)
  
  footprint <- footprint(follows_graph)
  
  occurring_transitions <- data.frame(from=character(0), to=character(0))
  for(row in dimnames(footprint)[[1]]) {
    
    for(col in dimnames(footprint)[[2]]) {
      
      if(footprint[row, col] == ">") {
        occurring_transitions <- rbind.data.frame(occurring_transitions,data.frame(from=row, to=col))
      }
    }
  }
  
  # look for conflicts w.r.t. transitivity
  successors_list <- lapply(rep(NA, length(dimnames(footprint)[[1]])), function(x) c())
  names(successors_list) <- dimnames(footprint)[[1]]
  
  # look for conflicts w.r.t. transitivity
  all_conflicts <- data.frame(from=c(), to=c())
  
  for(emitting_act in names(successors_list)) {
    
    # look up all possible level 1 successors
    possible_l1_successors <- names(footprint[emitting_act, ][which(footprint[emitting_act, ] == ">")])
    
    for(l1_successor in possible_l1_successors) {
      
      # look up all possible level 2 successors
      possible_l2_successors <- names(footprint[l1_successor, ][which(footprint[l1_successor, ] == ">")])
      
      for(l2_successor in possible_l2_successors) {
        
        # check if we have any connection to our initial activity
        if(emitting_act %in% names(footprint[l2_successor, ][which(footprint[l2_successor, ] %in% c(">", "||"))])) {
          
          all_conflicts <- rbind.data.frame(all_conflicts, data.frame(from=l1_successor, to=l2_successor))
        }
      }
    }
  }
  
  # exclude found conflicts from occurring_transitions
  for(conflict in 1:nrow(all_conflicts)) {
    
    if(nrow(occurring_transitions[occurring_transitions[, "from"]==all_conflicts[conflict, "from"]
                                  & occurring_transitions[, "to"]==all_conflicts[conflict, "to"], ]) > 0) {
      
      occurring_transitions <- occurring_transitions[-c(which(occurring_transitions[, "from"]==all_conflicts[conflict, "from"] 
                                                              & occurring_transitions[, "to"]==all_conflicts[conflict, "to"])), ]
    }
  }
  
  if(keep_non_occurring) {
    
    activities <- log[, unique(get(activity_colname))]
    sparse_network <- cpdag(empty.graph(nodes=activities))
    # adjust arcs of the network to the ones in occurring_transitions
    arcs(sparse_network) <- occurring_transitions
  } else {
    
    activities <- unique(c(occurring_transitions[, "from"], occurring_transitions[, "to"]))
    sparse_network <- cpdag(empty.graph(nodes=activities))
    # adjust arcs of the network to the ones in occurring_transitions
    arcs(sparse_network) <- occurring_transitions
  }
  
  # remove parents if parent count exceeds max_parent_count
  efg <- get_follows_graph(log=log, activity_colname=activity_colname, id_colname=id_colname, plot=FALSE, add_text=FALSE, sort=FALSE, horizon=max_trace_length, transition_threshold=transition_threshold)
  for(node in nodes(sparse_network)) {
    
    # check if parent count exceeds max_parent_count
    if(length(sparse_network$nodes[[node]]$parents) > max_parent_count) {
      
      parent_activities <- sparse_network$nodes[[node]]$parents
      
      parent_connections_to_node <- list()
      
      for(pa in parent_activities) {
        
        parent_connections_to_node[[pa]] <- efg[pa, node]
      }
      
      most_influential_parents <- sort(unlist(parent_connections_to_node), decreasing=TRUE)[1:max_parent_count]
      
      sparse_network$nodes[[node]]$parents <- names(most_influential_parents)
    }
  }
  
  
  return(sparse_network)
}

predict_case_activities_NO_only_above_50 <- function(full_case_data, network_structure, fitted_network, id_colname, activity_colname, n_samples) {
  
  case_id <- full_case_data[, unique(get(id_colname))]
  
  # queries
  evidence_sets <- lapply(1:nrow(full_case_data), function(row) full_case_data[row, lapply(.SD, function(x) if(any(is.na(x))){NULL} else {x}), .SDcols=nodes(network_structure)])
  
  predicted_distributions <- list()
  plot_nobs <- list()
  plot_dists_activities <- list()
  
  for(index in 1:length(evidence_sets)) {
    
    pred_place_dist <- tryCatch(cpdist(fitted_network, nodes=as.character(nodes(network_structure)), evidence=as.list(evidence_sets[[index]]), method='lw', n=n_samples),
                                error=function(cond) {return(NA)})
    
    suppressWarnings({if(!is.na(pred_place_dist)) {
      
      
      pred_place_dist_non_NA <- na.omit(pred_place_dist)
      
      activity_cols <- names(network_structure$nodes)
      long_dists_activities <- as.data.table(pivot_longer(pred_place_dist_non_NA[, activity_cols], cols=activity_cols, names_to="activity", values_to="place"))
      
      predicted_distributions[[index]] <- pred_place_dist_non_NA
      
      plot_nobs[[index]] <- nrow(pred_place_dist_non_NA)
      plot_dists_activities[[index]] <- long_dists_activities
    } else {
      
      predicted_distributions[[index]] <- NA
      plot_nobs[[index]] <- NA
      plot_dists_activities[[index]] <- NA
      
    }
    })
  }
  
  predicted_non_occurrence_probs <- suppressWarnings(lapply(predicted_distributions, function(x) {if(!is.na(x)) {apply(x, 2, function(x) sum(x=="NO")/length(x))} else {return(NA)}}))
  predicted_activities_without_NO <- suppressWarnings(lapply(predicted_distributions, function(x) {if(!is.na(x)) {apply(x, 2, function(x) Mode(x[x!="NO"])[1])} else {return(NA)}}))
  predicted_activities_interim <- predicted_activities_without_NO
  
  for(i in 1:length(predicted_activities_interim)) {
    
    predicted_activities_interim[[i]][predicted_non_occurrence_probs[[i]]>=0.5] <- "NO"
  }
  
  predicted_activities <- suppressWarnings(lapply(predicted_distributions, function(x) {if(!is.na(x)) {apply(x, 2, function(x) Mode(x)[1])} else {return(NA)}}))
  predicted_activities_without_NO_probs <- suppressWarnings(lapply(predicted_distributions, function(x) {if(!is.na(x)) {apply(x, 2, function(x) sum(x==Mode(x[x!="NO"])[1])/length(x))} else {return(NA)}}))
  predicted_activities_probs <- suppressWarnings(lapply(predicted_distributions, function(x) {if(!is.na(x)) {apply(x, 2, function(x) sum(x==Mode(x)[1])/length(x))} else {return(NA)}}))
  
  # carry over last feasible forecast to the NA forecasts if we get NA forecasts
  if(sum(is.na(predicted_distributions)) > 0) {
    
    last_feasible_pred_activities <- last(predicted_activities[!is.na(predicted_activities)])
    last_feasible_pred_activities_probs <- last(predicted_activities[!is.na(predicted_activities_probs)])
    
    if(!is.null(last_feasible_pred_activities)) {
      
      predicted_activities[is.na(predicted_activities)] <- lapply(1:sum(is.na(predicted_activities)), function(x) last_feasible_pred_activities)
      predicted_activities_probs[is.na(predicted_activities_probs)] <- lapply(1:sum(is.na(predicted_activities_probs)), function(x) last_feasible_pred_activities_probs)
    }
  }
  
  
  return(list(
    predicted_activities=predicted_activities,
    predicted_activities_probs=predicted_activities_probs,
    predicted_activities_without_NO=predicted_activities_without_NO,
    predicted_activities_without_NO_probs=predicted_activities_without_NO_probs))
}


# helper function for attaching count of activities inside a trace in case of loops
activity_count_helper <- function(activity_vector) {
  
  activity_counts <- table(activity_vector)
  
  current_increments <- as.list(rep(1, length(unique(activity_vector))))
  names(current_increments) <- unique(activity_vector)
  
  for(activity_index in 1:length(activity_vector)) {
    
    if(activity_counts[activity_vector[activity_index]] > 1) {
      
      current_activity <- activity_vector[activity_index]
      
      if(current_increments[[activity_vector[activity_index]]]!=1) {
        
        activity_vector[activity_index] <- paste0(activity_vector[activity_index], "_", current_increments[[activity_vector[activity_index]]])
      }
      
      current_increments[[current_activity]] <- current_increments[[current_activity]] + 1
    }
  }
  
  return(activity_vector)
}

manual_ndls <- function(actual_trace, predicted_trace) {
  
  # actual_trace and predicted_trace are the sequences as lists of character vectors
  
  if(is.null(predicted_trace[[1]])) {
    
    return(NA)
  } else {
    
    distance <- DamerauLevenshtein()(actual_trace, predicted_trace)
    normalized_distance <- distance/max(length(actual_trace[[1]]), length(predicted_trace[[1]]))
    nlds <- 1-normalized_distance
    return(nlds)
  }
}
