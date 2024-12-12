# Function to extract parameter names and values
extract_params <- function(params, group_name = NULL) {
  if (is.atomic(params)) {
    return(list(names = group_name, values = params))
  }
  
  names <- character(length(params))
  values <- vector("list", length(params))
  
  for (i in seq_along(params)) {
    names[i] <- if (is.null(group_name)) names(params)[i] else paste(group_name, names(params)[i], sep = " - ")
    values[[i]] <- params[[i]]
  }
  
  list(names = names, values = unlist(values))
}
