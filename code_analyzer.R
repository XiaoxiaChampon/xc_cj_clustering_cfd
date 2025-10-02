library(codetools) # Analyze function bodies
library(igraph)
library(visNetwork) # interactive HTML-based network graphs
library(htmltools) # To Escape HTML in tooltips to prevent rendering issues

get_all_expressions <- function(file_path)
{
  # read file as a list of lines and ignore file not having EOF
  code_lines <- readLines(file_path, warn = FALSE)
  
  # combine lines list into single string
  code_text <- paste(code_lines, collapse = "\n")
  
  # try to parse to identify all expressions
  exprs <- tryCatch(parse(text = code_text), error = function(e) NULL)
  
  exprs_list <- NULL
  # return expressions if its not null
  if (is.null(exprs) == FALSE)
  {
    exprs_list <- as.list(exprs)
  }
  
  return(exprs_list)
}

# === Step 1: Analyze internal dependencies across multiple files ===
analyze_internal_dependencies_multi <- function(file_paths) {
  all_code_lines <- character() # all code from all files
  all_exprs <- list() # parsed expressions
  env <- new.env()
  defined_functions <- character() # all functions
  function_file_map <- list() # maps function name to file it came from
  duplicates <- list() # list of duplicate function names
  
  for (file_path in file_paths) {
    code_lines <- readLines(file_path, warn = FALSE) # read file as a list of lines and ignore file not having EOF
    code_text <- paste(code_lines, collapse = "\n") # combine lines list into single string
    exprs <- tryCatch(parse(text = code_text), error = function(e) NULL)
    if (is.null(exprs)) next
    
    all_code_lines <- c(all_code_lines, code_lines)
    all_exprs <- c(all_exprs, as.list(exprs))
    
    for (expr in exprs) {
      # find all assignment calls where there is a LHS and a RHS
      if (is.call(expr) && (expr[[1]] == "<-" || expr[[1]] == "=")) {
        lhs <- expr[[2]]
        rhs <- expr[[3]]
        # A function definition should have following condition:
        # LHS is a symbol, RHS is a callable that starts with string "function"
        if (is.symbol(lhs) && is.call(rhs) && rhs[[1]] == "function") {
          fname <- as.character(lhs)
          ns_fname <- paste0(basename(file_path), "::", fname)
          
          if (!is.null(function_file_map[[fname]])) {
            duplicates[[fname]] <- unique(c(duplicates[[fname]], function_file_map[[fname]]))
          }
          defined_functions <- c(defined_functions, fname)
          function_file_map[[fname]] <- file_path
          assign(fname, eval(rhs), envir = env)
          assign(ns_fname, eval(rhs), envir = env)
        }
      }
    }
  }
  
  dependency_map <- list()
  unique_functions <- unique(defined_functions)
  for (fname in unique_functions) {
    f <- get(fname, envir = env)
    if (is.function(f)) {
      called <- codetools::findGlobals(f, merge = FALSE)$functions
      internal_calls <- intersect(called, unique_functions)
      dependency_map[[fname]] <- sort(unique(internal_calls))
    }
  }
  
  list(
    dependency_map = dependency_map,
    env = env,
    all_code_lines = all_code_lines,
    function_file_map = function_file_map,
    duplicates = duplicates
  )
}

# === Step 2: Plot using visNetwork with enhanced tooltips ===
plot_interactive_dependency_graph <- function(dep_info) {
  dep_map <- dep_info$dependency_map
  env <- dep_info$env
  all_code_lines <- dep_info$all_code_lines
  function_file_map <- dep_info$function_file_map
  
  functions <- names(dep_map)
  all_called <- unique(unlist(dep_map))
  all_nodes <- unique(c(functions, all_called))
  
  # Enhanced nodes with tooltip
  node_list <- lapply(all_nodes, function(fname) {
    src_file <- function_file_map[[fname]]
    ns_label <- paste0(basename(src_file), "::", fname)
    if (exists(fname, envir = env)) {
      fn <- get(fname, envir = env)
      args <- tryCatch(deparse(args(fn)), error = function(e) "N/A")
      args <- gsub("^function", "", args[1])
      
      body_exprs <- as.list(body(fn))
      ret_expr <- NULL
      return_call <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("return")), body_exprs)
      if (length(return_call) > 0) {
        ret_expr <- deparse(return_call[[1]][[2]])
      } else if (length(body_exprs) > 0) {
        ret_expr <- deparse(tail(body_exprs, 1)[[1]])
      }
      
      src_file <- function_file_map[[fname]]
      doc_line <- grep(paste0(fname, "\\s*<-\\s*function"), all_code_lines)[1] - 1
      doc <- if (!is.na(doc_line) && doc_line > 0 && grepl("^\\s*#'", all_code_lines[doc_line])) {
        sub("^\\s*#'\\s*", "", all_code_lines[doc_line])
      } else {
        ""
      }
      
      
      
      data.frame(
        id = fname,
        label = ns_label,
        shape = "box",
        color = "lightblue",
        title = paste0(
          "<b>Function:</b> ", fname, "<br>",
          "<b>Args:</b> ", htmlEscape(args), "<br>",
          "<b>Returns:</b> ", if (!is.null(ret_expr)) htmlEscape(ret_expr) else "?", "<br>",
          "<b>Source:</b> ", htmlEscape(src_file), "<br>",
          if (nzchar(doc)) paste0("<b>Description:</b> ", htmlEscape(doc), "<br>") else ""
        ),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        id = fname,
        label = ns_label,
        shape = "ellipse",
        color = "gray",
        title = paste0("<b>Function:</b> ", fname),
        stringsAsFactors = FALSE
      )
    }
  })
  
  nodes <- do.call(rbind, node_list)
  nodes <- nodes[!duplicated(nodes$id), ]
  
  edges <- do.call(rbind, lapply(names(dep_map), function(from) {
    to_list <- dep_map[[from]]
    if (length(to_list) == 0) return(NULL)
    data.frame(from = from, to = to_list, arrows = "to", width = 2, stringsAsFactors = FALSE)
  }))
  
  visNetwork(nodes, edges, width = "100%", height = "600px") %>%
    visEdges(smooth = list(enabled = TRUE, type = "dynamic")) %>%
    visNodes(font = list(size = 20)) %>%
    visInteraction(navigationButtons = TRUE, dragNodes = TRUE) %>%
    visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(
        gravitationalConstant = -30,
        centralGravity = 0.02,
        springLength = 120,
        springConstant = 0.05,
        damping = 0.8,
        avoidOverlap = 1
      ),
      stabilization = list(
        enabled = TRUE,
        iterations = 200,
        updateInterval = 20,
        onlyDynamicEdges = FALSE,
        fit = TRUE
      )
    ) %>%
    visEvents(
      dragStart = "function (params) {
    if (params.nodes.length > 0) {
      var id = params.nodes[0];
      this.body.data.nodes.update({id: id, fixed: {x: false, y: false}});
    }
  }",
      
      dragEnd = "function (params) {
    if (params.nodes.length > 0) {
      var id = params.nodes[0];
      this.body.data.nodes.update({id: id, fixed: {x: true, y: true}});
    }
  }"
    )
}

# === Step 3: Example usage ===
files <- c("./R/XCATFDA/R/catfda_experiments.R", "./R/XCATFDA/R/catfda_experiment_functions.R", "./R/XCATFDA/R/catfda_main.R")
dep_info <- analyze_internal_dependencies_multi(files)
print(dep_info$duplicates)  # shows duplicate function names and their original files
plot_interactive_dependency_graph(dep_info)
