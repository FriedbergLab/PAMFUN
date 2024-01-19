## ---------------------------
## Purpose of script: Helper functions to download kegg information
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------

# format KEGG output into list object.
format_kegg <- function(.x){
			attr_ind <- grep("^[A-Z]", .x)
			attrs <- sapply(.x[attr_ind], function(.y){
				strsplit(.y, split = " ")[[1]][1]
				#str_split_1(.y, "\\s{3,}")[[1]]
				}) %>% as.character()

			clean_attrs <- janitor::make_clean_names(attrs)
			res <- list()
			for(i in 1:(length(attrs)-1)){
				term <- attrs[[i]][1]
				replace <- paste0(rep(" ", nchar(term)), collapse = "")
				res[[ clean_attrs[i] ]] <- gsub(term, replace, .x[attr_ind[i]:(attr_ind[i+1]-1)])
			}
			return(res)
		}

# loop through EC queries and collect information.
multi_kegg <- function(.query, pattern = ""){
	is_null <- tryCatch({
		page <- read_html(.query) %>%
		html_text() %>%
		str_split("\n") %>%
		unlist()
	}, error = function(e) {
		if (grepl("404", e$message)) {
		return(NULL)
		} else {
		stop(e)
		}
	})
	if(is.null(is_null)){
		return(NULL)
	}

	start <- which(grepl("^ENTRY", page) == TRUE)
	end <- which(grepl("^///", page) == TRUE)

	page_info <- lapply(1:length(start), function(.x){
			temp <- page[start[.x]:end[.x]]
			format_kegg(temp)
		}) 
	page_names <- lapply(page_info, function(.x){
		regmatches(.x$entry[1], regexpr(pattern , .x$entry[1]))
	})
	names(page_info) <- unlist(page_names)
	return(page_info)
}

fetch_kegg_info <- function(.x){
	message(Sys.time(), " ", .x, " pathway query")
	
	# Helper function to fetch page and return data table or NULL
	fetch_page <- function(url, colnames) {
		page <- tryCatch({
			rvest::read_html(url) %>%
				rvest::html_text() %>%
				read.table(text = ., sep = "\t") %>%
				data.table::as.data.table()
		}, error = function(e) NULL)
		if(!is.null(page)) {
			colnames(page) <- colnames
		}
		return(page)
	}
	
	# get the genome information
	genome_page <- fetch_page(paste("https://rest.kegg.jp/link/genome/", .x, sep = ""), c("kegg_protein", "org_code"))

	# KEGG genome id to ncbi protein accessions
	protein_page <- fetch_page(paste("https://rest.kegg.jp/conv/ncbi-proteinid/", .x, sep = ""), c("kegg_protein", "ncbi_protein"))
	
	# list of modules 
	module_page <- fetch_page(paste("https://rest.kegg.jp/link/module/", .x, sep = ""), c("kegg_protein", "module_id"))
	
	# list of pathways
	pathway_page <- fetch_page(paste("https://rest.kegg.jp/link/pathway/", .x, sep = ""), c("kegg_protein", "pathway"))

	# link to protein orthologs
	ortholog_page <- fetch_page(paste("https://rest.kegg.jp/link/ko/", .x, sep = ""), c("kegg_protein", "ortholog"))

	# read organism entry page to only get complete modules present
	complete_page <- tryCatch({
		read_html(paste("https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_modules&org=",  .x, sep = "")) %>%
			html_text() %>%
			str_extract_all(pattern = "M[0-9]+")
	}, error = function(e) NULL)
	res <- list(genome_page, protein_page, module_page, pathway_page, ortholog_page, complete_page)
	names(res) <- c("genome", "protein", "module", "pathway", "ortholog", "complete")
	Sys.sleep(0.5)
	return(res)
}

retry_kegg <- function(input, max_attempts = 5, pause_time = 1) {
  for (i in 1:max_attempts) {
    #message(i)
	result <- tryCatch(
      expr = fetch_kegg_info(input),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      message("error: attempt ", i)
	  Sys.sleep(pause_time)
    } else {
		#message("not error")
		break
    }
  }
  Sys.sleep(pause_time)
  return(result)
}

# function to parse a module definition
remove_first_last_parentheses <- function(string) {
  if (substr(string, 1, 1) == "(" && substr(string, nchar(string), nchar(string)) == ")") {
    open_parentheses <- 0
    for (i in 2:(nchar(string) - 1)) {
      char <- substr(string, i, i)
      if (char == "(") {
        open_parentheses <- open_parentheses + 1
      } else if (char == ")") {
        open_parentheses <- open_parentheses - 1
        if (open_parentheses < 0) {
          return(string)
        }
      }
    }
    if (open_parentheses == 0) {
      return(substr(string, 2, nchar(string) - 1))
    }
  }
  return(string)
}

split_string <- function(string, sep) {
  split <- strsplit(string, "")[[1]]
  result <- c()
  temp <- ""
  open_parenthesis <- 0
  for (char in split) {
    if (char == "(") {
      open_parenthesis <- open_parenthesis + 1
    } else if (char == ")") {
      open_parenthesis <- open_parenthesis - 1
    }
    if (char == sep && open_parenthesis == 0) {
      result <- c(result, remove_first_last_parentheses(temp))
      temp <- ""
    } else {
      temp <- paste0(temp, char)
    }
  }
  result <- c(result, remove_first_last_parentheses(temp))
  return(result)
}

space_separate <- function(.x){split_string(.x, sep = ",") %>% as.list()}


# calculate the number of possible paths given a KEGG module definition
permute_steps <- function(.x){
	.y <- sapply(.x, split_string, sep = ",", simplify = FALSE) %>% unlist(recursive = FALSE) %>% unname()
	while(sum(.y != .x) > 0){
	.x <- .y
	.y <- sapply(.x, split_string, sep = ",", simplify = FALSE) %>% unlist(recursive = FALSE) %>% unname()
	}
	res <- sapply(.y, function(.z){if(grepl("\\(|\\)", .z)){
		.a <- split_string(.z, sep = " ") %>% as.list()
		.b <- lapply(.a, permute_steps)  %>% expand.grid() %>% apply(1, function(.x){paste(.x, collapse = " ")})
		return(.b)
		}else{return(.z)}}) %>% unlist(recursive = FALSE) %>% unname()
	return(res)
}

decipher <- function(.x){
	.a <- gsub("-\\(.*?\\)|-K\\d{5}", "", .x)
	.a <- gsub("\\+", " ", .a)
	a_steps <- split_string(.a, sep = " ") %>% as.list()
	a_combos <- lapply(a_steps, permute_steps) %>% expand.grid() %>% apply(1, function(.z){paste(.z, collapse = " ")})

  	.b <- gsub("-", " ", .x)
	.b <- gsub("\\+", " ", .b)
	b_steps <- split_string(.b, sep = " ") %>% as.list()
	b_combos <- lapply(b_steps, permute_steps) %>% expand.grid() %>% apply(1, function(.z){paste(.z, collapse = " ")})

	res <- unique(c(a_combos, b_combos))
	return(res)
}


separate_string <- function(x) {
  regex_ec <- "([0-9]+)\\.([0-9]+)\\.([0-9]+)\\.([0-9\\-]+)"
  regex_id <- "[K]\\d{5}"
  regex_rn <- "[R]\\d{5}"

  ec <- str_extract_all(x, regex_ec)[[1]]
  if(length(ec) == 0) {
	ec <- NA
  } else {
	ec[length(ec) == 0] <- NA
  }

  id <- str_extract_all(x, regex_id)[[1]]
  id <- gsub(","," ",id, perl = TRUE)
  res <- expand.grid(id, ec)
  colnames(res) <- c("kegg_id", "ec_number")
  return(res)
}

# to identify non-orthologous replacements
decipher2 <- function(.x){
        .a <- gsub("-\\(.*?\\)|-K\\d{5}", "", .x)
        .a <- gsub("\\+", " ", .a)
        a_steps <- split_string(.a, sep = " ") %>% as.list()
        a_combos <- lapply(a_steps, permute_steps)

        .b <- gsub("-", " ", .x)
        .b <- gsub("\\+", " ", .b)
        b_steps <- split_string(.b, sep = " ") %>% as.list()
        b_combos <- lapply(b_steps, permute_steps) 

        res <- unique(list(a_combos, b_combos))
        return(res)
}