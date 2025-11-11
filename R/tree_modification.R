# R/tree_modification.R

#' Adds LUCA-MRCA branch to APE tree
#'
#' @param intree ape (phylo) tree
#' @param branch_length length of the LUCA-MRCA branch (either this or total_depth required).
#' @param total_depth total final tree depth (used to calculate branch_length, required if branch_length is not specified)
#'
#' @returns phylo tree
#' @export

add_cenancestor_tree <- function(intree,branch_length=NULL,total_depth=NULL){
  if(!methods::is(intree,"phylo"))
    stop("ERROR: intree must be an ape phylogenetic tree (phylo class)")
  if(is.null(branch_length) && is.null(total_depth))
    stop("ERROR: either branch_length or total_depth required")

  if(!is.null(total_depth)){
    mrca_depth <- max(ape::node.depth.edgelength(intree))
    branch_length <- total_depth - mrca_depth
    if(branch_length < .phyfumr_env[['precision']])
      stop("ERROR: desired branch length impossible to achieve with a positive LUCA-MRCA branch length")
  }

  outtree <- list()
  outtree$tip.label <- intree$tip.label
  outtree$Nnode <- intree$Nnode + 1 #Add one internal node
  #Pre-pend the luca_branch length to the edge.lengths
  outtree$edge.length <- c(branch_length,intree$edge.length)
  #Add a first edge from LUCA to MRCA and increase the id of all internal nodes + 1
  n_leaves <- length(intree$tip.label)
  new_edge <- intree$edge
  new_edge[!new_edge %in% seq.int(to=n_leaves)] <- new_edge[!new_edge %in% seq.int(to=n_leaves)] + 1
  outtree$edge <- rbind(c(n_leaves+1,n_leaves+2), #LUCA - MRCA edge
                        new_edge ) #Add 1 to all except the leaves
  class(outtree) <- "phylo"
  return(outtree)
}

#' Adds LUCA-MRCA branch to NEXUS content
#'
#' @param trees_content NEXUS tree file content read using readLines
#' @param ptable .log file with a column with the luca_branch information
#' @param n_digits number of decimal digits used to print the branch length
#'
#' @returns NEXUS tree file content with the LUCA-MRCA branch length added
#' @keywords internal

add_cenancestor_nexus_content <- function(trees_content,ptable,n_digits=16){
  sprintf_format <- paste0(":%.",n_digits,"f);")

  #Find the lines with trees
  line_trees <- which(grepl("^\\s*tree\\s", trees_content, ignore.case = TRUE, perl = TRUE))
  offset <- line_trees[1]-1 #to translate tree_content line to ptable row
  n_trees <- length(line_trees)
  if(n_trees == 0)
    stop("ERROR: problem reading the .trees file, trees not found")

  #ptable sanity check
  if (nrow(ptable) < n_trees) {
    stop("ERROR: ptable has fewer rows than tree lines; check .log/.trees alignment")
  }

  #main loop to modify the trees, adding a ( at the start and :bl) before the ;
  final_trees_content <- trees_content
  for (tree_content_line in line_trees) {
    branch_length <- ptable$luca_branch[tree_content_line-offset]

    # Add the branch length at the end of the tree
    temp_line <- sub(";\\s*$", sprintf(sprintf_format, branch_length), trees_content[tree_content_line], perl = TRUE)

    # Add one extra '(' at the first opening parenthesis
    final_line <- sub(" (", " ((", temp_line, fixed = TRUE)

    final_trees_content[tree_content_line] <- final_line
  }

  return(final_trees_content)
}
