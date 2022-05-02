load_Hmat_SA_h5 <- function(h5_object, run=NULL) {
  
  root = "/H"
  
  if(!is.null(run)) root=paste0("/run", run, "/H")
  
  cn <- h5read(h5_object, paste(root, "axis1", sep="/"))
  rn <- h5read(h5_object,  paste(root, "block0_items", sep="/"))
  
  N.signatures <- length(rn)-2
  
  H <- h5read(h5_object,  paste(root, "block0_values", sep="/"))[1:N.signatures,]
  
  rownames(H) <- rn[1:N.signatures]
  colnames(H) <- cn
  
  htib <- as.data.frame(H) |> 
    rownames_to_column("Signature") |>
    pivot_longer(cols=-Signature, names_to="Tumor_Sample_Barcodes", values_to="Weight")
  
  return(htib)
}

load_Xmat_SA_h5 <- function(h5_object) {
  #from X=WH
  root = "/X"
  
  cn <- h5read(h5_object, paste(root, "axis1", sep="/"))
  rn <- h5read(h5_object,  paste(root, "axis0", sep="/"))
  
  X_0 <- h5read(h5_object,  paste(root, "block0_values", sep="/"))
  X_0_rn <- h5read(h5_object,  paste(root, "block0_items", sep="/"))
  rownames(X_0) <- X_0_rn
  
  X_1 <- h5read(h5_object,  paste(root, "block1_values", sep="/"))
  X_1_rn <- h5read(h5_object,  paste(root, "block1_items", sep="/"))
  rownames(X_1) <- X_1_rn
  
  # check for future debug
  all(rn %in% c(X_0_rn, X_1_rn))
  
  X <- rbind(X_0, X_1)
  colnames(X) <- cn
  
  # for X. tibble object
  xtib <- as.data.frame(X) |> 
    rownames_to_column("Pair") |>
    arrange(Pair) |>
    pivot_longer(cols=-Pair, names_to="Mutation", values_to="Weight")
  
  return(xtib)
}
