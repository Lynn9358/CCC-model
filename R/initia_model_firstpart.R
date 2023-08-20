## Null simulation data generating

df <- list()
df$meta <- list()
df$meta$cell <- c(1:500)
df$meta$x <- runif(500, min = 0, max = 100)
df$meta$y <- runif(500, min = 0, max = 100)
df$meta$gene1 <- runif(500, min = 0, max = 1)
df$meta$gene2 <- runif(500, min = 0, max = 1)
df$meta$gene3 <- runif(500, min = 0, max = 1)
df$meta$gene4 <- runif(500, min = 0, max = 1)
df$meta$locat <- data.frame(df$meta$x,df$meta$y)
df$para$dist <- as.matrix(dist(df$meta$locat))


##cell type decomposing
celltype_dec <- function(obj) {
  obj$meta$cell_type <- max.col(as.matrix(cbind(obj$meta$gene1, obj$meta$gene2, obj$meta$gene3, obj$meta$gene4)))
  return(obj)
}
df <- celltype_dec(df)

##plot different cell
plot(df$meta$locat, pch=19,cex = 0.5, col = ifelse(df$meta$cell_type %in% c(3,4), "grey", df$meta$cell_type))
table(df$meta$cell_type)


## Weighted K-NN for lr pairs
get_lr_wknn <- function(obj,sender,receiver, k){
  
  if(!all(c(sender, receiver)) %in% unique(obj$meta$cell_type)) {
    stop("Can not find sender/receiver in the data")}
  else {
    
    distance <- obj$para$dist[obj$meta$cell_type == sender, obj$meta$cell_type == receiver]
    receiver <- list()
    
    receiver$indices <- as.matrix(t(head(apply(distance, 1, function(x)  as.numeric(colnames(distance)[order(x,decreasing = FALSE)])),k)))
    
    receiver$distances <- t(head(apply(distance, 1, function(x) sort(x, decreasing = FALSE)), k))
    
    weighted_counts <- receiver$distances * matrix( c(1,sqrt(2),sqrt(3)) , nrow = nrow(receiver$distances), ncol = ncol(receiver$distances), byrow = TRUE)
    
    percentile_weight <- head(sort(weighted_counts), length(weighted_counts) * 0.5)
    index_wknn <- t(apply(weighted_counts, 1 , function(x) x %in% percentile_weight))
    
    lr_wknn  <- receiver$indices * index_wknn
    lr_wknn[lr_wknn == 0 ] <- NA
    
    return(lr_wknn)
  }
}


lr_example <- get_lr_wknn(df,
                          sender = 1,
                          receiver = 2,
                          k = 3)

lr_example

## create lr index
create_lr_index <- function(lrpair){
  A = rep(row.names(lrpair), each = ncol(lrpair))
  B = as.vector(t(lrpair))
  lr_index <- cbind(A,B)
  
  lr_index <- lr_index[complete.cases(lr_index),]
  return(lr_index)
}


lr_index_example <- create_lr_index(lr_example)


##  plot lr pairs 
find_row_list_number <- function(list_of_lists, lr_index) {
  for (i in 1: length(list_of_lists)) {
    if (all(as.numeric(lr_index) %in% list_of_lists[[i]])) {
      ind = i 
      break
    }
  }
  return(ind)
}



plot_lr_pair <- function(obj, lr_pair, arrows = T, cluster = list()){
  
  cell_index <- unique(as.numeric(c(rownames(lr_pair),as.numeric(lr_pair))))
  
  lr_index <- create_lr_index(lr_pair)
  
  cell_meta <- obj$meta$locat[unique(lr_index),]
  
  plot(df$meta$locat, pch=19,cex = 0.5, col = ifelse(obj$meta$cell_type %in% c(3,4), "grey", df$meta$cell_type))
  
  text(df$meta$x[cell_index],df$meta$y[cell_index], cell_index, pos = 3)
  
  
  for (i in 1:nrow(lr_index)) {
    A = lr_index[i,1]
    B = lr_index[i,2]
    x1 <- cell_meta[rownames(cell_meta) == A, 1]
    y1 <- cell_meta[rownames(cell_meta) == A, 2]
    x2 <- cell_meta[rownames(cell_meta) == B, 1]
    y2 <- cell_meta[rownames(cell_meta) == B, 2]
    
    
    if (length(cluster) == 0 ){
      if (arrows)  {
        arrows(x1, y1, x2, y2, length = 0.1)
      }
      else{
        lines(c(x1, x2), c(y1, y2), col = 1 )
      }
      
    }
    else{
      
      if (arrows)  {
        arrows(x1, y1, x2, y2, length = 0.1, col =as.numeric(find_row_list_number(cluster, lr_index[i, ])))
      }
      else{
        lines(c(x1, x2), c(y1, y2), col = as.numeric(find_row_list_number(cluster, lr_index[i, ])))
      }
    }
    
  }
}




#example  


plot_lr_pair(df, lr_example,)



## model for points 
find_next_connected_point <- function(this_point, point_column_index, lr_index) {
  opposite_column_index = 3 - point_column_index
  cluster_points = list(this_point)
  if (is.null(dim(lr_index))) return (unique(cluster_points))
  connected_points <- lr_index[lr_index[, point_column_index] == this_point, opposite_column_index]
  if (length(connected_points)>0) {
    for (connect_point in connected_points){
      if (is.null(dim(lr_index))) break
      lr_index = lr_index[! ( (lr_index[, point_column_index] == this_point) & (lr_index[, opposite_column_index] == connect_point) ), ]
      next_points = find_next_connected_point(connect_point, opposite_column_index, lr_index)
      cluster_points = rbind(cluster_points, next_points)
    }
  } 
  return (unique(cluster_points))
}

delete_points_from_df <- function(points, lr_index) {
  for (point in points){
    if (is.null(dim(lr_index))) break
    lr_index <- lr_index[(lr_index[, 1] != point) & (lr_index[, 2] != point), ]
  }
  return (lr_index)
}



find_point <- function (lr_index) {
  
  lr_index = rbind(lr_index, c(-1,-1))
  network = list()
  network_index= 1
  
  while (!is.null(dim((lr_index)))){
    current_point = lr_index[1,1]
    cluster_points = unique(find_next_connected_point(current_point, 1, lr_index))
    network[[network_index]] = cluster_points
    network_index = network_index + 1
    lr_index = delete_points_from_df(cluster_points, lr_index)
  }
  return(network)
}


#create network

create_lr_network <- function(obj, lr_pair, plot = TRUE) {
  
  
  lr_index <- create_lr_index(lr_pair)
  
  #To be con
  point_within <- find_point(lr_index)
  network = list()
  
  for( i in 1:length(point_within) ){
    network[[i]] <- lr_index[lr_index[,1] %in% point_within[[i]] & lr_index[,2] %in% point_within[[i]],]
    
  }
  
  if (plot == TRUE) {
    plot_lr_pair(df, lr_example, arrows = F, cluster = network )
  }
  
  
  return(network)
  
  
}


CLN <- create_lr_network(df, lr_example, plot = FALSE)


##dataset create ***
build_dataset <- function(df, lr) {
  lr_network <- create_lr_network(df, lr, plot = FALSE)
  dataset <- vector("list", length(lr_network))  
  
  # Initialize dataset as a list with the required number of elements
  
  for (j in 1:length(lr_network)) {
    single_network <- lr_network[[j]]
    
    for (i in 1:(length(single_network)/2)) {
      
      if (length(single_network)/2 < 2) {
        ligand <- as.numeric(single_network[[1]])
        receptor <- as.numeric(single_network[[2]])
      } else {
        row <- single_network[i, ]
        ligand <- as.numeric(row[1])
        receptor <- as.numeric(row[2])
      }
      
      gene_class_L <- df$meta$cell_type[ligand]
      gene_class_R <- df$meta$cell_type[receptor]
      
      expression <- cbind(df$meta$gene1, df$meta$gene2, df$meta$gene3, df$meta$gene4)
      x <- expression[ligand, gene_class_L]
      y <- expression[receptor, gene_class_R]
      
      new_expression_row <- c(x, y)
      dataset[[j]] <- rbind(dataset[[j]], new_expression_row)
    }
  }
  
  return(dataset)
}

df_example <- build_dataset(df, lr_example)

