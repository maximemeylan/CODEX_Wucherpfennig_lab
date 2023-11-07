## CODEX utils functions
## Maxime Meylan 2023
## Wucherpfennig Lab - Dana-Farber Cancer Institute

library(raster)
library(collapse)
library(matrixStats)
options(future.globals.maxSize = 15 * 1024 ^ 3)
myClip <- function(x, a, b){
  ifelse(x <= a,  a, ifelse(x >= b, b, x))
}
#wrapper function:
# Prepare dataframe from seurat object selected by a sample id (from orig.ident)
# Capture and parse output of shiny app
# Outputs list of selected cells 
interactive_points_keep <- function(seurat_object,sample_id){
  cell_ids <- names(which(seurat_object$orig.ident == sample_id))
  df<- as.data.frame(cbind(seurat_object$x[cell_ids], seurat_object$y[cell_ids]))
  colnames(df) <- c('x','y')
  cells_to_keep <- suppressWarnings(capture.output(interactive_points_selection(df, x = ~x, y = ~y,sample_id = sample_id) ))
  cells_to_keep <- paste0(cells_to_keep, collapse = '')
  cells_to_keep <- unlist(strsplit(cells_to_keep,','))
  return(cells_to_keep)
}

# Updated function for points selection 
interactive_points_keep_V2 <- function(seurat_object){
  df<- as.data.frame(cbind(seurat_object$x, seurat_object$y))
  colnames(df) <- c('x','y')
  cells_to_keep <- suppressWarnings(capture.output(interactive_points_selection(df, x = ~x, y = ~y,sample_id = unique(seurat_object$orig.ident)) ))
  cells_to_keep <- paste0(cells_to_keep, collapse = '')
  cells_to_keep <- unlist(strsplit(cells_to_keep,','))
  return(cells_to_keep)
}

# Updated function for points selection 
interactive_points_keep_downsample <- function(seurat_object,sample_id=NULL,n_cells=30000,pixel_to_mm=0.000325,return_type=c("seurat","ids")[1]){
  #Get the cell ids
  if(is.null(sample_id)){
  sample_id <- unique(seurat_object$orig.ident)
  }
  cell_ids <- colnames(seurat_object)[seurat_object$orig.ident == sample_id]
  # Downsample dataset
  n_cells <- ifelse(length(cell_ids) < n_cells, length(cell_ids),n_cells)
  cell_ids <- sample(cell_ids,n_cells,replace=F)
  df<- as.data.frame(cbind(seurat_object$x[cell_ids], seurat_object$y[cell_ids]))
  colnames(df) <- c('x','y')
  # Select boundaries of tissue using the downsampled dataset
  cells_to_keep <- suppressWarnings(capture.output(interactive_points_selection(df, x = ~x, y = ~y,sample_id = sample_id) ))
  cells_to_keep <- paste0(cells_to_keep, collapse = '')
  cells_to_keep <- unlist(strsplit(cells_to_keep,','))
  all_pts <- seurat_object@meta.data[,c("x",'y')]
  
  # Draw polygons of the boundaries
  points_d <- as.matrix(all_pts[cells_to_keep,])
  polygons <- concaveman(points_d)
  poly <- st_polygon(list(polygons))
  colnames(polygons) <- colnames(points_d)
  
  ## convert points to sf format
  points_sf <- st_as_sf(all_pts, coords = c("x", "y"))
  # Find which points are contained within the polygon
  contained <- st_contains(poly, points_sf)
  # Extract the indices of the points that are contained
  contained_indices <- unlist(contained)
  # Extract the points that are contained
  contained_points <- rownames(all_pts)[contained_indices]
  
  if(return_type=="seurat"){
    #compute surface of final shape
    area <- round(st_area(poly*pixel_to_mm),digits = 2)
    all_pts_mm <- all_pts*pixel_to_mm
    all_pts_mm$alpha <- 0.05
    all_pts_mm[contained_points,"alpha"] <- 0.5
    polygons_mm <- polygons*pixel_to_mm
    p <- ggplot(data = data.frame(all_pts_mm), aes(x=x, y=y)) +
      geom_point(size=0.01,alpha=all_pts_mm$alpha)+
      geom_polygon(data = as.data.frame(polygons_mm), aes(x = x, y = y), fill = "transparent", color = "red")+
      guides(alpha=F)+ coord_fixed() + ggtitle(bquote(.(sample_id) ~ " surface: " ~ .(area) ~ " mm"^2)) + theme_classic() + NoLegend() + coord_fixed()
    plot(p)
    subset_seurat <- subset(seurat_object,cells = unique(unlist(contained_points)))
    subset_seurat@misc$area <- area
    subset_seurat@misc$tissue_boundaries <-list(polygons)
    return(subset_seurat)
  }else{
    return(contained_points)
  }
}

#function for interactive selection
#take a dataframe containing x and y coordinates and id as rownames
#returns list of ids
interactive_points_selection <- function(dat, x, y, key = row.names(dat),sample_id) {
  title_ <- paste0("Interactive selection of points for sample: ",sample_id)
  ui <- miniPage(
    gadgetTitleBar(title_),
    plotlyOutput("plot1", height = "100%")
  )
  
  server <- function(input, output, session) {
    
    # mechanism for managing selected points
    keys <- reactiveVal()
    
    observeEvent(event_data("plotly_selected"), {
      key_new <- event_data("plotly_selected")$key
      key_old <- keys()
      keys(c(key_new, key_old))
    })
    
    output$plot1 <- renderPlotly({
      is_outlier <- key %in% keys()
      cols <- ifelse(is_outlier, "red", "grey90")
      dat %>%
        plot_ly(x = x, y = y) %>%
        add_markers(key = row.names(dat), color = I(cols), marker = list(size = 5))%>%
        plotly::layout(dragmode='lasso')%>%
        toWebGL()
    })
    # Return the most recent fitted model, when we press "done"
    observeEvent(input$done, {
      cat(keys(),fill = T,sep=',')
      stopApp()
    })
  }
  shinyApp(ui, server)
}



#Predict cluster labels using pre-computed NN between two umaps
# Need test umap object with nn between test and train umap
# method = 'vote' or 'dist'
predict_labels_nn <- function(umap_test,train_labels,method){
  ## code with barely any adaption from github.com/davpinto/fastknn
  k=ncol(umap_test$nn[[1]]$idx)
  #### Compute class membership probabilities
  label.mat <- matrix(train_labels[umap_test$nn[[1]]$idx], ncol = k)
  knn.prob <- switch(
    method,
    ## P(y_j | x_i) = sum(1/d(nn_i) * (y(nn_i) == y_j)) / sum(1/d(nn_i))
    'dist' = {
      sapply(levels(train_labels), function(cl, d, y) {
        rowSums(1/d * (y == cl)) / rowSums(1/d)
      }, d = pmax(umap_test$nn[[1]]$dist, 1e-15), y = label.mat, 
      simplify=FALSE, USE.NAMES=TRUE)
    },
    ## P(y_j | x_i) = sum(y(nn_i) == y_j) / k
    'vote' = {
      sapply(levels(train_labels), function(cl, y) {
        rowSums(y == cl) / ncol(y)
      }, y = label.mat, simplify=FALSE, USE.NAMES=TRUE)
    }
  )
  knn.prob <- as.matrix(do.call('cbind.data.frame', knn.prob))
  knn.prob <- sweep(knn.prob, 1, rowSums(knn.prob), "/")
  
  #### Assign class labels
  knn.label <- levels(train_labels)[max.col(knn.prob, ties.method = "first")]
  knn.label <- factor(knn.label, levels(train_labels))
  return(knn.label)
} 

#Refine labels to better match umap coordinates
# Takes a seurat object as input
# Returns a seurat object with a "label_refined" metadata
# For each cell and each cluster, measure distance to centroid
# Identify outliers 
# find non outiliers knn for each outlier
# Assign new label based on mode 
refine_labels <- function(sobj,label_name = "leiden",k_neighbors=200,n_threads=1,Kmean_k=3){
  labels <- setNames(sobj@meta.data[,label_name],colnames(sobj))
  #measure distance to centroid of each cluster
  for( x in levels(labels)){
    all_coords <- sobj@reductions$umap@cell.embeddings[which(labels == x),]
    #Find biggest spatial cluster
    km <- kmeans(all_coords,Kmean_k)
    big_clus <- which.max(table(km$cluster))
    centroid <- km$centers[big_clus,]
    #Compute distance to centroid
    dists <- raster::pointDistance(all_coords,centroid,lonlat = F)
    #Remove labels for cells that are far from the centroid
    ids <- dists > sd(dists)
    labels[names(which(ids))] <- NA
  }
  #define k closest neighbors between "outlier" cells and "non-outlier" cells
  knn_res <- dbscan::kNN(x = sobj@reductions$umap@cell.embeddings[which(!is.na(labels)),],
                         query = sobj@reductions$umap@cell.embeddings[which(is.na(labels)),],
                         k = k_neighbors)
  names_no_na <- rownames(sobj@reductions$umap@cell.embeddings[which(!is.na(labels)),])[knn_res$id]
  names_na <- rownames(sobj@reductions$umap@cell.embeddings[which(is.na(labels)),])
  #create a matrix with cluster labels instead of idx
  label_mat <- matrix(labels[names_no_na],nrow = nrow(knn_res$id),ncol = ncol(knn_res$id))
  #take the mode to assign refined labels
  new_labels <- setNames(fmode(t(label_mat),nthreads = n_threads),names_na)
  #Fill back refined labels
  labels[names(new_labels)] <- new_labels
  sobj<- AddMetaData(sobj,labels,col.name = paste0(label_name,"_refined"))
}

find_neighborhoods <- function(sobj,k,n_neighborhoods,cell_pop_var="cell_annot",group_var="cell_population"){
  # initiate neighborhoods metadata 
  sobj$neighborhoods <- NA
  t0<- Sys.time()
  region_ids <- unique(sobj@meta.data[,group_var])
  #For each region
  neihgbors_list <- lapply(region_ids,function(r_id){
    print(paste0('Processing: ', r_id))
    ids <- sobj@meta.data[,group_var]==r_id
    #Extract coordinates from seurat object
    coords<- data.frame(cbind(sobj$x[ids],
                              sobj$y[ids]))
    #Find k nearest neighbors
    knn_res <- kNN(coords,k = k,sort = T)
    #Add cell annotation (must be after knn)
    coords[,cell_pop_var] <- sobj@meta.data[ids,cell_pop_var]
    neighb <- data.frame(knn_res$id)
    neighb$id <- rownames(knn_res$id)
    #Melt to long format for vectorized annotation
    neighb_m <-  suppressMessages(melt(neighb,varnames = 'id'))
    neighb_m$annot <- coords[neighb_m$value,cell_pop_var]
    #Compute frequencies  
    freq <- table(neighb_m$id,neighb_m$annot)
    freq <- matrix(freq, ncol = ncol(freq), dimnames = dimnames(freq))
    #Account for absent cell population
    if(ncol(freq) != length(unique(sobj@meta.data[,cell_pop_var]))){
      missing_cat <- setdiff(unique(sobj@meta.data[,cell_pop_var]),colnames(freq))
      missing_mat <- matrix(rep(0,nrow(freq)), ncol = length(missing_cat))
      colnames(missing_mat) <- missing_cat
      freq <- cbind(freq,missing_mat)
    }
    return(freq)
  })
  #Aggregate regions
  neihgbors_agg<- do.call(rbind,neihgbors_list)
  #Add nearest neighbords to the seurat object
  sobj@misc <- list('nearest_neihgbors'=neihgbors_agg)
  print(paste0('Performing neighborhoods identification on ', length(region_ids),' regions'))
  #perform clutering to identify neighborhoods on frequency tables
  clust <- kmeans(x = neihgbors_agg,centers = n_neighborhoods)
  #Add neighborhoods information to seurat object
  sobj$neighborhoods[names(clust$cluster)] <- as.character(paste0('CN-',clust$cluster))
  t1<- Sys.time()
  print(t1-t0)
  return(sobj)
}

find_neighborhoods_V2 <- function(sobj,k,n_neighborhoods,cell_pop_var="cell_annot",pop_to_exclude="NULL",group_var="orig.ident"){
  # initiate neighborhoods metadata 
  sobj$neighborhoods <- NA
  t0<- Sys.time()
  region_ids <- unique(sobj@meta.data[,group_var])
  #For each region
  neihgbors_list <- lapply(region_ids,function(r_id){
    print(paste0('Processing: ', r_id))
    ids <- sobj@meta.data[,group_var]==r_id
    #Extract coordinates from seurat object
    coords<- data.frame(cbind(sobj$x[ids],
                              sobj$y[ids]))
    #Find k nearest neighbors
    knn_res <- kNN(coords,k = k,sort = T)
    #Add cell annotation (must be after knn)
    coords[,cell_pop_var] <- sobj@meta.data[ids,cell_pop_var]
    neighb <- data.frame(knn_res$id)
    neighb$id <- rownames(knn_res$id)
    #Melt to long format for vectorized annotation
    neighb_m <-  suppressMessages(melt(neighb,varnames = 'id'))
    neighb_m$annot <- coords[neighb_m$value,cell_pop_var]
    #Compute frequencies  
    freq <- table(neighb_m$id,neighb_m$annot)
    freq <- matrix(freq, ncol = ncol(freq), dimnames = dimnames(freq))
    #Account for absent cell population
    if(ncol(freq) != length(unique(sobj@meta.data[,cell_pop_var]))){
      missing_cat <- setdiff(unique(sobj@meta.data[,cell_pop_var]),colnames(freq))
      missing_mat <- matrix(0,nrow = nrow(freq),ncol= length(missing_cat))
      colnames(missing_mat) <- missing_cat
      freq <- cbind(freq,missing_mat)
    }
    #Remove unwanted population
    if(!is.null(pop_to_exclude)){
      pop_to_keep <- setdiff(colnames(freq),pop_to_exclude)
      cell_to_keep <- rownames(coords)[!coords[,3] %in% pop_to_exclude]
      freq <- (freq[cell_to_keep,c(pop_to_keep)] / rowSums(freq[cell_to_keep,c(pop_to_keep)],na.rm = T)*100)
      freq[is.na(freq)] <- 0
    }
    return(freq)
  })
  #Aggregate regions
  neihgbors_agg<- do.call(rbind,neihgbors_list)
  #Add nearest neighbords to the seurat object
  sobj@misc <- list('nearest_neihgbors'=neihgbors_agg)
  print(paste0('Performing neighborhoods identification on ', length(region_ids),' regions'))
  #perform clutering to identify neighborhoods on frequency tables
  clust <- kmeans(x = neihgbors_agg,centers = n_neighborhoods)
  #Add neighborhoods information to seurat object
  sobj$neighborhoods[names(clust$cluster)] <- as.character(paste0('CN-',clust$cluster))
  sobj$neighborhoods[is.na(sobj$neighborhoods)] <- "none"
  t1<- Sys.time()
  print(t1-t0)
  return(sobj)
}
# (Re-Annotation with XGBoost)
ReAnnoX <- function(seurat_object,
                    mode=c("train","predict")[1],
                    cells_to_classify=NULL, #cell labels, if NULL, uses all the cells
                    markers=NULL, #channel names, if NULL uses all the channels
                    assay="spatial", #Assay to fetch MFI values
                    label, #variable used to train the model, must be in the metadata of the seurat object
                    sampling_mode=c('dowsample','upsample','random')[1], 
                    n_random_sampling=1e5,
                    XGboost_model=NULL, # input pre-training model
                    training_data=NULL, # must be supplied if running predict
                    train_test_split=0.95, #percentage of total cells to assign to training (should be <1)
                    nrounds=500, #number of rounds to train xgboost
                    confidence_level=0.8, #probability over which the prediction is kept, otherwise it is assigned "low_confidence"
                    n_threads=1){
  set.seed(8)
  if(class(seurat_object) != "Seurat"){
    stop(paste0("seurat_object input is not a Seurat object"))
  }
  if(!assay %in% names(seurat_object@assays)){
    stop(paste0("Assay \"",assay,"\" was not found in the seurat object \n Please specify on the assay ", paste0(names(seurat_object@assays),collapse = " or ")))
  }
  # ids of all cells except cells to classify
  if(is.null(cells_to_classify)){
    ids <- colnames(seurat_object)
  }else{
    ids <- setdiff(colnames(seurat_object),cells_to_classify)
  }
  if(is.null(markers)){
    markers <- rownames(seurat_object@assays[[assay]])
  }
  switch(mode,
         train = {
           # get train and test datasets
           switch(sampling_mode,
                  downsample = {
                    print("performing downsampling")
                    downsampled_data <- downSample(x =  t(seurat_object@assays[[assay]]@counts[markers,ids]),
                                                   y = factor(seurat_object@meta.data[ids,label]),
                                                   yname = label, list = TRUE)
                    downsampled_data$x[,label] <- downsampled_data$y
                    # Convert factors to numeric indices for XGBoost
                    downsampled_data$x[,'label_num'] <- as.numeric(downsampled_data$y) - 1
                    original_levels <- levels(downsampled_data$y)
                    trainIndex <- createDataPartition(downsampled_data$y, p = train_test_split, list = FALSE, times = 1)
                    # Convert dataframes to DMatrix objects
                    dtrain <- xgb.DMatrix(data = as.matrix(downsampled_data$x[ trainIndex,markers]), label = downsampled_data$x[ trainIndex,"label_num"])
                    dtest <- xgb.DMatrix(data = as.matrix( downsampled_data$x[-trainIndex,markers]))
                    training_data <- list("train_labels"=factor(downsampled_data$x[ trainIndex,label]),
                                          "train_nums"= downsampled_data$x[ trainIndex,"label_num"],
                                          "test_labels"=factor(downsampled_data$x[ -trainIndex,label]),
                                          "test_nums"= downsampled_data$x[ -trainIndex,"label_num"])
                    rm(downsampled_data)
                    gc()
                  },
                  upsample = {
                    print("performing upsampling")
                    upsampled_data <- upSample(x =  t(seurat_object@assays[[assay]]@counts[markers,ids]), y = factor(seurat_object@meta.data[ids,label]), yname = label, list = TRUE)
                    upsampled_data$x[,label] <- upsampled_data$y
                    # Convert factors to numeric indices for XGBoost
                    upsampled_data$x[,'label_num'] <- as.numeric(upsampled_data$y) - 1
                    original_levels <- levels(upsampled_data$y)
                    trainIndex <- createDataPartition(upsampled_data$y, p = train_test_split, list = FALSE, times = 1)
                    # Convert dataframes to DMatrix objects
                    dtrain <- xgb.DMatrix(data = as.matrix(upsampled_data$x[ trainIndex,markers]), label = upsampled_data$x[ trainIndex,"label_num"])
                    dtest <- xgb.DMatrix(data = as.matrix( upsampled_data$x[-trainIndex,markers]))
                    training_data <- list("train_labels"=factor(upsampled_data$x[ trainIndex,label]),
                                          "train_nums"= upsampled_data$x[ trainIndex,"label_num"],
                                          "test_labels"=factor(upsampled_data$x[ -trainIndex,label]),
                                          "test_nums"= upsampled_data$x[ -trainIndex,"label_num"])
                    rm(upsampled_data)
                    gc()
                  },
                  random = {
                    print("performing random sampling")
                    sampled_ids <- sample(ids,n_random_sampling) 
                    train_ids <- sampled_ids %>% sample(n_random_sampling*train_test_split)
                    test_ids <- setdiff(sampled_ids,train_ids)
                    # Convert factors to numeric indices for XGBoost
                    labels_num <- as.numeric(factor(seurat_object@meta.data[train_ids,label])) - 1
                    original_levels <- levels(factor(seurat_object@meta.data[train_ids,label]))
                    # Convert dataframes to DMatrix objects
                    dtrain <- xgb.DMatrix(data = as.matrix(t(seurat_object@assays[[assay]]@counts[markers,train_ids])), label =labels_num)
                    dtest <- xgb.DMatrix(data = as.matrix(t(seurat_object@assays[[assay]]@counts[markers,test_ids])))
                    training_data <- list("train_labels"=factor(seurat_object@meta.data[train_ids,label]),
                                          "train_nums"=  as.numeric(factor(seurat_object@meta.data[train_ids,label])) - 1,
                                          "test_labels"=factor(seurat_object@meta.data[test_ids,label]),
                                          "test_nums"=  as.numeric(factor(seurat_object@meta.data[test_ids,label])) - 1)
                  }
           )
           if(is.null(XGboost_model)){
             # Set XGBoost parameters
             params <- list(
               objective = "multi:softprob",
               num_class = length(original_levels),
               "eval_metric" = "mlogloss",
               "eta" = 0.03,
               nthread = n_threads)
             # Train the model
             print("Cells sampled for training:")
             print(table(training_data$train_labels))
             print("Cells sampled for testing:")
             print(table(training_data$test_labels))
             print(paste0("training model with ",ncol(dtrain)," markers and ",nrow(dtrain), " cells on ",length(original_levels), " cell types" ))
             XGboost_model <- xgb.train(params, dtrain,nrounds = nrounds)
           }
           
           # Make predictions
           predictions_xgb <- predict(XGboost_model, dtest)
           
           # Extract predictions probs
           pred_prob <- matrix(predictions_xgb, ncol = length(original_levels), byrow = T) 
           
           # Get the class with max probability for each observation
           max_probs <-  rowMaxs(pred_prob)
           max_pred <- apply(pred_prob,1,which.max)
           pred_label <- factor(max_pred, levels = 1:length(original_levels), labels = original_levels)
           
           # If max probability is below threshold, replace prediction with NA
           pred_label[max_probs < confidence_level] <- NA
           
           # Calculate the confusion matrix
           cm_xg <- confusionMatrix(pred_label,training_data$test_labels)
           print(cm_xg$overall)
           print(cm_xg$byClass)
           return(list("XGboost_model"=XGboost_model,
                       "confusion_matrix"=cm_xg,
                       "training_data"=training_data))
         })
  predict = {
    #original_levels <- levels(factor(seurat_object@meta.data[cells_to_classify,label]))
    # Convert dataframes to DMatrix objects
    dtest <- xgb.DMatrix(data = as.matrix(t(seurat_object@assays[[assay]]@counts[markers,cells_to_classify])))
    
    # Make predictions
    predictions_xgb <- predict(XGboost_model, dtest)
    
    # Extract predictions probs
    pred_prob <- matrix(predictions_xgb, ncol = length(unique(training_data$train_labels)), byrow = T) 
    
    # Get the class with max probability for each observation
    max_probs <-  rowMaxs(pred_prob)
    max_pred <- apply(pred_prob,1,which.max)
    pred_label <- factor(max_pred, levels = 1:length(levels(training_data$train_labels)), labels = levels(training_data$train_labels))
    
    # If max probability is below threshold, replace prediction with NA
    pred_label[max_probs < confidence_level] <- NA
    return(list("predicted_labels"=pred_label))
  }
}
