# These packages must be installed once on each computer

#install.packages("ShinyImage")
#install.packages("image.dlib")
#install.packages("magick")
#install.packages("imager")
#install.packages("abind")
#install.packages("FNN")

# To get EBImage, required by ShinyImage
# install.packages("BiocManager")
# BiocManager::install("EBImage")

library(image.dlib)
library(magick)
library(imager)
library(abind)
library(FNN)
library(ShinyImage)

#######################################################################
# NB: File name format is xxx_y_z.png
# xxx is a unique three digit identifier for the frog
# y is a one character tag distinguishing multiple images of the same frog
# z can be any legal string of characters.
# _ are literal
# .png is the file suffix. This is the preferred file type
#
# For example, 003_a_field_2007_winter.png and 003_b_field_2008_spring.png 
# would be two images of frog 003.
#######################################################################

#----------------------------------------------------------------------
# match_targets
# Following Lowe 2004, the algorithm to generate match quality between two images is:
# 1) Using SIFT/SURF identify invariant features in each image. Features are
#    pixel areas that have clear contrast patterns -- for example edges and corners.
#    Features are described with a 64 item value vector
# 2) Using nearest neighbour search, find, for each feature in the target, the
#    most similar features from the standard. Simliarity is effectively the
#    Euclidian distance between features in 64-dimensional space
# 3) For each target feature, consider the 1st and 2nd best matching feature from
#    the standard. If the 1st and 2nd best matching are of similar distance, it
#    is unlikely that the "best" is a real match, as if it were, it should be
#    well distinguished from all other features in the standard. So a "good" match
#    is one where the 1st best matching feature is much better than the second best
#    matching feature. That is, the 1st is "special", and so is likely to be a real match.
#    The ratio argument defines how much closer the 1st match must be relative to the 2nd
#    in order to be considered special and real.
# 4) For each standard compute the proportion of total feature matches that are real. The
#    image with the highest proportion is the putative match.
#----------------------------------------------------------------------
# Standards is always a folder. Target can be a single file or a folder
match_targets<- function(targets, standards, ratio, best_only = TRUE)
{
  # Set up targets
  if (!file.exists(targets)){
    stop("Can't find target frog images")
  }
  
  # default is a single input file
  target_file_vector <- c(targets)
  
  # if it's a directory, reinitialise
  if (dir.exists(targets)){
    target_file_vector <- list.files(targets)
    target_file_vector <- paste(targets, "/", target_file_vector, sep = "")
  }
  
  # Prepare to gather up
  target_name_vector <- c()
  match_std_vector <- c()
  match_qual_vector <- c()
  std_img_vector <- c()
  
  # Compare all, gathering results into vectors
  match_df <- data.frame(TargetID = list(), StdId = list(), MatchQual = list(), StdImg = list())
  name_length <- 5 # 3 digit id, delimiter, 1 char tag. See above.
  
  for (target_file in target_file_vector)
  {
    # All the work happens here. See match_one_target, below
    matches <- match_one_target(target_file, standards, ratio)

    if (best_only) {
      # Matches comes out sorted. Best match is first
      target_name_vector <- c(target_name_vector, matches$TargetID[1])
      match_std_vector <- c(match_std_vector, matches$StdId[1])
      match_qual_vector <- c(match_qual_vector, matches$MatchQual[1])
      std_img_vector <- c(std_img_vector, matches$StdImg[1])
    } else {
      target_name_vector <- c(target_name_vector, matches$TargetID)
      match_std_vector <- c(match_std_vector, matches$StdId)
      match_qual_vector <- c(match_qual_vector, matches$MatchQual)
      std_img_vector <- c(std_img_vector, matches$StdImg)
    }

  }
  
  match_df <- data.frame(Target = target_name_vector, 
                         Standard = match_std_vector, 
                         MatchQual = match_qual_vector,
                         StdImg = std_img_vector)
  return(match_df)
} # end match_targets

#=============================================================


match_accuracy <- function(best_match_multiple_targets)
{
  target_ids <- best_match_multiple_targets$Target
  std_ids <- best_match_multiple_targets$Standard
  
  id_length = 3
  target_ids <- as.numeric(substr(target_ids, 1, id_length))
  std_ids <- as.numeric(substr(std_ids, 1, id_length))
  
  cross_tab <- table(target_ids, std_ids)
  perc_corr <- sum(target_ids == std_ids)/length(target_ids)
  
  accuracy_result <- list(PercCorrIds = perc_corr, ConfMatrix = cross_tab)
  
}
#=============================================================
#=============================================================



#----------------------------------------------------------------------
# Utility methods
#----------------------------------------------------------------------
# Using the python algorithm of "two closest"
# nn_list is the structure returned by get.knnx
good_match_pr <- function(nn_list, ratio = 0.6)
{
  good_points <- c()
  n_point_pairs <- nrow(nn_list$nn.dist)
  
  for (i in 1:n_point_pairs)
  {
    first_distance <- nn_list$nn.dist[i,1]
    second_distance <- nn_list$nn.dist[i,2]
    if (first_distance < ratio * second_distance){
      good_points <- c(good_points, nn_list$nn.index[i,1])
    }
  }
  
  good_match_proportion <- length(good_points)/n_point_pairs
  return(good_match_proportion)
}


#=======================================================================
# One target against a directory of standards.
match_one_target <- function(target_img_file_path, std_img_directory_path, ratio)
{

  name_length <- 5 # 3 digit id, delimiter, 1 char tag. See above.
  match_qual_vector <- c()
  std_id_vector <- c()
  std_img_vector <- c()
  
  # Get target features
  target_file_name <- basename(target_img_file_path)
  target_name <- substr(target_file_name, 1, name_length)
  target_img <- image_read(target_img_file_path)
  
  # Get names of standard image files. Will have to rebuild the path to read
  std_img_file_vector <- list.files(std_img_directory_path)
  
  # For each standard, get features, get best matches, use "first vs second" algorithm to assess quality
  for (std_img_file_name in std_img_file_vector)
  {
    std_name <- substr(std_img_file_name, 1, name_length)
    std_id_vector <- c(std_id_vector, std_name)
    
    full_std_file <- paste(std_img_directory_path, "/", std_img_file_name, sep = "")
    std_img <- image_read(full_std_file)
     
    match_qual <- img_match_quality(target_img, std_img, ratio)
    match_qual_vector <- c(match_qual_vector, match_qual)
    std_img_vector <- c(std_img_vector, full_std_file)
    
  } # end for each standard
  
  target_name_vector <- rep(target_name, length(std_id_vector))
  match_df <- data.frame(TargetID = target_name, StdId = std_id_vector, MatchQual = match_qual_vector, StdImg = std_img_vector)
  # high match proportion = more common features
  match_df <- match_df[order(match_df$MatchQual, decreasing = TRUE),] 
  
  return(match_df)
  
} # end match_target


#------------------------------------------------------------------------------------------
img_match_quality <- function(target_img, std_img, ratio)
{
  target_img_data <- image_data(target_img, channels = "rgb") # contains hex
  target_img_data <- as.integer(target_img_data, transpose = FALSE) # hex converted to decimal
  target_img_surf <- image_surf(target_img_data, max_points = 50)
  
  std_img_data <- image_data(std_img, channels = "rgb") # contains hex
  std_img_data <- as.integer(std_img_data, transpose = FALSE) # hex converted to decimal
  std_img_surf <- image_surf(std_img_data, max_points = 50, detection_threshold = 0)
  
  if (std_img_surf$points >=2 & target_img_surf$points >= 2) {
    
    # standard is data (arg1), target is query (arg2)
    k <- FNN::get.knnx(std_img_surf$surf, target_img_surf$surf, k = 2)
    good_pr <- good_match_pr(k, ratio)
  } else {
    good_pr <- -1 # to indicate the fail
  }
  return(good_pr)
}




