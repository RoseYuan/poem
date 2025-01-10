## code to prepare `sp_toys` dataset goes here
# Define the grid dimensions
n_rows <- 16
n_cols <- 15  

simulate_array <- function(n_rows=10, n_cols=10){
  # Create coordinates for the hexagonal grid
  array <- data.frame(row = integer(), col = integer())
  
  for (r in 0:(n_rows - 1)) {
    if (r %% 2 == 0) {
      # Even row: column indices are even (0, 2, 4, ...)
      cols <- seq(0, n_cols*2 - 1, by = 2)
    } else {
      # Odd row: column indices are odd (1, 3, 5, ...)
      cols <- seq(1, n_cols*2 - 1, by = 2)
    }
    # Add to coordinates data frame
    array <- rbind(array, data.frame(row = r, col = cols))
  }
  return(array)
}

simulate_hexagonal_grid <- function(n_rows=10, n_cols=10){
  hex_radius <- 1
  hex_height <- 1.5 * hex_radius
  
  # Initialize vectors to store the x and y coordinates
  x_coords <- numeric(n_rows * n_cols)
  y_coords <- numeric(n_rows * n_cols)
  
  index <- 1
  
  for (row in 1:n_rows) {
    for (col in 1:n_cols) {
      # Calculate x and y for the center of the hexagon
      x <- col * sqrt(3) * hex_radius
      y <- row * hex_height
      
      # Shift alternate rows
      if (row %% 2 == 0) {
        x <- x + sqrt(3)/2
      }
      
      # Store the coordinates
      x_coords[index] <- x
      y_coords[index] <- y
      index <- index + 1
    }
  }
  
  # Return the coordinates as a data frame
  coords <- data.frame(x = x_coords, y = y_coords)
  coords <- cbind(coords, simulate_array(n_rows, n_cols))
  return(coords)
}


sp_toys <- simulate_hexagonal_grid(n_rows, n_cols)
sp_toys$label <- NA
for (i in 1:dim(sp_toys)[1]) {
  if(sp_toys[i, "col"] > sp_toys[i, "row"] + 6){
    sp_toys[i, "label"] <- 1
  }else{
    sp_toys[i, "label"] <- 2
  }
}
sp_toys$label <- factor(sp_toys$label)

sp_toys$p2 <- sp_toys$label
sp_toys[sp_toys$row == 4 & (sp_toys$col %in% c(12,14,16,18,20)), "p2"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 5 & (sp_toys$col %in% c(13,15,17,19)), "p2"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 6 & (sp_toys$col %in% c(14,16,18,20)), "p2"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 7 & (sp_toys$col %in% c(15,17,19)), "p2"] <- factor(2, levels=c(1,2))

sp_toys$p1 <- sp_toys$label
for (i in unique(sp_toys$row)) {
  sp_toys[sp_toys$row==i & sp_toys$col==i+8, "p1"] <- factor(2, levels=c(1,2))
} 

sp_toys$p3 <- sp_toys$label
sp_toys[sp_toys$row == 3 & (sp_toys$col %in% c(21,23)), "p3"] <- factor(2, levels=c(1,2))

sp_toys[sp_toys$row == 4 & (sp_toys$col %in% c(12,14)), "p3"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 5 & (sp_toys$col %in% c(13,15)), "p3"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 6 & (sp_toys$col %in% c(14)), "p3"] <- factor(2, levels=c(1,2))

sp_toys[sp_toys$row == 5 & (sp_toys$col %in% c(21,23)), "p3"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 6 & (sp_toys$col %in% c(22)), "p3"] <- factor(2, levels=c(1,2))

sp_toys[sp_toys$row == 9 & (sp_toys$col %in% c(25)), "p3"] <- factor(2, levels=c(1,2))

sp_toys[sp_toys$row == 12 & (sp_toys$col %in% c(20,26)), "p3"] <- factor(2, levels=c(1,2))

sp_toys[sp_toys$row == 10 & (sp_toys$col %in% c(18,20)), "p3"] <- factor(2, levels=c(1,2))
sp_toys[sp_toys$row == 11 & (sp_toys$col %in% c(19)), "p3"] <- factor(2, levels=c(1,2))

sp_toys$label2 <- sp_toys$p2
sp_toys$p4 <- sp_toys$label


sp_toys$p5 <- sp_toys$p2
sp_toys[sp_toys$row == 4 & (sp_toys$col %in% c(20)), "p5"] <- factor(1, levels=c(1,2))
sp_toys[sp_toys$row == 5 & (sp_toys$col %in% c(19)), "p5"] <- factor(1, levels=c(1,2))
sp_toys[sp_toys$row == 6 & (sp_toys$col %in% c(20)), "p5"] <- factor(1, levels=c(1,2))
sp_toys[sp_toys$row == 7 & (sp_toys$col %in% c(19)), "p5"] <- factor(1, levels=c(1,2))

for (i in unique(sp_toys$row)) {
  if(i !=4 & i!=5 & i!=6 & i!=7){
    sp_toys[sp_toys$row==i & sp_toys$col==i+6, "p5"] <- factor(1, levels=c(1,2))
  }
} 

usethis::use_data(sp_toys, overwrite = TRUE)
