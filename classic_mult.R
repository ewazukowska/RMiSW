
classic_multiplication <- function(mtr1, mtr2) {
  # Number of operations: +-* /
  flops <- 0
  
  sum <- 0
  
  # The shape of the result matrix will be: <number of rows from the first matrix> x <number of columns from the second matrix>
  result_matrix <- matrix(0, nrow = nrow(mtr1), ncol = ncol(mtr2))
  
  # Number of columns in the first matrix = number of rows in the second matrix
  # (this variable will help us with the most inner iteration)
  size <- ncol(mtr1)
  
  # Iterate over rows of the first matrix
  for (i in 1:nrow(mtr1)) {
    
    # Iterate over columns of the second matrix
    for (j in 1:ncol(mtr2)) {
      
      # "k" denotes columns in the first matrix and rows in the second matrix
      for (k in 1:size) {
        
        sum <- sum + mtr1[i, k] * mtr2[k, j]
        flops <- flops + 2
        
      }
      
      result_matrix[i, j] <- sum
      sum <- 0
      
    }
  }
  
  # Return a list containing the result matrix and the number of operations
  return(list(result_matrix, flops))
}

# Example usage with two matrices
# Macierze M1 i M2
M1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
M2 <- matrix(c(1, 2, 3, 4, 5, 6, 1, 2, 1), nrow = 3, ncol = 3)

# Display matrices
print(M1)
print(M2)


result <- classic_multiplication(M1, M2)
print(result[[1]])  # Result matrix
print(result[[2]])  # Number of operations


