strassen_multiplication <- function(A, B) {
  n <- nrow(A)
  # Number of operations: +-* /
  flops <- 0
  
  # base condition
  if (n <= 2) {
    return(A %*% B)
  }
  
  # division into "sub" matrices
  new_size <- n %/% 2
  A11 <- A[1:new_size, 1:new_size]
  A12 <- A[1:new_size, (new_size + 1):n]
  A21 <- A[(new_size + 1):n, 1:new_size]
  A22 <- A[(new_size + 1):n, (new_size + 1):n]
  
  B11 <- B[1:new_size, 1:new_size]
  B12 <- B[1:new_size, (new_size + 1):n]
  B21 <- B[(new_size + 1):n, 1:new_size]
  B22 <- B[(new_size + 1):n, (new_size + 1):n]
  
  # Recursive multiplication of matrices
  P1 <- strassen_multiplication(A11 + A22, B11 + B22)
  P2 <- strassen_multiplication(A21 + A22, B11)
  P3 <- strassen_multiplication(A11, B12 - B22)
  P4 <- strassen_multiplication(A22, B21 - B11)
  P5 <- strassen_multiplication(A11 + A12, B22)
  P6 <- strassen_multiplication(A21 - A11, B11 + B12)
  P7 <- strassen_multiplication(A12 - A22, B21 + B22)
  
  # formula from the lecture
  C11 <- P1 + P4 - P5 + P7
  C12 <- P3 + P5
  C21 <- P2 + P4
  C22 <- P1 - P2 + P3 + P6
  
  # joining results into one matrix
  result <- matrix(0, n, n)
  result[1:new_size, 1:new_size] <- C11
  result[1:new_size, (new_size + 1):n] <- C12
  result[(new_size + 1):n, 1:new_size] <- C21
  result[(new_size + 1):n, (new_size + 1):n] <- C22
  
  return(list(result, flops))
}

# Example usage with two matrices
A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), nrow = 4)
B <- matrix(c(17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32), nrow = 4)

result <- strassen_multiplication(A, B)
print(result)
