gauss <- function(m, b) {
  m_copy <- m
  n <- nrow(m)
  
  m_copy <- cbind(m_copy, b)
  
  for (row in 1:(n - 1)) {
    coef_div <- m_copy[row, row]
    m_copy[row, ] <- m_copy[row, ] / coef_div
    
    for (j in (row + 1):n) {
      coef_sub <- m_copy[j, row]
      m_copy[j, ] <- m_copy[j, ] - m_copy[row, ] * coef_sub
    }
  }
  
  # return the matrix and b vector separately
  return(list(matrix = m_copy[, -ncol(m_copy)], vector = m_copy[, ncol(m_copy)]))
}





gauss_with_pivot <- function(m, b) {
  m_copy <- m
  n <- nrow(m)
  
  m_copy <- cbind(m_copy, b)
  
  for (i in 1:n) {
    if (i < (n - 1)) {
      pivot <- max(m_copy[(i + 1):n, i])
      max_row_index <- which(m_copy[(i + 1):n, i] == pivot) + i
      
      if (pivot > m_copy[i, i]) {
        temp_row <- m_copy[i, ]
        m_copy[i, ] <- m_copy[max_row_index, ]
        m_copy[max_row_index, ] <- temp_row
      }
    }
    
    for (j in (i + 1):n) {
      coef <- m_copy[j, i] / m_copy[i, i]
      m_copy[j, ] <- m_copy[j, ] - m_copy[i, ] * coef
    }
  }
  
  # return the matrix and b vector separately
  return(list(matrix = m_copy[, -ncol(m_copy)], vector = m_copy[, ncol(m_copy)]))
}




solve <- function(m, b, func) {
  # Wywołanie funkcji gauss_with_pivot lub gauss
  result <- func(m, b)
  coef_matrix <- result$matrix
  b_vector <- result$vector
  
  # Wsteczne podstawianie
  n <- length(b_vector)
  x <- numeric(n)
  
  for (i in (n - 1):1) {
    x[i] <- (b_vector[i] - sum(coef_matrix[i, (i + 1):n] * x[(i + 1):n])) / coef_matrix[i, i]
  }
  
  return(x)
}




m <- matrix(c(3, 2, 3, 4,
              4, 5, 6, 6,
              1, 2, 3, 7,
              1, 2, 3, 1), nrow = 4, byrow = TRUE)


new_column <- c(10, 11, 12, 3)




result <- gauss(m, new_column)
coef_matrix1 <- result$matrix
b_vector1 <- result$vector




# Przykładowa macierz 4x4
m <- matrix(c(3, 2, 3, 4,
              4, 5, 6, 6,
              1, 2, 3, 7,
              1, 2, 3, 1), nrow = 4, byrow = TRUE)

# Przykładowy wektor
b <- c(10, 11, 12, 3)

# Wywołanie funkcji gauss
result <- gauss(m, b)

# Wyświetlenie wyników
print("Macierz współczynników:")
print(result$matrix)

print("Wektor b:")
print(result$vector)