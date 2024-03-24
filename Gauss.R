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
  
  # Normalizacja ostatniego wiersza, jeśli to konieczne
  if (m_copy[n, n] != 1) {
    m_copy[n, ] <- m_copy[n, ] / m_copy[n, n]
  }
  
  # Return the matrix and b vector separately
  return(list(matrix = m_copy[, -ncol(m_copy)], vector = m_copy[, ncol(m_copy)]))
}



# Przykładowa macierz 4x4
m <- matrix(c(3, 2, 3, 4,
              4, 5, 6, 6,
              1, 2, 3, 7,
              1, 2, 3, 1), nrow = 4, byrow = TRUE)

# Przykładowy wektor

b <- matrix(c(10, 11, 12, 3), ncol = 1)

# Wywołanie funkcji gauss
result <- gauss(m, b)

# Wyświetlenie wyników
print("Macierz współczynników:")
print(result$matrix)

print("Wektor b:")
print(result$vector)



gauss_with_pivot <- function(m, b) {
  m_copy <- m
  n <- nrow(m)
  
  m_copy <- cbind(m_copy, b)
  
  for (i in 1:(n - 1)) {
    pivot <- max(m_copy[(i + 1):n, i])
    max_row_index <- which.max(m_copy[(i + 1):n, i] == pivot) + i
    
    if (pivot > m_copy[i, i]) {
      temp_row <- m_copy[i, ]
      m_copy[i, ] <- m_copy[max_row_index, ]
      m_copy[max_row_index, ] <- temp_row
    }
    
    for (j in (i + 1):n) {
      coef <- m_copy[j, i] / m_copy[i, i]
      m_copy[j, ] <- m_copy[j, ] - m_copy[i, ] * coef
    }
  }
  
  # Return the matrix and b vector separately
  return(list(matrix = m_copy[, -ncol(m_copy)], vector = m_copy[, ncol(m_copy)]))
}





# Przykładowa macierz 4x4
m <- matrix(c(3, 2, 3, 4,
              4, 5, 6, 6,
              1, 2, 3, 7,
              1, 2, 3, 1), nrow = 4, byrow = TRUE)

# Przykładowy wektor
b <- c(10, 11, 12, 3)

# Wywołanie funkcji gauss
result2 <- gauss_with_pivot(m, b)

# Wyświetlenie wyników
print("Macierz współczynników:")
print(result2$matrix)

print("Wektor b:")
print(result2$vector)





solve <- function(m, b, func) {
  # Gauss elimination
  result <- func(m, b)
  coef_matrix <- result$matrix
  b_vector <- result$vector
  
  n <- nrow(coef_matrix)
  x <- numeric(n)

  # Wsteczne podstawianie
  for (i in seq(n-1, 1, -1) {
    x[i] <- (b_vector[i] - sum(coef_matrix[i, (i + 1):n] * x[(i + 1):n])) / coef_matrix[i, i]
  }
  
  
  return(x)
}


#----------------------------

# Przykładowa macierz 4x4
m <- matrix(c(3, 2, 3, 4,
              4, 5, 6, 6,
              1, 2, 3, 7,
              1, 2, 3, 1), nrow = 4, byrow = TRUE)

# Przykładowy wektor
b <- matrix(c(10, 11, 12, 3), ncol = 1)

# Użycie funkcji solve
result3 <- solve(m, b, gauss_with_pivot)
print(result3)


