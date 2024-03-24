# RMiSW

Rozmiar macierzy: 20 + 5 = 25

```{r}
n <- 3
A <- matrix(runif(n^2), nrow = n)

```

Algorytm faktoryzacji LU:

```{r}
LU_factorization <- function(A) {
  n <- nrow(A)
  L <- matrix(0, nrow = n, ncol = n)
  U <- matrix(0, nrow = n, ncol = n)
  print(A)
  
  for (i in 1:n) {
    for (j in i:n) {
      U[i, j] <- A[i,j] - L[i,1:(i-1)] %*% U[1:(i-1),j]
    }
    for (j in i:n) {
      L[j,i] <- (A[j,i] - L[j,1:(i-1)] %*% U[1:(i-1),i]) / U[i,i]
    }
    L[i, i] <- 1
  }
  
  return(list(L = L, U = U))
}
```


```{r}
LU_decomposition <- LU_factorization(A)
L <- LU_decomposition$L
U <- LU_decomposition$U
```
```{r}
LU_factorization(A)
```
```{r}
L %*% U
```

```{r}
# Sprawdzenie poprawności LU faktoryzacji
is_allclose <- function(a, b, tol = 1e-10) {
  max_diff <- max(abs(a - b))
  return(max_diff < tol)
}

# Sprawdzenie, czy A = LU
if (is_allclose(A, L %*% U)) {
  print("LU faktoryzacja jest poprawna.")
} else {
  print("LU faktoryzacja jest niepoprawna.")
}
```


Algorytm faktoryzacji LU z pivotingiem:

```{r}
LU_factorization_pivot <- function(A, epsilon = 1e-10) {
  n <- nrow(A)
  md <- 1
  W <- 1:n
  L <- diag(1, nrow = n)
  U <- matrix(0, nrow = n, ncol = n)
  
  for (k in 1:(n-1)) {
    maxw <- k
    maxe <- abs(A[W[k], k])
    
    for (i in (k+1):n) {
      if (abs(A[W[i], k]) > maxe) {
        maxw <- i
        maxe <- abs(A[W[i], k])
      }
    }
    
    if (maxe <= epsilon) {
      return(NULL)  # Macierz jest zdegenerowana
    }
    
    if (maxw != W[k]) {
      md <- -md
      temp <- W[k]
      W[k] <- W[maxw]
      W[maxw] <- temp
    }
    
    U[k, ] <- A[W[k], ]
    
    for (i in (k+1):n) {
      L[i, k] <- A[W[i], k] / A[W[k], k]
      U[i, ] <- A[W[i], ] - L[i, k] * U[k, ]
    }
  }
  
  U[n, ] <- A[W[n], ]
  
  return(list(L = L, U = U))
}


```

```{r}
LU_decomposition_with_pivoting <- LU_factorization_pivot(A)
L <- LU_decomposition_with_pivoting$L
U <- LU_decomposition_with_pivoting$U
```

```{r}
LU_factorization_pivot(A)
```


```{r}
# Sprawdzenie poprawności LU faktoryzacji
is_allclose <- function(a, b, tol = 1e-10) {
  max_diff <- max(abs(a - b))
  return(max_diff < tol)
}

# Sprawdzenie, czy A = LU
if (is_allclose(A, L %*% U)) {
  print("LU faktoryzacja jest poprawna.")
} else {
  print("LU faktoryzacja jest niepoprawna.")
}
```

```{r}

LU_decomposition <- LU_factorization_with_pivoting(A)
L <- LU_decomposition$L
U <- LU_decomposition$U
P <- LU_decomposition$P

# Sprawdzenie, czy A = PLU
if (is_allclose(A, P %*% L %*% U)) {
  print("LU faktoryzacja z pivotingiem jest poprawna.")
} else {
  print("LU faktoryzacja z pivotingiem jest niepoprawna.")
}
```

```{r}
L %*% U
```

```{r}
A
```

# Eliminacja Gaussa

### Algorytm eliminacji Gaussa bez pivotingu generujący jedynki na przekątnej

```{r}
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
  
  # Zwrócenie macierzy i wektora b
  return(list(matrix = m_copy[, -ncol(m_copy)], vector = m_copy[, ncol(m_copy)]))
}
```

### Algorytm eliminacji Gaussa z pivotingiem

```{r}
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
```

```{r}
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
```

```{r}
date <- 28 + 6
birthday_matrix <- matrix(sample(1:9, date*date, replace=TRUE), nrow=34)
v <- t(matrix(runif(34, min=-100, max=100), nrow=1))
```

```
# calculate results for Gauss and Gauss with Pivot
x6 <- solve(birthday_matrix, v, gauss)
print(x6)

x7 <- solve(birthday_matrix, v, gauss_with_pivot)
print(x7)
```

```{r}
identical(x6, x7)
```


```{r}
# check if correct - Gauss
r1 <- birthday_matrix %*% x6
identical(r1, v)
v - r1
```

```{r}
# check if correct - Gauss with pivot
r2 <- birthday_matrix %*% x7
identical(r2, v)
v - r2
```