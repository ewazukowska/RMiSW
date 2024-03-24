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

