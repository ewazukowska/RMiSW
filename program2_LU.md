# 







# 



# 




# <center> Eliminacja Gaussa i Faktoryzacja LU </center>

---


<center> Rachunek Macierzowy i Statystyka Wielowymiarowa </center>

---
# 



# 


# 



# 



# 


# 


# 

# 


# 
# 

**Ewa Pelc**  
**Ewa Żukowska** 

***
# 
## Faktoryzacja LU




Rozmiar macierzy: 20 + 5 = 25

```{r}
n <- 25
A <- matrix(runif(n^2), nrow = n)

```
### Pseudokod algorytmu LU

1. Zainicjuj 𝐋 jako macierz identycznościową, 𝐈 o wymiarach 𝑛×𝑛 oraz 𝐔 = 𝐀.
2. Dla 𝑖 = 1 ,…, 𝑛 wykonaj krok 3.
3. Dla 𝑗 = 𝑖 + 1, …, 𝑛 wykonaj kroki 4-5.
4. Ustaw 𝑙𝑗𝑖=𝑢𝑗𝑖/𝑢𝑖𝑖.
5. Wykonaj 𝑈𝑗=(𝑈𝑗 − 𝑙𝑗𝑖 * 𝑈𝑖) (gdzie 𝑈𝑖, 𝑈𝑗 reprezentują odpowiednio wiersze 𝑖 oraz 𝑗 macierzy 𝐔).


### Algorytm faktoryzacji LU:

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
#### Sprawdzenie poprawności LU faktoryzacji

W celu sprawdzenia poprawności faktoryzacji zdefiniowałyśmy funkcję, która porównuje dwa obiekty dopuszczając ustalony błąd. 

```{r}
is_allclose <- function(a, b, tol = 1e-10) {
  max_diff <- max(abs(a - b))
  return(max_diff < tol)
}
```
#### Sprawdzenie, czy A = LU

```{r}
if (is_allclose(A, L %*% U)) {
  print("LU faktoryzacja jest poprawna.")
} else {
  print("LU faktoryzacja jest niepoprawna.")
}
```
### Pseudokod algorytmu LU z pivotingiem

1. Zainicjuj macierze 𝐋 = 𝐏 = 𝐈 o wymiarach 𝑛×𝑛 oraz 𝐔 = 𝐀.
2. Dla każdego 𝑖 = 1, …, 𝑛 wykonaj kroki 3 - 4, 8.
3. Niech 𝑘 = 𝑖. 
4. Dopóki 𝑢𝑖𝑖 = 0, wykonuj kroki 5-7.
5. Zamień miejscami wiersz 𝑈𝑖 z wierszem 𝑈𝑘+1.
6. Zamień miejscami wiersz 𝑃𝑖 z wierszem 𝑃𝑘+1.
7. Zwiększ wartość 𝑘 o 1.
8. Dla każdego 𝑗 = 𝑖+1, …, 𝑛 wykonaj kroki 9 - 10.
9. Ustaw wartość 𝑙𝑗𝑖 = 𝑢𝑗𝑖/𝑢𝑖𝑖.
Wykonaj operację 𝑈𝑗 = 𝑈𝑗 − 𝑙𝑗𝑖 * 𝑈𝑖 (gdzie 𝑈𝑖,𝑈𝑗 reprezentują odpowiednio 𝑖-ty i 𝑗-ty wiersz macierzy 𝐔).

![](diff1.png)

W trakcie obliczeń numerycznych może wystąpić utrata precyzji, co prowadzi do drobnych różnic w wynikach. 

Należy jednak zwrócić 
### Algorytm faktoryzacji LU z pivotingiem:

```{r}

LU_factorization_with_pivoting <- function(A) {
  n <- nrow(A)
  L <- diag(n)  # Macierz L inicjujemy jako macierz jednostkową
  P <- diag(n)  # Macierz permutacji P inicjujemy jako macierz jednostkową
  U <- A        # Macierz U inicjujemy jako macierz A
  
  for (i in (1:(n-1)) ){
    k <- i
    
    # Szukanie elementu maksymalnego w kolumnie i
    while (U[i, i] == 0) {
      if (k > n) {
        return("Macierz jest zdegenerowana.")
      }
      U[c(k, k+1), ] <- U[c(k+1, k), ]
      P[c(k, k+1), ] <- P[c(k+1, k), ]
      k <- k + 1
    }
    
    # Eliminate entries below i with row operations on U 
    # and reverse the row operations to manipulate L
    factor <- U[(i+1):n, i] / U[i, i]
    L[(i+1):n, i] <- factor
    U[(i+1):n, ] <- U[(i+1):n, ] - outer(factor, U[i, ])
  }
  
  return(list(P = P, L = L, U = U))
}

```

```{r}
LU_decomposition_with_pivoting <- LU_factorization_with_pivoting(A)
L <- LU_decomposition_with_pivoting$L
U <- LU_decomposition_with_pivoting$U
```

#### Sprawdzenie, czy A = PLU:

```{r}
if (is_allclose(A, P %*% L %*% U)) {
  print("LU faktoryzacja z pivotingiem jest poprawna.")
} else {
  print("LU faktoryzacja z pivotingiem jest niepoprawna.")
}

```

Poniższa macierz przedstawia różnicę, między macierzą A, a zaproponowanym przez nas rozkładem LU. Jak widać, na niektórych pozycjach wartości są różne, jednak różnica ta jest rzędu 10^16.

Funkcje wbudowane mogą stosować bardziej zaawansowane strategie wyboru elementu głównego w kroku faktoryzacji, co może prowadzić do dokładniejszych wyników.

![](diff.png)

# Eliminacja Gaussa


Macierz rozmiaru dzień urodzenia + miesiąc urodzenia

```{r}
date <- 28 + 6
birthday_matrix <- matrix(sample(1:9, date*date, replace=TRUE), nrow=date)
```

Wektor wyrazów wolnych do układu równań
```{r}
v <- t(matrix(runif(date, min=-10, max=20), nrow=1))
```

### Pseudokod algorytmu eliminacji Gaussa generującego jedynki

1. Dla każdego wiersza 'row' od 1 do (n - 1):
- a. coef_div = m_copy[row, row]
- b. Dla każdego wiersza 'j' od (row + 1) do n:
    - i. coef_sub = m_copy[j, row]
    - ii. m_copy[j, ] = m_copy[j, ] - m_copy[row, ] * coef_sub

2. Jeśli ostatni wiersz nie jest znormalizowany (m_copy[n, n] != 1):
  - a. m_copy[n, ] = m_copy[n, ] / m_copy[n, n] 

3. Zwróć macierz (bez ostatniej kolumny) i ostatnią kolumnę (wektor b) z m_copy.



### Algorytm eliminacji Gaussa generujący jedynki na przekątnej

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
![alt text](ep_images/gauss_elimination.png)

### Pseudokod algorytmu eliminacji Gaussa z pivotingiem

1. Dla każdej kolumny 'i' od 1 do (n - 1):
- a. Znajdź pivot jako maksymalną wartość w kolumnie 'i' dla wierszy od (i + 1) do n.
- b. Znajdź indeks wiersza z maksymalnym pivotem: max_row_index = indeks wiersza z maksymalnym pivotem + i.

    - Jeśli pivot > m_copy[i, i]:
          Zamień miejscami wiersze 'i' i 'max_row_index' w macierzy 'm_copy'.

    - Dla każdego wiersza 'j' od (i + 1) do n:
      - i. Oblicz współczynnik coef = m_copy[j, i] / m_copy[i, i].
      - ii. Odjęcie wiersza 'i' pomnożonego przez coef od wiersza 'j'.

2. Zwróć macierz (bez ostatniej kolumny) i ostatnią kolumnę (wektor b) z m_copy.


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

![alt text](ep_images/gauss_pivot.png)


### Pseudokod funkcji solve()

1. Wykonaj eliminację Gaussa na macierzy 'm' i wektorze 'b' za pomocą funkcji 'func', która zwraca macierz współczynników i wektor 'b' po eliminacji.
2. Przypisz macierz współczynników do 'coef_matrix' i wektor 'b' do 'b_vector'.
3. Określ liczbę równań 'n' na podstawie liczby wierszy macierzy 'coef_matrix'.
4. Zainicjuj wektor wynikowy 'x' jako pusty wektor numeryczny o długości 'n'.

5. Wsteczne podstawianie:
  - a. Dla każdego i od n do 1, malejąco:
    - i. Oblicz x[i] jako iloraz b_vector[i] i coef_matrix[i, i].
    - ii. Jeśli i < n, oblicz x[i] jako różnicę b_vector[i] i sumy iloczynów elementów coef_matrix[i, (i + 1):n] i x[(i + 1):n], podzieloną przez coef_matrix[i, i].

6. Zwróć wektor wynikowy 'x'.


### Funkcja solve() rozwiązująca układ równań
```{r}
solve <- function(m, b, func) {
  # Gauss elimination
  result <- func(m, b)
  coef_matrix <- result$matrix
  b_vector <- result$vector
  print(coef_matrix)
  
  n <- nrow(coef_matrix)
  x <- numeric(n)

  # Wsteczne podstawianie
  for (i in seq(n, 1, -1)) {
    
    x[i] <- b_vector[i]/coef_matrix[i, i]
    if (i < n) {
      x[i] <- (b_vector[i] - sum(coef_matrix[i, (i + 1):n] * x[(i + 1):n])) / coef_matrix[i, i]
    }
  }
  
  return(x)
}
```

Obliczenie wyników korzystając z obu algorytmów eliminacji Gaussa

```
# calculate results for Gauss and Gauss with Pivot
x6 <- solve(birthday_matrix, v, gauss)
x7 <- solve(birthday_matrix, v, gauss_with_pivot)
```

Sprawdzenie czy rozwiązania są identyczne dla obu metod eliminacji

```{r}
identical(x6, x7)
```

W celu sprawdzenia poprawności (czy oryginalna macierz przemnożona przez wektor X rozwiązań jest równa wektorowi b) ponownie wykorzystujemy funkcję *is_allclose()*, która dopuszcza ustalony błąd 

```{r}
# check if correct - Gauss
r1 <- birthday_matrix %*% x6
is_allclose(v,r1)
```

```{r}
# check if correct - Gauss with pivot
r2 <- birthday_matrix %*% x7
is_allclose(v,r2)
```

![alt text](ep_images/check_gauss.png)
![alt text](ep_images/check_pivot1.png)
