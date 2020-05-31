# sherali_adams
A small library of functions with few dependencies to run k rounds of Sherali-Adams relaxation on a system of the form Ax <= b represented by raw numpy matrix A and array b.
## notes
The constraints 0 <= y <= 1 will be added for all variables in the system returned. Send me a note if you'd like for this to be optional. 
## install
```
pip install sherali_adams
```
### Testing
```> nosetests``` 
```----------------------------------------------------------------------
Ran 5 tests in 0.423s

OK
```
### Examples

1. Run 1 round of SA on a system with 2 variables.  
  ```A = np.matrix([1,1],[1,1])```
   ```b = np.matrix([1,1])```
   ```(AA,bb) = run_SA(1,2,A,b)```
2. Find the original monomial corresponding to new variables. Here we find it for y_3 in the result above after 1 round of SA on a system with 2 variables.(0 based indexing)
  ``` monomial = invert(2,2,1)```
  ``` monomial == [0,1]```
  Thus y_3 = y_{0,1} which delinearizes to x_0x_1 
3. Supports dynamic programming/memoized mode and brute force. To use memoized:
  ```run_SA(k = 2,n = 5,A = A,b = b,memoize = True)```   
