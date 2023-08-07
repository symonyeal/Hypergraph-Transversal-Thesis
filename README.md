# ```fk_a``` + Visualization
This repository contains the source code for the implementation of [Fredman-Khachiyan Algorithm- A](https://www.semanticscholar.org/paper/On-the-Complexity-of-Dualization-of-Monotone-Normal-Fredman-Khachiyan/cf09761a7a863915f91346881df95484b5bee617?p2df) to detect if two given hypergraphs are transversals of each other or not. Furthermore, the code includes a visualization tool that charts the recursive instances of the algorithm.

## Usage
## Code Explanation
### Input
For the implementation of the FK-A algorithm, the main libraries used here are Numpy and Pandas. The inputs are two binary incidence matrices of two hypergraphs, where each row represents a hyper-edge and each column represents a variable.  
### Pre-processing Steps
Before the recursive step, the algorithm performs several pre-processing checks to verify that the given hypergraphs satify the conditions of two hypergraphs being transversal to each other in that recursion. 
* **Duplicate Rows Removal**: The function Unique(f, g) helps us to remove duplicate rows in the matrices. Numpy library has a build in function called _unique_, that will return the unique elements in the array.
Since our input is n-dimensional array, and we only want to remove the duplicate rows, we use axis = 0 to tell the function to look at the rows only: ```np.unique(np f, axis = 0)```.
* **Superset Removal** : For sets _A_ and _B_, _A_ is a superset of _B_ if all elements of _B_ are in set _A_. In FK_A, the function takes a matrix as an input where the row _A_ is a superset of row _B_ if all 1s in row _B_ are also 1s in row _A_. The function operates by sorting the array by descending order of the sums of the rows so that the first row will have the least number of 0s, then iterates through the array for each row, it identifies the columns that have 0 and checks if there exists another row with matching 0 columns. If such a row exists, then current row is identified as a superset and removed.
*  _The algorithm only records the columns with 0 value instead of 1 because if the current row has a column value 1, regardless of values on the 0 or 1 in the corresponding 1 column, the current row is considered a super set of the corresponding row if the algorithm only considered the 1 columns._
*  The input array is denoted as ```np arr```. First, it sorts the array with:  ```np_arr``` = ```sorted(np_arr, key= lambda x: sum(x), reverse=True)``` . ```Sorted``` sorts the rows in an array in a specific order, since the default order is ascending, set **reverse= True**, to return the result in descending order.
  ```key=lambda x: sum(x)``` means that the sums of each row are used as a key for the sort comparison. Next, ```deleted_rows``` stores the index of the superset rows that needs to be removed. Lastly, it deletes all rows based on the index in ```delete_rows``` and returs both ```np_arr``` and ```delete_rows```.
### Pre-conditions
After the duplicate and superset removals, FK-A checks the inputs against a set of conditions based on transversal relationships before decomposing and recursing. If the inputs fail any of the conditions, the algorithm stops and returns FALSE. The function ```check_pre``` tests the following pre-conditions:
* _f_ and _g_ must contain the same set of variables.
* Largest row sum in _f_ must be at most the number of rows in _g_, or vice versa.
* Let r<sub>s</sub> be the row sum of each row of _f_ and _g_. Then, $` \sum_{r_s \in f} 2^{-|r_s|} + \sum_{r_s \in g} 2^{-|r_s|}  \geq 1 `$. This guarantees the existence of a frequent variable.
* For any rows in _f_ and any rows in _g_, they must have at least one variable in common (must have one matching column with value 1). The function for this condition uses bit-wise operator  ```xor``` and iterates through each element in the array, and returns false if they are the same, then the loop pauses and returns true. Else, it checks the next element.

### Additional Transversal Condition Check
If the inputs passes ```check_pre```, the base case for the recursive algorithm is applied by the function ```check_base```. It takes _f_ and _g_ as inputs and checks for the following three conditions: 

**For simplicity, we will be addressing an (_f_,_g_) case, but they also apply for (_g_,_f_) cases as well.**

* If $` |f|= 1, |g|= 0 `$, then they are transversals if the single row in _f_ is a row of all 0s.
* If $`|f|= |g|= 1`$, then they are transversals if both hypergraphs have only one non-zero entry, and they are in the same variable.
* If $`|f|= 1, |g|= k, k \in \mathbb{N} `$, this is also a special condition that FK-A can terminate on. If _f_  is  _k_-_ary_ and _g_ is a $`k \cross k `$ array, then they are transversals when _g_  has all 1s on the diagonal and the row index of 1s in _g_ matches the column index of 1s in _f_.

### Splitting
* Once the algorithm finishes it's pre-conditions check and base case, it implements a recursive ```split``` function for FK-A. With each call, the algorithm chooses a splitting variable to decompose the input functions into two smaller instances of the problem and denotes it as Most Frequent Variable.
* A variable is identified as _Frequent_ if it's frequency exceeds $` \frac{1}{log_2 (|f|+ |g|)}`$. Due to the pre-conditions, all transversal hypergraphs are guaranteed to have frequent variables.
* A variable is **Most Frequent** if it attains maximal frequency in either _f_ or _g_, meaning for most frequent variable $`x`$, $` F_x = \max \{\frac{Count_E: E \in f, x \in E, \forall x \in f}{|f|}, \frac{Count_E: E \in g, x \in E, \forall x \in g}{|g|} \} `$. For variables with same frequency, the splitting variables are selected lexicographically.
* Using the splititng variable, the inputs _f_, _g_ are decomposed into $`f_0, g_0, f_1, g_1`$. Finally, ```split``` returns the most frequent variable, and the resultling $`f_0, g_0, f_1, g_1`$.

## ```fk(f,g)``` Function
* Input of the function is _f_, _g_ and a list variable called path that keeps track of the current position of the tree. ```def_fk(f, g, path)```. If the output is true, means that the input _f_ and _g_ are dual to each other; if is false, means they are not dual.
* At each recursion call, if the base case is not satisfied, it first prints the list variable path to show the current position of the tree, then run the split function on the current input _f_ and _g_.
* It takes the decomposed $`f_0, g_0, f_1, g_1`$ and concatenates and creates two new instances $` (f_0, g_0 \vee g_1), (g_0, f_0 \vee f_1) `$.
