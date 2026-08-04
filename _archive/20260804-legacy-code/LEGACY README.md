# FK-A Algorithm + Visualization

This repository implements the **Fredman-Khachiyan Algorithm-A (FK-A)** for solving the hypergraph dualization problem, also known as the transversal hypergraph problem (\texttt{Trans-Hyp}). The implementation includes a visualization tool for recursion trees and animations of the algorithm’s execution.

## Overview

The FK-A algorithm determines whether two hypergraphs $` G `$ and $` H `$, represented as sets of hyperedges, are transversals of each other. A hypergraph $` G `$ is a transversal of $` H `$ if every hyperedge in $` G `$ intersects at least one hyperedge in $` H `$, and vice versa. The algorithm works recursively, decomposing hypergraphs until base cases are reached. 


## Authors and Copyright

This implementation of the Fredman-Khachiyan Algorithm-A (FK-A) was jointly developed by **Saimon Islam** and **Lacey Liang**. The project includes enhancements and optimizations for superset removal, additional base case handling, and visual representation of the recursion tree.

This work is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this code under the terms of the GPL-3.0 license. Please ensure appropriate attribution to the authors when using this work.

**Copyright © 2024 Saimon Islam and Lacey Liang.**
---

## Preprocessing Steps

### Duplicate Hyperedges Removal

The algorithm removes duplicate hyperedges in the input hypergraphs, ensuring that all hyperedges are unique. This process enforces the Sperner property, a prerequisite for the algorithm.

### Superset Hyperedges Removal

The FK-A algorithm eliminates hyperedges that are supersets of other hyperedges. A hyperedge $` A `$ is considered a superset of hyperedge $` B `$ if $` B \subseteq A `$. By ensuring this property, the algorithm guarantees that the input hypergraphs remain Sperner throughout the process.

### Preconditions

Before initiating recursive decomposition, FK-A checks the following conditions:
1. The two hypergraphs must share the same set of variables (same number of vertices).
2. The largest hyperedge size in $` G `$ must not exceed the number of hyperedges in $` H `$, and vice versa.
3. Let $` r_s `$ be the size of each hyperedge in $` G `$ and $` H `$. Then, 
   $` \sum_{r_s \in G} 2^{-|r_s|} + \sum_{r_s \in H} 2^{-|r_s|} \geq 1 `$. 
   This guarantees the existence of a frequent variable.
4. Every hyperedge in $` G `$ must intersect with at least one hyperedge in $` H `$. This ensures transversality conditions are satisfied.

If any of these preconditions fail, the algorithm terminates and returns $` G \neq Tr(H) `$. 

---

## Recursive Decomposition

### Splitting Process

At each recursion step, the algorithm identifies the **Most Frequent Variable (MFV)**, which appears most often in the hyperedges of $` G `$ and $` H `$. This ensures an efficient reduction of the problem size.

#### Decomposition Rules

Using the MFV, the input hypergraphs are split as follows:
- $` G_0, G_1 `$: Hyperedges of $` G `$ where the MFV is absent or present, respectively.
- $` H_0, H_1 `$: Hyperedges of $` H `$ where the MFV is absent or present, respectively.

The algorithm recursively solves the problem for the pairs:
1. $` (G_1, H_0 \lor H_1) `$
2. $` (H_1, G_0 \lor G_1) `$

---

### Base Cases

The recursion terminates when:
1. $` |G| = 1, |H| = 0 `$: $` G `$ is a single hyperedge containing no variables, and $` H `$ is empty.
2. $` |G| = 1, |H| = 1 `$: $` G `$ and $` H `$ each contain one hyperedge with exactly one variable in common.
3. $` |G| = 1, |H| = k, k \in \mathbb{N} `$: $` G `$ contains a single $` k `$-ary hyperedge, and $` H `$ is a $` k \times k `$ collection of singleton hyperedges forming the identity structure.

These conditions allow the result to be computed directly without further recursion.

---

## Visualization

### Recursion Tree Structure

The recursion tree $` T `$ represents the hierarchical decomposition:
- **Nodes:** Represent recursive calls with subsets of $` G `$ and $` H `$. Each node stores:
  - The current hypergraphs $` G `$, $` H `$
  - Hyperedges removed during preprocessing
  - The splitting variable
- **Edges:** Correspond to the splitting process into $` G_0, G_1 `$, and $` H_0, H_1 `$. 

### Tree Properties

1. **Left Branch (\textbf{L}):** Represents the recursive call $` (G_1, H_0 \lor H_1) `$.
2. **Right Branch (\textbf{R}):** Represents the recursive call $` (H_1, G_0 \lor G_1) `$. 

For each node $` v `$, its depth $` d_v(T) `$ is defined as the length of the path from the root to $` v `$. Tracking this provides insight into the recursive complexity.

## Optimizations Over Native FK-A Implementation

### Superset Removal Optimization

In the native FK-A implementation, superset removal involves pairwise comparisons of hyperedges, resulting in a time complexity of $`O(n^2 \cdot m)`$ for an input hypergraph with $`n`$ hyperedges and $`m`$ variables. This naive approach compares each hyperedge against all others, checking whether one is a superset of another. 

Our optimized implementation improves this process by leveraging sorting and structured comparisons:
1. **Sorting by Hyperedge Size:** Hyperedges are sorted in descending order of size (sum of variables), ensuring larger hyperedges are processed first.
2. **Zero-Column Indexing:** Columns with zero values in each hyperedge are identified, significantly reducing unnecessary comparisons.
3. **Efficient Superset Detection:** By comparing the zero-indexed columns across hyperedges, supersets are detected and removed in a single pass, minimizing redundant operations.
4. **Deletion and Output:** The reduced hypergraph is produced alongside a record of removed hyperedges for transparency.

These steps reduce the computational overhead and improve the scalability of superset removal in sparse hypergraphs.

---

### Additional Base Case Optimization

The original FK-A implementation defines a base case when $`|G| \cdot |H| \leq 1`$, allowing for termination in $`O(1)`$. However, certain special cases are computationally expensive to reach under this condition in the native approach.

Our implementation introduces an additional pre-check for cases where:
- $`|G| = 1, |H| = k, k \in \mathbb{N}`$: If $`G`$ contains a single $`k`$-ary hyperedge and $`H`$ is a $`k \times k`$ collection of singleton hyperedges forming an identity structure, the result is computed directly without recursion. 

This enhancement avoids unnecessary recursive steps, simplifying the recursion tree and reducing its depth. For example:
- In the native implementation, $`k-1`$ additional recursive calls are required to terminate at trivial base cases like $`|G| = 1, |H| = 1`$. 
- In our optimized approach, these are resolved in $`O(k)`$, bypassing redundant operations.

By reducing the depth of recursion and addressing these edge cases directly, the overall efficiency of the algorithm is significantly improved.


---

## Summary

The FK-A algorithm balances preprocessing, recursive decomposition, and direct handling of base cases to ensure computational efficiency. This implementation achieves  clarity in visualizing hypergraph dualization and gives us additional insight into the algorithm's mechanism.

