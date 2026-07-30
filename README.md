# InterDependence Scores (IDS)

Code to compute IDS between sets of variables. Supports both numpy only users and PyTorch users.

## Prerequisites

Python version >= 3.7.

## Installation

Clone the repo and run the following command to install ids:

```bash
pip install -e .
```

## Example usage

The following two use cases are supported.

1. Given a data matrix, X, of size n x d (n samples, d variables), one can compute IDS between all pairs of variables using:

```python
ids.compute_IDS(X)
```

2.  Given data matrices, X and Y, of size n x d and n x c (n samples, d variables in X and c variables in Y), one can compute IDS between all d variables in X and all c variables of Y using:

```python
ids.compute_IDS(X, Y)
```

Notebooks demonstrating various use cases of IDS are provided in examples.

## Parameters

The `compute_IDS` function supports the following parameters:

1. `X` - a matrix of size n samples by d variables
2. `Y` (optional, default=None) - a matrix of size n samples by c variables
3. `num_terms` (optional, default=6) - number of terms used in the Taylor series approximation of the Gaussian kernel
4. `p_norm` (optional, default='max') - 'max' means using IDS-max, integer 1 means using IDS-1, integer 2 means using IDS-2
5. `p_val` (optional, default=False) - boolean indicating whether to return p-values or not
6. `num_tests` (optional, default=100) - number of permutation tests to run for computing p-values
7. `bandwidth_term` (optional, default=1/2) - constant multiplier in exponent of Gaussian kernel

## Handling sparse matrix and batching logic

TThe data is a big grid.

- rows are cells
- columns are genes

Call it n cells by d genes.

IDS wants to answer, for every pair of genes, "do these two move together in any way?" — including curved, non-straight-line relationships that plain correlation misses.

To catch those curved relationships, IDS does one trick before it measures anything.

It takes each gene's single column and expands it into 6 columns, which is called the "feature map". (ex. 1, x, x², … under a bell curve). This means a gene is not 1 column; it is 6 columns.

Then, IDS computes the correlation between every pair of those feature-columns, and finally squashes each gene-vs-gene 6×6 patch back down to a single dependence number.

> The end result is a d × d table: one IDS score per gene pair.

Now, where does the memory explode? Two separate places, and this is the key insight:

1. **The expanded data grid** — n rows × 6d columns. This grows with the number of _cells_. With 27 million cells, this grid is enormous.
2. **The correlation table** — 6d × 6d. This grows with the number of _genes_, squared. With 20,000 genes that's ~120,000 × 120,000 ≈ 14 billion numbers.

Either one alone can be too big for memory. That's when batching comes in — and there's one trick for each explosion.

Trick 1 (tames the genes² table):

- Do a few genes at a time (batch over genes d).
- Instead of expanding all 20,000 genes into one giant table, take a small handful of genes at a time.
- Expand just those, correlate just those.
- Then move to the next handful, and continue. Paste all the small blocks together into the full table.

Trick 2 (tames the cells grid):

- Do a few cells at a time (batch over cells n).
- Even one handful of genes still has a column that's 27-million-cells tall.
- So pour the cells through in slabs, and repeat the same block idea.

You may wonder... a correlation needs all the cells to compute — how can you look at them a slab at a time?

A correlation between two columns is built from exactly three ingredients:

1. each column's average
2. each column's spread (standard deviation)
3. how much they move together (covariance)

Every one of those three is just a sum — a sum of values, a sum of squares, a sum of products — and sums can be built up a slab at a time. So you never need all the cells in memory at once, and the final answer is exactly the same as if you had.

The paper is careful about the _order_, though: it does this in two passes.

- **Pass 1:** stream through all the cells once, just to compute each column's average, and save it.
- **Pass 2:** stream through again, subtract that saved average (this is "centering"), and only then measure how the columns move together.

Finding the true average first, then measuring around it, is the numerically safe order — and it's what keeps the answer accurate once the number of cells gets very large. (There's a faster single-pass version that skips Pass 1, but it can lose precision at huge cell counts.)
