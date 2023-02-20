# 'negacyclic_probabilistic_bound.py'

This python code computes for random nega-cyclic matrices their smallest and largest singular vaules for dimensions 2^i where i goes from 1 to 10.

We referred to this repository in an earlier version of our paper "Towards Classical Hardness of Module-LWE: The Linear Rank Case". The earlier version can be accessed via https://eprint.iacr.org/archive/2020/1020/1598495175.pdf. Note that since an update on the 16.03.2021, the lower bound on the smallest singular value is not needed anymore.


# `heuristics_bound.sage`

This sage code computes for general cyclotomic rings heuristic bounds on the smallest and largest singular values of the multiplication matrix.

We refer to this repository in our paper "Entropic Hardness of Module-LWE from Module-NTRU", accessible via https://eprint.iacr.org/2022/245

---

```C
usage: sage heuristics_bound.sage [-h] [--nu nu] [--d d] [--nb nb]
```

## Change the ordering 

First, we show that the ordering of the embedding does not have
any consequences on the lattice, nor the singular values. 
Indeed, changing the order corresponds to multiplying the
embedding by a permutation matrix, which is unimodular
(integer matrix whose inverse is also an integer matrix) and unitary
(orthonormal matrix).   

Assume that the first canonical embedding $\sigma^{(1)}$ is indexed so that for all $i \in [n/2]$, 
it holds that $\sigma^{(1)}_{i + n/2} = \overline{\sigma^{(1)}_i}$.
Consider the second canonical embedding $\sigma^{(2)}$ is indexed so that for all $i \in [n/2]$, 
it holds that $\sigma^{(2)}_{2i-1} = \overline{\sigma^{(2)}_{2i}}$.
Then, there is a permutation matrix $\mathbf{P}$ such that 
$\sigma^{(2)} = \mathbf{P}\sigma^{(1)}$. The multiplication matrices
are thus linked by

```math
M_{\sigma^{(2)}}(\cdot) = \mathbf{P}\cdot M_{\sigma^{(1)}}(\cdot) \cdot \mathbf{P}^T
```

Additionally, we have $M_{\sigma^{(2)}}(x) = \textrm{diag}(\sigma^{(1)}_1(x), \overline{\sigma^{(1)}_1(x)}, \ldots, \sigma^{(1)}_{n/2}(x), \overline{\sigma^{(1)}_{n/2}(x)})$.
In the paper we defined $\mathbf{U}_H^{(1)}$ as

```math
\mathbf{U}_H^{(1)} = \frac{1}{\sqrt{2}}\begin{bmatrix}\mathbf{I}_{n/2} & -i\mathbf{I}_{n/2} \\ \mathbf{I}_{n/2} & i\mathbf{I}_{n/2}\end{bmatrix}.
```

The multiplication matrix in $\sigma_H^{(1)}$ is thus $M_{\sigma_H^{(1)}}(\cdot) = (\mathbf{U}_H^{(1)})^\dagger \cdot M_{\sigma^{(1)}}(\cdot) \cdot \mathbf{U}_H^{(1)}$. By defining
$\mathbf{U}_H^{(2)}$ by

```math
\mathbf{U}_H^{(2)} = \textrm{diag}(\mathbf{R}, \ldots, \mathbf{R})\text{, with }\mathbf{R} = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 & i \\ 1 & -i \end{bmatrix},
```

we have $\mathbf{U}_H^{(2)} = \mathbf{P} \cdot \mathbf{U}_H^{(1)} \cdot \mathbf{P}^T$. The multiplication matrix in $\sigma_H^{(2)}$
is therefore

```math
M_{\sigma_H^{(2)}}(\cdot) = \mathbf{P}\cdot M_{\sigma_H^{(1)}}(\cdot) \cdot \mathbf{P}^T
```

The singular values are thus the same because $\mathbf{P}$ is unitary. From now on, we thus consider the second definition. We
omit the exponent $(2)$ for clarity.

## `R_lattice_basis` 

The generation of the basis of $R$ with respect to $\sigma_H$
is done by the Sage function $\texttt{minkowski\_embedding}$. It
gives the basis whose columns are the $\sigma_H(\zeta^0), \ldots, \sigma_H(\zeta^{n-1})$. The definition of $\sigma_H$ in the function matches the second order described above.

## `sample_rotation_matrix`

We need to sample $x \in R$ such that $\sigma_H(x)$ follows the discrete 
Gaussian distribution $\mathcal{D}_{\sigma_H(R ), \gamma}$. Then, we need 
to return the multiplication matrix $M_{\sigma_H}(x) = \mathbf{U}_H^\dagger \textrm{diag}(\sigma(x)) \mathbf{U}_H$. 
Let $x$ be in $R$. Because of the way we ordered the canonical embedding, we have

```math
\textrm{diag}(\sigma(x)) = \textrm{diag}(\mathbf{\Sigma}_1, \ldots, \mathbf{\Sigma}_{n/2}),
```

with $\mathbf{\Sigma}_i = \begin{bmatrix}\sigma_{i}(x) & 0 \\ 0 & \overline{\sigma_{i}(x)}\end{bmatrix}$. 
Computing $M_{\sigma_H}(x)$ leads to

```math
M_{\sigma_H}(x)= \textrm{diag}(\mathbf{\Sigma}_1', \ldots, \mathbf{\Sigma}_{n/2}'),
```

with $\mathbf{\Sigma}_i' = \mathbf{R}^\dagger \mathbf{\Sigma}_i \mathbf{R} = \begin{bmatrix}\mathfrak{R}(\sigma_{i}(x)) & -\mathfrak{I}(\sigma_i(x)) \\ \mathfrak{I}(\sigma_i(x)) & \mathfrak{R}(\sigma_i(x))\end{bmatrix}$. 

Since the output depends on the embedding of $x$, it is not necessary to 
get the ring element $x$. The output of $\mathcal{D}_{\sigma_H(R ), \gamma}$ is $\mathbf{v} = \sigma_H(x)$ 
for some $x \in R$. We then have

```math
\mathbf{\Sigma}_i' = \frac{1}{\sqrt{2}} \begin{bmatrix}\mathbf{v}_{2i-1} & -\mathbf{v}_{2i} \\ \mathbf{v}_{2i} & \mathbf{v}_{2i-1}\end{bmatrix}.
```

We can then express $M_{\sigma_H}(x)$ only in terms of $\mathbf{v}$.


# `heuristics_constants.sage`
---

```C
usage: sage heuristics_constants.sage [-h] [--nu nu] [--d d] [--nb nb]
```

This part is dedicated to assess the behavior of $C_\gamma$ and $c_\gamma$ with
respect to $\gamma$. Therefore, we can only consider $d \times d$ complex matrices
that correspond to one embedding of the matrix over $R$.

## Sampling 

As ($n/2$ of) the embeddings of a $d \times d$ matrix are independent,
we only need to generate $\texttt{nb} \cdot 2/n$ matrices $\mathbf{A} \in R^{d \times d}$.
It will give us $\texttt{nb}$ matrices over $\mathbb{C}^{d \times d}$ once
embedded with the first $n/2$ embeddings.  
As we now need to separate each embedding (we want $\sigma_1(\mathbf{A}), \ldots, \sigma_{n/2}(\mathbf{A})$),
we use a different sampling procedure to avoid unnecessary computations. So we
can sample $n/2$ independent $d \times d$ matrices over $\mathbb{C}$ by sampling $d^2$ 
vectors $(\mathbf{v}^{(k\ell)})_{k,\ell \in [d]}$. 

## Evaluation of constants

We have the proven relation 

```math
\mathbb{P}[s_{\min}(\mathbf{A}) > \varepsilon/\sqrt{d}] \leq C_\gamma \varepsilon + c_\gamma^d.
```

To evaluate $C_\gamma$ and $c_\gamma$ we sample many matrices, and we compare their smallest singular
value to $\varepsilon/\sqrt{d}$ for a range of values of $\varepsilon$. We thus obtain a 
statistical approximation of $\mathbb{P}[s_{\min}(\mathbf{A}) > \varepsilon/\sqrt{d}]$.
We then compute an affine function of $\varepsilon$ that upper-bounds these approximations
for all $\varepsilon$. This gives us approximations of $C_\gamma$ and $c_\gamma$.  

By repeating this process for several values of $\gamma$, we can determine the behavior
of $C_\gamma$ and $c_\gamma$ with respect to $\gamma$. This heuristics shows that
$C_\gamma$ behaves like $O(1/\gamma^2)$, and $c_\gamma$ seems negligible in $\gamma$.   
