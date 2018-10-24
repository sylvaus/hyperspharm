---
tab: Clenshaw-Curtis
title: Clenshaw-Curtis
math: true
---

The Clenshaw-Curtis quadrature [[1](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)] is a numerical integration allowing to compute exactly $\int_{-1}^{1} P(x) dx$ where $P(x)$ is a polynomial of order less or equal than $N$ with $N + 1$ points taken at $\cos(\frac{n\pi}{N})\;\text{for}\, n \in [0, N]$

The Clenshaw-Curtis quadrature can be written as:

$$\begin{eqnarray}\int_{-1}^{1} P(x)dx &=&  \int_{0}^{\pi} P(\cos\theta)\sin\theta d\theta 
\\&=& \sum_{n=0}^{N}P(\cos\left(\frac{n\pi}{N}\right)) w_n\end{eqnarray}$$

