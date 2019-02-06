---
tab: Clenshaw-Curtis
title: Clenshaw-Curtis
math: true
---

The Clenshaw-Curtis quadrature [[1](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)] is a numerical integration allowing to compute exactly $\int_{-1}^{1} P(x) dx$ where $P(x)$ is a polynomial of order less or equal than $N$ with $N + 1$ points taken at $\cos(\frac{n\pi}{N})\;\text{for}\, n \in [0, N]$

The Clenshaw-Curtis quadrature can be written as:

$$\begin{eqnarray}\int_{-1}^{1} P(x)dx &=&  \int_{0}^{\pi} P(\cos\theta)\sin\theta d\theta 
\\&=& \sum_{n=0}^{N}P(\cos\left(\frac{n\pi}{N}\right)) w_n\end{eqnarray}$$

$$\begin{eqnarray}\text{where} \quad w_n &=& \frac{1}{N}\sum_{k=0}^{N_{2max}} \frac{a_k b_n}{1 - (2 k)^2} \cos\left(\frac{2 n k \pi}{N}\right) 
\; \text{,} \; N_{2max} = \begin{cases} & \frac{N}{2} \; \text{if} \; N = 0 \; \text{even} \\ & \frac{N-1}{2} \; \text{otherwise}\end{cases}
\\ a_k &=& \begin{cases} & 1 \; \text{if} \; k = 0 \; \text{or} \; 2k= N \\ & 2 \; \text{otherwise}\end{cases}
\; \text{,} \; b_n = \begin{cases} & 1 \; \text{if} \; n = 0 \; \text{or} \; n= N \\ & 2 \; \text{otherwise}\end{cases}\end{eqnarray}$$ 

Demonstration (based on [[1](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)] ):

Since $P(x)$ is a polynomial of order maximum $N$, $P(\cos\theta)$ can be written as:<br>
$$\begin{eqnarray}P(\cos\theta) &=& \sum_{n=0}^{N}p_n\cos^n\theta
\\&=& \sum_{n=0}^{N}p'_n\cos\left(n\theta\right) \text{ using trigonometric formula}\href{http://mathworld.wolfram.com/TrigonometricPowerFormulas.html}{[2]}
\end{eqnarray}$$

This shows that the maximum frequency of the signal defined by $P(\cos\theta)$ is $N$ and thus it can 