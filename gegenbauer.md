---
tab: Gegenbauer
title: Gegenbauer
math: true
---

Gegenbauer Polynomials is defined by the two following initial values and the recurring formula [[1](https://en.wikipedia.org/wiki/Gegenbauer_polynomials)]:

$C^m_{0}(x) = 1$

$C^m_{1}(x) = 2mx$

$C^m_{l}(x) = \frac{1}{l} \left[ 2x{\left (m + l - 1 \right )}C^m_{l-1}(x) -\left(2 m + l - 2\right) C^m_{l-2}(x) \right]$

These polynomials can be normalized by the following factor[[1](https://en.wikipedia.org/wiki/Gegenbauer_polynomials)]: 

$Factor_l^m = \sqrt{\frac{2^{2m - 1}\left(m + l\right) l!}{\pi\left(2 m + l - 1\right)!}} \left(m - 1\right)! \quad \forall m > \text{-}\frac{1}{2}$

Thus:

$ \operatorname{N^m_l}{\left (x \right )} = \sqrt{\frac{2^{2m - 1}\left(m + l\right) l!}{\pi\left(2 m + l - 1\right)!}} \left(m - 1\right)!\operatorname{C^m_l}{\left (x \right )}$

$\operatorname{N^{1}_{0}}{\left (x \right )} = \sqrt{\frac{2}{\pi}} \quad \text{and} \quad \operatorname{N^{m+1}_0}{\left (x \right )} = \sqrt{\frac{2(m+1)}{2m+1}}\operatorname{N^{m}_0}{\left (x \right )}$

$\operatorname{N^{1}_{1}}{\left (x \right )} = 2 x \sqrt{\frac{2}{\pi}} \quad \text{and} \operatorname{N^{m+1}_1}{\left (x \right )} = \sqrt{\frac{2(m+2)}{2m+1}}\operatorname{N^{m}_1}{\left (x \right )}$

$ \operatorname{N^m_l}{\left (x \right )} = 2 x \sqrt{\frac{ (l + m) (l + m - 1)}{l (l + 2 m - 1)}} \operatorname{N^m_{l-1}}{\left (x \right )} - \sqrt{\frac{\left(l - 1\right) \left(l + m\right) \left(l + 2 m - 2\right)}{l \left(l + m - 2\right) \left(l + 2 m - 1\right)}} \operatorname{N^{m}_{l-2}}{\left (x \right )}$
