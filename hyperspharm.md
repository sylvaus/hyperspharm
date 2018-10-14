---
tab: HyperSpharm
title: HyperSpharm
math: true
---

Hyper Spherical Harmonics can be defined as [[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4033314/)]:

$$\begin{matrix} \operatorname{Z^n_{l,m}}{\left (\beta,\theta,\phi \right )} = \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )} \operatorname{NP^m_l}{\left (\cos\theta \right )} e^{i m \phi} \sin^{l}\theta \\ \text{where }\left\{\begin{matrix} \operatorname{NP^m_l}{\left (x \right )} &=& \frac{1}{\sqrt{4\pi}}\sqrt{\frac{\left(2 l + 1\right) \left(l - m\right)!}{\left(l + m\right)!}} &\operatorname{P^m_l}{\left (x \right )} \\ \operatorname{NG^m_l}{\left (x \right )} &=& \frac{2^{m - \frac{1}{2}} }{\sqrt{\pi}} \sqrt{\frac{\left(l + m\right) l!}{\left(l + 2 m - 1\right)!}} \left(m - 1\right)! &\operatorname{G^m_l}{\left (x \right )} \end{matrix}\right. \end{matrix}$$

And Hyper Spherical Transform as:

$$\operatorname{f^{m}_{l,n}} = \int_{0}^{2 \pi}\int_{0}^{\pi}\int_{0}^{\pi} \operatorname{Z^{m*}_{l,n}}{\left (\beta,\theta,\phi \right )} \operatorname{f}{\left (\beta,\theta,\phi \right )} \sin^{2}{\left (\beta \right )} \sin{\left (\theta \right )}\, d\beta\, d\theta\, d\phi $$