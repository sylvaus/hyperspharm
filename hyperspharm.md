---
tab: HyperSpharm
title: HyperSpharm
math: true
---

Hyper Spherical Harmonics can be defined as [[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4033314/)]:

$$\begin{matrix} \operatorname{Z^n_{l,m}}{\left (\beta,\theta,\phi \right )} = \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )} \operatorname{NP^m_l}{\left (\cos\theta \right )} e^{i m \phi} \sin^{l}\beta \\ \text{where }\left\{\begin{matrix} \operatorname{NP^m_l}{\left (x \right )} &=& \frac{1}{\sqrt{4\pi}}\sqrt{\frac{\left(2 l + 1\right) \left(l - m\right)!}{\left(l + m\right)!}} &\operatorname{P^m_l}{\left (x \right )} \\ \operatorname{NG^m_l}{\left (x \right )} &=& \frac{2^{m - \frac{1}{2}} }{\sqrt{\pi}} \sqrt{\frac{\left(l + m\right) l!}{\left(l + 2 m - 1\right)!}} \left(m - 1\right)! &\operatorname{G^m_l}{\left (x \right )} \end{matrix}\right. \end{matrix}$$

And Hyper Spherical Transform as:

$$\begin{eqnarray} \operatorname{f^{m}_{l,n}} &=& \int_{0}^{2 \pi}\int_{0}^{\pi}\int_{0}^{\pi} \operatorname{Z^{m*}_{l,n}}{\left (\beta,\theta,\phi \right )} \operatorname{f}{\left (\beta,\theta,\phi \right )} \sin^{2}{\beta} \sin{\theta}\, d\beta\, d\theta\, d\phi  
\\ &=&  \int_{0}^{\pi}\left(\int_{0}^{\pi} \left(\int_{0}^{2 \pi} \operatorname{f}{\left (\beta,\theta,\phi \right )} e^{-i m \phi} \, d\phi\right)\, \operatorname{NP^m_l}{\left (\cos\theta \right )}\sin{\theta}d\theta\right)\, \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )}\sin^{l+2}{\beta}d\beta
\\ &=&  \int_{0}^{\pi}\left(\int_{0}^{\pi} \operatorname{f}^{m}{\left (\beta,\theta \right )} \, \operatorname{NP^m_l}{\left (\cos\theta \right )}\sin{\theta}d\theta\right)\, \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )}\sin^{l+2}{\beta}d\beta
\\ &=&  \int_{0}^{\pi} \operatorname{f}^{m}_{l}{\left (\beta \right )} \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )}\sin^{l+2}{\beta}d\beta\end{eqnarray}$$

We can numerically calculate $\operatorname{f}^{m}{\left (\beta,\theta \right )}$ using an FFT since 
$\operatorname{FFT}(\operatorname{f}) \simeq [\frac{N}{P}\int_{0}^{P} \operatorname{f}{\left(\phi \right )} e^{i m \phi} \, d\phi \quad\text{for}\, m \in [0, N-1] ]$

Thus:

$$\left[\operatorname{f}^{m}{\left (\beta_j,\theta_k \right )}\right] = \frac{N}{P}\operatorname{FFT}(\operatorname{f}{\left (\beta_j,\theta_k,\phi \right )}) $$

Knowing the $\operatorname{f}^{m}{\left (\beta,\theta \right )}$, we can calulate $\operatorname{f}^{m}_{l}{\left (\beta \right )}$ using [Clenshaw Curtis quadrature](clenshaw_curtis.md):

$$\begin{eqnarray}\operatorname{f}^{m}_{l}{\left (\beta_j \right )} 
&=& \int_{0}^{\pi} \operatorname{f}^{m}{\left (\beta_j,\theta \right )} \, \operatorname{NP^m_l}{\left (\cos\theta \right )}\sin{\theta}d\theta 
\\ &\simeq& \sum_{k=0}^{k_{max}} \operatorname{f}^{m}{\left (\beta_j,\theta_k \right )} \, \operatorname{NP^m_l}{\left (\cos\theta_k \right )}w_k\end{eqnarray}$$

where $\theta_k = \frac{k\pi}{k_{max}}$ and the $w_k$ are the  Clenshaw Curtis quadrature weights (see [link](clenshaw_curtis.md))

Then to compute $$\int_{0}^{\pi} \operatorname{f}^{m}_{l}{\left (\beta \right )} \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )}\sin^{l+2}{\beta}d\beta$$: <br>
&nbsp;&nbsp;&nbsp;&nbsp; If $l$ is even, the [Chebyshev–Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) can bes used:<br>
$$\begin{eqnarray}\int_{0}^{\pi} \operatorname{f}^{m}_{l}{\left (\beta \right )} \operatorname{NG^{l+1}_{n-l}}{\left (\cos\beta \right )}\sin^{l+2}{\beta}d\beta
&=&  \int_{0}^{\pi} \operatorname{f}^{m}_{l}{\left (\beta \right )} \operatorname{NG^{l+1}_{n-l}}{\left(\cos\beta \right)}\left(1 - \cos^2{\beta} \right)^{\frac{l}{2}}\sin^2{\beta}d\beta
\\ &=& \int_{-1}^{1} \operatorname{f}^{m}_{l}{\left (\arccos{x} \right )} \operatorname{NG^{l+1}_{n-l}}{\left(x \right)}\left(1 - x^2 \right)^{\frac{l}{2}}\sqrt{1 - x^2}dx
\\ &\simeq& \sum_{j=0}^{j_{max}} \operatorname{f}^{m}_{l}{\left (\arccos{x_j} \right )} \operatorname{NG^{l+1}_{n-l}}{\left(x_j \right)}\left(1 - x_j^2 \right)^{\frac{l}{2}} w_j
\\ &\simeq& \sum_{j=0}^{j_{max}} \operatorname{f}^{m}_{l}{\left (\beta_j \right )} \operatorname{NG^{l+1}_{n-l}}{\left(\cos{\beta_j} \right)}\left(1 - \cos^2{\beta_j} \right)^{\frac{l}{2}} w_j
\\ &\simeq& \sum_{j=0}^{j_{max}} \operatorname{f}^{m}_{l}{\left (\beta_j \right )} \operatorname{NG^{l+1}_{n-l}}{\left(\cos{\beta_j} \right)}\sin^l{\beta_j} w_j
\end{eqnarray}$$