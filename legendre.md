---
tab: Legendre
title: Legendre
math: true
---
# Legendre Polynomial
The mathematical formula of the Legendre Polynomials is[[1](https://en.wikipedia.org/wiki/Legendre_polynomials)]: 

$$P_n(x) = \frac{1}{2^n n!}\frac{\partial^n }{\partial x^n}(x^2-1)^n$$

and can be obtained recursively using the formula[[1](https://en.wikipedia.org/wiki/Legendre_polynomials)]:

$$P_0(x) = 1,\; P_1(x) = x,\; (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)$$

# Associated Legendre Polynomial
The Legendre Polynomials can be generalized and this generalization is called Associated Legendre Polynomials. 
Their definition can be given as a derivative[[2](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Definition_for_non-negative_integer_parameters_%E2%84%93_and_m)]:    

$$P_l^m(x) = \frac{(-1)^m}{2^ll!}(1-x^2)^{m/2}\frac{\partial^{l+m} }{\partial x^{l+m}}(x^2-1)^l$$

or as recursion[[2](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Definition_for_non-negative_integer_parameters_%E2%84%93_and_m)] :   

$$(l-m+1)P_{l+1}^m(x) = (2l+1)xP_l^m(x) - (l+m)P_{l-1}^m(x)$$

However, those two formulas are not suited for numerical computation. For numerical computation, the following formula can be use (demonstrated in [[3](http://www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf)]):

1). Compute $P_l^{-l}(x)$ with:  $P_l^{-l}(x)=\frac{(-1)^l}{2^ll!}(1-x^2)^{l/2}$

2). Compute $P_l^{1-l}(x)$ with:  $P_l^{1-l}(x) = \frac{-2lx}{\sqrt{1-x^2}} P_l^{-l}(x)$

3). Compute $P_l^{-\left\|m\right\|}(x)$ with the recursive formula:    

$$P_l^{m+2}(x) = a_l^m P_l^{m+1}(x) -b_l^m P_l^{m}(x)\text{ where } \begin{cases} & a_l^m = 2(m+1)\frac{x}{\sqrt{1-x^2}}\\ & b_l^m = (l-m)(l+m+1) \end{cases}$$ 

4). If m is positive then obtain $P_l^{m}(x)$ using:    

$$P_l^{\left\|m\right\|}(x) = (-1)^m \frac{(l+\left\|m\right\|)!}{(l-\left\|m\right\|)!}P_l^{-\left\|m\right\|}(x)$$

# Fully Normalized Legendre Polynomial

Fully Normalized Legendre Polynomials have different definitions. In this section, the following one will be used[[3](https://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html)]:

$$N_l^{m}(x) = \sqrt{(l + 1/2)\frac{(l-m)!}{(l+m)!}}P_l^{m}(x)$$

They can be computed using the following two steps[[4](http://mitgcm.org/~mlosch/geoidcookbook/node11.html)]:

1). Compute Sectoral (l=m) $N_m^{m}(x)$ using:

$$N_m^{m}(x) = \left(\sqrt{1-x^2}\right)^m \sqrt{\frac{3}{2}}\;\prod_{m}^{i=2}\sqrt{\frac{2i+1}{2i}}$$

2). Compute non-Sectoral $N_l^{m}(x)$ using the recursive formula:

$$N_l^{m}(x) = a_l^{m}N_{l-1}^{m}(x) - b_l^{m}N_{l-2}^{m}(x)\begin{cases} & a_l^{m} = x\sqrt{\frac{(2l-1)(2l+1)}{(l-m)(l+m)}}\\ & b_l^{m} = \sqrt{\frac{(2l+1)(l+m-1)(l-m-1)}{(l-m)(l+m)(2l-3)}} \end{cases}$$

### Identities formulas

$N_l^{m}(x) = (-1)^mN_l^{-m}(x)$

$N_{m+1}^{m+1}(x) = \sqrt{1 - x^2}\sqrt{\frac{2l+3}{2l+2}}\;N_m^{m}(x)$

