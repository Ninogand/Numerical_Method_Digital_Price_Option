# Implicit Euler and Theta Method

## Introduction
The project focuses on applying numerical methods to solve the Black-Scholes differential equations for option pricing.  
In this case study, the presented algorithms were employed to price a **European digital put option**:  

$$
\phi(S) = \mathbf{1}_{\{S < K\}}
$$

The dynamics of the underlying asset $$S$$ follow a Geometric Brownian Motion (**GBM**):  

$$
dS_t = \mu S_t \, dt + \sigma S_t \, dW_t
$$

Where:  
- $$S_t$$ is the asset price at time $$t$$.  
- $$\mu$$ is the drift rate (mean growth).  
- $$\sigma$$ is the volatility.  
- $$W_t$$ is a standard Brownian motion.  

The GBM dynamics combine two contributions: a deterministic term from the exponential drift $$\mu$$ and a stochastic term from $$W_t$$. This implies that the underlying's dynamics are stochastic and proportional to its volatility $$\sigma$$. Assuming constant volatility simplifies the model, making qualitative observations easier.  

An option is a contract whose value depends on the underlying financial asset. These contracts grant the right (but not the obligation) to buy (**Call**) or sell (**Put**) a specific quantity of the underlying at a predetermined **strike price** by/at expiration. Options are also called *derivatives* since their value derives from an underlying asset.  

### By exercise style:
- **European options**: Exercisable only at expiration.  
- **American options**: Exercisable at any time before expiration.  

### By payoff type:
- **Vanilla options**: Standard buy/sell rights at market value.  
- **Barrier options**: Activation/termination conditions based on the underlying's price.  
- **Asian options**: Payoff depends on the average underlying price during the option's life.  
- **Lookback options**: Payoff depends on the maximum/minimum underlying price during the option's life.  
- **Digital options**: Fixed payoff if the underlying meets a specific condition.  

---

## Objectives and Learning Outcomes
The project aims to analyze the performance and limitations of two numerical methods (**Implicit Euler** and **Theta Method**) for solving the PDE of the digital put option: 

$$
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0
$$

Where:  
- $$V = V(S, t) = \phi(S)$$ is the digital option value.  
- $$\phi(S) = \mathbf{1}_{\{S < K\}}$$ is the indicator function (1 if $$S < K$$, else 0).  

### Implicit Euler
The **Implicit Euler method** is a numerical technique for solving differential equations:  

$$
y_{n+1} = y_n + h \cdot f(t_n, y_n)
$$

For pricing a digital option governed by Equation, we apply the Implicit Euler method to both the underlying price $$S \in [S_L, S_U]$$ and time $$t \in [0, T]$$.  

**Discretization of $$S$$**   

$$
S_i = i \cdot \Delta S, \quad \text{for } i = 0, 1, \dots, N
$$  

**Discretization of $$t$$**  

$$
t_j = j \cdot \Delta t, \quad \text{for } j = 0, 1, \dots, M
$$  

The PDE (Equation \ref{eq:PDE_dig_put}) is approximated as follows:  

**Time derivative**    

$$
\frac{\partial V}{\partial t} \approx \frac{V_i^{j+1} - V_i^j}{\Delta t}
$$  

**First spatial derivative**  

$$
\frac{\partial V}{\partial S} \approx \frac{V_{i+1}^j - V_{i-1}^j}{2 \Delta S}
$$  

**Second spatial derivative**  

$$
\frac{\partial^2 V}{\partial S^2} \approx \frac{V_{i+1}^j - 2 V_i^j + V_{i-1}^j}{(\Delta S)^2}
$$  

Combining these, the discretized PDE becomes:    
**Discretized Equation**  

$$
\frac{V_i^{j+1} - V_i^j}{\Delta t} + \frac{1}{2} \sigma^2 S_i^2 \frac{V_{i+1}^{j+1} - 2 V_i^{j+1} + V_{i-1}^{j+1}}{(\Delta S)^2} + r S_i \frac{V_{i+1}^{j+1} - V_{i-1}^{j+1}}{2 \Delta S} - r V_i^{j+1} = 0
$$  

In matrix form:    
**Matrix Formulation**  

$$
\mathbf{A} \cdot \mathbf{V}^{j+1} = \mathbf{V}^j
$$

where $$\mathbf{A}$$ is a diagonal matrix defined as:    

$$
\begin{aligned}
A_{i,i} &= 1 + \Delta t \left( \frac{\sigma^2 S_i^2}{(\Delta S)^2} + r \right), \\
A_{i,i+1} &= -\Delta t \left( \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} + \frac{r S_i}{2 \Delta S} \right), \\
A_{i,i-1} &= -\Delta t \left( \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} - \frac{r S_i}{2 \Delta S} \right).
\end{aligned}
$$  

**Boundary Conditions**  
For a digital put:    

$$
A_{0,0} = 0, \quad A_{N,N} = e^{-r(T-t)}.
$$  

### Theta Method
The **Theta Method** combines Implicit and Explicit Euler methods. The parameter $$ \Theta $$ balances stability (Implicit) and accuracy (Explicit). Here, $$\Theta = \frac{1}{2}$$ (Crank-Nicolson).  

**Discretized Equation**    

$$
\begin{aligned}
\frac{V_i^{j+1} - V_i^j}{\Delta t} &+ \theta \left( \frac{1}{2} \sigma^2 S_i^2 \frac{V_{i+1}^{j+1} - 2 V_i^{j+1} + V_{i-1}^{j+1}}{(\Delta S)^2} + r S_i \frac{V_{i+1}^{j+1} - V_{i-1}^{j+1}}{2 \Delta S} - r V_i^{j+1} \right) \\
&+ (1 - \theta) \left( \frac{1}{2} \sigma^2 S_i^2 \frac{V_{i+1}^j - 2 V_i^j + V_{i-1}^j}{(\Delta S)^2} + r S_i \frac{V_{i+1}^j - V_{i-1}^j}{2 \Delta S} - r V_i^j \right) = 0
\end{aligned}
$$  

In matrix form:    

$$
\begin{aligned}
&\mathbf{A} \cdot \mathbf{V}^{j+1} = \mathbf{B} \cdot \mathbf{V}^j, \\
&\mathbf{V}^{j+1} = [V_0^{j+1}, V_1^{j+1}, \dots, V_N^{j+1}]^T, \\
&\mathbf{V}^j = [V_0^j, V_1^j, \dots, V_N^j]^T.
\end{aligned}
$$  

**Matrix $$\mathbf{A}$$**      

$$
\begin{aligned}
A_{i,i-1} &= \theta \left( -\frac{r S_i}{2 \Delta S} + \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} \right), \\
A_{i,i} &= -\frac{1}{\Delta t} + \theta \left( -\frac{\sigma^2 S_i^2}{(\Delta S)^2} - r \right), \\
A_{i,i+1} &= \theta \left( \frac{r S_i}{2 \Delta S} + \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} \right).
\end{aligned}
$$  

**Matrix $$\mathbf{B}$$**    

$$
\begin{aligned}
B_{i,i-1} &= (1 - \theta) \left( -\frac{r S_i}{2 \Delta S} + \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} \right), \\
B_{i,i} &= -\frac{1}{\Delta t} + (1 - \theta) \left( -\frac{\sigma^2 S_i^2}{(\Delta S)^2} - r \right), \\
B_{i,i+1} &= (1 - \theta) \left( \frac{r S_i}{2 \Delta S} + \frac{\sigma^2 S_i^2}{2 (\Delta S)^2} \right).
\end{aligned}
$$  

**Boundary Conditions**  

$$
A_{0,0} = 0, \quad A_{N,N} = e^{-r(T-t)}, \quad B_{0,0} = 0, \quad B_{N,N} = e^{-r(T-t)}
$$  

---

## Methodology  
After implementing the numerical methods (see attached file `Equity_Metodinumerici.ipynb`), we tested them for varying grid sizes $$N$$ (spatial) and $$M$$ (temporal) to minimize the **RMSE** (Root Mean Square Error) against the closed-form Black-Scholes solution:  
$$V_{\text{digital}} = e^{-rT} \cdot \Phi(d_2)$$  
where: 

$$
d_2 = \frac{\ln\left(\frac{S}{K}\right) + \left(r - \frac{\sigma^2}{2}\right)T}{\sigma \sqrt{T}}
$$

and $$\Phi(d_2)$$ is the CDF of the standard normal distribution.  

**RMSE Formula**  

$$
\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^n (V_{\text{analyt}} - V_{\text{num}})^2}
$$

---

## Results

### General Results
The table below compares RMSE for different grid sizes $$N$$ (spatial) and $$M$$ (temporal):  

| **N (Spatial)** | **M (Temporal)** | **Implicit Euler** | **Theta Method** | **Difference** |  
|------------------|-------------------|---------------------|------------------|----------------|  
| 50               | 50                | $$1.386 \times 10^{-2}$$ | $$1.395 \times 10^{-2}$$ | $$-8.600 \times 10^{-5}$$ |  
| 100              | 150               | $$1.261 \times 10^{-2}$$ | $$1.235 \times 10^{-2}$$ | $$2.620 \times 10^{-4}$$ |  
| 200              | 300               | $$4.980 \times 10^{-4}$$ | $$4.570 \times 10^{-4}$$ | $$4.100 \times 10^{-5}$$ |  
| 800              | 800               | $$2.776 \times 10^{-3}$$ | $$2.702 \times 10^{-3}$$ | $$7.400 \times 10^{-5}$$ |  

*Table 1: Comparison of numerical results (RMSE) between Implicit Euler and Theta Method.*  

The best results occur at $$N = 200$$ and $$M = 300$$.  

### Temporal Grid Contribution
Testing the impact of temporal grid refinement ($$N = 200$$):  

| **N** | **M** | **Implicit Euler** | **Theta Method** | **Difference** |  
|-------|-------|---------------------|------------------|----------------|  
| 200   | 200   | $$5.170 \times 10^{-4}$$ | $$4.680 \times 10^{-4}$$ | $$4.900 \times 10^{-5}$$ |  
| 200   | 500   | $$4.850 \times 10^{-4}$$ | $$4.550 \times 10^{-4}$$ | $$3.000 \times 10^{-5}$$ |  
| 200   | 1000  | $$4.750 \times 10^{-4}$$ | $$4.570 \times 10^{-4}$$ | $$1.800 \times 10^{-5}$$ |  
| 200   | 2000  | $$4.700 \times 10^{-4}$$ | $$4.600 \times 10^{-4}$$ | $$1.000 \times 10^{-5}$$ |  

*Table 2: Impact of temporal grid refinement on RMSE.*  

### Spatial Grid Contribution
Testing the impact of spatial grid refinement ($$M = 800$$):  

| **N** | **M** | **Implicit Euler** | **Theta Method** | **Difference** |  
|-------|-------|---------------------|------------------|----------------|  
| 10    | 800   | $$3.620 \times 10^{-2}$$ | $$3.621 \times 10^{-2}$$ | $$-1.000 \times 10^{-5}$$ |  
| 100   | 800   | $$1.256 \times 10^{-2}$$ | $$1.251 \times 10^{-2}$$ | $$4.700 \times 10^{-5}$$ |  
| 200   | 800   | $$4.770 \times 10^{-4}$$ | $$4.560 \times 10^{-4}$$ | $$2.100 \times 10^{-5}$$ |  
| 800   | 800   | $$2.776 \times 10^{-3}$$ | $$2.702 \times 10^{-3}$$ | $$7.400 \times 10^{-5}$$ |  

*Table 3: Impact of spatial grid refinement on RMSE.*  

---

## Discussion & Conclusion

### General Results  
Finer grids generally yield better results. However, the Theta Method underperforms Implicit Euler for coarse grids due to its second-order temporal accuracy, where discretization error scales with $$(\Delta t)^2$$.  

### Temporal Grid Analysis  
Table 2 confirms that the Theta Method generally outperforms Implicit Euler as $$M$$ increases. However, excessively fine temporal grids may reduce performance.  

### Spatial Grid Analysis  
Table 3 shows that spatial refinement improves accuracy up to a threshold ($$N = 200$$). Beyond this, errors accumulate comparably for both methods.  

**Conclusion**: Both methods benefit from grid refinement, but the Theta Method’s advantage diminishes with very fine grids. The choice between methods depends on computational resources and desired accuracy.  
