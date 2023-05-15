import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math

def generate_gamma_rv(n, mean):
    randomExp = rd.gammavariate(n, mean)
    return randomExp

def generate_sqrt_gamma_rv(n, mean):
    return np.sqrt(generate_gamma_rv(n, mean)/n)

# def sqrt_gamma_pdf(x, n, mean):
#     if x < 0 or n <= 0 or mean <= 0:
#         return 0

#     term1 = 2 * x * n
#     arg = n * x * x
#     critical_part = 1
#     for i in range(1, n):
#             critical_part *= arg / i * math.exp(-arg / mean / (n - 1))
#     pdf = term1 * (critical_part) / (mean ** n)
#     return pdf

def sqrt_gamma_pdf(x, n, mean):
    if x < 0 or n <= 0 or mean <= 0:
        return 0

    term1 = 2 * x * n
    arg = n * x * x
    
    log_critical_part = 0
    for i in range(1, n):
        log_critical_part += math.log(arg / i) - math.log(mean) - (arg / mean / (n - 1))

    pdf = term1 * math.exp(log_critical_part)
    return pdf

# Parameters
n = 5000
mean = 1
num_samples = 10000
x_values = np.linspace(0.75, 1.25, 1000)

# Generate samples and calculate empirical distribution
samples = [generate_sqrt_gamma_rv(n, mean) for _ in range(num_samples)]
hist, bins = np.histogram(samples, bins=50, density=True)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Calculate theoretical distribution
pdf_values = [sqrt_gamma_pdf(x, n, mean) for x in x_values]

# Plot empirical and theoretical distributions
plt.hist(samples, bins=50, density=True, alpha=0.6, label='Empirical')
plt.plot(x_values, pdf_values, label='Theoretical', linewidth=2)
plt.legend()
plt.xlabel('x')
plt.ylabel('Probability density')
plt.title('Empirical vs Theoretical Sqrt-Gamma Distribution')
plt.show()