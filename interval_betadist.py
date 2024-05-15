import scipy.stats as stats

# Define the parameters for each cell
params = {
    'Cell 1': (11, 11),
    'Cell 2': (31, 21),
    'Cell 3': (6, 6)
}

# Define the parameters for each cell with empirical Bayes prior
params_empirical = {
    'Cell 1': (35.31, 29.69),
    'Cell 2': (55.31, 39.69),
    'Cell 3': (30.31, 24.69)
}

# Calculate and print the 2.5th and 97.5th percentiles for each Beta distribution
for cell, (alpha, beta) in params.items():
    lower_bound = stats.beta.ppf(0.025, alpha, beta)
    upper_bound = stats.beta.ppf(0.975, alpha, beta)
    print(f"{cell} 95% credible interval: ({lower_bound:.4f}, {upper_bound:.4f})")


# Define the parameters for each cell with empirical Bayes prior
params_empirical = {
    'Cell 1': (35.31, 29.69),
    'Cell 2': (55.31, 39.69),
    'Cell 3': (30.31, 24.69)
}

# Calculate and print the 2.5th and 97.5th percentiles for each Beta distribution
for cell, (alpha, beta) in params_empirical.items():
    lower_bound = stats.beta.ppf(0.025, alpha, beta)
    upper_bound = stats.beta.ppf(0.975, alpha, beta)
    print(f"{cell} 95% credible interval: ({lower_bound:.4f}, {upper_bound:.4f})")



