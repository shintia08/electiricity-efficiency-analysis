# electiricity-efficiency-analysis

Technical efficiency of electricity distribution units in Indonesia is analyzed using Stochastic Frontier Analysis (SFA) with a Bayesian estimation approach to obtain more robust efficiency estimates under assumption violations and the presence of outliers.

The dataset is derived from PLN Statistics 2024 (Unaudited) published by PT PLN (Persero), covering 33 electricity distribution units across Indonesia. Input variables include the length of medium-voltage distribution lines, average electricity tariff, and connected customer load, while electricity sales is used as the output variable.

Results show that the Cobb-Douglas (CD) frontier with a Half-Normal inefficiency specification provides the most appropriate model. However, diagnostic tests indicate a violation of the normality assumption in the error term, leading to Bayesian estimation with a Student’s t disturbance distribution for more stable and robust parameter estimation.

Technical efficiency scores range from 0.6306 to 0.9622 with an average of 0.8719, indicating relatively high efficiency across most units. Efficiency scores are further combined with connected customer load and analyzed using K-Medoids clustering, resulting in three distinct clusters with good validity measures, including a Silhouette coefficient of 0.657, Davies-Bouldin Index of 0.399, and R² of 0.815.
