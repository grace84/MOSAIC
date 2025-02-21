# MOSAIC: Integrating Genomic and Clinical Data for Mutation Clustering in SARS-CoV-2

### Clustering of Logistic Regression Model (CLRM)
```math
y_{n,m} | Z_{m,k}=1 \sim \text{Bernoulli}(\pi_{n, k})  \text{ with } \pi_{n, k} = \frac{\exp( \mathbf{x}_{n}^{\top} \boldsymbol{\beta}_{k})}{1 + \exp( \mathbf{x}_{n}^{\top} \boldsymbol{\beta}_{k})}, 
```
```math
\mathbf{Z}_m \sim \text{Multinomial}(\boldsymbol{\lambda}) \text{ with } \boldsymbol{\lambda} = [\lambda_{1}, \ldots, \lambda_{K}]^{\top}, \lambda_{k} \ge 0 \text{ and } \sum_{k=1}^K \lambda_{k} = 1.
```
- $y_{n,m}$ denote whether mutation m is present (1) or absent (0) in the genomic sequence of the n-th individual;
- $Z_{m,k}$ indicates if mutation m belongs to cluster k;
- $\pi_{n, k}$ indicates the probability that mutation m exists in the n-th genome sequence, given that we know the clinical features of individual n and the mutation m belongs to
cluster k;
- $\mathbf{x}_{n}$ is a $(D+1) \times 1$ vector representing the intercept and covariates associated with n-th patient;
- $\boldsymbol{\beta}_k$ is a $(D+1) \times 1$ vector storing the logistic regression intercept and coefficients associated with mutations in cluster $k$.
