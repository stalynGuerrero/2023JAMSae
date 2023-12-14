# Session 4 - Area model for labor market statistics




The National Labour Force Survey (NLFS) is a key survey in Jamaica conducted by the Statistical Institute of Jamaica (STATIN). This survey provides detailed and up-to-date information on the dynamics of the country's labor market. Some highlights of the NLFS include:

1. **Employment and Unemployment Measurement:** The survey gathers comprehensive data on the employment status of working-age individuals, identifying the employed, unemployed populations, and the unemployment rate across different demographic segments and regions of the country.

2. **Underemployment and Labor Conditions:** In addition to measuring unemployment, the survey assesses employment conditions, including underemployment, involuntary part-time work, and other forms of inadequate employment.

3. **Demographic Variables:** Important demographic data such as age, gender, education, and geographic location are collected, enabling the identification of specific labor patterns within different population groups.

4. **Frequency and Scope:** The survey is conducted periodically to capture changes in the labor market over time and covers a broad, representative sample of the population to ensure precise and reliable results.


## Definition of the Multinomial Model

-   Let $K$ be the number of categories of the variable of interest $Y \sim multinomial\left(\boldsymbol{\theta}\right)$, with $\boldsymbol{\theta}=\left(p_{1},p_{2},\dots ,p_{k}\right)$ and $\sum_{k=1}^{K}p_{k}=1$.

-   Let $N_i$ be the number of elements in the i-th domain and $N_{ik}$ be the number of elements in the k-th category. Note that $\sum_{k=1}^{K}N_{ik}=N_{i}$ and $p_{ik}=\frac{N_{ik}}{N_{i}}$.

-   Let $\hat{p}_{ik}$ be the direct estimation of $p_{ik}$ and $v_{ik}=Var\left(\hat{p}_{ik}\right)$, and denote the estimator of the variance by $\hat{v}_{ik}=\widehat{Var}\left(\hat{p}_{ik}\right)$

Note that the design effect changes between categories; therefore, the first step is to define the effective sample size per category. This is:

The estimation of $\tilde{n}$ is given by $\tilde{n}_{ik} = \frac{(\tilde{p}_{ik}\times(1-\tilde{p}_{ik}))}{\hat{v}_{ik}},$

$\tilde{y}_{ik}=\tilde{n}_{ik}\times\hat{p}_{ik}$

Then, $\hat{n}_{i} = \sum_{k=1}^{K}\tilde{y}_{ik}$

From where it follows that $\hat{y}_{ik} = \hat{n}_i\times \hat{p}_{ik}$

Let $\boldsymbol{\theta}=\left(p_{1},p_{2}, p_{3}\right)^{T}=\left(\frac{N_{i1}}{N_{i}},\frac{N_{i2}}{N_{i}}\frac{N_{i3}}{N_{i}}\right)^{T}$, then the multinomial model for the i-th domain would be:

$$
\left(\tilde{y}_{i1},\tilde{y}_{i2},\tilde{y}_{i3}\right)\mid\hat{n}_{i},\boldsymbol{\theta}_{i}\sim multinomial\left(\hat{n}_{i},\boldsymbol{\theta}_{i}\right)
$$ 
Now, you can write $p_{ik}$ as follows:

$\ln\left(\frac{p_{i2}}{p_{i1}}\right)=\boldsymbol{X}_{i}^{T}\beta_{2} + u_{i2}$ and
$\ln\left(\frac{p_{i3}}{p_{i1}}\right)=\boldsymbol{X}_{i}^{T}\beta_{3}+ u_{i3}$



Given the restriction $1 = p_{i1} + p_{i2} + p_{i3}$ then 
$$p_{i1} + p_{i1}(e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta_{2}}}+  u_{i2})+p_{i1}(e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta}_{3}} + u_{i3})$$ from where it follows that 

$$
p_{i1}=\frac{1}{1+e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta_{2}}}+ u_{i2}+e^{\boldsymbol{X_{i}}^{T}\boldsymbol{\beta_{3}}}+ u_{i3}}
$$

The expressions for $p_{i2}$ and $p_{i3}$ would be:

$$
p_{i2}=\frac{e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta}_{2}} + u_{i2}}{1+e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta_{2}}}+ u_{i2}+e^{\boldsymbol{X_{i}}^{T}\boldsymbol{\beta_{3}}}+ u_{i3}}
$$

$$
p_{i3}=\frac{e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta}_{3}}+ u_{i3}}{1+e^{\boldsymbol{X}_{i}^{T}\boldsymbol{\beta_{2}}}+ u_{i2}+e^{\boldsymbol{X_{i}}^{T}\boldsymbol{\beta_{3}}}+ u_{i3}}
$$

## Loading Libraries

  -   The `survey` library is a statistical analysis tool in R that allows working with complex survey data, such as stratified, multistage, or weighted surveys. It provides functions for parameter estimation, sample design, analysis of variance and regression, and calculation of standard errors.

  -   The `tidyverse` library is a collection of R packages used for data manipulation and visualization. It includes `dplyr`, `ggplot2`, `tidyr`, and others, characterized by its focus on 'tidy' or organized programming, making data exploration and analysis easier.

  -   The `srvyr` library is an extension of the `survey` library that integrates `survey` functions with `dplyr` syntax, facilitating the manipulation of complex survey data. It includes functions for grouping, filtering, and summarizing survey data using 'tidy' syntax.

  -   The `TeachingSampling` library is an R tool used for teaching statistical sampling methods. It includes functions for simulating different types of samples, estimating parameters, calculating standard errors, and constructing confidence intervals, among others.

  -   The `haven` library is an R tool that allows importing and exporting data in different formats, including SPSS, Stata, and SAS. It works with survey data files, providing functions for labeling variables, coding missing data, and converting data across formats.

  -   The `bayesplot` library is an R tool used for visualization and diagnostics of Bayesian models. It includes functions for plotting posterior distributions, convergence diagnostics, residual diagnostic plots, and other graphics related to Bayesian analysis.

  -   The `patchwork` library is an R tool that allows simple and flexible combination of plots. This library makes creating complex plots easier by allowing the combination of multiple plots into a single visualization, especially useful in data analysis and modeling.

  -   The `stringr` library is an R tool used for string manipulation. It includes functions for extracting, manipulating, and modifying text strings, particularly useful in data cleaning and preparation before analysis.

  -   The `rstan` library is an R tool used for Bayesian model estimation using the Markov Chain Monte Carlo (MCMC) method. This library allows specifying and estimating complex models using a simple and flexible language, offering various tools for diagnosis and visualization of results.
  
  

```r
library(survey)
library(tidyverse)
library(srvyr)
library(TeachingSampling)
library(haven)
library(bayesplot)
library(patchwork)
library(stringr)
library(rstan)
```

## Reading the survey and direct estimates

This code performs several operations on a labor survey in Jamaica, represented by the `encuesta` object, which is read from a file in RDS format. Here's the breakdown:

1. **Data Reading:** The code reads data from the Jamaican labor survey from an RDS file located at 'Resources/05_Employment/01_data_JAM.rds'. The data is stored in the `encuesta` object.

2. **Data Transformation:** Through a sequence of operations using the `%>%` pipe and the `transmute()` function from the `dplyr` package, the following transformations are performed on the survey:
    - Specific columns `dam2`, `RFACT`, `PAR_COD`, `CONST_NUMBER`, `ED_NUMBER`, `STRATA`, and `EMPSTATUS` are selected from the survey.
    - More descriptive names are assigned to some columns such as `fep` for `RFACT`, `upm` by combining `PAR_COD`, `CONST_NUMBER`, and `ED_NUMBER`, `estrato` using conditions to define the value based on `STRATA`, and `empleo_label` and `empleo` representing specific categories derived from `EMPSTATUS` with labeled levels and categorical values.

In summary, the code reads a labor survey in Jamaica and performs a series of transformations on selected columns, renaming and reorganizing them for future analyses or processing.



```r
encuesta <- readRDS('Recursos/05_Empleo/01_data_JAM.rds')
## 
id_dominio <- "dam2"

encuesta <-
  encuesta %>%
  transmute(
    dam2,
    fep = RFACT,
    upm = paste0(PAR_COD , CONST_NUMBER, ED_NUMBER),
    estrato = ifelse(is.na(STRATA) ,strata,STRATA),
    empleo_label = as_factor(EMPSTATUS ,levels  = "labels"),
    empleo = as_factor(EMPSTATUS ,levels  = "values") 
  )
```

The presented code defines the sampling design for the analysis of the "survey" in R. The first line sets an option for handling singleton PSU (primary sampling units), indicating that adjustments need to be applied in standard error calculations. The second line uses the "as_survey_design" function from the "survey" library to define the sampling design. The function takes "encuesta" as an argument and the following parameters:

  -   `strata`: The variable defining the strata in the survey, in this case, the "estrato" variable.
  
  -   `ids`: The variable identifying the PSUs in the survey, here, the "upm" variable.
  
  -   `weights`: The variable indicating the survey weights of each observation, in this case, the "fep" variable.
  
  -   `nest`: A logical parameter indicating whether the survey data is nested or not. In this case, it's set to "TRUE" because the data is nested by domain.

Together, these steps allow defining a sampling design that takes into account the sampling characteristics and the weights assigned to each observation in the survey. This is necessary to obtain precise and representative estimations of the parameters of interest.



```r
options(survey.lonely.psu= 'adjust' )
diseno <- encuesta %>%
  as_survey_design(
    strata = estrato,
    ids = upm,
    weights = fep,
    nest=T
  )
```

The following code conducts a descriptive analysis based on a survey design represented by the object `diseno`.

1. **Grouping and Filtering:** It uses the `%>%` function to chain operations. Initially, it groups the data by the domain identifier (`id_dominio`) using `group_by_at()` and subsequently filters observations where the variable `empleo` falls within the range of 3 to 5.

2. **Variable Summary:** With the `summarise()` function, it computes various summaries for different categories of the variable `empleo`. These summaries include the weighted count for employed, unemployed, and inactive individuals (`n_ocupado`, `n_desocupado`, `n_inactivo`). Furthermore, it utilizes the `survey_mean()` function to obtain weighted mean estimates for each category of `empleo`, considering the variable type (`vartype`) and design effect (`deff`).



```r
indicador_dam <-
  diseno %>% group_by_at(id_dominio) %>% 
  filter(empleo %in% c(3:5)) %>%
  summarise(
    n_ocupado = unweighted(sum(empleo == 3)),
    n_desocupado = unweighted(sum(empleo == 4)),
    n_inactivo = unweighted(sum(empleo == 5)),

    Ocupado = survey_mean(empleo == 3,
      vartype = c("se",  "var"),
      deff = T
    ),
    Desocupado = survey_mean(empleo == 4,
                          vartype = c("se",  "var"),
                          deff = T
    ),
    Inactivo = survey_mean(empleo == 5,
                          vartype = c("se",  "var"),
                          deff = T
    )
  )
```

3. **Upms counts by domains:**  This code performs operations on the survey data. First, it selects the columns id_dominio and upm, removes duplicate rows, and then counts the number of unique upm values for each id_dominio. Subsequently, it performs an inner join of these results with an existing object indicador_dam based on the id_dominio column, thus consolidating information about the quantity of unique upm values per identified domain in the survey.


```r
indicador_dam <- encuesta %>% select(id_dominio, upm) %>%
  distinct() %>% 
  group_by_at(id_dominio) %>% 
  tally(name = "n_upm") %>% 
  inner_join(indicador_dam, by = id_dominio)
#Save data----------------------------
saveRDS(indicador_dam,'Recursos/05_Empleo/indicador_dam.Rds' )
```


## Domain Selection

After conducting the necessary validations, the rule is set to include in the study the domains that have:

- Two or more PSUs per domain.
- An estimated design effect greater than 1 in all categories.



```r
indicador_dam1 <- indicador_dam %>%
  filter(n_upm >= 2,
         Desocupado_deff > 1,
         Ocupado_deff > 1,
         Inactivo_deff > 1)  %>%
  mutate(id_orden = 1:n())

saveRDS(object = indicador_dam1, "Recursos/05_Empleo/02_base_modelo.Rds")
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrr}
\toprule
dam2 & n\_upm & n\_ocupado & n\_desocupado & n\_inactivo & Ocupado & Ocupado\_se\\
\midrule
0101 & 17 & 784 & 40 & 487 & 0.6000 & 0.0188\\
0201 & 11 & 733 & 40 & 420 & 0.6033 & 0.0213\\
0202 & 11 & 617 & 25 & 332 & 0.6448 & 0.0470\\
0203 & 12 & 645 & 17 & 363 & 0.6409 & 0.0184\\
0204 & 11 & 638 & 69 & 402 & 0.5574 & 0.0258\\
\addlinespace
0207 & 11 & 473 & 22 & 374 & 0.5350 & 0.0193\\
0209 & 8 & 434 & 4 & 307 & 0.5666 & 0.0405\\
0211 & 13 & 713 & 11 & 396 & 0.6232 & 0.0314\\
0212 & 9 & 639 & 49 & 301 & 0.6179 & 0.0274\\
0302 & 10 & 676 & 91 & 327 & 0.5907 & 0.0285\\
\bottomrule
\end{tabular}
\end{table}


## Modeling in `STAN`

The code presents the implementation of a multinomial logistic response area model using the programming language `STAN`. In this model, it is assumed that the response variable in each domain follows a multinomial distribution. The parameters governing the relationship between the predictor variables and the response variable are assumed to be different in each domain and are modeled as random effects.

The *functions* section defines an auxiliary function called `pred_theta()`, used to predict the values of the response variable in the unobserved domains. The `data` section contains the model's input variables, including the number of domains, the number of response variable categories, direct estimates of the response variable in each domain, observed covariates in each domain, and covariates corresponding to the unobserved domains.

The *parameters* section defines the model's unknown parameters, including the *beta* parameter matrix, containing coefficients relating covariates to the response variable in each category. Standard deviations of the random effects are also included.

The *transformed parameters* section defines the `theta` parameter vector, containing the probabilities of belonging to each category of the response variable in each domain. Random effects are used to adjust the values of `theta` in each domain.

The *model* section defines the model structure and includes prior distributions for the unknown parameters. Particularly, a normal distribution is used for the coefficients of the beta matrix. Finally, it calculates the likelihood function of the multinomial distribution for the direct estimates of the response variable in each domain.

The *generated quantities* section is used to compute predictions of the response variable in the unobserved domains using the previously defined auxiliary function.


```
functions {
  matrix pred_theta(matrix Xp, int p, matrix beta){
  int D1 = rows(Xp);
  real num1[D1, p];
  real den1[D1];
  matrix[D1,p] theta_p;
  matrix[D1,p] tasa_pred;
  
  for(d in 1:D1){
    num1[d, 1] = 1;
    num1[d, 2] = exp(Xp[d, ] * beta[1, ]' ) ;
    num1[d, 3] = exp(Xp[d, ] * beta[2, ]' ) ;
    
    den1[d] = sum(num1[d, ]);
  }
  
  for(d in 1:D1){
    for(i in 2:p){
    theta_p[d, i] = num1[d, i]/den1[d];
    }
    theta_p[d, 1] = 1/den1[d];
   }

for(d in 1:D1){
    tasa_pred[d, 1] = theta_p[d,2]/(theta_p[d,1] + theta_p[d,2]);// TD
    tasa_pred[d, 2] = theta_p[d,1];                              // TO
    tasa_pred[d, 3] = theta_p[d,1] + theta_p[d,2];               // TP
    }

  return tasa_pred  ;
  }
  
}

data {
  int<lower=1> D; // número de dominios 
  int<lower=1> P; // categorías
  int<lower=1> K; // cantidad de regresores
  int hat_y[D, P]; // matriz de datos
  matrix[D, K] X_obs; // matriz de covariables
  int<lower=1> D1; // número de dominios 
  matrix[D1, K] X_pred; // matriz de covariables
}
  

parameters {
  matrix[P-1, K] beta;// matriz de parámetros 
  vector<lower=0>[P-1] sigma_u;       // random effects standard deviations
  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[P-1] L_u;
  matrix[P-1, D] z_u;                  
}

transformed parameters {
  simplex[P] theta[D];// vector de parámetros;
  real num[D, P];
  real den[D];
  matrix[D,P] tasa_obs;
  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[P-1, D] u; // random effect matrix
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
  
  for(d in 1:D){
    num[d, 1] = 1;
    num[d, 2] = exp(X_obs[d, ] * beta[1, ]' + u[1, d]) ;
    num[d, 3] = exp(X_obs[d, ] * beta[2, ]' + u[2, d]) ;
    
    den[d] = sum(num[d, ]);
  }
  
  for(d in 1:D){
    for(p in 2:P){
    theta[d, p] = num[d, p]/den[d];
    }
    theta[d, 1] = 1/den[d];
  }
  
  for(d in 1:D){
    tasa_obs[d, 1] = theta[d,2]/(theta[d,1] + theta[d,2]);// TD
    tasa_obs[d, 2] = theta[d,1];                                // TO
    tasa_obs[d, 3] = theta[d,1] + theta[d,2];               // TP
    }

}

model {
  L_u ~ lkj_corr_cholesky(1); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 10000);
  // sigma_u ~ cauchy(0, 50);
  sigma_u ~ inv_gamma(0.0001, 0.0001);
  
  for(p in 2:P){
    for(k in 1:K){
      beta[p-1, k] ~ normal(0, 10000);
    }
    }
  
  for(d in 1:D){
    target += multinomial_lpmf(hat_y[d, ] | theta[d, ]); 
  }
}

  
generated quantities {
  matrix[D1,P] tasa_pred;
  matrix[2, 2] Omega;
  Omega = L_u * L_u'; // so that it return the correlation matrix
  
 tasa_pred = pred_theta(X_pred, P, beta);
}

```

## Preparing supplies for `STAN`

   1. Reading and adaptation of covariates
  

```r
statelevel_predictors_df <- 
  readRDS('Recursos/05_Empleo/03_statelevel_predictors_dam.rds') %>% 
  mutate(id_orden =1:n())

head(statelevel_predictors_df[,1:9],10) %>% tba()
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrrrr}
\toprule
dam2 & area1 & sex2 & age2 & age3 & age4 & age5 & tiene\_sanitario & tiene\_electricidad\\
\midrule
0101 & 1.0000 & 0.5087 & 0.2694 & 0.2297 & 0.1689 & 0.0672 & 0.0019 & 0.7596\\
0102 & 1.0000 & 0.4754 & 0.2857 & 0.2261 & 0.1527 & 0.0683 & 0.0011 & 0.9064\\
0103 & 1.0000 & 0.5037 & 0.3095 & 0.2015 & 0.1312 & 0.0449 & 0.0152 & 0.6930\\
0201 & 0.5147 & 0.5060 & 0.2962 & 0.2090 & 0.1844 & 0.0711 & 0.0138 & 0.2342\\
0202 & 0.9986 & 0.5376 & 0.2625 & 0.2226 & 0.2238 & 0.0961 & 0.0028 & 0.3852\\
\addlinespace
0203 & 0.9754 & 0.5432 & 0.2454 & 0.2254 & 0.2388 & 0.1160 & 0.0015 & 0.3326\\
0204 & 1.0000 & 0.5300 & 0.3151 & 0.2022 & 0.2034 & 0.0776 & 0.0042 & 0.5720\\
0205 & 1.0000 & 0.5182 & 0.3057 & 0.2286 & 0.1981 & 0.0768 & 0.0013 & 0.8060\\
0206 & 1.0000 & 0.5157 & 0.3192 & 0.1959 & 0.1552 & 0.0496 & 0.0290 & 0.0285\\
0207 & 1.0000 & 0.5097 & 0.3099 & 0.1966 & 0.1691 & 0.0538 & 0.0465 & 0.1581\\
\bottomrule
\end{tabular}
\end{table}
  
2. Select the model variables and create a covariate matrix.



```r
names_cov <-
  c(
   "dam2",
    "ODDJOB","WORKED",
    "stable_lights_mean",
    "accessibility_mean",
    "urban.coverfraction_sum",
    "id_orden"
  )
X_pred <-
  anti_join(statelevel_predictors_df %>% select(all_of(names_cov)),
            indicador_dam1 %>% select(dam2))
```

The code block identifies which domains will be the predicted ones.


```r
X_pred %>% select(dam2, id_orden) %>% 
  saveRDS(file = "Recursos/05_Empleo/dam_pred.rds")
```

Creating the covariate matrix for the unobserved (`X_pred`) and observed (`X_obs`) domains
  

```r
## Obteniendo la matrix 
X_pred %<>%
  data.frame() %>%
  select(-dam2,-id_orden)  %>%  as.matrix()

## Identificando los dominios para realizar estimación del modelo

X_obs <- inner_join(indicador_dam1 %>% select(dam2),
                    statelevel_predictors_df[,names_cov]) %>%
  arrange(id_orden) %>%
  data.frame() %>%
  select(-dam2, -id_orden)  %>%  as.matrix()
```
  
3. Calculating the n_cash and the $\tilde{y}$
  

```r
D <- nrow(indicador_dam1)
P <- 3 # Ocupado, desocupado, inactivo.
Y_tilde <- matrix(NA, D, P)
n_tilde <- matrix(NA, D, P)
Y_hat <- matrix(NA, D, P)

# n efectivos ocupado
n_tilde[,1] <- (indicador_dam1$Ocupado*(1 - indicador_dam1$Ocupado))/
  indicador_dam1$Ocupado_var
Y_tilde[,1] <- n_tilde[,1]* indicador_dam1$Ocupado


# n efectivos desocupado
n_tilde[,2] <- (indicador_dam1$Desocupado*(1 - indicador_dam1$Desocupado))/
  indicador_dam1$Desocupado_var
Y_tilde[,2] <- n_tilde[,2]* indicador_dam1$Desocupado

# n efectivos Inactivo
n_tilde[,3] <- (indicador_dam1$Inactivo*(1 - indicador_dam1$Inactivo))/
  indicador_dam1$Inactivo_var
Y_tilde[,3] <- n_tilde[,3]* indicador_dam1$Inactivo
```

 Now, we validate the consistency of the calculations carried out
  

```r
ni_hat = rowSums(Y_tilde)
Y_hat[,1] <- ni_hat* indicador_dam1$Ocupado
Y_hat[,2] <- ni_hat* indicador_dam1$Desocupado
Y_hat[,3] <- ni_hat* indicador_dam1$Inactivo
Y_hat <- round(Y_hat)

hat_p <- Y_hat/rowSums(Y_hat)
par(mfrow = c(1,3))
plot(hat_p[,1],indicador_dam1$Ocupado)
abline(a = 0,b=1,col = "red")
plot(hat_p[,2],indicador_dam1$Desocupado)
abline(a = 0,b=1,col = "red")
plot(hat_p[,3],indicador_dam1$Inactivo)
abline(a = 0,b=1,col = "red")
```
  
![](Recursos/05_Empleo/04_validando.png)

4. Compiling the model



```r
X1_obs <- cbind(matrix(1,nrow = D,ncol = 1),X_obs)
K = ncol(X1_obs)
D1 <- nrow(X_pred)
X1_pred <- cbind(matrix(1,nrow = D1,ncol = 1),X_pred)

sample_data <- list(D = D,
                    P = P,
                    K = K,
                    y_tilde = Y_hat,
                    X_obs = X1_obs,
                    X_pred = X1_pred,
                    D1 = D1)


library(rstan)
model_bayes_mcmc2 <- stan(
   # Stan program
  file = "Recursos/05_Empleo/modelosStan/00 Multinomial_simple_no_cor.stan", 
  data = sample_data,    # named list of data
  verbose = TRUE,
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
)

saveRDS(model_bayes_mcmc2,
        "Recursos/05_Empleo/05_model_bayes_multinomial_cor.Rds")
```


## Model validation

Model validation is essential to assess a model's ability to accurately and reliably predict future outcomes. In the case of a multinomial response area model, validation focuses on measuring the model's accuracy in predicting different response categories. The main objective of validation is to determine if the model can generalize well to unseen data and provide accurate predictions. This involves comparing the model's predictions to observed data and using evaluation metrics to measure model performance. Model validation is crucial to ensure prediction quality and the model's reliability for use in future applications.



```r
infile <- paste0("Recursos/05_Empleo/05_model_bayes_multinomial_cor.Rds")
model_bayes <- readRDS(infile)

#--- Exporting Bayesian Multilevel Model Results ---#

paramtros <- summary(model_bayes)$summary %>% data.frame()

tbla_rhat <- mcmc_rhat_data(paramtros$Rhat) %>% 
  group_by(description) %>% 
  tally() %>% mutate(Porcen = n/sum(n)*100)

tbla_rhat %>% tba()
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrr}
\toprule
description & n & Porcen\\
\midrule
hat(R) <= 1.05 & 599 & 100\\
\bottomrule
\end{tabular}
\end{table}

### Fixed effects {-}


```r
efecto_fijo <-  grep(pattern = "beta", 
                     x = rownames(paramtros),
                     value = TRUE)

p_fijo <- traceplot(model_bayes, pars = efecto_fijo)
p_fijo
```
![](Recursos/05_Empleo/07_betas_fijos.png)

### Random effects {-}


```r
efecto_aleatorio <-  grep(pattern = "z_u", 
                          x = rownames(paramtros),
                          value = TRUE)

p_alea1 <- traceplot(model_bayes, pars = efecto_aleatorio[1:26])
p_alea2 <- traceplot(model_bayes, pars = efecto_aleatorio[27:(26*2)])
p_alea3 <- traceplot(model_bayes, pars = efecto_aleatorio[(26*2):78])
```

![](Recursos/05_Empleo/08a_trece_alea1.png)
![](Recursos/05_Empleo/08b_trece_alea2.png)
![](Recursos/05_Empleo/08c_trece_alea3.png)
### Estimated values for the correlation matrix {-}


```r
## Valores estimados para la matriz de correlación
omega12 <- summary(model_bayes, pars = "Omega[1,2]")$summary

plot_omega <- Plot_dens_draws(model_bayes, pars = "Omega[1,2]")
```

![](Recursos/05_Empleo/09_omega.png)


### Posterior predictive distribution


```r
theta_dir <- indicador_dam1 %>%  
  transmute(dam2,
    n = n_desocupado + n_ocupado + n_inactivo,
        Ocupado, Desocupado, Inactivo) 

color_scheme_set("brightblue")
theme_set(theme_bw(base_size = 15))
y_pred_B <- as.array(model_bayes, pars = "theta") %>%
  as_draws_matrix()
  
rowsrandom <- sample(nrow(y_pred_B), 100)

theta_1<-  grep(pattern = "1]",x = colnames(y_pred_B),value = TRUE)
theta_2<-  grep(pattern = "2]",x = colnames(y_pred_B),value = TRUE)
theta_3<-  grep(pattern = "3]",x = colnames(y_pred_B),value = TRUE)
y_pred1 <- y_pred_B[rowsrandom,theta_1 ]
y_pred2 <- y_pred_B[rowsrandom,theta_2 ]
y_pred3 <- y_pred_B[rowsrandom,theta_3 ]

p1 <- ppc_dens_overlay(y = as.numeric(theta_dir$Ocupado), y_pred1)/
  ppc_dens_overlay(y = as.numeric(theta_dir$Desocupado), y_pred2)/
  ppc_dens_overlay(y = as.numeric(theta_dir$Inactivo), y_pred3)
```


![](Recursos/05_Empleo/10_ppc.png)
  

## Parameter estimation.

This code block starts by importing an `dam_pred.rds` file and sets some variables for the number of parameters and the number of domains. Then, it extracts summaries from a model named `model_bayes` for observed and predicted rates. Next, it organizes these rates into ordered matrices, assigning appropriate column names and converting them into data frames. These data frames are concatenated with the `dam2` column from the original `indicador_dam1` and `dam_pred` data, respectively, creating two data frames (`tasa_obs_ordered` and `tasa_pred_ordered`) containing organized parameter estimates ready for further analysis.




```r
dam_pred <- readRDS("Recursos/05_Empleo/dam_pred.rds")
P <- 3 
D <- nrow(indicador_dam1)
D1 <- nrow(dam_pred)
## Estimación del modelo. 
theta_dir <- indicador_dam1 
tasa_obs <- summary(model_bayes,pars = "tasa_obs")$summary
tasa_pred <-  summary(model_bayes,pars = "tasa_pred")$summary

## Ordenando la matrix de theta 
tasa_obs_ordenado <- matrix(tasa_obs[,"mean"], 
                            nrow = D,
                            ncol = P,byrow = TRUE) 

colnames(tasa_obs_ordenado) <- c("TD_mod", "TO_mod", "TP_mod")
tasa_obs_ordenado%<>% as.data.frame()
tasa_obs_ordenado <- cbind(dam2 = indicador_dam1$dam2,
                           tasa_obs_ordenado)



tasa_pred_ordenado <- matrix(tasa_pred[,"mean"], 
                             nrow = D1,
                             ncol = P,byrow = TRUE)

colnames(tasa_pred_ordenado) <- c("TD_mod", "TO_mod", "TP_mod")
tasa_pred_ordenado%<>% as.data.frame()
tasa_pred_ordenado <- cbind(dam2 = dam_pred$dam2, tasa_pred_ordenado)

estimaciones_obs <- full_join(theta_dir,
                              bind_rows(tasa_obs_ordenado, tasa_pred_ordenado))
```

## Estimation of Standard Deviation and Coefficient of Variation

This code block computes the standard deviations (sd) and coefficients of variation (cv) for the `theta` parameters, both observed and predicted. Initially, the `summary()` function from the `rstan` package is used to extract the `sd` values for observed and predicted `theta` parameters from the Bayesian estimation model (`model_bayes`). Subsequently, the `sd` values are organized into matrices ordered by `dam2` and given corresponding names. Using these matrices, another matrix is generated containing the coefficients of variation for the observed `theta` parameters (`theta_obs_ordenado_cv`). Similarly, ordered matrices are constructed by `dam2` for both the `sd` and `cv` values of the predicted `theta` parameters (`theta_pred_ordenado_sd` and `theta_pred_ordenado_cv`, respectively).



```r
tasa_obs_ordenado_sd <- matrix(tasa_obs[,"sd"], 
                               nrow = D,
                               ncol = P,byrow = TRUE) 

colnames(tasa_obs_ordenado_sd) <- c("TD_mod_sd", "TO_mod_sd", "TP_mod_sd")
tasa_obs_ordenado_sd%<>% as.data.frame()
tasa_obs_ordenado_sd <- cbind(dam2 = indicador_dam1$dam2,
                              tasa_obs_ordenado_sd)
tasa_obs_ordenado_cv <- tasa_obs_ordenado_sd[,-1]/tasa_obs_ordenado[,-1]

colnames(tasa_obs_ordenado_cv) <- c("TD_mod_cv", "TO_mod_cv", "TP_mod_cv")

tasa_obs_ordenado_cv <- cbind(dam2 = indicador_dam1$dam2,
                              tasa_obs_ordenado_cv)

tasa_pred_ordenado_sd <- matrix(tasa_pred[,"sd"], 
                                nrow = D1,
                                ncol = P,byrow = TRUE)

colnames(tasa_pred_ordenado_sd) <- c("TD_mod_sd", "TO_mod_sd", "TP_mod_sd")
tasa_pred_ordenado_sd%<>% as.data.frame()
tasa_pred_ordenado_sd <- cbind(dam2 = dam_pred$dam2, tasa_pred_ordenado_sd)

tasa_pred_ordenado_cv <- tasa_pred_ordenado_sd[,-1]/tasa_pred_ordenado[,-1]

colnames(tasa_pred_ordenado_cv) <- c("TD_mod_cv", "TO_mod_cv", "TP_mod_cv")

tasa_pred_ordenado_cv <- cbind(dam2 = dam_pred$dam2, tasa_pred_ordenado_cv)
```

The last step is to consolidate the bases obtained for the point estimate, standard deviation and coefficient of variation.


```r
tasa_obs_ordenado <-
  full_join(tasa_obs_ordenado, tasa_obs_ordenado_sd) %>%
  full_join(tasa_obs_ordenado_cv)

tasa_pred_ordenado <-
  full_join(tasa_pred_ordenado, tasa_pred_ordenado_sd) %>%
  full_join(tasa_pred_ordenado_cv)


estimaciones_obs <- full_join(indicador_dam1,
                              bind_rows(tasa_obs_ordenado, tasa_pred_ordenado))



saveRDS(object = estimaciones_obs, file = "Recursos/05_Empleo/11_estimaciones.rds")
tba(head(estimaciones_obs[,1:7],10))
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrr}
\toprule
dam2 & n\_upm & n\_ocupado & n\_desocupado & n\_inactivo & Ocupado & Ocupado\_se\\
\midrule
0101 & 17 & 784 & 40 & 487 & 0.6000 & 0.0188\\
0201 & 11 & 733 & 40 & 420 & 0.6033 & 0.0213\\
0202 & 11 & 617 & 25 & 332 & 0.6448 & 0.0470\\
0203 & 12 & 645 & 17 & 363 & 0.6409 & 0.0184\\
0204 & 11 & 638 & 69 & 402 & 0.5574 & 0.0258\\
\addlinespace
0207 & 11 & 473 & 22 & 374 & 0.5350 & 0.0193\\
0209 & 8 & 434 & 4 & 307 & 0.5666 & 0.0405\\
0211 & 13 & 713 & 11 & 396 & 0.6232 & 0.0314\\
0212 & 9 & 639 & 49 & 301 & 0.6179 & 0.0274\\
0302 & 10 & 676 & 91 & 327 & 0.5907 & 0.0285\\
\bottomrule
\end{tabular}
\end{table}

## Metodología de Benchmarking 

1. People counts aggregated by dam2, people over 15 years of age.
  

```r
conteo_pp_dam <- readRDS("Recursos/05_Empleo/12_censo_mrp.rds") %>% 
   filter(age > 1)  %>% 
   mutate(dam = str_sub(dam2,1,2)) %>% 
  group_by(dam, dam2) %>% 
  summarise(pp_dam2 = sum(n),.groups = "drop") %>% 
  group_by(dam) %>% 
  mutate(pp_dam = sum(pp_dam2))
head(conteo_pp_dam) %>% tba()
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{llrr}
\toprule
dam & dam2 & pp\_dam2 & pp\_dam\\
\midrule
01 & 0101 & 26282 & 64195\\
01 & 0102 & 19385 & 64195\\
01 & 0103 & 18272 & 64195\\
01 & 0198 & 256 & 64195\\
02 & 0201 & 44833 & 443954\\
\addlinespace
02 & 0202 & 35177 & 443954\\
\bottomrule
\end{tabular}
\end{table}

2. Estimation of the `theta` parameter at the level that the survey is representative.
  

```r
indicador_agregado <-
  diseno %>% 
   mutate(dam = str_sub(dam2,1,2)) %>% 
  group_by(dam) %>% 
  filter(empleo %in% c(3:5)) %>%
  summarise(
    TO = survey_mean(empleo == 3,
                     vartype = c("se", "cv")),
    TD = survey_ratio(empleo == 4,
                      empleo != 5,
                      vartype = c("se",  "cv")),
    TP = survey_mean(empleo != 5,
                     vartype = c("se",  "cv"))
  ) %>% select(dam,TO,TD, TP)


tba(indicador_agregado)
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrr}
\toprule
dam & TO & TD & TP\\
\midrule
01 & 0.5990 & 0.0609 & 0.6378\\
02 & 0.6065 & 0.0462 & 0.6359\\
03 & 0.5722 & 0.1573 & 0.6790\\
04 & 0.5459 & 0.1137 & 0.6160\\
05 & 0.5132 & 0.1029 & 0.5720\\
\addlinespace
06 & 0.5328 & 0.1499 & 0.6267\\
07 & 0.5003 & 0.0683 & 0.5369\\
08 & 0.5741 & 0.1177 & 0.6506\\
09 & 0.5506 & 0.0726 & 0.5938\\
10 & 0.5815 & 0.0666 & 0.6230\\
\addlinespace
11 & 0.6059 & 0.0400 & 0.6311\\
12 & 0.5746 & 0.0909 & 0.6320\\
13 & 0.5753 & 0.0951 & 0.6358\\
14 & 0.5946 & 0.0852 & 0.6500\\
\bottomrule
\end{tabular}
\end{table}

Organizing the output as a vector.


```r
temp <-
  gather(indicador_agregado, key = "agregado",
         value = "estimacion", -dam) %>%
  mutate(nombre = paste0("dam_", dam,"_", agregado))

Razon_empleo <- setNames(temp$estimacion, temp$nombre)
```
  
3. Define the weights by domains.
  

```r
names_cov <-  "dam"
estimaciones_mod <- estimaciones %>% 
  transmute(
    dam = str_sub(dam2,1,2),
    dam2,
    TO_mod,TD_mod,TP_mod) %>% 
  inner_join(conteo_pp_dam ) %>% 
  mutate(wi = pp_dam2/pp_dam)
```
  
 4. Create dummy variables
  

```r
estimaciones_mod %<>%
  fastDummies::dummy_cols(select_columns = names_cov,
                          remove_selected_columns = FALSE)

Xdummy <- estimaciones_mod %>% select(matches("dam_")) %>% 
  mutate_at(vars(matches("_\\d")) ,
            list(TO = function(x) x*estimaciones_mod$TO_mod,
                 TD = function(x) x*estimaciones_mod$TD_mod,
                 TP = function(x) x*estimaciones_mod$TP_mod)) %>% 
  select((matches("TO|TD|TP"))) 

# head(Xdummy) %>% tba()
```

Some validations carried out


```r
colnames(Xdummy) == names(Razon_empleo)
data.frame(Modelo = colSums(Xdummy*estimaciones_mod$wi),
Estimacion_encuesta = Razon_empleo)
```
  
  
 5. Calculate the weight for each level of the variable.
  
#### Occupancy Rate {-}
    

```r
library(sampling)
names_ocupado <- grep(pattern = "TO", x = colnames(Xdummy),value = TRUE)

gk_TO <- calib(Xs = Xdummy[,names_ocupado], 
            d =  estimaciones_mod$wi,
            total = Razon_empleo[names_ocupado],
            method="logit",max_iter = 5000,) 

checkcalibration(Xs = Xdummy[,names_ocupado], 
                 d =estimaciones_mod$wi,
                 total = Razon_empleo[names_ocupado],
                 g = gk_TO)
```

#### Unemployment Rate {-} 
    

```r
names_descupados <- grep(pattern = "TD", x = colnames(Xdummy),value = TRUE)

gk_TD <- calib(Xs = Xdummy[,names_descupados], 
                    d =  estimaciones_mod$wi,
                    total = Razon_empleo[names_descupados],
                    method="logit",max_iter = 5000,) 

checkcalibration(Xs = Xdummy[,names_descupados], 
                 d =estimaciones_mod$wi,
                 total = Razon_empleo[names_descupados],
                 g = gk_TD)
```

#### Participation Rate {-}


```r
names_inactivo <- grep(pattern = "TP", x = colnames(Xdummy),value = TRUE)

gk_TP <- calib(Xs = Xdummy[,names_inactivo], 
                    d =  estimaciones_mod$wi,
                    total = Razon_empleo[names_inactivo],
                    method="logit",max_iter = 5000,) 

checkcalibration(Xs = Xdummy[,names_inactivo], 
                 d =estimaciones_mod$wi,
                 total = Razon_empleo[names_inactivo],
                 g = gk_TP)
```
  
6. Validate the results obtained.
  

```r
par(mfrow = c(1,3))
hist(gk_TO)
hist(gk_TD)
hist(gk_TP)
```


![](Recursos/05_Empleo/13_plot_gks.png)  


7. Estimates adjusted by the weighter
  

```r
estimacionesBench <- estimaciones_mod %>%
  mutate(gk_TO, gk_TD, gk_TP) %>%
  transmute(
    dam,dam2,
    wi,gk_TO, gk_TD, gk_TP,
    TO_Bench = TO_mod*gk_TO,
    TD_Bench = TD_mod*gk_TD,
    TP_Bench = TP_mod*gk_TP
  ) 
```

8. Validation of results.
  

```r
tabla_validar <- estimacionesBench %>%
  group_by(dam) %>% 
  summarise(TO_Bench = sum(wi*TO_Bench),
            TD_Bench = sum(wi*TD_Bench),
            TP_Bench = sum(wi*TP_Bench)) %>% 
  inner_join(indicador_agregado)
tabla_validar %>% tba()
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrr}
\toprule
dam & TO\_Bench & TD\_Bench & TP\_Bench & TO & TD & TP\\
\midrule
01 & 0.5990 & 0.0609 & 0.6378 & 0.5990 & 0.0609 & 0.6378\\
02 & 0.6065 & 0.0462 & 0.6359 & 0.6065 & 0.0462 & 0.6359\\
03 & 0.5722 & 0.1573 & 0.6790 & 0.5722 & 0.1573 & 0.6790\\
04 & 0.5459 & 0.1137 & 0.6160 & 0.5459 & 0.1137 & 0.6160\\
05 & 0.5132 & 0.1029 & 0.5720 & 0.5132 & 0.1029 & 0.5720\\
\addlinespace
06 & 0.5328 & 0.1499 & 0.6267 & 0.5328 & 0.1499 & 0.6267\\
07 & 0.5003 & 0.0683 & 0.5369 & 0.5003 & 0.0683 & 0.5369\\
08 & 0.5741 & 0.1177 & 0.6506 & 0.5741 & 0.1177 & 0.6506\\
09 & 0.5506 & 0.0726 & 0.5938 & 0.5506 & 0.0726 & 0.5938\\
10 & 0.5815 & 0.0666 & 0.6230 & 0.5815 & 0.0666 & 0.6230\\
\addlinespace
11 & 0.6059 & 0.0400 & 0.6311 & 0.6059 & 0.0400 & 0.6311\\
12 & 0.5746 & 0.0909 & 0.6320 & 0.5746 & 0.0909 & 0.6320\\
13 & 0.5753 & 0.0951 & 0.6358 & 0.5753 & 0.0951 & 0.6358\\
14 & 0.5946 & 0.0852 & 0.6500 & 0.5946 & 0.0852 & 0.6500\\
\bottomrule
\end{tabular}
\end{table}


9. Save results
  

```r
estimaciones <- inner_join(estimaciones,estimacionesBench)
saveRDS(object = estimaciones, file = "Recursos/05_Empleo/15_estimaciones_Bench.rds")
```



### Grafico de validación del Benchmarking

1. Perform model estimates before and after Benchmarking


```r
estimaciones_agregada <- estimaciones %>%
  group_by(dam) %>% 
  summarise(
    TO_mod = sum(wi * TO_mod),
    TD_mod = sum(wi * TD_mod),
    TP_mod = sum(wi * TP_mod),
    TO_bench = sum(wi * TO_Bench),
    TD_bench = sum(wi * TD_Bench),
    TP_bench = sum(wi * TP_Bench))
```

2. Obtain the confidence intervals for the direct estimates


```r
indicador_agregado <-
  diseno %>% 
   mutate(dam = str_sub(dam2,1,2)) %>% 
  group_by(dam) %>% 
  filter(empleo %in% c(3:5)) %>%
  summarise(
     nd = unweighted(n()),
    TO = survey_mean(empleo == 3,
                     vartype = c("ci")),
    TD = survey_ratio(empleo == 4,
                      empleo != 5,
                      vartype = c("ci")),
    TP = survey_mean(empleo != 5,
                     vartype = c("ci"))
  ) 

data_plot <- left_join(estimaciones_agregada, indicador_agregado)
```

3. Select the results for an indicator (Occupancy rate)


```r
temp_TO <- data_plot %>% select(dam,nd, starts_with("TO"))


temp_TO_1 <- temp_TO %>% select(-TO_low, -TO_upp) %>%
  gather(key = "Estimacion",value = "value", -nd,-dam) %>% 
  mutate(Estimacion = case_when(Estimacion == "TO_mod" ~ "Area model",
                                Estimacion == "TO_bench" ~ "Area model (bench)",
                                Estimacion == "TO"~ "Direct estimator"))

lims_IC_ocupado <-  temp_TO %>%
  select(dam,nd,value = TO,TO_low, TO_upp) %>% 
  mutate(Estimacion = "Direct estimator")
```


4. Make the graph


```r
p_TO <- ggplot(temp_TO_1,
                    aes(
                      x = fct_reorder2(dam, dam, nd),
                      y = value,
                      shape = Estimacion,
                      color = Estimacion
                    )) +
  geom_errorbar(
    data = lims_IC_ocupado,
    aes(ymin = TO_low ,
        ymax = TO_upp, x = dam),
    width = 0.2,
    linewidth = 1
  )  +
  geom_jitter(size = 3) +
  labs(x = "Dam", title = "Occupancy Rate", y = "", 
       color= "Estimation", shape = "Estimation")
```

![](Recursos/05_Empleo/16_plot_uni_TO.png) 

5. Repeat the process with the other indicators.

![](Recursos/05_Empleo/16_plot_uni_TD.png) 
![](Recursos/05_Empleo/16_plot_uni_TP.png) 

## Labor market maps.

The code loads the libraries `sf` and `tmap`. Subsequently, it reads a shapefile containing geographic information and uses the 'inner_join' function to merge it with previously calculated survey estimates.



```r
library(sf)
library(tmap)
ShapeSAE <- read_sf("Shapefile/JAM2_cons.shp")

P1_empleo <- tm_shape(ShapeSAE %>%
                           inner_join(estimaciones))
```

### Occupancy Rate {-}

The following code creates the map of the occupancy rate using the tmap library. It begins by setting some global mapping options with tmap_options. Then, it uses the P1_employment variable as the base and adds polygon layers using tm_polygons. In this layer, the variable "TO_Bench" is defined as representing the occupancy rate, a title is assigned to the map, a color palette ("Blues") is chosen, and the data classification style ("quantile") is set. Subsequently, with tm_layout, various aspects of the map's presentation are defined, such as the position and size of the legend, the aspect ratio, and the text size of the legend. Finally, the resulting map is displayed.


```r
tmap_options(check.and.fix = TRUE)
Mapa_TO <-
  P1_empleo +
  tm_polygons("TO_Bench",
          title = "Occupancy Rate",
          palette = "-Blues",
          colorNA = "white",
          style= "quantile") +
  tm_layout( 
    legend.only = FALSE,
    legend.height = -0.3,
    legend.width = -0.5,
    asp = 1.5,
    legend.text.size = 3,
    legend.title.size = 3)

Mapa_TO
```


![](Recursos/05_Empleo/17_map_TO.png) 

### Unemployment Rate {-} 



```r
Mapa_TD <-
  P1_empleo + tm_polygons(
    "TD_Bench",
    title =  "Unemployment Rate",
    palette = "YlOrRd",
    colorNA = "white",
    style= "quantile"
  ) + tm_layout( 
    legend.only = FALSE,
    legend.height = -0.3,
    legend.width = -0.5,
    asp = 1.5,
    legend.text.size = 3,
    legend.title.size = 3)


Mapa_TD
```

![](Recursos/05_Empleo/17_map_TD.png) 

### Inactivo {-} 



```r
Mapa_TP <-
  P1_empleo + tm_polygons(
    "TP_Bench",
    title =  "Participation Rate",
    colorNA = "white",
    palette = "YlGn",
    style= "quantile"
  ) + tm_layout( 
    legend.only = FALSE,
    legend.height = -0.3,
    legend.width = -0.5,
    asp = 1.5,
    legend.text.size = 3,
    legend.title.size = 3)

Mapa_TP
```

![](Recursos/05_Empleo/17_map_TP.png) 
    
  
