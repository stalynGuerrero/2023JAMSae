# Session 2- Generalized Variance Function



\justify 

One of the most important inputs in the area model is the variance of the direct estimator at the domain level, which cannot be calculated directly. Accordingly, this value must be estimated from the data collected in each domain. However, in domains with very small sample sizes, these estimations will not perform well. Hence, it is very useful to use a **smoothing** model for variances to eliminate noise and volatility from these estimations and extract the true signal of the process.

Hidiroglou (2019) states that $E_{{MP}}\left(\hat{\theta}^{dir}_d\right)$
$=\boldsymbol{x}^{T}_{d}\boldsymbol{\beta}$ and $V_{\mathscr{MP}}\left(\hat{\theta}^{dir}_d\right)=\sigma_{u}^2+\tilde{\sigma}^2_{d}$, where the subscript $\mathscr{MP}$ refers to the double inference that must be taken into account in these adjustments and defines the joint probability measure between the model and the sampling design.

-   $\mathscr{M}$ refers to the probability measure induced by modeling and the inclusion of auxiliary covariates ($\boldsymbol{x}_{d}$).

-   $\mathscr{MP}$ refers to the probability measure induced by the complex sampling design that yields direct estimations.

The solution proposed here is known as the Generalized Variance Function, which involves fitting a log-linear model to the estimated direct variance. Starting from the fact that an unbiased estimator of $\sigma^2$ denoted by $\hat{\sigma}^2$ is available, it follows that:
$$
E_{\mathscr{MP}}\left(\hat{\sigma}_{d}^{2}\right)=E_{\mathscr{M}}\left(E_{\mathscr{P}}\left(\hat{\sigma}_{d}^{2}\right)\right)=E_{\mathscr{M}}\left(\sigma_{d}^{2}\right)=\tilde{\sigma}_{d}^{2}
$$

The above equality can be interpreted as an unbiased and simple estimator of $\tilde{\sigma}_{d}^{2}$, denoted as $\hat{\sigma}_{d}^{2}$. However, this sampling estimator is unstable when the sample size is small, which is precisely the dominant paradigm in small area estimation. Rivest and Belmonte (2000) consider smoothing models for the estimation of direct variances defined as follows:

$$
\log\left(\hat{\sigma}_{d}^{2}\right)=\boldsymbol{z}_{d}^{T}\boldsymbol{\alpha}+\boldsymbol{\varepsilon}_{d}
$$

Where $\boldsymbol{z}_{d}$ is a vector of explanatory covariates that are functions of $\boldsymbol{x}_{d}$, $\boldsymbol{\alpha}$ is a vector of parameters to be estimated, $\boldsymbol{\varepsilon}_{d}$ are random errors with zero mean and constant variance, assumed to be identically distributed conditionally on $\boldsymbol{z}_{d}$. From the above model, the smoothed estimation of the sampling variance is given by:

$$
\tilde{\sigma}_{d}^{2}=E_{\mathscr{MP}}\left(\sigma_{d}^{2}\right)=\exp\left(\boldsymbol{z}_{d}^{T}\boldsymbol{\alpha}\right)\times\Delta
$$

Where $E_{\mathscr{MP}}\left(\varepsilon_{d}\right)=\Delta$. There's no need to specify a parametric distribution for the errors of this model. Using the method of moments, the following unbiased estimator for $\Delta$ is obtained:

$$
\hat{\Delta}=\frac{\sum_{d=1}^{D}\hat{\sigma}_{d}^{2}}{\sum_{d=1}^{D}\exp\left(\boldsymbol{z}_{d}^{T}\boldsymbol{\alpha}\right)}
$$


Similarly, using standard procedures in linear regression, the estimation of the regression parameter coefficients is given by the following expression:

$$
\hat{\boldsymbol{\alpha}}=\left(\sum_{d=1}^{D}\boldsymbol{z}_{d}\boldsymbol{z}_{d}^{T}\right)^{-1}\sum_{d=1}^{D}\boldsymbol{z}_{d}\log\left(\hat{\sigma}_{d}^{2}\right)
$$

Finally, the smoothed estimator of the sampling variance is defined as:

$$
\hat{\tilde{\sigma}}_{d}^{2}=\exp\left(\boldsymbol{z}_{d}^{T}\hat{\boldsymbol{\alpha}}\right)\hat{\Delta}
$$

**Survey Data:**

The following code processes data using various R packages such as `survey`, `tidyverse`, `srvyr`, `TeachingSampling`, and `haven`.

1. **Library Loading:** The code loads the necessary libraries (`survey`, `tidyverse`, `srvyr`, `TeachingSampling`, `haven`) required for data manipulation and analysis.

2. **Survey Data Set Reading:** The code reads the dataset named 'data_JAM.rds' using the `readRDS` function.

3. **Data Manipulation:**
   - It creates new variables (`dam2`, `fep`, `upm`, `estrato`, `ingreso`, `pobreza`) using the `mutate` function from the `dplyr` package.
   - Filters the dataset based on a specific condition using the `filter` function.



```r
library(survey)
library(tidyverse)
library(srvyr)
library(TeachingSampling)
library(haven)
 
#read in data set

# Q518A: Gross average income From Employment
  
encuesta <- read_rds("Recursos/02_FGV/01_data_JAM.rds") 

encuesta <-
  encuesta %>%
  mutate(
    dam2,
    fep = RFACT/4,
    upm = paste0(PAR_COD , CONST_NUMBER, ED_NUMBER),
    estrato = ifelse(is.na(STRATA) ,strata,STRATA),
    ingreso = case_when(!is.na(Q518A) & !is.na(Q518B) ~ Q518A + Q518B, 
                        !is.na(Q518A) & !is.na(Q518B) ~ NA_real_,
                        !is.na(Q518A) & is.na(Q518B) ~ Q518A,
                        is.na(Q518A) & !is.na(Q518B) ~ Q518B),
    pobreza = ifelse(ingreso < 50,1,0)
  ) %>% filter(ingreso > 11)
```

*dam2*: Corresponds to the code assigned to the country's second administrative division.

**The income definition is structured in this manner to illustrate the process utilized in the small area estimation methodology for poverty estimation.**


\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrlll}
\toprule
ID & QUARTER & EMPSTATUS & STRATA & CLUSTER & CONST\_NUMBER\\
\midrule
02010105902170101 & 3 & 5 & 014 & 014011 & 01\\
02020207900130103 & 3 & 5 & 235 & 235472 & 02\\
02020404901820105 & 3 & 5 & 239 & 239488 & 04\\
04010400700370102 & 3 & 5 & 036 & 036073 & 04\\
04030406400160106 & 3 & 5 & 221 & 221453 & 04\\
\addlinespace
05020506600550101 & 3 & 5 & 164 & 164352 & 05\\
06030106300920101 & 3 & 5 & 122 & 122227 & 01\\
06030200300140102 & 3 & 5 & 228 & 228464 & 02\\
06030204100960101 & 3 & 5 & 084 & 084180 & 02\\
06030208501550102 & 3 & 5 & 234 & 234474 & 02\\
\bottomrule
\end{tabular}
\end{table}

In the following code block, the libraries `survey` and `srvyr` are used to create a sampling design from a survey database. The sampling design encompasses information about primary sampling units (PSUs), sampling weights (wkx), and strata (estrato) utilized in the sampling. Additionally, the "survey.lonely.psu" option is employed to adjust sample sizes within groups of primary sampling units that lack other primary sampling units within the same group.


```r
library(survey)
library(srvyr)
options(survey.lonely.psu = "adjust")
id_dominio <- "dam2"

diseno <-
  as_survey_design(
    ids = upm,
    weights = fep,
    strata = estrato,
    nest = TRUE,
    .data = encuesta
  )

#summary(diseno)
```
´

1. **Indicator Calculation:**
   - Groups the sampling design by domain ID and calculates indicators related to poverty (`n_pobreza` and `pobreza`). `n_pobreza` counts the number of instances where `pobreza` equals 1 (indicating poverty), while `pobreza` computes the survey mean of the `pobreza` variable with specific variance estimations.



```r
# Calculating indicators related to poverty and counts of UPMS per domain

# Indicator Calculation:
# Grouping the sampling design by domain ID and summarizing variables 
# related to poverty.
indicador_dam <-
  diseno %>% group_by_at(id_dominio) %>% 
  summarise(
    n_pobreza = unweighted(sum(pobreza == 1)),
    pobreza = survey_mean(pobreza,
                          vartype = c("se",  "var"),
                          deff = T
    )
  )
```

2. **Counts of UPMS per Domain:**
   - Extracts domain ID and UPMs from the survey dataset, obtaining unique UPMs per domain and counting them.
   - Joins the count of unique UPMs per domain with the previously calculated indicators related to poverty.


```r
# Counts of UPMS per domain:
# Selecting domain ID and UPMs, obtaining unique UPMs per domain, 
# and counting them.
# Joining the count of unique UPMs per domain with previously calculated 
# indicators related to poverty.
indicador_dam <- encuesta %>% select(id_dominio, upm) %>%
  distinct() %>% 
  group_by_at(id_dominio) %>% 
  tally(name = "n_upm") %>% 
  inner_join(indicador_dam, by = id_dominio)
# saveRDS(directodam2, "Recursos/02_FGV/indicador_dam.Rds")
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrr}
\toprule
dam2 & n\_upm & n\_pobreza & pobreza & pobreza\_se & pobreza\_var & pobreza\_deff\\
\midrule
0101 & 17 & 237 & 0.9980 & 0.0021 & 0.0000 & 5.228000e-01\\
0102 & 12 & 197 & 0.9836 & 0.0089 & 0.0001 & 1.053700e+00\\
0103 & 9 & 65 & 1.0000 & 0.0000 & 0.0000 & 3.065785e+32\\
0201 & 11 & 244 & 1.0000 & 0.0000 & 0.0000 & NaN\\
0202 & 11 & 141 & 0.9391 & 0.0292 & 0.0009 & 2.317700e+00\\
\addlinespace
0203 & 12 & 101 & 0.8117 & 0.0775 & 0.0060 & 4.795200e+00\\
0204 & 11 & 224 & 1.0000 & 0.0000 & 0.0000 & NaN\\
0205 & 9 & 85 & 0.9646 & 0.0331 & 0.0011 & 2.932600e+00\\
0206 & 8 & 102 & 1.0000 & 0.0000 & 0.0000 & NaN\\
0207 & 11 & 149 & 1.0000 & 0.0000 & 0.0000 & NaN\\
\bottomrule
\end{tabular}
\end{table}

Domains with 5 or more UPMs and all those with a deff greater than 1 are now filtered.



```r
base_sae <- indicador_dam %>%
  filter(n_upm >= 5,
         pobreza_deff > 1) 
```


Next, the transformation $\log(\hat{\sigma}^2_d)$ is performed. Additionally, the selection of the municipality identifier columns (`dam2`), the direct estimation (`pobreza`), the number of people in the domain (`nd`), and the estimated variance for the direct estimation (`vardir`) is carried out, the latter being transformed using the `log()` function.


```r
baseFGV <-  base_sae %>% 
  select(dam2, pobreza, nd = n_pobreza, vardir = pobreza_var) %>%
  mutate(ln_sigma2 = log(vardir))
```

## Graphical Analysis

The first graph, `p1`, displays a scatter plot of the variable `ln_sigma2` against the variable `pobreza`, with a smooth line representing a trend estimation. The x-axis is labeled as _pobreza_.

The second graph, `p2`, exhibits a scatter plot of the variable `ln_sigma2` against the variable `nd`, with a smooth line indicating a trend estimation. The x-axis is labeled as _Tamaño de muestra_ (Sample Size).

The third graph, `p3`, demonstrates a scatter plot of the variable `ln_sigma2` in relation to the product of `pobreza` and `nd`, with a smooth line representing a trend estimation. The x-axis is labeled as _Número de pobres_ (Number of Poor).

The fourth graph, `p4`, shows a scatter plot of the variable `ln_sigma2` against the square root of the variable `pobreza`, with a smooth line representing a trend estimation. The x-axis is labeled as _Raiz cuadrada de pobreza_ (Square Root of Poverty).

Overall, these graphs are designed to explore the relationship between `ln_sigma2` and different independent variables such as `pobreza`, `nd`, and the square root of poverty. Choosing to use the "loess" function to smooth the lines instead of a straight line aids in visualizing general trends in the data more effectively.



```r
theme_set(theme_bw())

# pobreza vs Ln_sigma2 #
p1 <- ggplot(baseFGV, aes(x = pobreza, y = ln_sigma2)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("poverty")

# Sample Size vs Ln_sigma2 #
p2 <- ggplot(baseFGV, aes(x = nd, y = ln_sigma2)) + 
  geom_point() +
  geom_smooth(method = "loess") + 
  xlab("Sample size")

# Number of Poor vs Ln_sigma2 #
p3 <- ggplot(baseFGV, aes(x = pobreza * nd, y = ln_sigma2)) + 
  geom_point() +
  geom_smooth(method = "loess") + 
  xlab("Number of poor")

# Square Root of Poverty vs Ln_sigma2 #
p4 <- ggplot(baseFGV, aes(x = sqrt(pobreza), y = ln_sigma2)) + 
  geom_point() +
  geom_smooth(method = "loess") + 
  xlab("Square root of poverty")

library(patchwork)
(p1 | p2) / (p3 | p4)
```

![](02_S1_FGV_files/figure-latex/unnamed-chunk-9-1.pdf)<!-- --> 


## Variance Model

The code fits a multiple linear regression model (using the `lm()` function), where `ln_sigma2` is the response variable and the predictor variables include `pobreza`, `nd`, and various transformations of these variables. The goal of this model is to estimate the generalized variance function (FGV) for the observed domains.


```r
library(gtsummary)
FGV1 <- lm(ln_sigma2 ~ -1 +  pobreza +
             I(pobreza*nd),
     data = baseFGV)

tbl_regression(FGV1) %>% 
  add_glance_table(include = c(r.squared, adj.r.squared))
```


\begin{tabular}{l|c|c|c}
\hline
**Characteristic** & **Beta** & **95\% CI** & **p-value**\\
\hline
pobreza & -12 & -20, -4.5 & 0.003\\
\hline
I(pobreza * nd) & 0.02 & -0.04, 0.07 & 0.5\\
\hline
R² & 0.345 &  & \\
\hline
Adjusted R² & 0.311 &  & \\
\hline
\end{tabular}

After obtaining the model estimation, the value of the constant $\Delta$ must be obtained, for which the following code is used.


```r
delta.hat = sum(baseFGV$vardir) / 
  sum(exp(fitted.values(FGV1)))
```



From which it is derived that $\Delta = 0.110434$. Finally, it is possible to obtain the smoothed variance by executing the following command.


```r
hat.sigma <- 
  data.frame(dam2 = baseFGV$dam2,
             hat_var = delta.hat * exp(fitted.values(FGV1)))

baseFGV <- left_join(baseFGV, hat.sigma)
tba(head(baseFGV, 10))
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrr}
\toprule
dam2 & pobreza & nd & vardir & ln\_sigma2 & hat\_var\\
\midrule
0102 & 0.9836 & 197 & 0.0001 & -9.4372 & 0e+00\\
0103 & 1.0000 & 65 & 0.0000 & -76.3620 & 0e+00\\
0202 & 0.9391 & 141 & 0.0009 & -7.0701 & 0e+00\\
0203 & 0.8117 & 101 & 0.0060 & -5.1159 & 0e+00\\
0205 & 0.9646 & 85 & 0.0011 & -6.8149 & 0e+00\\
\addlinespace
0212 & 0.8304 & 59 & 0.0101 & -4.5936 & 0e+00\\
0301 & 0.9419 & 45 & 0.0013 & -6.6569 & 0e+00\\
0302 & 0.9109 & 159 & 0.0017 & -6.3681 & 0e+00\\
0401 & 0.8063 & 319 & 0.0006 & -7.4369 & 4e-04\\
0402 & 0.8311 & 259 & 0.0012 & -6.7172 & 1e-04\\
\bottomrule
\end{tabular}
\end{table}

Model validation for the FGV


```r
par(mfrow = c(2, 2))
plot(FGV1)
```

![](02_S1_FGV_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 

Smoothed variance prediction


```r
base_sae <- left_join(indicador_dam,
                      baseFGV %>% select(id_dominio, hat_var),
                      by = id_dominio) %>%
  mutate(
    pobreza_var = ifelse(is.na(hat_var), NA_real_, pobreza_var),
    pobreza_deff = ifelse(is.na(hat_var), NA_real_, pobreza_deff)
  )
```

Now, we make a line graph to see the volatility and the estimates of the variances.


```{.r .fold-hide}
nDom <- sum(!is.na(base_sae$hat_var))
temp_FH <- base_sae %>% filter(!is.na(hat_var))
ggplot(temp_FH %>%
         arrange(n_pobreza), aes(x = 1:nDom)) +
  geom_line(aes(y = pobreza_var, color = "VarDirEst")) +
  geom_line(aes(y = hat_var, color = "FGV")) +
  labs(y = "Varianzas", x = "Sample size", color = " ") +
  scale_x_continuous(breaks = seq(1, nDom, by = 10),
labels = temp_FH$n_pobreza[order(temp_FH$n_pobreza)][seq(1, nDom, by = 10)]) +
  scale_color_manual(values = c("FGV" = "Blue", "VarDirEst" = "Red"))
```

![](02_S1_FGV_files/figure-latex/unnamed-chunk-15-1.pdf)<!-- --> 


This code performs several transformations on the dataset `base_sae`:

1. **Creation of new variables:**
   - `pobreza_deff`: Replaces NaN values with 1 if they exist; otherwise, it keeps the original value.
   - `deff_FGV`: Computes a new Design Effect (DEFF) using the formula `hat_var / (pobreza_var / pobreza_deff)` when `pobreza_var` is not equal to 0.
   - `n_eff_FGV`: Calculates the effective number of surveyed individuals as `n_pobreza / deff_FGV`.

2. **Modification of the variable `pobreza`:**
   - If `hat_var` is NA, it replaces `pobreza` values with NA; otherwise, it retains the original value.




```r
base_FH <- base_sae %>%
  mutate(
    pobreza_deff = ifelse(is.nan(pobreza_deff), 1, pobreza_deff),
    deff_FGV = ifelse(pobreza_var == 0 ,
      1,
      hat_var / (pobreza_var / pobreza_deff) #Fórmula del nuevo DEFF
    ),
    # Criterio MDS para regularizar el DeffFGV
    n_eff_FGV = n_pobreza / deff_FGV, #Número efectivo de personas encuestadas

     pobreza = ifelse(is.na(hat_var), NA_real_, pobreza) 
  )


#saveRDS(object = base_FH, "Recursos/02_FGV/base_FH_2020.rds")
```


