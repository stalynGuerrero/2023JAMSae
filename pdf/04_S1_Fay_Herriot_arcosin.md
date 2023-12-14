## Area models - ArcSin transformation.



In its most basic conception, the **Fay-Herriot** model is a linear combination of covariates. However, the result of this combination can take values that fall outside the acceptable range for a proportion; that is, generally, the Fay-Herriot estimator $\theta \in \mathbb{R}$, whereas the direct estimator $\theta \in (0,1)$. The arcsine transformation is given by:

$$
\hat{z}_d = \arcsin\left( \sqrt{ \hat{\theta}_d} \right)
$$ 

where 

$$
Var\left( \hat{z}_d \right) = \frac{\widehat{DEFF}_d}{4\times n_d} = \frac{1}{4\times n_{d, \text{effective}} }
$$

The Fay-Herriot model is defined as follows:


\begin{align*}
Z_d \mid \mu_d,\sigma^2_d &  \sim  N(\mu_d, \sigma^2_d)\\
\mu_d & = \boldsymbol{x}^{T}_{d}\boldsymbol{\beta} + u_d \\
\theta_d & =  \left(\sin(\mu_d)\right)^2
\end{align*}


where $u_d \sim N(0 , \sigma^2)$.

Let the prior distributions for $\boldsymbol{\beta}$ and $\sigma_{u}^{2}$ be given by:


\begin{align*}
\boldsymbol{\beta}	\sim	N\left(0,1000 \right)\\
\sigma_{u}^{2}	\sim	\text{IG}\left(0.0001,0.0001\right)
\end{align*}


### Estimation procedure

Reading the database that resulted in the previous step and selecting the columns of interest


```r
library(tidyverse)
library(magrittr)

base_FH <- readRDS('Recursos/04_FH_Arcosin/01_base_FH.Rds') %>% 
transmute(dam2,                            ## id dominios
          pobreza,
          T_pobreza = asin(sqrt(pobreza)),  ## creando zd
          n_effec = n_eff_FGV,              ## n efectivo
          varhat = 1/(4*n_effec)            ## varianza para zd
)
```

Joining the two databases.


```r
statelevel_predictors_df <-
  readRDS('Recursos/03_FH_normal/02_statelevel_predictors_dam.rds') %>%
  mutate(id_order = 1:n())
base_FH <-
  full_join(base_FH, statelevel_predictors_df, by = "dam2")
tba(base_FH[, 1:8] %>% head(10))
```

\begin{table}[H]
\centering
\centering
\begin{tabular}[t]{lrrrrrrr}
\toprule
dam2 & pobreza & T\_pobreza & n\_effec & varhat & area1 & sex2 & age2\\
\midrule
0101 & NA & NA & NA & NA & 1.0000 & 0.5087 & 0.2694\\
0102 & 0.9836 & 1.4422 & 1029.570 & 2.000000e-04 & 1.0000 & 0.4754 & 0.2857\\
0103 & 1.0000 & 1.5708 & 0.000 & 2.226882e+57 & 1.0000 & 0.5037 & 0.3095\\
0201 & NA & NA & NA & NA & 0.5147 & 0.5060 & 0.2962\\
0202 & 0.9391 & 1.3215 & 5881.354 & 0.000000e+00 & 0.9986 & 0.5376 & 0.2625\\
\addlinespace
0203 & 0.8117 & 1.1219 & 6965.029 & 0.000000e+00 & 0.9754 & 0.5432 & 0.2454\\
0204 & NA & NA & NA & NA & 1.0000 & 0.5300 & 0.3151\\
0205 & 0.9646 & 1.3814 & 11790.358 & 0.000000e+00 & 1.0000 & 0.5182 & 0.3057\\
0206 & NA & NA & NA & NA & 1.0000 & 0.5157 & 0.3192\\
0207 & NA & NA & NA & NA & 1.0000 & 0.5097 & 0.3099\\
\bottomrule
\end{tabular}
\end{table}

Selecting the covariates for the model.


```r
names_cov <-
    c(
       "ODDJOB","WORKED",
      "stable_lights_mean",
      "accessibility_mean",
      "urban.coverfraction_sum"
    )
```

### Preparing Inputs for `STAN`

1. Splitting the database into observed and unobserved domains

Observed domains.

```r
data_dir <- base_FH %>% filter(!is.na(T_pobreza))
```

Unobserved domains.

```r
data_syn <-
  base_FH %>% anti_join(data_dir %>% select(dam2))
```

2. Defining the fixed-effects matrix.


```r
## Observed domains
Xdat <- cbind(inter = 1,data_dir[,names_cov])

## Unobserved domains
Xs <-  cbind(inter = 1,data_syn[,names_cov])
```

3. Creating a parameter list for `STAN`.


```r
sample_data <- list(
  N1 = nrow(Xdat),       # Observed.
  N2 = nrow(Xs),         # Unobserved.
  p  = ncol(Xdat),       # Number of regressors.
  X  = as.matrix(Xdat),  # Observed Covariates.
  Xs = as.matrix(Xs),    # Unobserved Covariates
  y  = as.numeric(data_dir$T_pobreza),
  sigma_e = sqrt(data_dir$varhat)
)
```

4. Compiling the model in `STAN`.
  

```r
library(rstan)
fit_FH_arcoseno <- "Recursos/04_FH_Arcosin/modelosStan/15FH_arcsin_normal.stan"
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE) # speed up running time 
model_FH_arcoseno <- stan(
  file = fit_FH_arcoseno,  
  data = sample_data,   
  verbose = FALSE,
  warmup = 2500,         
  iter = 3000,            
  cores = 4              
)
saveRDS(model_FH_arcoseno,
        "Recursos/04_FH_Arcosin/02_model_FH_arcoseno.rds")
```


```r
model_FH_arcoseno <- readRDS("Recursos/04_FH_Arcosin/02_model_FH_arcoseno.rds")
```


#### Model results for the observed domains. 



```r
library(bayesplot)
library(patchwork)
library(posterior)

y_pred_B <- as.array(model_FH_arcoseno, pars = "theta") %>% 
  as_draws_matrix()
rowsrandom <- sample(nrow(y_pred_B), 100)

y_pred2 <- y_pred_B[rowsrandom, ]
ppc_dens_overlay(y = as.numeric(data_dir$pobreza), y_pred2)
```
![ppc arcosin](Recursos/04_FH_Arcosin/03_ppc_Arcosin.png)

Graphical analysis of the convergence of $\sigma^2_u$ chains.


```r
posterior_sigma2_u <- as.array(model_FH_arcoseno, pars = "sigma2_u")
(mcmc_dens_chains(posterior_sigma2_u) +
    mcmc_areas(posterior_sigma2_u) ) / 
  mcmc_trace(posterior_sigma2_u)
```

![Posterior sigma2](Recursos/04_FH_Arcosin/04_sigma2.png)



To validate the convergence of all chains, the *R-hat* is used.


```r
parametros <- summary(model_FH_arcoseno, 
                      pars =  c("theta", "theta_pred")
                      )$summary %>%   data.frame()
p1 <- mcmc_rhat(parametros$Rhat)
p1
```

![Rhat](Recursos/04_FH_Arcosin/05_Rhat.png)


Estimation of the FH of poverty in the observed domains.


```r
theta_FH <-   summary(model_FH_arcoseno,pars =  "theta")$summary %>%
  data.frame()
data_dir %<>% mutate(pred_arcoseno = theta_FH$mean, 
                     pred_arcoseno_EE = theta_FH$sd,
                     Cv_pred = pred_arcoseno_EE/pred_arcoseno)
```

Estimation of the FH of poverty in the NOT observed domains.


```r
theta_FH_pred <- summary(model_FH_arcoseno,pars =  "theta_pred")$summary %>%
  data.frame()
data_syn <- data_syn %>% 
  mutate(pred_arcoseno = theta_FH_pred$mean,
         pred_arcoseno_EE = theta_FH_pred$sd,
         Cv_pred = pred_arcoseno_EE/pred_arcoseno)
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-15)Estimation}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
dam2 & pobreza & pred\_arcoseno & pred\_arcoseno\_EE & Cv\_pred\\
\midrule
0101 & NA & 0.8630 & 0.1463 & 0.1696\\
0201 & NA & 0.8797 & 0.1426 & 0.1621\\
0204 & NA & 0.8412 & 0.1592 & 0.1893\\
0206 & NA & 0.8980 & 0.1296 & 0.1443\\
0207 & NA & 0.8958 & 0.1321 & 0.1475\\
\addlinespace
0208 & NA & 0.7818 & 0.1893 & 0.2421\\
0209 & NA & 0.8910 & 0.1313 & 0.1474\\
0210 & NA & 0.8990 & 0.1290 & 0.1435\\
0211 & NA & 0.8872 & 0.1389 & 0.1566\\
0502 & NA & 0.7649 & 0.1889 & 0.2470\\
\bottomrule
\end{tabular}
\end{table}



consolidating the bases of estimates for observed and UNobserved domains.


```r
estimacionesPre <- bind_rows(data_dir, data_syn) %>% 
  select(dam2, theta_pred = pred_arcoseno) %>% 
  mutate(dam = substr(dam2,1,2))
```


## Benchmark Process

1. From the census extract the total number of people by DAM2



```r
total_pp <- readRDS(file = "Recursos/04_FH_Arcosin/06_censo_mrp.rds") %>% 
   mutate(dam = substr(dam2,1,2))


N_dam_pp <- total_pp %>%   group_by(dam,dam2) %>%  
            summarise(total_pp = sum(n) ) %>% 
  group_by(dam) %>% mutate(dam_pp = sum(total_pp))

tba(N_dam_pp %>% data.frame() %>% slice(1:10),
    cap = "Number of people by DAM2")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-17)Number of people by DAM2}
\centering
\begin{tabular}[t]{llrr}
\toprule
dam & dam2 & total\_pp & dam\_pp\\
\midrule
01 & 0101 & 35750 & 88799\\
01 & 0102 & 26458 & 88799\\
01 & 0103 & 26591 & 88799\\
02 & 0201 & 58937 & 573065\\
02 & 0202 & 43700 & 573065\\
\addlinespace
02 & 0203 & 35309 & 573065\\
02 & 0204 & 47230 & 573065\\
02 & 0205 & 35694 & 573065\\
02 & 0206 & 38832 & 573065\\
02 & 0207 & 42103 & 573065\\
\bottomrule
\end{tabular}
\end{table}


2. Obtaining direct estimates by DAM or the level of aggregation at which the survey is representative.

In this code, an RDS file of a survey (`07_data_JAM.rds`) is read, and the `transmute()` function is used to select and transform the variables of interest.


```r
encuesta <- readRDS("Recursos/04_FH_Arcosin/07_encuesta.rds") %>% 
   mutate(dam = substr(dam2,1,2))
```

The code is conducting survey data analysis using the `survey` package in R. Initially, an object `design` is created as a survey design using the `as_survey_design()` function from the `srvyr` package. This design includes primary sampling unit identifiers (`upm`), weights (`fep`), strata (`estrato`), and survey data (`encuesta`). Subsequently, the `design` object is grouped by the variable "Aggregate," and the mean of the variable "pobreza" with a confidence interval for the entire population is calculated using the `survey_mean()` function. The result is stored in the `directoDam` object and displayed in a table.



```r
library(survey)
library(srvyr)
options(survey.lonely.psu = "adjust")

diseno <-
  as_survey_design(
    ids = upm,
    weights = fep,
    strata = estrato,
    nest = TRUE,
    .data = encuesta
  )
directoDam <- diseno %>% 
    group_by(dam) %>% 
  summarise(
    theta_dir = survey_mean(pobreza, vartype = c("ci"))
    )
tba(directoDam %>% slice(1:10), cap = "Direct estimation")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-19)Direct estimation}
\centering
\begin{tabular}[t]{lrrr}
\toprule
dam & theta\_dir & theta\_dir\_low & theta\_dir\_upp\\
\midrule
01 & 0.9926 & 0.9848 & 1.0004\\
02 & 0.9760 & 0.9638 & 0.9881\\
03 & 0.9191 & 0.8566 & 0.9816\\
04 & 0.8160 & 0.7742 & 0.8578\\
05 & 0.7486 & 0.6943 & 0.8029\\
\addlinespace
06 & 0.9419 & 0.8993 & 0.9844\\
07 & 0.8019 & 0.7076 & 0.8962\\
08 & 0.5933 & 0.4904 & 0.6961\\
09 & 0.8628 & 0.8167 & 0.9090\\
10 & 0.7379 & 0.6348 & 0.8410\\
\bottomrule
\end{tabular}
\end{table}


3. Carry out the consolidation of information obtained in *1* and *2*.


```r
temp <- estimacionesPre %>%
  inner_join(N_dam_pp ) %>% 
  inner_join(directoDam )

tba(temp %>% slice(1:10), cap = "Join datas")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-20)Join datas}
\centering
\begin{tabular}[t]{lrlrrrrr}
\toprule
dam2 & theta\_pred & dam & total\_pp & dam\_pp & theta\_dir & theta\_dir\_low & theta\_dir\_upp\\
\midrule
0102 & 0.9832 & 01 & 26458 & 88799 & 0.9926 & 0.9848 & 1.0004\\
0103 & 0.8759 & 01 & 26591 & 88799 & 0.9926 & 0.9848 & 1.0004\\
0202 & 0.9391 & 02 & 43700 & 573065 & 0.9760 & 0.9638 & 0.9881\\
0203 & 0.8117 & 02 & 35309 & 573065 & 0.9760 & 0.9638 & 0.9881\\
0205 & 0.9645 & 02 & 35694 & 573065 & 0.9760 & 0.9638 & 0.9881\\
\addlinespace
0212 & 0.8303 & 02 & 58509 & 573065 & 0.9760 & 0.9638 & 0.9881\\
0301 & 0.9419 & 03 & 41762 & 93896 & 0.9191 & 0.8566 & 0.9816\\
0302 & 0.9109 & 03 & 52134 & 93896 & 0.9191 & 0.8566 & 0.9816\\
0401 & 0.8079 & 04 & 49909 & 81732 & 0.8160 & 0.7742 & 0.8578\\
0402 & 0.8307 & 04 & 31823 & 81732 & 0.8160 & 0.7742 & 0.8578\\
\bottomrule
\end{tabular}
\end{table}

4. With the organized information, calculate the weights for the Benchmark



```r
R_dam2 <- temp %>% group_by(dam) %>% 
  summarise(
  R_dam_RB = unique(theta_dir) / sum((total_pp  / dam_pp) * theta_pred)
) %>%
  left_join(directoDam, by = "dam")

tba(R_dam2 %>% arrange(desc(R_dam_RB)), cap = "Weights for the Benchmark")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-21)Weights for the Benchmark}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
dam & R\_dam\_RB & theta\_dir & theta\_dir\_low & theta\_dir\_upp\\
\midrule
13 & 1.1793 & 0.6461 & 0.5226 & 0.7696\\
02 & 1.1158 & 0.9760 & 0.9638 & 0.9881\\
01 & 1.0996 & 0.9926 & 0.9848 & 1.0004\\
09 & 1.0669 & 0.8628 & 0.8167 & 0.9090\\
14 & 1.0550 & 0.8792 & 0.8376 & 0.9207\\
\addlinespace
10 & 1.0465 & 0.7379 & 0.6348 & 0.8410\\
05 & 1.0098 & 0.7486 & 0.6943 & 0.8029\\
06 & 1.0031 & 0.9419 & 0.8993 & 0.9844\\
04 & 0.9991 & 0.8160 & 0.7742 & 0.8578\\
08 & 0.9973 & 0.5933 & 0.4904 & 0.6961\\
\addlinespace
03 & 0.9940 & 0.9191 & 0.8566 & 0.9816\\
07 & 0.9474 & 0.8019 & 0.7076 & 0.8962\\
11 & 0.8631 & 0.4277 & 0.3056 & 0.5498\\
12 & 0.0480 & 0.0292 & 0.0120 & 0.0463\\
\bottomrule
\end{tabular}
\end{table}
calculating the weights for each domain.


```r
pesos <- temp %>% 
  mutate(W_i = total_pp / dam_pp) %>% 
  select(dam2, W_i)
tba(pesos %>% slice(1:10), cap = "Weights")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-22)Weights}
\centering
\begin{tabular}[t]{lr}
\toprule
dam2 & W\_i\\
\midrule
0102 & 0.2980\\
0103 & 0.2995\\
0202 & 0.0763\\
0203 & 0.0616\\
0205 & 0.0623\\
\addlinespace
0212 & 0.1021\\
0301 & 0.4448\\
0302 & 0.5552\\
0401 & 0.6106\\
0402 & 0.3894\\
\bottomrule
\end{tabular}
\end{table}


5. Perform FH Benchmark Estimation


```r
estimacionesBench <- estimacionesPre %>%
  left_join(R_dam2, by = c("dam")) %>%
  mutate(theta_pred_RBench = R_dam_RB * theta_pred) %>%
  left_join(pesos) %>% 
  select(dam, dam2, W_i, theta_pred, theta_pred_RBench)  

  tba(estimacionesBench %>% slice(1:10), cap = "Estimation Benchmark")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-23)Estimation Benchmark}
\centering
\begin{tabular}[t]{llrrr}
\toprule
dam & dam2 & W\_i & theta\_pred & theta\_pred\_RBench\\
\midrule
01 & 0102 & 0.2980 & 0.9832 & 1.0811\\
01 & 0103 & 0.2995 & 0.8759 & 0.9632\\
02 & 0202 & 0.0763 & 0.9391 & 1.0478\\
02 & 0203 & 0.0616 & 0.8117 & 0.9057\\
02 & 0205 & 0.0623 & 0.9645 & 1.0762\\
\addlinespace
02 & 0212 & 0.1021 & 0.8303 & 0.9264\\
03 & 0301 & 0.4448 & 0.9419 & 0.9363\\
03 & 0302 & 0.5552 & 0.9109 & 0.9054\\
04 & 0401 & 0.6106 & 0.8079 & 0.8072\\
04 & 0402 & 0.3894 & 0.8307 & 0.8299\\
\bottomrule
\end{tabular}
\end{table}

6. Validation: FH Estimation with Benchmark


```r
estimacionesBench %>% group_by(dam) %>%
  summarise(theta_reg_RB = sum(W_i * theta_pred_RBench)) %>%
  left_join(directoDam, by = "dam") %>% 
  tba(cap = "FH Estimation with Benchmark")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-24)FH Estimation with Benchmark}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
dam & theta\_reg\_RB & theta\_dir & theta\_dir\_low & theta\_dir\_upp\\
\midrule
01 & 0.9926 & 0.9926 & 0.9848 & 1.0004\\
02 & 0.9760 & 0.9760 & 0.9638 & 0.9881\\
03 & 0.9191 & 0.9191 & 0.8566 & 0.9816\\
04 & 0.8160 & 0.8160 & 0.7742 & 0.8578\\
05 & 0.7486 & 0.7486 & 0.6943 & 0.8029\\
\addlinespace
06 & 0.9419 & 0.9419 & 0.8993 & 0.9844\\
07 & 0.8019 & 0.8019 & 0.7076 & 0.8962\\
08 & 0.5933 & 0.5933 & 0.4904 & 0.6961\\
09 & 0.8628 & 0.8628 & 0.8167 & 0.9090\\
10 & 0.7379 & 0.7379 & 0.6348 & 0.8410\\
\addlinespace
11 & 0.4277 & 0.4277 & 0.3056 & 0.5498\\
12 & 0.0292 & 0.0292 & 0.0120 & 0.0463\\
13 & 0.6461 & 0.6461 & 0.5226 & 0.7696\\
14 & 0.8792 & 0.8792 & 0.8376 & 0.9207\\
\bottomrule
\end{tabular}
\end{table}

### Results Validation

This code conducts data analysis and visualization using the `ggplot2` library. Specifically, it merges two `data frames` using the `left_join()` function, groups the data by the `dam` variable, and performs some operations to transform the `thetaFH` and `theta_pred_RBench` variables. Afterwards, it utilizes the `gather()` function to organize the data in a long format and visualizes it with `ggplot()`.

The resulting visualization displays points in different shapes and colors, representing various estimation methods. Additionally, it includes two dashed lines that depict the upper and lower confidence intervals for the observed values in the `theta_dir` variable.



```r
temp <- estimacionesBench %>% left_join(
bind_rows(
data_dir %>% select(dam2, thetaFH = pred_arcoseno),
data_syn %>% select(dam2, thetaFH = pred_arcoseno))) %>% 
group_by(dam) %>% 
summarise(thetaFH = sum(W_i * theta_pred),
          theta_RBench = sum(W_i * theta_pred_RBench)
          ) %>%   
left_join(directoDam, by = "dam")  %>% 
mutate(id = 1:n())

temp %<>% gather(key = "Metodo",value = "Estimacion",
                -id, -dam, -theta_dir_upp, -theta_dir_low)

p1 <- ggplot(data = temp, aes(x = id, y = Estimacion, shape = Metodo)) +
  geom_point(aes(color = Metodo), size = 2) +
  geom_line(aes(y = theta_dir_low), linetype  = 2) +
  geom_line(aes(y = theta_dir_upp),  linetype  = 2) +
  theme_bw(20) + 
  scale_x_continuous(breaks = temp$id,
    labels =  temp$dam) +
  labs(y = "", x = "")

# ggsave(plot = p1,
#        filename = "Recursos/04_FH_Arcosin/08_validar_bench.png",width = 16,height = 12)
p1 
```

![](Recursos/04_FH_Arcosin/08_validar_bench.png)



## Poverty Map

This code block loads various packages (`sf`, `tmap`) and performs several operations. Initially, it conducts a `left_join` between the benchmark-adjusted estimates (`estimacionesBench`) and the model estimates (`data_dir`, `data_syn`), utilizing the `dam2` variable as the key for the join. Subsequently, it reads a `Shapefile` containing geospatial information for the country. Then, it creates a thematic map (`tmap`) using the `tm_shape()` function and adds layers using the `tm_polygons()` function. The map represents a variable `theta_pred_RBench` utilizing a color palette named "YlOrRd" and sets the intervals' breaks for the variable with the variable `brks_lp`. Finally, the `tm_layout()` function sets some design parameters for the map, such as the aspect ratio (asp).



```r
library(sf)
library(tmap)

estimacionesBench %<>% left_join(
bind_rows(
data_dir %>% select(dam2, pobreza, pred_arcoseno_EE , Cv_pred),
data_syn %>% select(dam2,pobreza, pred_arcoseno_EE , Cv_pred)))

## Leer Shapefile del país
ShapeSAE <- read_sf("Shapefile/JAM2_cons.shp")


mapa <- tm_shape(ShapeSAE %>%
                   left_join(estimacionesBench,  by = "dam2"))

tmap_options(check.and.fix = TRUE)
Mapa_lp <-
  mapa + tm_polygons(
    c("pobreza","theta_pred_RBench"),
    title = "Poverty map",
    palette = "YlOrRd",
    colorNA = "white"
  ) + tm_layout(asp = 1.5)

tmap_save(Mapa_lp, 
          filename = "Recursos/04_FH_Arcosin/09_map.png",
           width = 2500,
  height = 2000,
  asp = 0)
Mapa_lp
```

![](Recursos/04_FH_Arcosin/09_map.png)



