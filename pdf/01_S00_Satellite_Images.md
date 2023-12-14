# Session 1- Census and satellite information



## Use of Satellite Imagery and SAE

One of the pioneering articles in small area estimation was the paper by Singh, R, et al. (2002), which addressed crop yield estimation for the tehsils (sub-administrative units) of Rohtak district in Haryana, India.

Raster images represent the world through a set of contiguous equally spaced cells known as pixels. These images contain information like a geographic information system and a coordinate reference system. Images store an identifier, a value in each pixel (or a vector with different values), and each cell is associated with a color scale.

Images can be obtained in raw and processed forms. The former contains only color layers, while the latter also contains values that have been processed in each cell (vegetation indices, light intensity, type of vegetation).

Raw information can be used to train desired features (roads, crop types, forest/non-forest). Fortunately, in Google Earth Engine, we find many processed indicators associated with a pixel. These indicators can be aggregated at a geographical area level.



## Satellite Image Data Sources

Some of the main sources of satellite images include:

  * [USGS Earth Explorer](http://earthexplorer.usgs.gov/)
  
  * [Land Processes Distributed Active Archive Center (LP DAAC)](https://lpdaacsvc.cr.usgs.gov/appeears/)
  
  * [NASA Earthdata Search](https://search.earthdata.nasa.gov/search)
  
  * [Copernicus Open Access Hub](https://scihub.copernicus.eu/)
  
  * [AWS Public Dataset - Landsat](https://aws.amazon.com/public-data-sets/landsat/)

However, most of these sources are centralized within **Google Earth Engine**, which allows searching for satellite image data sources. GEE can be managed through APIs in different programming languages: JavaScript (by default), Python, and R (rgee package).


## Google Earth Engine

Create an account at [this link](https://earthengine.google.com/). Once logged in, you can search for datasets of interest:

![Night Lights Image](Recursos/01_Session1/01_lights.png)

* Upon searching for the dataset, you can open a code editor provided by Google in JavaScript.

* Copy and paste the syntax provided by the dataset search to visualize the raster image and obtain statements allowing for the retrieval of the dataset of interest later in R.

![Syntax in JavaScript](Recursos/01_Session1/02_query.png)

## Installing rgee

* Download and install Anaconda or Conda from [here](https://www.anaconda.com/products/individual).

* Open **Anaconda Prompt** and set up a working environment (Python environment JAM2023) using the following commands:

```
conda env list
conda create -n JAM2023 python=3.9
activate JAM2023
pip install google-api-python-client
pip install earthengine-api
pip install numpy
```

* List available Python environments in Anaconda Prompt:

```
conda env list
```

* Once you've identified the path of the JAM2023 environment, set it in R (**remember to change \\ to /**).
* Install `reticulate` and `rgee`, load packages for spatial processing, and set up the working environment as follows:



```r
library(reticulate) # Connection with Python
library(rgee) # Connection with Google Earth Engine
library(sf) # Package for handling geographic data
library(dplyr) # Package for data processing
library(magrittr)

rgee_environment_dir = "C:/Users/gnieto/Anaconda3/envs/JAM2023/python.exe"

# Set up Python (Sometimes not detected and R needs to be restarted)
reticulate::use_python(rgee_environment_dir, required=T)

rgee::ee_install_set_pyenv(py_path = rgee_environment_dir, py_env = "JAM2023")

Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)
```

* Once the environment is configured, you can initialize a Google Earth Engine session as follows:



```r
rgee::ee_Initialize(drive = T)
```

![Session started successfully](Recursos/01_Session1/03_Figura.PNG)


**Notes:**

- Each session must be initialized with the command `rgee::ee_Initialize(drive = T)`.

- JavaScript commands invoking methods with "." are replaced by the dollar sign ($), for example:


```r
ee.ImageCollection().filterDate()  # JavaScript
ee$ImageCollection()$filterDate()  # R
```


### Downloading Satellite Information

* **Step 1**: Have the shapefiles ready.


```r
shape <- read_sf("Shapefile/JAM2_cons.shp") 
plot(shape["geometry"])
```

![Shapefile](Recursos/01_Session1/04_JAM.PNG)

* **Step 2**: Select the image file you want to process, for example, **night lights**.


```r
lights <- ee$ImageCollection("NOAA/DMSP-OLS/NIGHTTIME_LIGHTS") %>%
  ee$ImageCollection$filterDate("2013-01-01", "2014-01-01") %>%
  ee$ImageCollection$map(function(x) x$select("stable_lights")) %>%
  ee$ImageCollection$toBands()
```

* **Step 3**: Download the information.


```r
## Takes about 10 minutes 
lights_shape <- map(unique(shape$dam2),
                 ~tryCatch(ee_extract(
                   x = lights,
                   y = shape["dam2"] %>% filter(dam2 == .x),
                   ee$Reducer$mean(),
                   sf = FALSE
                 ) %>% mutate(dam2 = .x),
                 error = function(e)data.frame(dam2 = .x)))

lights_shape %<>% bind_rows()

tba(lights_shape, cap = "Average of night lights")
```

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-7)Average standardized night lights}
\centering
\begin{tabu} to \linewidth {>{\raggedright}X>{\raggedleft}X}
\toprule
dam2 & stable\_lights\_mean\\
\midrule
0101 & 0.9393\\
0102 & 1.6857\\
0103 & 1.6900\\
0201 & -0.4380\\
0202 & 1.4627\\
\addlinespace
0203 & 1.3519\\
0204 & 1.6333\\
0205 & 1.7522\\
0206 & 1.7522\\
0207 & 1.7444\\
\bottomrule
\end{tabu}
\end{table}

Repeat the routine for:

- Soil type: **crops-coverfraction** (Percentage of crop cover) and **urban-coverfraction** (Percentage of urban cover) available at <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global#description>

- Travel time to the nearest hospital or clinic (**accessibility**) and travel time to the nearest hospital or clinic using non-motorized transport (**accessibility_walking_only**) information available at <https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_accessibility_to_healthcare_2019>

- Human modification, considering human settlements, agriculture, transportation, mining, energy production, and electrical infrastructure. You can find satellite information at the following link: <https://developers.google.com/earth-engine/datasets/catalog/CSP_HM_GlobalHumanModification#description>


* **Paso 4**: Consolidate the information.

\begin{table}[H]

\caption{(\#tab:unnamed-chunk-8)Standardized satellite predictors}
\centering
\begin{tabu} to \linewidth {>{\raggedright}X>{\raggedleft}X>{\raggedleft}X}
\toprule
dam2 & stable\_lights\_mean & crops.coverfraction\_mean\\
\midrule
0101 & 0.9393 & -0.5459\\
0102 & 1.6857 & -0.7090\\
0103 & 1.6900 & -0.3571\\
0201 & -0.4380 & -0.0874\\
0202 & 1.4627 & -0.6237\\
\addlinespace
0203 & 1.3519 & -0.6402\\
0204 & 1.6333 & -0.5050\\
0205 & 1.7522 & -0.6844\\
0206 & 1.7522 & -0.4289\\
0207 & 1.7444 & -0.4662\\
\bottomrule
\end{tabu}
\end{table}

![Night Lights Sum](Recursos/01_Session1/06_stable_lights_sum.png)

![Night Lights Satellite](Recursos/01_Session1/07_Luces_nocturnas.PNG)

![Crop Cover](Recursos/01_Session1/08_crops.coverfraction_sum.png)

![Crop Cover Satellite](Recursos/01_Session1/09_Suelo_cultivo.PNG)



![Urban Cover Sum](Recursos/01_Session1/10_urban.coverfraction_sum.png)

![Urban Cover Satellite](Recursos/01_Session1/11_Suelo_urbano.PNG)


![Human Modification Sum](Recursos/01_Session1/12_gHM_sum.png)

![Human Modification Satellite](Recursos/01_Session1/13_Modificación_humana.PNG)




![Average Travel Time to Hospital Sum](Recursos/01_Session1/14_accessibility_sum.png)

![Average Travel Time to Hospital Satellite](Recursos/01_Session1/15_Distancia_Hospitales.PNG)



![Average Travel Time to Hospital by Non-Motorized Vehicle Sum](Recursos/01_Session1/16_accessibility_walking_only_mean.png)

![Average Travel Time to Hospital by Non-Motorized Vehicle Satellite](Recursos/01_Session1/17_Distancia_Hospitales_caminando.PNG)


## Population and Housing Censuses

It's necessary to define the variables for the country you want to work with. As a first step, access to the country's census data is required. You can access it from the following link: <https://redatam.org/en/microdata>, where you'll find a *.zip* file with the microdata for the country. To read this dataset, you'll need to use the *redatam.open* function from the `redatam` library. This function directly depends on the census dictionary from REDATAM software, which is a file with a .dicx extension and should be located in the same folder as the data being read. This is how an object is created within R that merges the dictionary with the microdata from the census database. After performing a process in R using REDATAM syntax, we have the following table:


\centering
\begin{tabu} to \linewidth {>{\raggedright}X>{\raggedleft}X>{\raggedleft}X>{\raggedleft}X>{\raggedleft}X>{\raggedleft}X}
\toprule
dam2 & area1 & sex2 & age2 & age3 & age4\\
\midrule
0101 & 1.0000 & 0.5087 & 0.2694 & 0.2297 & 0.1689\\
0102 & 1.0000 & 0.4754 & 0.2857 & 0.2261 & 0.1527\\
0103 & 1.0000 & 0.5037 & 0.3095 & 0.2015 & 0.1312\\
0201 & 0.5147 & 0.5060 & 0.2962 & 0.2090 & 0.1844\\
0202 & 0.9986 & 0.5376 & 0.2625 & 0.2226 & 0.2238\\
\addlinespace
0203 & 0.9754 & 0.5432 & 0.2454 & 0.2254 & 0.2388\\
0204 & 1.0000 & 0.5300 & 0.3151 & 0.2022 & 0.2034\\
0205 & 1.0000 & 0.5182 & 0.3057 & 0.2286 & 0.1981\\
0206 & 1.0000 & 0.5157 & 0.3192 & 0.1959 & 0.1552\\
0207 & 1.0000 & 0.5097 & 0.3099 & 0.1966 & 0.1691\\
\bottomrule
\end{tabu}

### Mapas de las variables con información censal. 



\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/area1} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/sex2} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/tiene_sanitario} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/tiene_electricidad} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/tiene_acueducto} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/tiene_gas} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/eliminar_basura} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/tiene_internet} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/material_paredes} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/material_techo} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/TRANMODE_PRIVATE_CAR} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/ODDJOB} 
\includegraphics[width=4\linewidth]{Recursos/01_Session1/18_plot_Censo/WORKED} 

