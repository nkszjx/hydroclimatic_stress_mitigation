# hydroclimatic_geographic_stress_mitigation
These are the codes that were used to quantify hydroclimatic and geographic stress mitigation by urbanization for global urban green spaces

The data sources for this study are available at [A global dataset of urban vegetation structure and ecosystem functioning](https://figshare.com/articles/dataset/_b_A_global_dataset_of_urban_vegetation_structure_and_ecosystem_functioning_b_/33919303).

## Results
![](/fig.jpg)


## Package pre-requisites
The R codes running environment are required. 

```
library(parallel)
library(lmtest)
library(DescTools)
library(foreign)
library(Matrix)
library(lfe)  
library(magrittr)
library(margins)
library(naniar)
library(dplyr)
library(plotly)
library(zoo)
library(readxl)
library(mice)
library(rio)
library(orca)
library(DMwR2)
library(car)
library(AER)
library(lme4)
library(ggplot2)
library(brms)
library(mgcv)
library(future)
library(plm)
library(dplyr)
library(data.table)
library(fixest)
library(AER)
library(clubSandwich)
library(lmtest)
library(ggpubr)
library(marginaleffects)
library(webshot)
library(htmlwidgets)
library(scales)
library(scatterplot3d)
library(extrafont)
library(showtext)
library(future)
library(future.apply)
library(stargazer)
library(tidyr)

```
## The datasets utilized in this study are as follows:
[FVC and NPP datasets](https://www.glass.hku.hk/archive/). 
[City boundary dataset](https://data-starcloud.pcl.ac.cn/iearthdata/14). 
[Impervious surface dataset](http://irsip.whu.edu.cn/resv2/dataweb2.php). 
[TerraClimate dataset](https://www.climatologylab.org/terraclimate.html). 
[Digital elevation map](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM). 
[Soil dataset](https://www.fao.org/soils-portal/data-hub/). 
[GDP dataset](https://databank.worldbank.org/). 
[HDI dataset](https://hdr.undp.org/data-center/). 
[Population density dataset](https://hub.worldpop.org/project/categories?id=18). 
[Human settlement dataset](https://human-settlement.emergency.copernicus.eu). 

## Acknowledgement

Our ideas may partly come from the papers:

[Surrounding greenness is associated with lower risk and burden of low birth weight in Iran](https://www.nature.com/articles/s41467-023-43425-6#:~:text=By%20involving%20~4%20million%20Iranian,risks%20of%20LBW%20and%20TLBW)

[Air pollution lowers Chinese urbanites’ expressed happiness on social media](https://www.nature.com/articles/s41562-018-0521-2)




## Citation

Jianhua Guo†,*, Bowei Yu†,*, Jianghao Wang, Zhongchang Sun, Dong Liang, and Huadong Guo*. "Urban environments buffer vegetation structure but not ecosystem functioning under hydroclimatic stress." Submitted to Nature XXX, Sep 2026




