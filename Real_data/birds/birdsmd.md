BIRDS
================
EB
2025-10-20

    knitr::opts_chunk$set(echo = TRUE)

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
setwd("C:/Users/Asus/Desktop/FactorModels/birds/data")
da = read.csv("data.csv", stringsAsFactors=TRUE)
#source("sampler.R")
```

## Birds co-occurrence

This script reproduces the analysis of the “birds” dataset (Section
4.1). The dataset contains information on the co-occurrence patterns of
$p=50$ bird species in Finland; data were collected over nine years
(2006-2014) across $S=200$ locations.

In this data analysis exercise, we are interested in fitting a factor
model to analyze the bird species occurrence in a unsupervised manner,
including only the “group” information (observation site), as the APAFA
model foresees.

Then, we will investigate the relation between environmental covariates
available and latent factors found by the model a posteriori

\##Prepare the data

Y is the (nxp) matrix containing the counts (as presence/absence) of
species (p=50) in the S=200 location (group) over time (average group
dimension 4.5).

``` r
da = read.csv( ( "data/data.csv"), stringsAsFactors=TRUE)

colnames(da)[-c(1:9)]
```

    ##  [1] "Phylloscopus_trochilus"  "Turdus_iliacus"         
    ##  [3] "Fringilla_coelebs"       "Turdus_philomelos"      
    ##  [5] "Carduelis_spinus"        "Parus_major"            
    ##  [7] "Anthus_trivialis"        "Erithacus_rubecula"     
    ##  [9] "Cuculus_canorus"         "Muscicapa_striata"      
    ## [11] "Ficedula_hypoleuca"      "Turdus_pilaris"         
    ## [13] "Loxia_curvirostra"       "Columba_palumbus"       
    ## [15] "Sylvia_borin"            "Corvus_corone"          
    ## [17] "Dendrocopos_major"       "Emberiza_citrinella"    
    ## [19] "Turdus_merula"           "Sylvia_curruca"         
    ## [21] "Regulus_regulus"         "Prunella_modularis"     
    ## [23] "Phoenicurus_phoenicurus" "Parus_caeruleus"        
    ## [25] "Parus_montanus"          "Motacilla_alba"         
    ## [27] "Carduelis_chloris"       "Phylloscopus_collybita" 
    ## [29] "Tringa_ochropus"         "Parus_cristatus"        
    ## [31] "Sylvia_communis"         "Grus_grus"              
    ## [33] "Larus_canus"             "Pica_pica"              
    ## [35] "Phylloscopus_sibilatrix" "Turdus_viscivorus"      
    ## [37] "Hirundo_rustica"         "Sylvia_atricapilla"     
    ## [39] "Tetrao_tetrix"           "Pyrrhula_pyrrhula"      
    ## [41] "Dryocopus_martius"       "Certhia_familiaris"     
    ## [43] "Alauda_arvensis"         "Corvus_corax"           
    ## [45] "Saxicola_rubetra"        "Garrulus_glandarius"    
    ## [47] "Numenius_arquata"        "Carpodacus_erythrinus"  
    ## [49] "Gallinago_gallinago"     "Corvus_monedula"

``` r
Y = matrix((da[, -c(1:9)]) > 0, nrow = nrow(da))
Y = apply(Y, MARGIN = 2, FUN = as.numeric)
head(Y)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]    1    1    1    1    1    1    1    1    0     1     1     1     0     0
    ## [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
    ## [3,]    1    1    1    1    1    1    1    1    1     1     1     0     0     1
    ## [4,]    1    0    1    1    1    1    1    1    0     1     1     0     1     1
    ## [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1     0
    ## [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
    ##      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
    ## [1,]     1     1     1     1     1     1     0     0     1     1     0     1
    ## [2,]     1     1     0     1     1     1     0     0     1     1     0     1
    ## [3,]     1     1     0     1     1     1     1     0     1     0     1     1
    ## [4,]     1     1     0     0     1     1     0     0     1     1     0     1
    ## [5,]     1     1     1     1     1     1     1     1     1     1     0     1
    ## [6,]     1     1     1     1     1     1     1     0     1     1     1     1
    ##      [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
    ## [1,]     1     0     0     1     1     0     1     0     1     0     1     1
    ## [2,]     1     0     0     1     1     1     1     1     0     0     1     1
    ## [3,]     1     0     0     1     1     0     1     1     1     0     1     1
    ## [4,]     1     0     0     1     1     0     1     0     0     0     1     1
    ## [5,]     1     1     1     1     0     1     1     0     1     1     0     1
    ## [6,]     1     0     0     1     1     1     1     0     1     1     1     1
    ##      [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50]
    ## [1,]     0     0     0     0     1     1     0     0     0     1     0     0
    ## [2,]     0     0     1     0     1     0     0     0     0     1     0     0
    ## [3,]     0     0     0     0     0     0     0     0     0     1     0     0
    ## [4,]     0     0     0     1     1     1     0     0     0     1     0     0
    ## [5,]     0     0     0     1     1     0     0     0     0     0     1     1
    ## [6,]     1     1     1     1     0     1     0     1     0     0     0     1

``` r
colnames(Y) = colnames(da)[-c(1:9)]
colnames(Y)
```

    ##  [1] "Phylloscopus_trochilus"  "Turdus_iliacus"         
    ##  [3] "Fringilla_coelebs"       "Turdus_philomelos"      
    ##  [5] "Carduelis_spinus"        "Parus_major"            
    ##  [7] "Anthus_trivialis"        "Erithacus_rubecula"     
    ##  [9] "Cuculus_canorus"         "Muscicapa_striata"      
    ## [11] "Ficedula_hypoleuca"      "Turdus_pilaris"         
    ## [13] "Loxia_curvirostra"       "Columba_palumbus"       
    ## [15] "Sylvia_borin"            "Corvus_corone"          
    ## [17] "Dendrocopos_major"       "Emberiza_citrinella"    
    ## [19] "Turdus_merula"           "Sylvia_curruca"         
    ## [21] "Regulus_regulus"         "Prunella_modularis"     
    ## [23] "Phoenicurus_phoenicurus" "Parus_caeruleus"        
    ## [25] "Parus_montanus"          "Motacilla_alba"         
    ## [27] "Carduelis_chloris"       "Phylloscopus_collybita" 
    ## [29] "Tringa_ochropus"         "Parus_cristatus"        
    ## [31] "Sylvia_communis"         "Grus_grus"              
    ## [33] "Larus_canus"             "Pica_pica"              
    ## [35] "Phylloscopus_sibilatrix" "Turdus_viscivorus"      
    ## [37] "Hirundo_rustica"         "Sylvia_atricapilla"     
    ## [39] "Tetrao_tetrix"           "Pyrrhula_pyrrhula"      
    ## [41] "Dryocopus_martius"       "Certhia_familiaris"     
    ## [43] "Alauda_arvensis"         "Corvus_corax"           
    ## [45] "Saxicola_rubetra"        "Garrulus_glandarius"    
    ## [47] "Numenius_arquata"        "Carpodacus_erythrinus"  
    ## [49] "Gallinago_gallinago"     "Corvus_monedula"

``` r
# We next read the datafile containing species traits, and include in the TrData dataframe data on migratory strategy and body mass
alltraits = read.csv(file.path("data/traits.csv"), stringsAsFactors = TRUE)
cbind(as.character(alltraits[,1]), (colnames(da)[-c(1:9)]))
```

    ##       [,1]                      [,2]                     
    ##  [1,] "Phylloscopus_trochilus"  "Phylloscopus_trochilus" 
    ##  [2,] "Turdus_iliacus"          "Turdus_iliacus"         
    ##  [3,] "Fringilla_coelebs"       "Fringilla_coelebs"      
    ##  [4,] "Turdus_philomelos"       "Turdus_philomelos"      
    ##  [5,] "Carduelis_spinus"        "Carduelis_spinus"       
    ##  [6,] "Parus_major"             "Parus_major"            
    ##  [7,] "Anthus_trivialis"        "Anthus_trivialis"       
    ##  [8,] "Erithacus_rubecula"      "Erithacus_rubecula"     
    ##  [9,] "Cuculus_canorus"         "Cuculus_canorus"        
    ## [10,] "Muscicapa_striata"       "Muscicapa_striata"      
    ## [11,] "Ficedula_hypoleuca"      "Ficedula_hypoleuca"     
    ## [12,] "Turdus_pilaris"          "Turdus_pilaris"         
    ## [13,] "Loxia_curvirostra"       "Loxia_curvirostra"      
    ## [14,] "Columba_palumbus"        "Columba_palumbus"       
    ## [15,] "Sylvia_borin"            "Sylvia_borin"           
    ## [16,] "Corvus_corone"           "Corvus_corone"          
    ## [17,] "Dendrocopos_major"       "Dendrocopos_major"      
    ## [18,] "Emberiza_citrinella"     "Emberiza_citrinella"    
    ## [19,] "Turdus_merula"           "Turdus_merula"          
    ## [20,] "Sylvia_curruca"          "Sylvia_curruca"         
    ## [21,] "Regulus_regulus"         "Regulus_regulus"        
    ## [22,] "Prunella_modularis"      "Prunella_modularis"     
    ## [23,] "Phoenicurus_phoenicurus" "Phoenicurus_phoenicurus"
    ## [24,] "Parus_caeruleus"         "Parus_caeruleus"        
    ## [25,] "Parus_montanus"          "Parus_montanus"         
    ## [26,] "Motacilla_alba"          "Motacilla_alba"         
    ## [27,] "Carduelis_chloris"       "Carduelis_chloris"      
    ## [28,] "Phylloscopus_collybita"  "Phylloscopus_collybita" 
    ## [29,] "Tringa_ochropus"         "Tringa_ochropus"        
    ## [30,] "Parus_cristatus"         "Parus_cristatus"        
    ## [31,] "Sylvia_communis"         "Sylvia_communis"        
    ## [32,] "Grus_grus"               "Grus_grus"              
    ## [33,] "Larus_canus"             "Larus_canus"            
    ## [34,] "Pica_pica"               "Pica_pica"              
    ## [35,] "Phylloscopus_sibilatrix" "Phylloscopus_sibilatrix"
    ## [36,] "Turdus_viscivorus"       "Turdus_viscivorus"      
    ## [37,] "Hirundo_rustica"         "Hirundo_rustica"        
    ## [38,] "Sylvia_atricapilla"      "Sylvia_atricapilla"     
    ## [39,] "Tetrao_tetrix"           "Tetrao_tetrix"          
    ## [40,] "Pyrrhula_pyrrhula"       "Pyrrhula_pyrrhula"      
    ## [41,] "Dryocopus_martius"       "Dryocopus_martius"      
    ## [42,] "Certhia_familiaris"      "Certhia_familiaris"     
    ## [43,] "Alauda_arvensis"         "Alauda_arvensis"        
    ## [44,] "Corvus_corax"            "Corvus_corax"           
    ## [45,] "Saxicola_rubetra"        "Saxicola_rubetra"       
    ## [46,] "Garrulus_glandarius"     "Garrulus_glandarius"    
    ## [47,] "Numenius_arquata"        "Numenius_arquata"       
    ## [48,] "Carpodacus_erythrinus"   "Carpodacus_erythrinus"  
    ## [49,] "Gallinago_gallinago"     "Gallinago_gallinago"    
    ## [50,] "Corvus_monedula"         "Corvus_monedula"

``` r
rownames(alltraits) # 50 species: matches the p=50 columns of the Y matrix
```

    ##  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
    ## [16] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30"
    ## [31] "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45"
    ## [46] "46" "47" "48" "49" "50"

``` r
TrData = data.frame(Species = alltraits$Species, Migration=alltraits$Migration, LogMass = log(alltraits$Mass))
# scale the numeric data
TrData$LogMass = scale(TrData$LogMass)


# X and Traits formulae
XFormula = ~ hab + poly(clim, degree = 2,raw = TRUE) -1
TrFormula = ~ Migration + LogMass -1


# we now should find groups of birds defined by the philogenetic tree
#install.packages("ape")
library(ape)
phyloTree = read.tree("data/CTree.tre")
plot(phyloTree)
```

![](birdsmd_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# from pica pica to prunella modularis (40) we have the same order: passeriformes
# dendrocopos_major and dryocopus_martius has order piciformes
#  Cuculiformes
#  gruiformes
#  Charadriiformes
#  Columbiformes
#  tetrao tetrix
# 7 different orders with one of them that has 40 elements: not so beautiful.

# to separate the big group we can consider the family or superfamily
# 1-5 corvidae
# 6-18 Sylvioidea
# 19 Certhioidea
# 20-29 muscicapoidea
# 30 rguloidea
# 31-40 passeroidea
# an alternative is using
# corvidae (5) - Sylvioidea (13) - muscicapoidea (10) - passeroidea (10) as further variables.
# we can decide to ignore or consider also picidae (2) - Scolopacidae (3)

# Finally we identified 4 or 6 group variables (depending if we ignore the smallest groups or not).

# New TrData

corvoidea =  sylvioidea = muscicapoidea = passeroidea = picidae = scolopacidae = rep(0, 50)
corvoidea[which(colnames(Y) %in% phyloTree$tip.label[1:5])] = 1
sylvioidea[which(colnames(Y) %in% phyloTree$tip.label[6:18])] = 1
muscicapoidea[which(colnames(Y) %in% phyloTree$tip.label[20:29])] = 1
passeroidea[which(colnames(Y) %in% phyloTree$tip.label[31:40])] = 1
picidae[which(colnames(Y) %in% phyloTree$tip.label[41:42])] = 1
scolopacidae[which(colnames(Y) %in% phyloTree$tip.label[45:47])] = 1

TrData_taxonomy = data.frame(TrData,
                             corvoidea,
                             sylvioidea,
                             muscicapoidea,
                             passeroidea,
                             picidae,
                             scolopacidae)
TrFormula_taxonomy = ~ Migration + LogMass + corvoidea +
  sylvioidea + muscicapoidea + passeroidea + picidae + scolopacidae - 1
TrData_taxonomy[, 4]
```

    ##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
    ## [39] 0 0 0 0 0 1 0 1 0 0 0 1

## Prepare data to fit the model

``` r
# maximum number of factors
d = 10
k = 10
p = ncol(Y)

ns = table(da$Route) # units in each study
S = length(unique(da$Route)) # studies
n = sum(ns) # total units in all studies

#shared latent loadings and factors
shared_loadings = rnorm(p * d, sd = 1)
Lambda = matrix(shared_loadings, ncol = d, nrow = p)
Lambda
```

    ##              [,1]         [,2]        [,3]         [,4]        [,5]        [,6]
    ##  [1,]  0.20589449  0.611195653  0.12010808 -0.441184581  0.14515944 -0.28035366
    ##  [2,]  1.96056990  1.689903915 -0.82941873 -0.193474759 -1.19395843  0.59018113
    ##  [3,]  2.03999321 -1.326404288  1.44284575 -1.077562030 -0.53375750  0.76033383
    ##  [4,]  1.01545160  0.413659113  0.53384278  1.147093754  1.00564512 -0.49410208
    ##  [5,] -1.41792769 -0.008187162  1.78710055  0.943229252  0.39305628  0.07022077
    ##  [6,] -1.12765950  0.083935192  0.50510860  0.173388136  0.77408813  1.36463087
    ##  [7,]  2.49006955 -0.118692320 -1.43749979 -1.894970324 -1.99748862  1.07779210
    ##  [8,]  1.04743997 -0.778330233 -0.52880565 -0.464263505 -0.72736915  1.55075965
    ##  [9,]  0.49615037  1.811470477 -1.08430087  0.939859005  0.82489460 -1.93147026
    ## [10,] -0.03776979 -2.100636147 -0.99076331 -0.107025433 -1.49983413  0.76620785
    ## [11,] -1.63973722  0.932243383  0.44435664 -0.395431954 -1.60983574 -1.14515861
    ## [12,] -0.57705836 -0.247235182 -0.23389690  1.932151836  0.68929357 -0.12749856
    ## [13,]  0.16146851  1.011487412  1.66133115  0.189449680  0.29934269  0.76725718
    ## [14,]  0.41521516  0.459630860 -1.20234075  0.669797458 -1.28742442 -1.73212586
    ## [15,]  0.76958363 -0.283103754  0.07710205 -1.206590000  1.18218212 -0.40856465
    ## [16,] -0.13142282  0.443827421  1.27572813  0.002913505  0.70977617  0.22210313
    ## [17,]  1.43960024  1.232614963 -0.14053197  0.708228262  1.29645853  1.19223645
    ## [18,] -0.31563172 -2.086894344  0.28027614 -0.483622474 -1.05634983  0.91737669
    ## [19,]  1.55948009  2.011426639 -0.58057360 -0.015156173 -0.13968256 -0.63051578
    ## [20,] -0.20414269 -0.122368443 -1.17052030  0.436091108  0.10185875 -0.98905340
    ## [21,] -0.46565370 -1.135712305 -0.63084153  0.089014585  0.55856725 -0.39512898
    ## [22,]  0.36747021  1.355721059  0.84212588  0.032330031  1.67772084  0.41713113
    ## [23,] -0.18704450  2.378504681 -0.60291141  0.131970781  1.42720208  0.99483949
    ## [24,] -1.74877922 -0.507664396 -0.15377378 -0.070553736  1.65463663  0.67290100
    ## [25,] -0.59247356  0.734714780  1.78381298  0.028805839  1.71849977 -2.47776399
    ## [26,] -0.89801368  1.218539422 -0.44185057 -0.498309025 -0.79143727 -0.65406677
    ## [27,]  0.12653630 -1.292552832  0.28238344  0.697925487 -0.96808554  0.87690651
    ## [28,]  0.95070864 -0.919127008 -0.86340101 -0.608543174 -0.03403764 -0.27903424
    ## [29,]  1.34011154  1.002112749 -0.45412212 -0.086551349  1.23796525  1.21561172
    ## [30,]  1.07939845 -0.613448599  0.69622356 -0.338555545 -0.74002668  0.23195681
    ## [31,]  0.58508254 -0.718550065 -0.44108642  0.488811243  0.41904273 -0.37545492
    ## [32,]  1.10996817 -0.704007274 -0.70853036  1.085767606 -0.47298198  0.92378131
    ## [33,]  0.62552699 -0.510021797 -1.32331518  2.016392133  0.71667091 -0.61855356
    ## [34,] -0.73228366  0.388629825  0.34411263 -0.342430979  1.10258535 -0.89574903
    ## [35,] -0.58041097 -1.115099513 -1.01608539 -1.405087334 -0.51775612 -0.12055221
    ## [36,]  0.13307615 -1.296743256 -0.62097575  0.638492092  1.43827527 -0.17729853
    ## [37,] -0.76934269 -0.974017902  3.59634838  0.792184564 -2.03397732  0.34355641
    ## [38,]  0.63327489  0.780283776 -2.26878939  0.681770949 -0.84750647  0.60133300
    ## [39,]  0.26226114  1.225675503 -0.43872826  1.161908406 -1.12868564 -0.68105628
    ## [40,]  0.25140159 -1.199787624 -0.52357892  1.681062205 -0.86510584  0.48223305
    ## [41,] -1.16698827  0.027539188 -0.71581376  0.672688413 -0.60136068 -1.08836661
    ## [42,]  0.62351046  0.168928725  0.75135133  0.647691651  1.25876763  1.73468549
    ## [43,]  0.76499765  0.486761344 -0.48020012  1.582295884 -0.08099774  1.15108199
    ## [44,] -1.48617357 -1.535743682  1.02896178  1.710403671  1.05355250  1.20920582
    ## [45,] -0.87246919  2.035955221  0.40938971  0.691175073  0.29095480 -0.13734825
    ## [46,]  0.16025028  0.897193440 -0.90433682  1.574551877 -1.38532130 -0.19223705
    ## [47,]  0.42317784  0.246088313 -1.54699539 -1.058763763 -0.90843942  1.55043879
    ## [48,] -0.23500381  1.044497219 -0.62899893 -1.377677194  0.66459613  1.42697247
    ## [49,] -0.02698086  1.809713264 -0.60422172  0.891875914  0.30564714 -0.41928567
    ## [50,] -0.01785824 -0.009194083  0.56122593 -0.973056133 -0.82959106 -0.33890954
    ##               [,7]       [,8]        [,9]       [,10]
    ##  [1,] -1.308767803  0.7023497  1.29555440 -0.82563012
    ##  [2,] -0.398411756  1.6116980 -1.22282910  0.72527297
    ##  [3,]  0.304595138 -0.2765673 -0.95402480 -0.91527694
    ##  [4,]  0.198901495  0.6317935 -0.78777575 -0.83390507
    ##  [5,]  0.402492188 -0.6258058  0.64686684 -2.46393441
    ##  [6,] -0.586424434  1.5086162 -0.34032012 -0.15820062
    ##  [7,]  0.340474118 -0.3338879  0.66480808  0.40960789
    ##  [8,] -0.220836402  1.1172513  1.96944860 -0.64691660
    ##  [9,]  1.123171581 -1.3774250  0.78445483  0.28475548
    ## [10,] -0.518756551 -0.5894279  1.00609145 -0.52921523
    ## [11,]  0.198422005  1.0722930 -0.78971042  0.10933345
    ## [12,] -1.927602221 -0.6220650  1.55786654  1.08789222
    ## [13,]  1.356157981 -0.4917111  0.05501983  0.07424398
    ## [14,]  0.100865829  0.7676607  0.96567040  1.59182244
    ## [15,]  0.047160446 -0.2286881 -0.29255994 -1.93369600
    ## [16,]  0.232824497 -0.1802512 -0.13299220 -0.18094098
    ## [17,]  0.982456265  0.7842741 -1.24591252  0.45047686
    ## [18,] -1.315186034  0.1151185 -1.23298510 -2.19255958
    ## [19,] -0.932066414  1.1255542 -0.26810200 -0.81615118
    ## [20,] -0.655200995 -1.3140100  0.04300514 -0.54828425
    ## [21,]  1.687845047 -1.4513013 -0.63523236 -0.33813167
    ## [22,]  0.315858897 -0.9225566 -0.31956252  1.53537679
    ## [23,]  0.299383903 -1.2649507 -0.06558706  0.75127055
    ## [24,]  2.295500828  0.2565084 -0.33105370 -0.40441677
    ## [25,]  0.624757493 -0.9532196 -1.24870176 -1.95523431
    ## [26,] -0.631886370  1.0552689 -0.08559175 -0.77581622
    ## [27,]  0.542149410  1.2564608 -1.02228858  1.14393577
    ## [28,]  0.021706160  0.3844282 -1.18143550  0.90656580
    ## [29,]  0.537107121 -0.7975399  1.65513542 -1.31256998
    ## [30,]  0.434836049 -0.9001719  1.45732398 -1.13501203
    ## [31,]  0.406952008 -1.0819109 -1.52127054  0.51938883
    ## [32,]  1.032373538 -0.3861501  1.16137757 -1.69338313
    ## [33,] -0.443724725  0.4294528  0.30346211  1.92380151
    ## [34,]  0.309908205  1.4477822 -0.01557730 -0.94838724
    ## [35,]  0.609745657 -1.0250061  0.26583747 -0.87944996
    ## [36,]  0.650384659 -0.1234489 -0.06758851 -1.00015698
    ## [37,]  2.510164271  0.2003431  0.65557069 -0.47557236
    ## [38,]  0.002591751  1.0128225  1.01865382  0.48240148
    ## [39,] -0.857577678  1.8928244  0.17315216  2.86982149
    ## [40,] -0.171436724 -0.2570032 -1.00544207 -1.05627348
    ## [41,]  0.401780591  0.2071195 -0.33613998 -0.35476521
    ## [42,] -0.878832196  0.5025037  0.80314360 -0.55526810
    ## [43,] -0.656748638 -0.7507851  0.81439499 -0.95488037
    ## [44,]  0.833446901 -2.3594701  0.04274434  0.34587869
    ## [45,] -0.839805108  0.9469267 -0.40543082 -0.33696420
    ## [46,] -1.597265210 -0.4615725  1.41098174  1.25919035
    ## [47,]  0.515531620  1.6919139 -0.73723291  0.65705782
    ## [48,]  0.688659144  1.6120012  0.44988251 -1.04470256
    ## [49,] -0.313976029 -0.5262762  0.67847664  0.12930614
    ## [50,]  0.099477738 -0.3480572 -0.91388041 -1.77015428

``` r
eta = matrix(NA, ncol = d, nrow = n)
for (h in 1:d) {
  eta[, h] = rnorm(n)
}
# specific active factors
ks = c(rep(1, S))#specific factors study (important just for the ground truth in simulations)
ks[cumsum(ks) > k] = 0 # set some factors to "inactive"
sum(ks)
```

    ## [1] 10

``` r
#assign group labels (location labels)
group = rep(NA, n)
for (s in 1:S) {
  nscumpre = ifelse(s > 1, sum(ns[1:(s - 1)]) + 1, 1)
  nscum = sum(ns[1:s])
  group[nscumpre:nscum] = s
}
group
```

    ##   [1]   1   1   1   1   2   2   2   2   2   3   3   3   3   3   3   3   4   4
    ##  [19]   4   4   4   4   4   5   5   5   5   5   6   6   6   6   7   7   7   7
    ##  [37]   8   8   8   8   9   9   9   9   9   9  10  10  10  10  10  11  11  11
    ##  [55]  11  11  11  11  12  12  12  12  12  12  13  13  13  13  13  13  13  13
    ##  [73]  13  14  14  14  14  14  14  14  15  15  15  15  16  16  16  16  16  16
    ##  [91]  16  17  17  17  17  18  18  18  18  18  19  19  19  19  19  20  20  20
    ## [109]  20  20  21  21  21  21  21  22  22  22  22  23  23  23  23  23  24  24
    ## [127]  24  24  25  25  25  25  25  25  25  26  26  26  26  26  27  27  27  27
    ## [145]  27  28  28  28  28  28  29  29  29  29  29  30  30  30  30  31  31  31
    ## [163]  31  32  32  32  32  33  33  33  33  33  34  34  34  34  35  35  35  35
    ## [181]  35  36  36  36  36  37  37  37  37  37  37  37  38  38  38  38  38  39
    ## [199]  39  39  39  40  40  40  40  40  40  41  41  41  41  42  42  42  42  42
    ## [217]  42  42  43  43  43  43  43  44  44  44  44  44  44  45  45  45  45  45
    ## [235]  46  46  46  46  46  47  47  47  47  48  48  48  48  49  49  49  49  49
    ## [253]  49  50  50  50  50  50  50  50  51  51  51  51  51  51  51  52  52  52
    ## [271]  52  52  52  53  53  53  53  53  54  54  54  54  55  55  55  55  56  56
    ## [289]  56  56  56  56  56  56  56  57  57  57  57  57  58  58  58  58  58  59
    ## [307]  59  59  59  59  59  60  60  60  60  60  60  61  61  61  61  61  61  62
    ## [325]  62  62  62  62  62  62  62  62  63  63  63  63  64  64  64  64  64  65
    ## [343]  65  65  65  66  66  66  66  67  67  67  67  68  68  68  68  69  69  69
    ## [361]  69  70  70  70  70  71  71  71  71  72  72  72  72  73  73  73  73  74
    ## [379]  74  74  74  75  75  75  75  76  76  76  76  76  76  77  77  77  77  77
    ## [397]  78  78  78  78  78  78  79  79  79  79  79  79  79  80  80  80  80  81
    ## [415]  81  81  81  81  81  81  81  82  82  82  82  82  82  82  82  82  83  83
    ## [433]  83  83  84  84  84  84  84  84  84  85  85  85  85  85  86  86  86  86
    ## [451]  87  87  87  87  87  87  87  88  88  88  88  89  89  89  89  90  90  90
    ## [469]  90  90  91  91  91  91  92  92  92  92  92  93  93  93  93  93  94  94
    ## [487]  94  94  94  95  95  95  95  95  96  96  96  96  96  96  96  97  97  97
    ## [505]  97  98  98  98  98  99  99  99  99  99 100 100 100 100 100 101 101 101
    ## [523] 101 101 102 102 102 102 103 103 103 103 103 104 104 104 104 105 105 105
    ## [541] 105 105 105 105 106 106 106 106 106 107 107 107 107 107 108 108 108 108
    ## [559] 109 109 109 109 110 110 110 110 110 110 110 110 111 111 111 111 111 111
    ## [577] 112 112 112 112 112 112 113 113 113 113 113 114 114 114 114 115 115 115
    ## [595] 115 115 115 116 116 116 116 116 116 117 117 117 117 118 118 118 118 119
    ## [613] 119 119 119 119 119 119 119 119 120 120 120 120 120 120 121 121 121 121
    ## [631] 121 121 121 121 121 122 122 122 122 122 122 122 122 123 123 123 123 124
    ## [649] 124 124 125 125 125 126 126 126 126 127 127 127 128 128 128 129 129 129
    ## [667] 130 130 130 131 131 131 131 132 132 132 133 133 133 134 134 134 135 135
    ## [685] 135 136 136 136 136 136 137 137 137 137 137 138 138 138 139 139 139 140
    ## [703] 140 140 141 141 141 142 142 142 143 143 143 144 144 144 145 145 145 146
    ## [721] 146 146 147 147 147 148 148 148 149 149 149 149 150 150 150 150 151 151
    ## [739] 151 151 152 152 152 153 153 153 154 154 154 155 155 155 156 156 156 157
    ## [757] 157 157 158 158 158 158 158 159 159 159 159 160 160 160 161 161 161 162
    ## [775] 162 162 163 163 163 164 164 164 165 165 165 166 166 166 167 167 167 167
    ## [793] 168 168 168 169 169 169 170 170 170 171 171 171 172 172 172 173 173 173
    ## [811] 173 174 174 174 174 175 175 175 175 176 176 176 177 177 177 178 178 178
    ## [829] 179 179 179 180 180 180 181 181 181 182 182 182 183 183 183 183 184 184
    ## [847] 184 184 185 185 185 186 186 186 186 187 187 187 188 188 188 188 188 189
    ## [865] 189 189 190 190 190 191 191 191 192 192 192 193 193 193 194 194 194 195
    ## [883] 195 195 195 196 196 196 196 196 196 197 197 197 197 198 198 198 198 198
    ## [901] 198 199 199 199 199 199 199 199 200 200 200 200 200 200

``` r
table(group)
```

    ## group
    ##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
    ##   4   5   7   7   5   4   4   4   6   5   7   6   9   7   4   7   4   5   5   5 
    ##  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40 
    ##   5   4   5   4   7   5   5   5   5   4   4   4   5   4   5   4   7   5   4   6 
    ##  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60 
    ##   4   7   5   6   5   5   4   4   6   7   7   6   5   4   4   9   5   5   6   6 
    ##  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80 
    ##   6   9   4   5   4   4   4   4   4   4   4   4   4   4   4   6   5   6   7   4 
    ##  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
    ##   8   9   4   7   5   4   7   4   4   5   4   5   5   5   5   7   4   4   5   5 
    ## 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 
    ##   5   4   5   4   7   5   5   4   4   8   6   6   5   4   6   6   4   4   9   6 
    ## 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 
    ##   9   8   4   3   3   4   3   3   3   3   4   3   3   3   3   5   5   3   3   3 
    ## 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 
    ##   3   3   3   3   3   3   3   3   4   4   4   3   3   3   3   3   3   5   4   3 
    ## 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 
    ##   3   3   3   3   3   3   4   3   3   3   3   3   4   4   4   3   3   3   3   3 
    ## 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 
    ##   3   3   4   4   3   4   3   5   3   3   3   3   3   3   4   6   4   6   7   6

``` r
# convert into dummy variables the group labels
X = matrix(NA, ncol = S, nrow = n)
X[, 1:S] = model.matrix(rep(1, n) ~ -1 + as.factor(group))

# initialize specific loadings and factors
specific_loadings = rnorm(p * k, sd = 1)
Gamma = matrix(specific_loadings, ncol = k, nrow = p)
phi_ = phi = matrix(NA, ncol = k, nrow = n)
for (h in 1:k) {
  phi[, h] = phi_[, h] = rnorm(n)
}
for (i in 1:n){
for (h in 1:k) {
  phi[i, h] = phi_[i, h] * (group[i] == h) # set 0s blocks of different groups
}
}
head(phi,9)
```

    ##              [,1]       [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,] -1.11550728  0.0000000    0    0    0    0    0    0    0     0
    ##  [2,]  0.69465200  0.0000000    0    0    0    0    0    0    0     0
    ##  [3,]  0.02872015  0.0000000    0    0    0    0    0    0    0     0
    ##  [4,]  1.02850162  0.0000000    0    0    0    0    0    0    0     0
    ##  [5,]  0.00000000 -0.8470869    0    0    0    0    0    0    0     0
    ##  [6,]  0.00000000 -0.6122160    0    0    0    0    0    0    0     0
    ##  [7,]  0.00000000 -0.6076651    0    0    0    0    0    0    0     0
    ##  [8,]  0.00000000 -0.4008893    0    0    0    0    0    0    0     0
    ##  [9,]  0.00000000  0.8520305    0    0    0    0    0    0    0     0

## Set prior hyperparameters

``` r
alpha_eta = 5 # number of active shared factors a priori
v_eta = c(rbeta(d - 1, shape1 = 1, shape2 = alpha_eta), 1)
w_eta = v_eta * c(1, cumprod(1 - v_eta[-d]))                     # weights
z_eta = rep(d, d)
alpha_phi = 5 # number of active specific factors a priori
v_phi = c(rbeta(k - 1, shape1 = 1, shape2 = alpha_phi), 1)
w_phi = v_phi * c(1, cumprod(1 - v_phi[-k]))                     # weights
z_phi = rep(k, k)
whichgroup = unique(sapply(1:n , function (k)
  (which(X[k, ] == 1))))

# initailize regression coefficients for modelling the specific factors' activation pattern
# based on group labels
betas = matrix(0, ncol = k, nrow = S)
```

## define the state of the MCMC chain

``` r
state = list(
  y = Y,
  # response
  Lambda = Lambda + rnorm(prod(dim(Lambda))),
  Lambda_ = Lambda + rnorm(prod(dim(Lambda)), sd =0.1),
  # shared loadings
  eta = (eta),
  # shared factors
  Gamma = Gamma + rnorm(prod(dim(Gamma)), sd = 0.1),
  #specific loadings
  phi = (phi) ,
  phi_ = (phi_) + rnorm(prod(dim(phi_)), sd = 0.1),
  #sparse and non sparse specific factors
  # list of specific covariance matrices
  n = n,
  ns = ns,
  X = X,
  S = S,
  d = d,
  k = k,
  p = p,
  #prior
  a_lambda = 1,
  b_lambda = 2,
  a_gamma = 1,
  b_gamma = 2,
  tau_eta = c(c(rep(1, 10)), c(rep(0, d - 10))),
  tau_phi = c(c(rep(1, 10)), c(rep(0, k - 10))),
  z_eta = z_eta,
  z_phi = z_phi,
  w_eta = w_eta,
  w_phi = w_phi,
  v_eta = v_eta,
  v_phi = v_phi,
  alpha_eta = 5,
  alpha_phi = 5,
  # equal tO expected n.of active factors
  betas = betas,
  ps = matrix(rbinom(n * k, 1, 0.1), ncol = k)
)

copy_state = state
```

## Run the Gibbs sampler

``` r
Gaussian = F
#pre-allocate memory for results
maxiter = 1
ris_phi = array(dim = c(maxiter, dim(state$phi)))
ris_eta = array(dim = c(maxiter, dim(state$eta)))
ris_th = array(dim = c(maxiter, dim(state$th)))
ris_ps = array(dim = c(maxiter, dim(state$ps)))

ris_beta = array(dim = c(maxiter, dim(state$betas)))
ris_eta = array(dim = c(maxiter, dim(state$eta)))
ris_lambda = array(dim = c(maxiter, dim(state$Lambda)))
ris_gamma = array(dim = c(maxiter, dim(state$Gamma)))
ris_tau_eta = array(dim = c(maxiter, length(state$tau_eta)))
ris_tau_phi = array(dim = c(maxiter, length(state$tau_phi)))

set.seed(111)

maxiter=10000
#maxiter=10000
iter=1
run=F
if (run==T){
for (iter in 1:maxiter) {
  cat(iter)
  state = Gibbs_Kernel_non_gauss(state)
  #shrinkage tau
  ris_tau_eta[iter, ] = state$tau_eta
  ris_tau_phi[iter, ] = state$tau_phi
  #factors
  ris_phi[iter, , ] = state$phi
  ris_eta[iter, , ] = state$eta
 
  #beta
  ris_beta[iter, , ] = state$betas

  #local activation psi
  ris_ps[iter, , ] = state$ps
 
  #loadings
  ris_lambda[iter, , ] = state$Lambda
  ris_gamma[iter, , ] = state$Gamma
 
  #print iteration and number of active factors
  if(iter%%100==0){
    print(state$tau_eta)
    print(state$tau_phi)
    }
}
}
#save.image("birds_results.RData")
```

# Upload results

``` r
load("birds_results.RData")
```

``` r
par(mfrow=c(1,1))
#specific factors vs type of habitat
levels(da[,6])=c("Broadleaved","Conifer","Open","Urban","Wetlands")
for (i in 6:6) plot(colMeans(ris_phi[5000:10000,,1])~da[,i],
                    xlab=colnames(da)[i],ylab=expression(varphi[1]))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Posterior mean of the factors ordering units by latitude

``` r
# Define the matrix to be visualized 
# Define the color palette with red and blue
custom_colors <- c("red","white", "blue")

# Create a custom color function that maps values to colors
color_function <- colorRampPalette(custom_colors)
colors <- color_function(0.3*length(GG1))
 
#head(colMeans((ris_phi[5000:10000,,])*ris_ps[5000:10000,,]))

par(mfrow=c(1,2))
image(t(colMeans((ris_phi[5000:10000,,1:5])*
                    ris_ps[5000:10000,,1:5])),
ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```

![](birdsmd_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Posterior mean ordering by latitude and habitat

``` r
oor=order(da$Habitat)
cuts=table(da$Habitat[oor])/sum(table(da$Habitat[oor]))
image(t(colMeans((ris_phi[5000:10000,oor,1:5])*
                   ris_ps[5000:10000,oor,1:5])),
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Upper and lower bounds of credible intervals for specific factors

``` r
colQuant90=function (m) apply(m,c(2,3), function(x) quantile(x,0.9))
colQuant10=function (m)apply(m,c(2,3),function(x) quantile(x,0.1))

par(mfrow=c(1,2))

image(t(colQuant10((ris_phi[5000:10000,,1:5])*
                   ris_ps[5000:10000,,1:5]))
      ,  main="q=0.1",
      ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```

![](birdsmd_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
image(t(colQuant90((ris_phi[5000:10000,,1:5])*
                     ris_ps[5000:10000,,1:5])),
      main="q=0.9",
      ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```

![](birdsmd_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
colQuant90=function (m) apply(m,c(2,3), function(x) quantile(x,0.9))
colQuant10=function (m)apply(m,c(2,3),function(x) quantile(x,0.1))


 

image(t(colQuant10((ris_phi[5000:10000,oor,1:5])*
                     ris_ps[5000:10000,oor,1:5]))
      ,  main="q=0.1",
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
image(t(colQuant90((ris_phi[5000:10000,oor,1:5])*
                     ris_ps[5000:10000,oor,1:5])),
      main="q=0.9",
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# covariance matrix
G=colMeans(colMeans(ris_tau_phi[5000:10000, ])*(ris_gamma[5000:10000,, ]))
G
```

    ##              [,1]          [,2]         [,3]         [,4]        [,5]
    ##  [1,]  1.42339442 -1.3615038239  0.004461156 -1.800490362 -0.73783587
    ##  [2,]  0.67965061 -1.1538778315  0.292678587 -0.869366216 -0.35406064
    ##  [3,]  1.42703049 -1.2678069733 -0.082215852 -0.560396485  1.07491704
    ##  [4,]  0.80751092 -1.3399529564 -0.028331256 -1.221098934  0.45734470
    ##  [5,]  0.82003146 -1.1664052105 -0.071927281 -1.294006618  0.73455631
    ##  [6,]  1.61826991 -0.9197783581 -0.016520694 -0.154234617  0.77875839
    ##  [7,]  0.54367477 -1.2050684256 -0.002772315 -1.058337386  0.95207151
    ##  [8,]  1.02891489 -1.2092645875 -0.140548517 -0.183039893  1.00452862
    ##  [9,]  0.20470908 -0.6286026529  0.199068051 -0.731936801 -0.03120222
    ## [10,]  0.41123101 -0.3002943906 -0.038243170 -0.504429306  0.55346326
    ## [11,]  0.78793931 -0.2516283596 -0.159776292 -0.005046139  0.50240056
    ## [12,]  0.85124263 -0.8868492931  0.350578553  0.187879903  0.05059663
    ## [13,]  0.11922388 -0.4091561613 -0.031201161 -0.543685693  0.31281854
    ## [14,]  0.97701100 -0.7544672025  0.405283328  0.633476279  1.06582226
    ## [15,]  1.26027135 -1.0393196574  0.090006644  1.245083773  0.88800502
    ## [16,]  0.92462999 -0.3371668064  0.251953232  0.172802566  0.19024929
    ## [17,]  0.41563084 -0.2278742340  0.103886871 -0.206391836  0.73253987
    ## [18,]  0.78416389 -0.7309783284  0.126945553  1.000625165  1.25570909
    ## [19,]  1.11096012 -0.5414836835  0.150705674  1.334269621  1.19146461
    ## [20,]  0.35922916 -0.2786138388 -0.074073137  0.511716768  0.85945503
    ## [21,]  0.24603430 -0.6099112732  0.213566153  0.417184672  1.19742218
    ## [22,]  0.02282406 -0.7138952929  0.077203378  0.186797503  0.55939699
    ## [23,]  0.06804494  0.0230388452 -0.366704030 -1.043255916 -0.05526636
    ## [24,]  0.83904429 -0.1904379383  0.125010153  1.217885880  0.85736896
    ## [25,] -0.21734943 -0.3930393970 -0.002573745 -0.113401379  0.41605576
    ## [26,]  0.67718668 -0.0479038386  0.019659426  0.329468658  0.42751404
    ## [27,]  0.67123998 -0.1396522849  0.180103112  0.504862733  0.56103903
    ## [28,] -0.08085610 -0.3605854078  0.109778389  0.925911014  1.07340284
    ## [29,] -0.30203543 -0.3469521104 -0.113286698  0.245270004  0.65268452
    ## [30,] -0.07422771 -0.0810788608 -0.166592992  0.958302798  0.81765456
    ## [31,]  0.88475088 -0.0786859133  0.181200878  1.675580450  1.28268695
    ## [32,] -0.15751268  0.0872761479  0.144786628 -0.054494315  0.49658081
    ## [33,]  0.38046595 -0.1712232492 -0.097246519  0.606136773  0.19000197
    ## [34,]  0.63630962 -0.1364705180  0.413303675  0.730061212  0.51040947
    ## [35,]  0.08290974 -0.1230980405 -0.132633185  0.736638612  0.59215997
    ## [36,] -0.25735883  0.0942393917 -0.043942023 -0.068885558  0.41010745
    ## [37,]  0.34679376 -0.1034432675  0.005766582  0.718418442  0.57584538
    ## [38,]  0.60632314  0.2184077877 -0.005208729  1.574722808  1.07264573
    ## [39,] -0.36191943  0.0006376525  0.004464052 -0.089476971  0.83502912
    ## [40,] -0.31852971 -0.1513989067  0.089690292  0.005033480  0.49530705
    ## [41,] -0.18068674 -0.0192790406  0.092978392  0.257575116  0.75816706
    ## [42,] -0.09085583 -0.0196577065  0.060786755  1.073797461  0.85790184
    ## [43,]  0.55530505  0.3985647125  0.393215330  1.713527713  0.98061622
    ## [44,] -0.10277652  0.0549816686  0.046194442  0.074925480  0.21467934
    ## [45,] -0.14191029 -0.0129209391  0.252846390  0.045856232  0.66942116
    ## [46,] -0.19138870 -0.0527336559  0.168226600  0.819794023  0.58903840
    ## [47,] -0.23470313 -0.2694085297  0.358288078  0.470938482  0.88954633
    ## [48,] -0.08955774 -0.4457514685 -0.311458867  0.894203952  0.73484576
    ## [49,] -0.61864419 -0.1146549691  0.103735966 -0.150571164  0.12194357
    ## [50,]  0.48166652  0.2645390712  0.364982722  1.625617952  0.82199875
    ##                [,6]         [,7]          [,8]          [,9]        [,10]
    ##  [1,]  0.0008392635 -0.016824175 -0.0011334395  0.0168253868  0.095224161
    ##  [2,]  0.0139618034  0.039177710  0.0119688418  0.0224069633 -0.008779649
    ##  [3,] -0.0446555522 -0.003687957  0.0202403275  0.0166958441  0.007096530
    ##  [4,]  0.0119475660  0.013146413 -0.0005738770 -0.0340153498 -0.019664456
    ##  [5,] -0.0387340361  0.001464265  0.0037101217  0.0050125845  0.067353773
    ##  [6,] -0.0017164341  0.044563354 -0.0173929359  0.0335920136 -0.054996720
    ##  [7,] -0.0012263008  0.027096433  0.0592413865 -0.0021036661 -0.024870420
    ##  [8,]  0.0171756632  0.048721958  0.0281796428  0.0017418423 -0.031271394
    ##  [9,] -0.0032473188  0.029221141 -0.0140690574 -0.0345378208  0.076047726
    ## [10,] -0.0320942299 -0.012220057  0.0187572125 -0.0149315409  0.052161466
    ## [11,] -0.0086301454 -0.026485210 -0.0078433552  0.0118313409  0.067814937
    ## [12,]  0.0032292633 -0.048239584  0.0040309357  0.0112368986 -0.032897483
    ## [13,]  0.0046873718  0.029475735 -0.0007447232  0.0057026827 -0.037692734
    ## [14,] -0.0034300487  0.025562003 -0.0015134767 -0.0185334919  0.075735521
    ## [15,]  0.0554390698  0.016599942  0.0034041718  0.0198952306  0.033255975
    ## [16,] -0.0271785111  0.022581327  0.0013101020 -0.0282220368  0.170511881
    ## [17,] -0.0098243028 -0.025325167 -0.0098306022 -0.0026520743 -0.029131953
    ## [18,] -0.0124616429  0.024833284 -0.0067991375  0.0431519474  0.101520674
    ## [19,]  0.0679373312 -0.039326102 -0.0355056292 -0.0084481919  0.061073072
    ## [20,] -0.0050627160 -0.072435778 -0.0091099095 -0.0390601413  0.004093658
    ## [21,] -0.0095103076 -0.044865002 -0.0016834764 -0.0223432438  0.010742340
    ## [22,] -0.0052185813 -0.015631971 -0.0109451665 -0.0221897175  0.031026545
    ## [23,] -0.0340453394 -0.015974732  0.0261021342  0.0330409755 -0.063758571
    ## [24,] -0.0536630154  0.015447591 -0.0209305942  0.0038243624 -0.002071017
    ## [25,]  0.0124349759 -0.047335732  0.0245308753  0.0332748469  0.041126782
    ## [26,]  0.0029793297 -0.075307024 -0.0412027669  0.0334773841 -0.105578520
    ## [27,] -0.0217264004 -0.006448082  0.0454183475 -0.0026292359 -0.034213252
    ## [28,] -0.0552938165  0.025716199  0.0065710241  0.0412353460  0.039599196
    ## [29,] -0.0058612692  0.049593357  0.0130355063 -0.0305895923 -0.001343884
    ## [30,]  0.0699603775 -0.021388374 -0.0372013705 -0.0035613423  0.001587935
    ## [31,] -0.0188552206  0.003767734  0.0291956536  0.0291150067 -0.060455928
    ## [32,] -0.0446345005  0.011395649 -0.0048342047  0.0264944103 -0.018461671
    ## [33,]  0.0493487083  0.021045399  0.0145482807  0.0154922328  0.002356927
    ## [34,]  0.0352761170 -0.004976012 -0.0129006096  0.0009067167 -0.066001866
    ## [35,]  0.0652124760  0.003679355  0.0563902995  0.0325026104  0.114378251
    ## [36,] -0.0124510404 -0.080751771 -0.0072530283 -0.0063316596  0.084757306
    ## [37,]  0.0748714506 -0.006031496  0.0101968389  0.0134455957  0.058669090
    ## [38,] -0.0015988642  0.014495913  0.0156495604  0.0114076394  0.009905974
    ## [39,] -0.0269434962  0.035253737 -0.0137538740  0.0024117308 -0.044162913
    ## [40,] -0.0515581215 -0.055784950  0.0482264767 -0.0194545297 -0.014504762
    ## [41,]  0.0052048383 -0.061742489  0.0179333838 -0.0107983815 -0.062618395
    ## [42,] -0.0202418264 -0.010730240 -0.0417296274  0.0157424087  0.054349033
    ## [43,] -0.0250623364  0.061986766 -0.0037777808 -0.0222739252  0.182920101
    ## [44,]  0.0110645713 -0.031768538  0.0221769048 -0.0255729927  0.007567578
    ## [45,]  0.0102390750  0.004921886 -0.0194977122  0.0231484955  0.013371552
    ## [46,]  0.0263247123 -0.058906324 -0.0099687766 -0.0098577671  0.124497254
    ## [47,]  0.0028188679  0.013205405  0.0307979754 -0.0153341552 -0.055421775
    ## [48,] -0.0035858843  0.004981030 -0.0072611719  0.0081454227  0.101664092
    ## [49,]  0.0097575263 -0.042198340  0.0196404395 -0.0200316922  0.015831011
    ## [50,]  0.0026305246 -0.013146595  0.0326801325 -0.0276357112  0.045485189

``` r
L=colMeans(colMeans(ris_tau_eta[5000:10000, ])*
             (ris_lambda[5000:10000,, ]))
L
```

    ##                [,1]         [,2]         [,3]         [,4]          [,5] [,6]
    ##  [1,]  3.271290e-01 -0.203280953  0.269712611 -0.030921149 -0.0069605111    0
    ##  [2,]  1.464034e-01 -0.033748143  0.170645717 -0.118301748 -0.1103905525    0
    ##  [3,]  4.090885e-01 -0.237096810  0.301220304  0.003082606  0.0226637175    0
    ##  [4,]  3.175262e-01 -0.295249364  0.355387881  0.017732554  0.0120651572    0
    ##  [5,]  3.022904e-01 -0.232184122  0.287777642  0.084107991  0.0279528883    0
    ##  [6,]  2.886941e-01 -0.136578882  0.262234951 -0.017330104  0.0032615398    0
    ##  [7,]  4.523378e-01 -0.246089539  0.384581029 -0.029721010  0.0084848874    0
    ##  [8,]  3.345256e-01 -0.312644338  0.378929689  0.048947056  0.0373257604    0
    ##  [9,]  2.276799e-01 -0.092333917  0.186062219  0.098921374  0.0039599303    0
    ## [10,]  7.348724e-02 -0.030506058  0.141013068  0.002940930  0.0486578816    0
    ## [11,]  4.347305e-02 -0.022627311  0.059487339  0.054169472  0.0162304475    0
    ## [12,]  2.537137e-02  0.097316369 -0.017617735  0.015788268  0.0291828636    0
    ## [13,]  9.884554e-02 -0.025826133  0.153352185 -0.071935212  0.0134983857    0
    ## [14,]  2.254772e-01 -0.165036831  0.199754179  0.000192909 -0.0301370203    0
    ## [15,]  2.034938e-01 -0.090731406  0.240932360 -0.116611750 -0.0731300399    0
    ## [16,]  4.389830e-02  0.087975468 -0.043863843  0.013006943  0.0327450278    0
    ## [17,]  6.837058e-02 -0.033829405  0.098245752 -0.027643480 -0.0334977205    0
    ## [18,]  2.187096e-01 -0.010437037  0.113535262 -0.167140929  0.0098209077    0
    ## [19,]  2.024593e-01 -0.218551287  0.270212023  0.018168491 -0.0684436542    0
    ## [20,]  1.247934e-01 -0.099909406  0.135964123 -0.035455446 -0.0006077091    0
    ## [21,]  1.455496e-01 -0.191128578  0.247470783  0.015284056  0.0479804612    0
    ## [22,]  1.220931e-01 -0.216043091  0.211818012  0.019840405  0.0042700740    0
    ## [23,] -2.485792e-02  0.021548768 -0.035743679  0.118380179  0.0618354434    0
    ## [24,]  6.477862e-02 -0.087918202  0.092854105  0.065904438  0.0609839677    0
    ## [25,]  1.109308e-01 -0.107879758  0.093043977  0.043536936  0.0437089402    0
    ## [26,] -4.607957e-02  0.107616984 -0.038398718  0.006902376  0.0449429551    0
    ## [27,] -6.896702e-02  0.132108065 -0.050045083 -0.119536768  0.0013785646    0
    ## [28,]  1.201853e-01 -0.117358585  0.147398685 -0.049251488 -0.0034417681    0
    ## [29,]  1.155586e-01 -0.121758318  0.120630085  0.014608417  0.0262137705    0
    ## [30,]  1.007331e-01 -0.063341397  0.171805062  0.018816922  0.0759786774    0
    ## [31,] -6.619958e-02  0.088305808  0.022354803 -0.050207811 -0.0418049638    0
    ## [32,]  8.420517e-02  0.012881153  0.028322789  0.039590953  0.0069168591    0
    ## [33,] -1.064372e-01  0.119171111 -0.118517265  0.142855293  0.0735810222    0
    ## [34,] -1.146141e-01  0.151239049 -0.170645039 -0.042927267 -0.0049927463    0
    ## [35,] -1.452676e-02 -0.044019925  0.095290777  0.055761106 -0.0639714748    0
    ## [36,]  1.678566e-02 -0.004474891  0.098692370  0.059489271 -0.0142583403    0
    ## [37,] -6.707216e-02  0.142990643 -0.073862100 -0.058119822  0.0021174829    0
    ## [38,] -7.254449e-02 -0.061849120  0.115584236  0.006560162 -0.0764365160    0
    ## [39,]  4.652047e-02 -0.058013252  0.030680000  0.074573793  0.0090752267    0
    ## [40,] -2.501403e-02 -0.104170592  0.015508027  0.052174865  0.0199496009    0
    ## [41,]  3.070377e-03  0.041619032  0.027784011 -0.028943288  0.0094111048    0
    ## [42,]  1.712485e-02 -0.101923848  0.175361735  0.022916934  0.0185511705    0
    ## [43,]  2.486515e-02  0.258465758 -0.123906360 -0.235835488  0.0124349761    0
    ## [44,] -1.375685e-02  0.014859143 -0.005765718  0.056003066  0.0345924657    0
    ## [45,] -6.554746e-05  0.199079678 -0.147914375 -0.188318292  0.0196290481    0
    ## [46,]  1.563637e-02 -0.006506866  0.061110667  0.020918895  0.0006643288    0
    ## [47,] -1.208315e-02  0.207793506 -0.284911479  0.006259155  0.0809434645    0
    ## [48,] -1.184321e-01  0.138421139 -0.113622864 -0.093905138 -0.0202269722    0
    ## [49,] -3.246849e-02  0.051219864 -0.046910987  0.021428308 -0.0201635517    0
    ## [50,] -1.514338e-01  0.190734545 -0.218210651  0.080217376  0.0806449485    0
    ##       [,7] [,8] [,9] [,10]
    ##  [1,]    0    0    0     0
    ##  [2,]    0    0    0     0
    ##  [3,]    0    0    0     0
    ##  [4,]    0    0    0     0
    ##  [5,]    0    0    0     0
    ##  [6,]    0    0    0     0
    ##  [7,]    0    0    0     0
    ##  [8,]    0    0    0     0
    ##  [9,]    0    0    0     0
    ## [10,]    0    0    0     0
    ## [11,]    0    0    0     0
    ## [12,]    0    0    0     0
    ## [13,]    0    0    0     0
    ## [14,]    0    0    0     0
    ## [15,]    0    0    0     0
    ## [16,]    0    0    0     0
    ## [17,]    0    0    0     0
    ## [18,]    0    0    0     0
    ## [19,]    0    0    0     0
    ## [20,]    0    0    0     0
    ## [21,]    0    0    0     0
    ## [22,]    0    0    0     0
    ## [23,]    0    0    0     0
    ## [24,]    0    0    0     0
    ## [25,]    0    0    0     0
    ## [26,]    0    0    0     0
    ## [27,]    0    0    0     0
    ## [28,]    0    0    0     0
    ## [29,]    0    0    0     0
    ## [30,]    0    0    0     0
    ## [31,]    0    0    0     0
    ## [32,]    0    0    0     0
    ## [33,]    0    0    0     0
    ## [34,]    0    0    0     0
    ## [35,]    0    0    0     0
    ## [36,]    0    0    0     0
    ## [37,]    0    0    0     0
    ## [38,]    0    0    0     0
    ## [39,]    0    0    0     0
    ## [40,]    0    0    0     0
    ## [41,]    0    0    0     0
    ## [42,]    0    0    0     0
    ## [43,]    0    0    0     0
    ## [44,]    0    0    0     0
    ## [45,]    0    0    0     0
    ## [46,]    0    0    0     0
    ## [47,]    0    0    0     0
    ## [48,]    0    0    0     0
    ## [49,]    0    0    0     0
    ## [50,]    0    0    0     0

``` r
#### cov
LL <- tcrossprod(L, L)

ordine <- 50*as.numeric(as.factor(TrData_taxonomy$Migration=="R"))+
  100*as.numeric(as.factor(TrData_taxonomy$Migration)=="S")+
  150*as.numeric(as.factor(TrData_taxonomy$Migration)=="L")+
  as.numeric((TrData_taxonomy$LogMass))
ordine=order(ordine)
```
