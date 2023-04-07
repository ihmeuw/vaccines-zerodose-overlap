---
title: "Vaccines Zero-Dose Overlap"
output: html_document
date: "2023-04-07"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This repo contains cleaned code used to:

1) Generate estimates of routine first-dose diphtheria-tetanus-pertussis-containing-vaccine (DTP1) coverage.
2) Assess the potential for integrated health service delivery opportunities by analysing the degree of overlap in prevalence of children who have never received a dose of DTP (no-DTP) with other important health indicators in select countries. 
  
This code was used in producing the analyses described in the following citation:

  Haeuser, E.; Nguyen, J.Q.; Rolfe, S.; Nesbit, O.; Fullman, N.; Mosser, J.F. Assessing Geographic Overlap between Zero-Dose Diphtheria–Tetanus–Pertussis Vaccination Prevalence and Other Health Indicators. Vaccines 2023, 11, 802. https://doi.org/10.3390/vaccines11040802

This repo contains three directories:

1) *lbd_core* contains general code for model based geostatistics used across Local Burden of Disease projects at IHME.
2) *vaccine* contains code that is specific to data processing and estimating vaccination coverage.
3) *vaccines* contains wrapper scripts, supporting data and code for the overlap analyses.

The paper can be found here: https://doi.org/10.3390/vaccines11040802. The Global Health Data Exchange (GHDx) record, with tabular estimates and secondary analytical results, can be found here: http://ghdx.healthdata.org/record/ihme-data/dtp-vaccines-zerodose-overlap.