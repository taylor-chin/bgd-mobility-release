# Code for "Bias in mobility datasets drives divergence in modeled outbreak dynamics" by Chin et al. (2023)

Data
-----------------------

CDRs used in the analysis are available upon request from a2i (see Data Availability statement in manuscript for more detail). Meta’s Data For Good is available to researchers from Meta upon request. Descriptive data sources, including shapefiles and data on population, mean household income, and nightlight levels, are publicly available (see Methods). 

Code
-----------------------
Two scripts have functions used in multiple scripts:

* `mat_utils.R` -  matrix functions used in `utils.R` and multiple scripts
* `utils.R` - loads and processes shapefiles and demographic data sources; includes functions for processing CDR data for origin-destination matrices

Scripts are run in numerical order from `01_process_cdr.R` to `06_descrip_compare.R`:

* `01_process_cdr.R` - processes raw CDR data to get origin-destination matrices and matrices with beta distribution parameters
* `02_impute_cdr_2017.R` - processes Operator 1's 2017 data for imputing missing Dhaka values for Operators 1 and 3 in 2020
* `03_fb_movement.R` - processes Meta's Data For Good into origin-destination matrices
* `04_gravity_model.R` - creates upazila- and district-level gravity models using parameter estimates from the literature
* `05_spatial_seir.R` - code for the metapopulation model
* `06_descrip_compare.R` - descriptive comparisons (results Figs. 1-2) of mobility data sources using matrices created in `01_process_cdr.R`, `03_fb_movement.R`, and `04_gravity_model.R`.

We then used `cluster_metapop.R` to run the functions from `05_spatial_seir.R` on a cluster. Main simulation results figures and supplement figures are found in the following:

* `07_spatial_seir_results.R` - code for Figs. 3-5
* `08_supp_figs.R` - code for Supplemental Figures





