# roessler_netzl_et_al2022
Code repository for "BA.2 and BA.5 omicron differ immunologically from both BA.1 omicron and pre-omicron variants" by Roessler, A., Netzl, A., et al. 2022

The raw data for figure 1, 2 and S1 is in data/titer_data/221121_Raw data.xlsx.

The raw titer data for making the maps is data/titer_data/ in the excel sheets 220309_Omicron I_raw data.xlsx, 220309_Omicron II_raw data.xlsx, 220308_further data for Omicron III_raw data.xlsx. Running the script code/excel_to_titertable.R will create the titer_table.csv in data/titer_data used for map creation. The script code/make_map.R creates the initial map, the code/reactivity_adjusted_maps.R creates the manuscript map with P.1.1 reactivity adjusted by -1 2-fold.

Running the ablandscape_fit.R followed by ablandscape_plot.R in the code/ directory creates the antibody landscapes in the figures/ directory. The labelled map figure is generated by running the figures/labelled_map/plot_labelled_map.R script.

All code for SOM diagnostics are in the som/ directory, additional analyses made upon revision are in the revision/ directory.


All code was run under R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics":
R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing,
  Vienna, Austria. URL https://www.R-project.org/.
  
  
Antigenic maps were constructed using the Racmacs package, Version 1.1.35:
Wilks S (2022). _Racmacs: R Antigenic Cartography Macros_. https://acorg.github.io/Racmacs,
  https://github.com/acorg/Racmacs.
  
  
Antibody landscapes were constructed using the ablandscapes package, Version 1.1.0: https://github.com/acorg/ablandscapes





