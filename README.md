# Heliconius_forewing_band

Analysis and figure code for


### Perfect mimicry between <i>Heliconius</i> butterflies is constrained by genetics and development

Steven M. Van Belleghem, Paola A. Alicea Roman, Heriberto Carbia Gutierrez, Brian A. Counterman and Riccardo Papa

For questions, please email Steven M. Van Belleghem (vanbelleghemsteven@hotmail.com) and/or Paola A. Alicea Roman (paola.alicea4@upr.edu)


### Requirements:
````
# Requirements to run patternize analysis and plotting (Best install patternize through github):
install.packages("devtools")
library(devtools)
install_github("StevenVB12/patternize")
# Morpho has committed a change that affects computeTransform used in
# the patternize functions patLanRGB and patLanK (not available in CRAN)
install_github("zarquon42b/Morpho")

# Requirements for landmark analysis and mutant analysis:
install.packages("Momocs", "raster", "viridis")

# Requirements to plot distributions:
install.packages("ggplot2", "maps", "mapdata", "plyr", "alphahull", "G1DBN", "geosphere", "RColorBrewer")
````

### Analysis and plotting scripts:

<i>plot_landmarks.R</i> - Script for analysis of landmarks <i>H. erato/H. melpomene</i>.

<i>Plot_distributions.R</i> - Script to plot distributions of <i>H. erato/H. melpomene</i>.

<i>patternize_erato.R</i> - Script to extract color patterns using all 18 landmarks <i>H. erato</i>.\n
<i>patternize_erato_sub.R</i> - Script to extract color patterns using subset landmarks <i>H. erato</i>.\n
<i>patternize_melpomene.R</i> - Script to extract color patterns using all 18 landmarks <i>H. melpomene</i>.\n
<i>patternize_melpomene_sub.R</i> - Script to extract color patterns using subset landmarks <i>H. melpomene</i>.\n

<i>plot_Heatmaps_comimics.R</i> - Script to plot comparisons color patterns <i>H. erato/H. melpomene</i> using all 18 landmarks.\n
<i>plot_Heatmaps_comimics_sub.R</i> - Script to plot comparisons color patterns <i>H. erato/H. melpomene</i> using subset landmarks.\n
<i>plot_Heatmaps_comimics_sub_male.R</i> - Script to plot comparisons color patterns <i>H. erato/H. melpomene</i> using subset landmarks and only male samples.\n

<i>plot_PCA.R</i> - Script to plot PCA of color patterns <i>H. erato/H. melpomene</i> using all 18 landmarks.\n
<i>plot_PCA_sub.R</i> - Script to plot PCA of color patterns <i>H. erato/H. melpomene</i> using subset landmarks.\n

<i>plot_mutant_comparison.R</i> - Script to compare color patterns <i>H. e. demophoon/H. m. rosina</i> to WntA mutant phenotypes (using subset landmarks).


