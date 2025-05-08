
Aitchison distance matrices for one-way permutational analyses of variance (PERMANOVA) with “sampling site” as the factor


Rarefaction wins?

https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.708716/full


Compositional nature of metabarcoding 

Aitchison distance is a compositional euclidean distance

The Aitchison distance matrices were used to conduct one-way permutational analyses of variance (PERMANOVA) with “sampling site” as the factor, using vegan v. 2.5-7 (999 permutations; Oksanen et al., 2020). Post hoc pairwise comparisons (999 permutations) including a Benjamini-Hochberg correction (Benjamini and Hochberg, 1995) were performed using RVAideMemoire v. 0.9–78 (Hervé, 2020). Tests for the homogeneity of multivariate dispersions (PERMDISP; Anderson, 2006) were performed with vegan (999 permutations; Oksanen et al., 2020).

To explore correlations between community structure and sediment, sediment porewater and water column geochemical variables (Table 3), linear dependencies between variables were first identified by a principal component analysis (PCA) using centered and scaled data (Figure 2), and by inspecting variance inflation factors (VIFs). Two sets of up to five explanatory variables (Table 4) were then selected for distance-based redundancy analyses (db-RDA; Legendre and Anderson, 1999; Ramette, 2007) performed with vegan (Oksanen et al., 2020). Distance-based redundancy analyses were selected as the ordination method based on the inspection of primary axis lengths using detrended correspondence analyses (Lepš and Šmilauer, 2003). Variables in the first model (Model 1) were selected due to their minimal multicollinearity (VIFs of ≤6), while the second model (Model 2) was limited to variables similar to those often measured during environmental monitoring surveys (see e.g., the HELCOM Monitoring Manual5) and had VIFs of ≤13.7 (Table 4). Both models were run using CLR-transformed 16S or 18S rRNA gene sequence data, resulting in a total of four model runs. Environmental variables were projected onto the ordination space using the envfit() function in vegan (Oksanen et al., 2020), with fitting performed on linear scores of the ordination axes. Following global significance tests for the ordinations, the variance explained by each variable was compared using permutation tests with the remaining variables included as covariates, using the anova() function in vegan (999 permutations; Oksanen et al., 2020).




https://forum.qiime2.org/t/why-is-there-no-principal-component-and-redundancy-analysis-in-qiime2/23670/2

EICODE solves both of these issues and has the advantage of being a PCA implemented for QIIME 2. Basically, DEICODE takes your data, performs a partial CLR transform, and then uses mathematical magic (sparse matrix closure) to solve an ordination while dealing with the zeros. The underlying math is beyond me, but when I tried to get it through my linear algebra-less brain, I found the github tutorials super helpful.

However, there is a way to get that pseudo feature ranking. Adonis will provide an effect size ranking for covariates of interest with both continuous and categorical data. You can even adjust for other factors (as long as you set up the equation correctly.) The current qiime2 adonis implementation does this process one at a time, so you'd have to harvest the data yourself to turn it into a visualization. You can see one example of that adonis work in panel B of figure 1 from He et al, showing their goegraphic effect is 5x any other factor:


" Adonis will provide an effect size ranking for covariates of interest with both continuous and categorical data."


Quarter + Region

or 

Quarter + Nerr 


RDA on 


 
β-diversity distances (e.g. PERMANOVA)