

gpsfb = subset_taxa(GPfr, Phylum == "Bacteroidetes")
gpsfbg = tax_glom(gpsfb, "Family")
plot_tree(gpsfbg, color="SampleType", shape="Class", size="abundance")