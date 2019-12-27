require(pals)

colors_redblu10 = c("#3360a4","#1a8af9","#06b8ff","#b3defb","#ffffff","#fdc82f","#f67028","#e9302a","firebrick3","firebrick4")

colors_hmm6=c("Enh"="gold","Others"="mediumaquamarine","Quies"="bisque","ReprPC"="slategray2","Tss"="firebrick1","Tx"="green4")
colors_hmm18 = c("TssA"="#ff0000","TssFlnk"="#ff4500","TssFlnkU"="#ff4500","TssFlnkD"="#ff4500","Tx"="#008000","TxWk"="#006400","EnhG1"="#c2e105","EnhG2"="#c2e105","EnhA1"="#ffc34d","EnhA2"="#ffc34d","EnhWk"="#ffff00","ZNF_Rpts"="#66cdaa","Het"="#8a91d0","TssBiv"="#cd5c5c","EnhBiv"="#bdb76b","ReprPC"="#808080","ReprPCWk"="#c0c0c0","Quies"="bisque")

colors_ct <- setNames(alphabet2()[c("caramel","orange","zinnia","jade","turquoise","amethyst","navy","violet","ultraviolet","lavender","lavender","pink","kingcrab","turquoise","blue","lavender","damson","pink","amethyst","orange","lavender","orange","zinnia","turquoise","lavender","jade","blue")],c("PBMC","CD14","CD56","CD19","CD3","CD4T","CD8T","CD4_naive","CD4_memory","Th","Th1","Th17","Treg","CD8_naive","CD8_memory","naive_T","memory_T","effector_T","Lymphoid","Myeloid","Adaptive","Monocytes","NK_cells","Tcells","Naive_Tcells","Bcells","acCD8_Tcells"))
colors_ct <- c(colors_ct,"Innate"="#F9DB60","B_naive"="#a4e5b8","naive_B"="#a4e5b8","B_memory_nonclass_switched"="#138537","B_memory_class_switched"="#169880","B_memory"="#169880","memory_B"="#169880","CD15"="#e39012","HSC"="#FBEF7E","CD34"="#FBEF7E","CD34_mobilized"="#FBDD7E","Erythroblasts"="tomato","Erythrocytes"="tomato","Macrophages_alt_activated"="#da3835","Macrophages_inflammatory"="#e0566a","Macrophages"="#e0566a","Megakaryocytes"="firebrick2","Plasma_cells"="#c3f960","DCs"="#8dc53f","HSCs"="#FBEF7E","pDCs"="#8ee3bf","CD16_monocytes"="#e0566a","Treg_naive"="#efbe8e","Th2"="#ca66fc","Tfh"="#9c05e8")
colors_ct <- colors_ct[order(names(colors_ct))]

shape_sex <- c("Females"=16,"F"=16,"Males"=17,"M"=17)
colors_sex <- c("Females"="#ED6A6F","F"="#ED6A6F","Males"="#829CBD","M"="#829CBD")
colors_sex_darken1 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"black"))(10)[3])
colors_sex_darken2 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"black"))(10)[6])
colors_sex_lighten1 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"white"))(10)[3])
colors_sex_lighten2 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"white"))(10)[6])

colors_age <- c("HO"="darkorange","Age3"="darkorange","Opening"="darkorange","Opening peaks"="darkorange","Opening peaks in"="darkorange","Opening peaks with aging"="darkorange","Opening with aging"="darkorange","Up"="darkorange","Up with aging"="darkorange","Upregulated"="darkorange","Upregulated with aging"="darkorange",
                "HM"="#49c06c","Age2"="#49c06c",
                "HY"="darkmagenta","Age1"="darkmagenta","Closing"="darkmagenta","Closing peaks"="darkmagenta","Closing peaks in"="darkmagenta","Closing peaks with aging"="darkmagenta","Closing with aging"="darkmagenta","Down"="darkmagenta","Down with aging"="darkmagenta","Downregulated"="darkmagenta","Downregulated with aging"="darkmagenta")
colors_age_darken1 <- sapply(colors_age,function(h) colorRampPalette(c(h,"black"))(10)[3])
colors_age_darken2 <- sapply(colors_age,function(h) colorRampPalette(c(h,"black"))(10)[6])
colors_age_lighten1 <- sapply(colors_age,function(h) colorRampPalette(c(h,"white"))(10)[3])
colors_age_lighten2 <- sapply(colors_age,function(h) colorRampPalette(c(h,"white"))(10)[6])

lighter <- function(clr,step=6) {
  colorRampPalette(c(clr,"white"))(10)[step]
}
darker <- function(clr,step=3) {
  colorRampPalette(c(clr,"black"))(10)[step]
}