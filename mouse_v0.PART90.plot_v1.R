#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Atrial Asymmetry
####  2024-02-22 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#### ResMP MoMP markers -- K. Lavine's Nature Medicine Paper ####
Mac_marker_1 <- list(
    TRMP_gene = c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA',
                  'PDGFC', 'WLS', 'DAB2', 'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3',
                  'EGFL7', 'ECM1', 'SDC3'),
    MoMP_gene = c('CCL17', 'CCR5', 'CCR2', 'CXCL9', 'CXCL3', 'CXCL2', 'CXCL10', 'CSF2RA', 'CCL5', 'CCL20',
                  'TREM1', 'TET2', 'NFKBIE', 'IL1R1', 'IL1R2', 'NFKB1', 'REL', 'MAP2K3', 'IL1A', 'NLRP3',
                  'TRAF1', 'MAPK6', 'SRC', 'IL1B', 'NOD2', 'MAP3K8', 'RELB', 'SOCS3', 'MYD88',
                  'MMP9', 'IL20', 'AREG', 'PTX3', 'IL23A', 'IL10', 'OSM', 'TIMP1', 'EREG', 'IL27')
)
for(i in 1:L(Mac_marker_1)){Mac_marker_1[[i]] <- ConvertGeneSpecies(Mac_marker_1[[i]], 'human', 'mouse')}

#### ResMP MoMP markers -- S. Epelman's Science Immunology Paper ####
genes <- read_excel(paste0(db_dir, 'sciimmunol.abf7777_table_s7.xlsx'),
                    sheet = 3, col_names = T, skip = 1) ## These genes are from the E. Slava's Science Immun Paper
genes <- genes[genes$avg_logFC > 0.75 & genes$p_val_adj < 0.001, ]
Mac_marker_2 <- split(genes$gene, genes$cluster)
for(i in 1:L(Mac_marker_2)){Mac_marker_2[[i]] <- ConvertGeneSpecies(Mac_marker_2[[i]], 'human', 'mouse')}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~ Library Info ####
lib_meta.df <- read.csv(paste0(Docu_dir, 'mouse_sample_meta.csv'))

##~~  Full RNA data with ambiguous cells ####
full_amb.srt <- readRDS('integrated/PART19.consolidated.srt.rds')
Idents(full_amb.srt) <- 'Cell_type'
full_amb.srt$Cell_state <- revalue(full_amb.srt$Cell_state, replace = c(
        'FB1' = 'FB-L',
        'FB2' = 'FB-R',
        'EpiC1' = 'EpiC-L',
        'EpiC2' = 'EpiC-R',
        'EndoC1' = 'EndoC-L',
        'EndoC2' = 'EndoC-R',
        'MP1' = 'MP-R',
        'MP2' = 'MP-L')
)

Table(full_amb.srt$Cell_state, full_amb.srt$Cell_type)

##~~  Clean RNA data without ambiguous cells ####
Table(full_amb.srt$Cell_state, full_amb.srt$Non_ambiguous)
full.srt <- DropMetaLevels(full_amb.srt[, full_amb.srt$Non_ambiguous])
Table(full.srt$Cell_state, full.srt$Cell_type)

##~~  Markers and DEGs ####
deg <- readRDS('analysis/PART20.markers.srt_misc.rds')

##~~  Clean Teichmann Human Heart ####
human.srt <- readRDS('integrated/PART90.human_la_ra.srt.rds')
human.srt$age <- as.numeric(str_split(human.srt$age_group, '-', simplify = T)[, 1])
human.srt$sex <- human.srt$gender
human.srt$dataset <- '2020_Nature_STeichmann'
human.srt$Cell_type <- as.vector(human.srt$cell_type)
human.srt$Cell_state <- as.vector(human.srt$cell_states)


##~~  ST Data - Whole heart ####
st.srt <- readRDS('integrated/PART90.whole_heart_st_for_plotting.srt.rds')
st.srt$Zone2 <- factor(st.srt$Zone2, levels = c('LA', 'RA', 'LV', 'RV', 'IVS', 'Other'))
pt.size.wh <- 1

##~~  Global Milo Result ####
milo_cont <- readRDS('analysis/PART90.global_ra_vs_la.milo_result.rds')


##~~  Cell Type CellChat ####
full_cell_type.list <- readRDS('analysis/PART90.full_cell_type.cellchat_list.rds')
ra_cell_type.cch <- full_cell_type.list[[1]]
la_cell_type.cch <- full_cell_type.list[[2]]


##~~  CCR2KO mouse CD45+ scRNA-seq  ####
ccr2null.srt <- readRDS('integrated/PART30.ccr2null_immune_data_annotated_clean.srt.rds')
ccr2null.srt <- DropMetaLevels(ccr2null.srt[, ccr2null.srt$Group2 %in% c('WT', 'Ccr2Null')])
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  For Publications 2025-03  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~ Plot 1 Milo - Mouse RA vs LA Cell Abundance ####
cutoff <- 4
milo_cont[[2]]$logFC[milo_cont[[2]]$logFC > cutoff] <- cutoff
milo_cont[[2]]$logFC[milo_cont[[2]]$logFC < -cutoff] <- -cutoff

p1.1 <- plotNhoodGraphDA.my(milo_cont[[1]], milo_res = milo_cont[[2]],
                             alpha = 0.05, edge_alpha = 0.1, edge_color = 'grey90') +
        labs(title = 'Control RA vs LA') +
        scale_fill_distiller(palette = 'RdBu', limits = c(-cutoff, cutoff)) +
        theme(aspect.ratio = 1,
              axis.line = element_line(color = 'black'))
p1.2 <- DimPlot2(tmp.srt, group.by = 'Cell_state', cols = Color_cell_state,
                  pt.size = 0.1, alpha = 0.75, reduction = 'clean_RNA_umap')
p1.1 + p1.2
PlotPDF('1.1.milo.ra_vs_la', 10, 5)
p1.1 + p1.2
dev.off()


##~~ Plot 2 Bar - RA Cell state go enrichment ####
FB_R <- deg$Cell_state_rna_marker$FB$gene[deg$Cell_state_rna_marker$FB$cluster == 'FB2' &
                                                  deg$Cell_state_rna_marker$FB$avg_log2FC > 0.5]
EpiC_R <- deg$Cell_state_rna_marker$EpiC$gene[deg$Cell_state_rna_marker$EpiC$cluster == 'EpiC2' &
                                                      deg$Cell_state_rna_marker$EpiC$avg_log2FC > 0.5]
EndoC_R <- deg$Cell_state_rna_marker$EC$gene[deg$Cell_state_rna_marker$EC$cluster == 'EndoC2' &
                                                      deg$Cell_state_rna_marker$EC$avg_log2FC > 0.5]
M1_R <- deg$Cell_state_rna_marker$Myeloid$gene[deg$Cell_state_rna_marker$Myeloid$cluster == 'MP1' &
                                                     deg$Cell_state_rna_marker$Myeloid$avg_log2FC > 0.5]
gl <- list(Reduce(intersect, list(FB_R, EpiC_R, EndoC_R)),
           Reduce(intersect, list(FB_R, EpiC_R, EndoC_R, M1_R)))
common <- ModuleEnrichment(gl, human_or_mouse = 'mouse')
Reactome <- common$Reactome[[1]]
Reactome <- Reactome[order(Reactome$p.adjust), ]
Reactome$Description <- factor(Reactome$Description, levels = Reactome$Description)
GO <- common$GO[[1]]
GO <- GO[order(GO$p.adjust), ]
GO$Description <- factor(GO$Description, levels = GO$Description)

p2.1 <- ggplot(Reactome[Reactome$Description %in% c(
    'Respiratory electron transport',
    'Cellular responses to stimuli',
    'Cellular response to chemical stress',
    'Antigen processing-Cross presentation',
    'The NLRP3 inflammasome'), ],
    aes(x = -log10(p.adjust), y = Description))+
    geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
    scale_y_discrete(limit = rev) +
    labs(x = '-Log10 p value', y = 'Enriched Reactome Pathways') +
    theme_Publication(aspect.ratio = 0.5)
p2.2 <- ggplot(GO[GO$Description %in% c(
    'response to oxidative stress',
    'antigen processing and presentation',
    'positive regulation of inflammatory response',
    'regulation of innate immune response',
    'myeloid leukocyte migration'), ],
    aes(x = -log10(p.adjust), y = Description))+
    geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
    scale_y_discrete(limit = rev) +
    labs(x = '-Log10 p value', y = 'Enriched GO Biological Processes') +
    theme_Publication(aspect.ratio = 0.5)
p2.1/p2.2
PlotPDF('1.2.bar.ra_signature_enriched', 5, 5)
p2.1 / p2.2
dev.off()


##~~ Plot 3 Box - Mouse M1 vs M2 scoring ####
gl <- list(
    GOBP_In = read.table(paste0('external/GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE',
                                '_TO_ANTIGENIC_STIMULUS.v2023.2.Hs.grp'), header = T)[,1],
    RC_Anti = read.table(paste0('external/REACTOME_CD163_MEDIATING_AN_ANTI_',
                                'INFLAMMATORY_RESPONSE.v2023.2.Hs.grp'), header = T)[,1]
)

for(i in 1:L(gl)){gl[[i]] <- ConvertGeneSpecies(gl[[i]], 'human', 'mouse')}
tmp.srt <- full.srt[, full.srt$Cell_state %in% c('MP-L', 'MP-R')]
tmp.srt <- AddModuleScore2(tmp.srt, assay = 'RNA',
                           features = c(Mac_marker_combine,
                                        gl,
                                        mhc_genes),
                           names = c('MoMP', 'TRMP', 'GOBP_In', 'RC_Anti', 'MHC-I', 'MHC-II'),
                           return_z = T)

a <- split(tmp.srt$MoMP, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE
a <- split(tmp.srt$TRMP, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE
a <- split(tmp.srt$GOBP_In, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE
a <- split(tmp.srt$RC_Anti, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE
a <- split(tmp.srt$`MHC-I`, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE
a <- split(tmp.srt$`MHC-II`, tmp.srt$Group2)
t.test(a$RA, a$LA)$p.val < 1e-8 ## TRUE

b <- split(tmp.srt$MoMP, tmp.srt$Group1)
any(c(t.test(b$LA_M, b$LA_F)$p.val < 1e-8, t.test(b$RA_M, b$RA_F)$p.val < 1e-8)) ## FALSE
b <- split(tmp.srt$TRMP, tmp.srt$Group1)
all(c(t.test(b$LA_M, b$LA_F)$p.val < 1e-8, t.test(b$RA_M, b$RA_F)$p.val < 1e-8)) ## FALSE

p3.list <- list()
features <- c('MoMP', 'TRMP', 'GOBP_In', 'RC_Anti', 'MHC-I', 'MHC-II')
for(i in 1:6){
    p3.list[[i]] <- BoxPlot2(tmp.srt, feature = features[i], group.by = 'Group2', min = 0.05, max = 0.95)
}
p3 <- wrap_plots(p3.list, ncol = 2) &
    scale_y_continuous(limits = c(-2, 2.4)) &
    theme_Publication(aspect.ratio = 1.5) &
    NoLegend()
p3
PlotPDF('1.3.box.macrophage_signature', 4, 6)
p3
dev.off()


##~~ Plot 4 ST - M1 vs M2 scoring ####
tmp.srt <- AddModuleScore2(st.srt,
                           features = c(Mac_marker_combine['MoMP'], mhc_genes[2]),
                           names = c('M1', 'MHC-II'),
                           assay = 'ST',
                           return_z = T)

p4.1 <- DimPlotST(tmp.srt, group.by = 'Zone2', pt.sizes = pt.size.wh, cols = Color_st, ncol = 1, legend = 1) &
    theme(aspect.ratio = 1) &
    RestoreLegend()

p4.2 <- FeaturePlotST_Dark(tmp.srt,
                   features = c('M1', 'MHC-II'),
                   pt.sizes = pt.size.wh,
                   minvals = rep(0, 2),
                   maxvals = rep(1.5, 2),
                   ncol = 1,
                   asp = 1)
p4.3 <- BoxPlot2(tmp.srt[, tmp.srt$Zone %in% c('LA', 'LV', 'RA', 'RV', 'IVS')], cols = Color_st,
                feature = 'M1', group.by = 'Zone2', min = 0.01, max = 0.99) +
    BoxPlot2(tmp.srt[, tmp.srt$Zone %in% c('LA', 'LV', 'RA', 'RV', 'IVS')], cols = Color_st,
            feature = 'MHC-II', group.by = 'Zone2', min = 0.01, max = 0.99) &
    scale_y_continuous(limits = c(-2.3, 4.3))
p4 <- p4.1 / (p4.2[[1]] | p4.2[[2]]) / p4.3
p4
PlotPDF('1.4.st.macrophage_signature', 8, 8)
p4
dev.off()


##~~ Plot 5 Box - Human M1 vs M2 scoring ####
tmp.srt <- AddModuleScore2(human.srt, features = c(Mac_marker_human), names = c('M1', 'M2'), return_z = T)
tmp.srt <- tmp.srt[, tmp.srt$cell_type %in% c('Myeloid') & tmp.srt$cell_states != 'Mast']
a <- split(tmp.srt$M1, tmp.srt$region)
all(c(t.test(a$RA, a$LA)$p.val < 1e-8, t.test(a$RA, a$LV)$p.val < 1e-8, t.test(a$RA, a$RV)$p.val < 1e-8)) ## TRUE
b <- split(tmp.srt$M2, tmp.srt$region)
all(c(t.test(b$LA, b$RA)$p.val < 1e-8, t.test(b$LA, b$LV)$p.val < 1e-8, t.test(b$LA, b$RV)$p.val < 1e-8)) ## TRUE

p5 <- BoxPlot2(tmp.srt, feature = 'M1', group.by = 'region', min = 0.05, max = 0.95) /
    BoxPlot2(tmp.srt, feature = 'M2', group.by = 'region', min = 0.05, max = 0.95) &
    theme_Publication(aspect.ratio = 1.5)
p5
PlotPDF('1.5.box.macrophage_signature_in_human', 4, 8)
p5
dev.off()


##~~ Plot 6 CellChat - M1 to All cells ####
tmp.srt <- DropMetaLevels(full.srt[, full.srt$Cell_state %in% c('MP-R', 'MP-L')])
Idents(tmp.srt) <- 'Cell_state'
mye_deg <- FindAllMarkers(tmp.srt, only.pos = T)
mye_deg <- mye_deg[mye_deg$p_val_adj < 0.05, ]
intersect(la_cell_type.cch@LR$LRsig$ligand, mye_deg$gene[mye_deg$cluster == 'MP-L'])

r_net <- ra_cell_type.cch@LR$LRsig[ra_cell_type.cch@LR$LRsig$ligand %in% mye_deg$gene[mye_deg$cluster == 'MP-R'], ]
l_net <- la_cell_type.cch@LR$LRsig[la_cell_type.cch@LR$LRsig$ligand %in% mye_deg$gene[mye_deg$cluster == 'MP-L'], ]

PlotPDF('1.6.chord.mp_to_all_signal', 8, 8)
netVisual_chord_gene(ra_cell_type.cch,
                     sources.use = 'Myeloid', targets.use = c('CM', 'FB', 'EpiC', 'EC', 'Mural', 'Lymphoid'),
                     slot.name = "net",
                     signaling = r_net$pathway_name,
                     title.name = 'RA - MP to All Cells',
                     scale = T,
                     link.border = T,
                     big.gap = 35, small.gap = 5,
                     lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30, annotationTrackHeight = 0.05)
netVisual_chord_gene(la_cell_type.cch,
                     sources.use = 'Myeloid', targets.use = c('CM', 'FB', 'EpiC', 'EC', 'Mural', 'Lymphoid'),
                     slot.name = "net",
                     signaling = l_net$pathway_name,
                     title.name = 'LA - MP to All Cells',
                     scale = T,
                     link.border = T,
                     big.gap = 25, small.gap = 3,
                     lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30, annotationTrackHeight = 0.05)
dev.off()


##~~ Plot 8 EP on CCR2-null  - CCR2 null has low PR interval ####
tmp.srt <- DropMetaLevels(ccr2null.srt[, ccr2null.srt$Cell_state %in% c('Mono1', 'Mono2', 'MP1', 'MP2', 'MP3','Nphil')])
tmp.srt <- DropMetaLevels(tmp.srt[, tmp.srt$Group1 %in% c('WT_2_F', 'Ccr2Null_2_F')])
tmp.srt$tmp <- revalue(tmp.srt$Cell_state, replace = c(
    'Mono1' = 'MoMP',
    'Mono2' = 'TRMP',
    'MP1' = 'TRMP',
    'MP2' = 'TRMP',
    'MP3' = 'TRMP'
))
tmp.srt <- DownsampleByMeta2(tmp.srt, meta_var = 'Group2')
tmp.srt <- RunUMAP(tmp.srt, dims = 1:30, reduction = 'mnn', min.dist = 2, n.neighbors = 100)
FeaturePlot2(tmp.srt, features = c('Ccr2', 'mRFP1'), split.by = 'Group2', 
             min.cutoff = 0.5, max.cutoff = 2, reduction = 'umap') &
    scale_color_distiller(palette = 'RdYlBu', limits = c(0.5, 2))
p0 <- DimPlot2(tmp.srt, group.by = 'tmp', cols = mycol_20s,  reduction = 'umap', pt.size = 2)
p0


data <- DropMetaLevels(tmp.srt[, tmp.srt$tmp %in% c('MoMP', 'TRMP', 'Nphil')])
total <- Table(data$Group1)
data <- as.data.frame(Table(data$tmp, data$Group1))
data$Frac <- data$Freq/total[data$Var2]
p1 <- ggplot(data) +
    geom_bar(aes(x = Var2, y = Frac, fill = Var2), stat = 'identity') +
    scale_fill_manual(values = c('grey50', 'red3')) +
    facet_grid(~Var1, scales = 'free_y') +
    labs(x = '', y = 'Fraction of myeloid cells', fill = '') +
    theme_Publication(aspect.ratio = 1/2) +
    NoLegend() 
p1

pr <- read_csv('../../pacing/PR.csv', )
t.test(pr$`Group A (Control)`, pr$`Group B (CCR2^KO)`)$p.value

p2 <- ggplot(melt(pr), aes(x = variable, y = value)) +
    geom_boxplot(outliers = F) +
    geom_beeswarm(aes(color = variable), cex = 3) +
    scale_color_manual(values = c('grey30', 'red3')) +
    scale_y_continuous(limits = c(32.5, 43.5), breaks = seq(31, 43, 2)) +
    labs(x = '', y = 'PR interval (ms)', color = '', caption = 'p < 0.021') +
    theme_Publication(aspect.ratio = 1) +
    NoLegend()
p2
PlotPDF('1.8.ccr2null_scrna.pacing_pr_interval', 8, 4)
p0 + p1 + p2
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----

