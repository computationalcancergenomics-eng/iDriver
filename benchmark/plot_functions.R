library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

AUPR_LinePlot <- function(preprocessed_pRes, n, ann_PCAWG_ID, element,
                          tissue, based_on, include_legend = F, 
                          selected_methods_plots = NULL,
                          include_lable = F, font_size = c(26, 26, 26),
                          by_nTop = 100){
  
  METHODs <- names(preprocessed_pRes)
  sigCriteria <- define_significant_criteria("fixedNumberOfElems", n)
  
  AUPR_tab <- prepare_AUPR(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on)
  n <- min(nrow(preprocessed_pRes[[1]]), n)
  is_TP <- AUPR_tab$label
  
  
  #if(!is.na(selected_methods_plots)){
   # is_TP <- is_TP[which(names(is_TP) %in% selected_methods_plots)]
  #} 
  recall_fraction <- lapply(is_TP, function(s){
    X=c()
    for (i in 1:n) {
      x <- sum(s[1:i])
      X=c(X,x)
    }
    X
  })
  
  x <- data.frame(t(do.call(cbind, recall_fraction)))
  colnames(x) <- paste0("Top ", 1:n)
  x$method <- names(recall_fraction)
  
  dat.m <- reshape2::melt(x,id.vars="method",
                          measure.vars=paste0("Top ", 1:n))
  
  
  colnames(dat.m) <- c("method", "top_n", "#True Positives")
  dat.m$method <- ifelse(dat.m$method == 'oncodriveFML_cadd', 'oncodriveFML (CADD)', dat.m$method)
  dat.m$method <- ifelse(dat.m$method == 'ncDriver_combined', 'ncDriver', dat.m$method)
  dat.m$method <- ifelse(dat.m$method == 'oncodriveFML_vest3', 'oncodriveFML (vest3)', dat.m$method)
  
  # p <- ggplot(dat.m, aes(x=top_n, y=`#True Positives`, group=method))+
  #   geom_line(aes(color=method),linewidth=1) +
  #   # scale_color_manual(values=create_method_colours(METHODs))
  #   scale_color_manual(values=create_method_colours())
  # p + scale_x_discrete(breaks = paste0("Top ", seq(0, n, by = 100)))+
  #   theme_bw()+
  #   # labs(title = paste0(element, " measures in ", tissue))+
  #   theme(axis.text.x = element_text(angle = 45,hjust = 1),
  #         panel.grid.major = element_blank(),
  #         axis.title.x = element_blank(),
  #         panel.grid.minor = element_blank())+
  #   theme(legend.text = element_text(size = 12),
  #         
  #         axis.title.y = element_text(size = 14))
  
  p <- ggplot(dat.m, aes(x = top_n, y = `#True Positives`, group = method)) +
    geom_line(aes(color = method), linewidth = 1) +
    scale_color_manual(values = create_method_colours()) +
    scale_x_discrete(breaks = paste0("Top ", seq(0, n, by = by_nTop))) 
    
                                       # Remove legend
  if (TRUE) { #(multiPanel == T
    p <- p+ theme_bw() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = font_size[1]),  # x-axis text size
      axis.text.y = element_text(size = font_size[2]),                         # y-axis text size
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),                       # y-axis title size
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  } else {
    p <- p+ theme(axis.text.x = element_text(angle = 45, hjust = 1, size = font_size[1]),
                  axis.text.y = element_text(size = font_size[2]),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),           # remove box border
                  panel.background = element_blank(),       # ensure no background override
                  axis.line = element_line(colour = "black")
                  
        )          
  }
  if (include_legend) {
    p <- p+ guides(color = guide_legend(title = NULL))  # Show legend without title
  } else {
    p <- p + theme(legend.position = "none" )     
  }
  
  if (include_lable) {
    p <- p + labs(title = paste0( tissue))+
      theme(plot.title = element_text(size = 20))
  }
  
  
}



generate_paths <- function(cohorts, base_paths) {
  
  
  all_paths <- lapply(cohorts, function(cohort) {
    gsub("\\{cohort\\}", cohort, base_paths)
  })
  names(all_paths) <- cohorts
  return(all_paths)
}

generate_singlepanel_AUPR_plots <- function(preprocessed_pRes, n, ann_PCAWG_ID, element,
                                            path_save_plot, tissue, based_on,
                                            include_legend = F,
                                            selected_methods_plots = NULL){
  
  p = AUPR_LinePlot(preprocessed_pRes, n, ann_PCAWG_ID, element,
                                tissue, based_on, include_legend, selected_methods_plots)
    
    
    dir.create(paste0(path_save_plot, '_AUPR_linePlots/'),
               showWarnings = F,
               recursive = T)
  
  ggsave(paste0(element,"_" , tissue, "_", based_on, "_", n, "_nTPs_Plot.png"), device = "png",
         width = 10,
         #height = 9,
         path = paste0(path_save_plot, '_AUPR_linePlots/'))
}


# Function to generate multi-panel plots
generate_multipanel_AUPR_plots <- function(all_paths, cohorts, n, ann_PCAWG_ID, newRESULTS,
                                           path_save_plot, based_on, selected_methods_plots,
                                           base_method = 'iDriver',
                                           include_legend = F,
                                           selected_methods_tables = NULL, rows = 4) {
  
  plot_list <- list()
  
  # Generate individual plots for each cohort
  for (tissue in cohorts) {
    print(tissue)
    
    cohort <- ifelse(tissue == 'Pan_Cancer', 'Pancan', tissue)
    path_pRes <- paste0('../extdata/procInput/PCAWG_results/', cohort, '.RData')
    PATHs_newResults <- all_paths[[tissue]]
    
    
    
    dir.create(paste0(path_save_plot), showWarnings = F, recursive = T)
    
    tmp_plot_list <- list()
    for (element in c('CDS', 'NC')) {
      preprocessed_pRes <- preprocess_pRes(path_pRes, PATHs_newResults, newRESULTS,
                                           element, base_method)
      
      plot <- AUPR_LinePlot(
        preprocessed_pRes = preprocessed_pRes,
        n = n,
        ann_PCAWG_ID = ann_PCAWG_ID,
        element = element,
        tissue = tissue,
        based_on = based_on,
        include_legend = F,
        selected_methods_plots = selected_methods_plots,
        include_lable = T, font_size = c(15, 15, 18),
        by_nTop = 20
      )
      tmp_plot_list[[paste0(element, '_', tissue)]] <- plot
    }
    plot_list <- append(plot_list, tmp_plot_list)
  }
  
  # Define layout
  ncol <- 6
  nrow <- 5
  plots_per_page <- ncol * nrow
  
  # Loop over plot chunks and save as PNGs
  for (i in seq(1, length(plot_list), by = plots_per_page)) {
    page_num <- ceiling(i / plots_per_page)
    png(file = file.path(path_save_plot, 
                         paste0("AUPR_MultiPanel_", n, "_page", page_num, ".png")),
        width = 17.5, height = 19, units = "in", res = 350)  # adjust size/res
    
    grid_plot <- do.call(grid.arrange, 
                         c(plot_list[i:min(i + plots_per_page - 1, length(plot_list))], 
                           ncol = ncol, nrow = nrow))
    
    dev.off()  # Close PNG device after each page
  }
  
  return(invisible(plot_list))
}


save_AUPR_singleCohort <- function(n, PATH_newResults, newRESULT,
                                   element, path_save_plot, based_on, tissue,
                                   selected_methods_plots = NA,
                                   base_method = 'iDriver',
                                   path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                                   include_legend = F,
                                   selected_methods_tables = NA
){
  
  if (tissue == 'Pan_Cancer') {
    cohort = 'Pancan'
  } else {
    cohort = tissue
  }
  
  path_pRes <- paste0('../extdata/procInput/PCAWG_results/', cohort, '.RData')
  
  preprocessed_pRes <- preprocess_pRes(path_pRes, PATH_newResults, newRESULT,
                                       element, base_method, selected_methods = selected_methods_tables)
  
  
  dir.create(paste0(path_save_plot),
             showWarnings = F,
             recursive = T)
  
  
  
  if(!is.na(selected_methods_plots)){
    if (!is.null(selected_methods_plots)) {
      selected_methods_plots <- c(selected_methods_plots, newRESULT)
    }
  }
  
  ann_PCAWG_ID <- fread(path_ann_PCAWG_ID)
  
  
  generate_singlepanel_AUPR_plots(preprocessed_pRes, n, ann_PCAWG_ID, element,
                                  path_save_plot, tissue, based_on, include_legend,
                                  selected_methods_plots)
  
}

save_AUPR_allCohorts <- function(n, PATHs_newResults_allCohorts, newRESULTS,
                                 path_save_plot, based_on,
                                 selected_methods_plots, cohorts, 
                                 base_method = 'iDriver',
                                 path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement.csv',
                                 selected_methods_tables = NULL
){
  
  # generate the paths from a pattern of paths
  all_paths <- generate_paths(cohorts, PATHs_newResults_allCohorts)
  
  # if(!is.na(selected_methods_plots)){
  #   if (!is.null(selected_methods_plots)) {
      # selected_methods_plots <- c(selected_methods_plots, newRESULTS)
  #   }
  # }
  
  selected_methods_plots <- newRESULTS
  ann_PCAWG_ID <- fread(path_ann_PCAWG_ID)
  
  
  multi_panel_plot <- generate_multipanel_AUPR_plots(
    all_paths = all_paths,
    cohorts,
    n,
    ann_PCAWG_ID,
    newRESULTS,
    path_save_plot,
    based_on,
    selected_methods_plots = selected_methods_plots,
    base_method,
    selected_methods_tables
  )
  
}

create_method_colours <- function(){
  
  METHODs <- c("ActiveDriverWGS", "dNdScv", "iDriver", 
               "ExInAtor", "LARVA", "MutSig", 
               "NBR", "ncdDetect", "ncDriver",
               "oncodriveFML (CADD)", "oncodriveFML (vest3)",
               "regDriver", "DriverPower",
               "Population-level",
               "Only-recurrence",
               "iDriver (no-FI 'Beta')",
               "Tippett", "Fisher",
               'Synonymous', "Non-synonymous") 
  
  
  all_colors <- c('#b2df8a', '#247740', 'darkred',
                  'grey80',  '#a1cbd2', '#E2A9B2',
                  "darkblue", '#d7c154', '#e7d998', 
                  '#653780', '#cab2d6', '#F89F5B',
                  '#317d89',  
                  'grey','darkblue','skyblue', 'gold',  'black', '#fde0dd', '#54278f'
  )
  
  
  METHODs <- factor(METHODs, levels = unique(METHODs))
  method_colors <- all_colors[1:length(METHODs)]
  names(method_colors) <- levels(METHODs)
  
  return(method_colors)
}

save_legend_methodColors <- function(path_save){
  # Generate colors
  method_colors <- create_method_colours()
  
  # Dummy data for legend generation
  df <- data.frame(METHOD = factor(names(method_colors), levels = names(method_colors)),
                   x = 1, y = seq_along(method_colors))
  
  # Dummy plot
  p <- ggplot(df, aes(x, y, color = METHOD)) +
    geom_point(size = 4) +
    scale_color_manual(values = method_colors, name = NULL) +  # Remove legend title
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14)  # Adjust font size here
    )
  
  # Extract legend
  legend <- cowplot::get_legend(p)
  
  # Save legend to file
  ggsave(path_save, legend, width = 4, height = 6)
  
}



extract_topK_methods <- function(cancer_types, measure, k, 
                                 path_save_tables, element_type){
  
  # Loop through each cancer type and process the corresponding table
  all_data <- list()
  for (cancer in cancer_types) {
    
    file_path <- gsub("\\{cancer\\}", cancer, path_save_tables)
    file_path <- gsub("\\{elem\\}", element_type, file_path)
    
    df <- read.csv(file_path)
    
    df$method <- ifelse(df$method == 'oncodriveFML_cadd', 'oncodriveFML (CADD)', df$method)
    df$method <- ifelse(df$method == 'ncDriver_combined', 'ncDriver', df$method)
    df$method <- ifelse(df$method == 'oncodriveFML_vest3', 'oncodriveFML (vest3)', df$method)
    
    # # Extract top-5 methods based on the dynamic measure
    # top5_methods <- df %>%
    #   arrange(desc(!!sym(measure))) %>%
    #   head(k) %>%
    #   select(method, all_of(measure))  # Use the dynamic measure here
    
    top5_methods <- df %>%
      arrange(desc(!!sym(measure))) %>%
      head(k) %>%
      select(method, all_of(measure))
    
    
    # Add a column for cancer type to distinguish later
    top5_methods$cancer_type <- cancer
    
    # Append to the list
    all_data[[cancer]] <- top5_methods
  }
  
  # Combine all data frames into one
  final_df <- bind_rows(all_data)
  
  # Sort the data frame in decreasing order of the dynamic measure within each cancer type
  final_df <- final_df %>%
    group_by(cancer_type) %>%
    arrange(desc(!!sym(measure)), .by_group = TRUE)
  
  final_df
}

create_methods_tables_allCohorts <- function(hyperCancers, cancer_types, 
                                             path_save_tables_withHypers_eMETrun,
                                             path_save_tables_rmHypers, measure, k,
                                             element_type){
  
  if (length(hyperCancers) == 0) {
    
    methods <- extract_topK_methods(cancer_types, measure, k,
                                    path_save_tables_rmHypers, element_type)
    
  } else {
    k_methods_withHypers_eMETrun <- extract_topK_methods(hyperCancers, measure, k,
                                                         path_save_tables_withHypers_eMETrun,
                                                         element_type)
    
    k_methods_woHypers <- extract_topK_methods(cancer_types, measure, k,
                                               path_save_tables_rmHypers, element_type)
    
    
    methods <- rbind(k_methods_withHypers_eMETrun, k_methods_woHypers)
    
  }
  
  methods_tables <- data.frame(table(methods$method))
  
  colnames(methods_tables) <- c('method', 'Freq')
  
  methods_tables <- methods_tables[order(methods_tables$Freq, decreasing = T),]
  
  list('Table' = methods_tables, 'allCancers' = methods)
}

save_grouped_barPlotAllCohorts <- function(save_name, path_save, 
                                           hyperCancers, cancer_types,
                                           path_save_tables_withHypers_eMETrun,
                                           path_save_tables_rmHypers, measure, k,
                                           element_type){
  
  dir.create(path_save, recursive = T, showWarnings = F)
  
  
  LIST = create_methods_tables_allCohorts(hyperCancers, cancer_types,
                                          path_save_tables_withHypers_eMETrun,
                                          path_save_tables_rmHypers, measure, k,
                                          element_type)
  
  data <- LIST$allCancers
  data$cancer_type <- ifelse(data$cancer_type == 'Pan_Cancer', 'Pancancer', data$cancer_type)
  data$cancer_type <- ifelse(data$cancer_type == 'Pancan-no-skin-melanoma-lymph', 'Pancancer*', data$cancer_type)
  
  # Sort the data frame in decreasing order of the dynamic measure within each cancer type
  final_df <- data %>%
    group_by(cancer_type) %>%
    arrange(desc(!!sym(measure)), .by_group = TRUE)
  
  # Create custom color mapping using the create_method_colours function
  method_colors <- create_method_colours()
  
  # Plot the grouped bar plot with decreasing order
  ggplot(final_df, aes(x = reorder(cancer_type, -!!sym(measure)), y = !!sym(measure), fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = .6), width = 0.6) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 21),
          panel.grid = element_blank(), 
          legend.title = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black")) +
    labs(#title = "Top-5 Methods Based on F1 Score Across Cancer Types", 
      x = "",
      y = paste(measure, "")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = method_colors)  # Use custom colors
  
  
  
  ggsave(paste0(measure, "grouped_barPlot_", element_type, '_', save_name,".png"),
         device = "png", width = 14, height = 6,
         bg = 'white',
         path = path_save)
  
  
}


save_kTop_methods <- function(save_name, path_save, 
                              hyperCancers, cancer_types,
                              path_save_tables_withHypers_eMETrun,
                              path_save_tables_rmHypers, measure, k_top,
                              element_type){
  
  save_name <- paste0(measure, save_name, '_', '_ktopMethods', element_type)
  dir.create(path_save, recursive = T, showWarnings = F)
  
  plot_list <- list()
  
  for (k in 1:k_top) {
    LIST = create_methods_tables_allCohorts(hyperCancers, cancer_types, 
                                            path_save_tables_withHypers_eMETrun, 
                                            path_save_tables_rmHypers, measure, k, 
                                            element_type)
    data <- LIST$Table
    data$Freq <- as.integer(data$Freq)
    print(data)
    
    # Get colors for methods in the dataset
    method_colors <- create_method_colours()
    
    # Base plot
    p <- ggplot(data, aes(x = Freq, y = reorder(method, Freq), fill = method)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = method_colors) +
      scale_x_continuous(breaks = scales::pretty_breaks())+
      labs(x = "Frequency", y = "Method", fill = "Method") +
      theme_minimal() +
      labs(
        title = paste0("k top methods = ", k),
        x = NULL,
        y = NULL
      )+
      theme(text = element_text(size = 16), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.title = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none")  # Hide legend in main plots
    
    plot_list[[k]] <- p
  }
  
  
  # Arrange plots without legend
  ncol <- 2
  nrow <- 2
  plots_per_page <- ncol * nrow
  
  plot_png_path <- file.path(path_save, paste0(save_name, ".png"))
  png(plot_png_path, width = 10, height = 4, units = "in", res = 500)
  
  for (i in seq(1, length(plot_list), by = plots_per_page)) {
    grid_plot <- do.call(grid.arrange, c(plot_list[i:min(i + plots_per_page - 1, length(plot_list))], 
                                         ncol = ncol, nrow = nrow))
  }
  
  dev.off()
}



extract_cohorts_tableDir <- function(path_save_tables){
  
  tables_dir <- dirname(path_save_tables)
  
  unique(unlist(lapply(strsplit(list.files(tables_dir), '_'),
                       function(s){
                         s[1]
                       })))
}


exreact_saved_hypr_nonhypers <- function(path_save_tables_nonHypers, path_save_tables_Hypers){
  
  
  nonHyperCancers <- c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph', "Liver-HCC",
                       "Breast-AdenoCa", "Panc-AdenoCA", "Kidney-RCC", "Ovary-AdenoCA",
                       "Prost-AdenoCA", "Bladder-TCC",
                       "Panc-Endocrine",  "CNS-Medullo", "Bone-Leiomyo", 
                       "Lymph-CLL", "Bone-Osteosarc", "Cervix-SCC",
                       "Breast-LobularCa", "Kidney-ChRCC", "Thy-AdenoCA",
                       "CNS-Oligo", "Myeloid-MPN",
                       "Myeloid-AML", "Bone-Epith", "CNS-PiloAstro", "Cervix-AdenoCA",
                       "Breast-DCIS", "Bone-Cart"
  )
  
  hyperCancers <- c(
    'Skin-Melanoma',  'Stomach-AdenoCA',   'Uterus-AdenoCA', # 'Lymph-BNHL' 
    'ColoRect-AdenoCA',  'Eso-AdenoCa', 'Head-SCC', 'Lung-AdenoCA', 'Lung-SCC',
    'Biliary-AdenoCA',  'CNS-GBM'
  )
  
  
  saved_hypers <- extract_cohorts_tableDir(path_save_tables_Hypers)
  saved_nonHypers <- extract_cohorts_tableDir(path_save_tables_nonHypers)
  
  print(saved_hypers)
  print('****')
  print((saved_nonHypers))
  print('----1-----')
  saved_nonHypers <- saved_nonHypers[which(saved_nonHypers %in% nonHyperCancers)]
  saved_hypers <- saved_hypers[which(saved_hypers %in% hyperCancers)]
  
  print((saved_hypers))
  print('****')
  print((saved_nonHypers))
  print('----2-----')
  list('withHyper' = saved_hypers, 'nonHypers' = saved_nonHypers)
  
  
}

###################### lolipop plots ############################

map_muts <- function(gr, path_to_GEs, element,callable_GE_train, 
                     count_blocks=TRUE, test=NULL, train=NULL){
  
  if (element == "test") {
    
    load(path_to_GEs)
    
    
    testGE_gr <- unlist(testGE)
    
    g <- names(testGE_gr)
    s <- strsplit(g, "[::]")
    GenomicElement <- unlist(lapply(s, function(x){x[1]}))
    
    lenElement <- width(testGE_gr)
    
    mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
    
    callable_GE <- relist(testGE_gr, testGE)
    lenElement_new <- sum(width(callable_GE))
    
  } else if (element == "train") {
    
    #callable_GE <- make_callable_trainGE(path_to_GEs, path_to_callable)
    callable_GE <- callable_GE_train
    lenElement_new <- sum(width(callable_GE))
  }
  
  
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  # if we run this function for each element, length(unique(queryHits(ov))) is nMut and
  # length(unique(subjectHits(ov))) is number of mutated blocks. >>> ov <- findOverlaps(gr, testGE$`gc19_pc.cds::gencode::ARID1A::ENSG00000117713.13`, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_Mutated <- callable_GE[idx_callable_GE]  
  
  if (element == "test") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    GenomicElement=unique(mcols(GE_Mutated, level = "within")[,"GenomicElement"]), 
                                    elemenntLength = sum(mcols(GE_Mutated, level = "within")[,"lenElement"]),
                                    name = names(GE_Mutated))
  } else if (element == "train") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    name = names(GE_Mutated))
  }
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  if (count_blocks) {
    
    GE_blocks <- unlist(callable_GE)
    
    ov_block <- findOverlaps(gr, GE_blocks, ignore.strand = TRUE)
    idx_gr_block <- queryHits(ov_block)
    idx_GE_block <- subjectHits(ov_block)
    
    gr_blocks_mappedGE <- gr[idx_gr_block]
    GE_blocks_mutated <- GE_blocks[idx_GE_block]
    
    mcols(gr_blocks_mappedGE) <- DataFrame(mcols(gr_blocks_mappedGE), 
                                           binID=names(GE_blocks_mutated),
                                           block_ID=paste0(names(GE_blocks_mutated), "**",1:length(gr_blocks_mappedGE)),
                                           block_length=width(GE_blocks_mutated))
    
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new, gr_blocks = gr_blocks_mappedGE)
  } else {
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new)
  }
  
}

generate_mut_df <- function(gr, testGEs, sigHits){
  
  all_elements <- map_muts(gr, path_to_GEs, path_to_callable, element = 'test',
                           callable_GE_train = '', count_blocks = T)
  gr_annotated <- all_elements$MAF_GE_mapped
  
  gr_annotated <- gr_annotated[which(gr_annotated$name %in% sigHits)]
  
  mutation <- as.data.frame(cbind(gr_annotated$name,
                                  as.character(gr_annotated$var_type),
                                  as.numeric(start(gr_annotated))))
  colnames(mutation) <- c("binID", "var_type", "X_position")
  
  mutation
}


lollipop_full_length <- function(gene, mutation, testGE, main = TRUE){
  
  
  GeneMutation <- mutation[mutation$binID == gene, ]
  
  coordinats <- subset(as.data.frame(table(as.character(GeneMutation$var_type),
                                           GeneMutation$X_position)), Freq != 0)
  
  colnames(coordinats) <- c("Variant_type", "X", "freq")
  coordinats$Variant_type <- gsub("SNP", "SNV", coordinats$Variant_type)
  
  y_ylim <- if(max(coordinats$freq) > 5) {
    y_ylim = max(coordinats$freq)
  }else{
    y_ylim = 5}
  
  element_coordinates <- testGE[[gene]]
  L <- sum(width(element_coordinates))
  
  S <- start(element_coordinates) # - (start(element_coordinates)[1])    #exon start
  E <- end(element_coordinates) #- start(element_coordinates)[1]        #exon End
  c <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1")) #pallet
  
  gene_info = define_element_type2(gene)
  
  # Set the font size for axis labels, axis numbers, and main title
  par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  
  # Define the plot
  plot(x = 1,                 
       xlab = "Position of Mutations", 
       ylab = "#Mutations",
       xlim = c(min(as.numeric(as.character(coordinats$X))), 
                max(as.numeric(as.character(coordinats$X)))), 
       ylim = c(-1, y_ylim),
       frame = FALSE,
       type = "n",
       if (main) {
         main = paste0( gene_info[[2]][1], " ", gene_info[[1]][1])
       }
  )
  
  # Add a gray rectangle at the bottom of the plot
  rect(min(as.numeric(as.character(coordinats$X))) - 100, 
       -.3, 
       max(as.numeric(as.character(coordinats$X))),
       0, density = NA, angle = 45, col = "grey", border = NA)
  
  # Add blue rectangles for element coordinates
  for(i in 1:length(element_coordinates)){
    rect(S[i], -.4, E[i], 0.2, density = NA, angle = 45, col = '#386cb0', border = NA)
  }
  
  # Add segments for mutation positions
  for(i in 1:nrow(coordinats)){
    segments(as.numeric(as.character(coordinats$X[i])),
             0, as.numeric(as.character(coordinats$X[i])),
             coordinats$freq[i], lwd = .2)
  }
  
  # Add points for mutations
  for (i in 1:nrow(coordinats)) {
    if (coordinats$Variant_type[i] == 'SNV') {
      point_type = 16
      col = "#4dac26" 
    } else if (coordinats$Variant_type[i] == 'INS' | coordinats$Variant_type[i] == 'DEL') {
      point_type = 17
      col = 'red' 
    } else if (coordinats$Variant_type[i] == 'DNP' | coordinats$Variant_type[i] == 'TNP' | coordinats$Variant_type[i] == 'ONP') {
      point_type = 8
      col = '#984ea3'
    }
    
    points(x = as.numeric(as.character(coordinats$X[i])), 
           y = coordinats$freq[i], 
           pch = point_type, 
           col = col)
  }
  
  nMut <- nrow(GeneMutation)
  
  # Define legend labels, point types, and colors
  legend_labels <- c("SNV", "INS/DEL") #, "DNV/TNV/ONV")
  legend_pch <- c(16, 17, 8)
  legend_colors <- c("#4dac26", "red") #, "#984ea3")
  
  # Add the legend
  legend(x = "topright",
         legend = legend_labels,
         pch = legend_pch,
         col = legend_colors,
         title = paste0("Element length = ", L),
         ncol = 1, horiz = FALSE, cex = 1.2, pt.cex = 1.3)
  
}


define_element_type2 <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElements <- unlist(lapply(strsplit(GenomicElement, "[.]"), function(x){x[length(x)]}))
  gene_name <- unlist(lapply(s, function(x){x[5]}))
  
  GEs <- c()
  for (GenomicElement in GenomicElements) {
    if (GenomicElement == 'enhancers') {
      GenomicElement = 'enhancer'
    } else if (GenomicElement == '3utr') {
      GenomicElement = "3'UTR"
    } else if (GenomicElement == '5utr') {
      GenomicElement = "3'UTR"
    } else if (GenomicElement == 'cds') {
      GenomicElement = 'CDS'
    } else if (GenomicElement == 'promCore') {
      GenomicElement = 'Core Promoter'
    } else if (GenomicElement == 'ss') {
      GenomicElement = 'Splice site'
    } 
    
    GEs <- c(GEs, GenomicElement)
  }
  
  
  list(GEs, gene_name)
  
}


save_lolipops <- function(mutation, testGE, path_out, save_name = ''){
  selected_genes <- unique(mutation$binID)
  dir.create(path_out, showWarnings = F, recursive = T)
  
  info <- define_element_type2(selected_genes)
  
  elemType <- info[[1]]
  gene_name <- info[[2]]
  
  for (i in 1:length(selected_genes)) {
    print(i)
    file_name <- paste0(gene_name[i], '_',elemType[i])
    png( paste0(path_out, file_name, save_name, ".png"))
    
    lollipop_full_length(selected_genes[i], mutation, testGE)
    dev.off()
  }
  
}


library(gridExtra)
library(gridGraphics)

save_multiple_lolipops <- function(novelHits, testGE, path_out, save_name = 'allSigHits') {
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  plot_list <- list()
  
  for (cohort in unique(novelHits$cohort)) {
    
    exclude_lymph_melanoma <- ifelse(cohort ==  "Pancan-no-skin-melanoma-lymph", T, F)
    
    donorInfo <- select_cohort(path_donorInfo, cohort, exclude_lymph_melanoma,
                               exclude_hyper_mutated)
    
    print(cohort)
    gr_cohort = gr[which(gr$D_id %in% donorInfo$donor_id)]
    cohort_hits <- novelHits[which(novelHits$cohort == cohort),]
    mutation <- generate_mut_df(gr = gr_cohort, testGEs = testGE, sigHits = cohort_hits$PCAWG_ID)
    
    
    selected_genes <- unique(mutation$binID)
    
    info <- define_element_type2(selected_genes)
    elemType <- info[[1]]
    gene_name <- info[[2]]
    
    plot_list_cohort <- list()
    for (i in seq_along(selected_genes)) {
      message("Generating plot ", i, ": ", gene_name[i], " (", elemType[i], ")")
      
      # Capture base plot into a grid object
      grid_plot <- tryCatch({
        grid.newpage()
        lollipop_full_length(selected_genes[i], mutation, testGE, main = F)
        grid.echo()
        grid.grab()
      }, error = function(e) {
        message("Error capturing plot for ", selected_genes[i], ": ", e$message)
        NULL
      })
      
      if (!is.null(grid_plot)) {
        title_grob <- grid::textGrob(
          paste0(gene_name[i], " (", elemType[i], ")", " in ", cohort),
          gp = grid::gpar(fontface = "bold"),
          x = 0.5, y = unit(1, "npc") - unit(2, "mm"),
          just = "top"
        )
        plot_with_title <- gridExtra::arrangeGrob(grid_plot, top = title_grob)
        plot_list_cohort[[length(plot_list_cohort) + 1]] <- plot_with_title
      }
    }
    plot_list <- append(plot_list, plot_list_cohort)
  }
  
  # Save plots in 2x3 layout per page
  pdf(file = file.path(path_out, paste0(save_name, ".pdf")), width = 12, height = 9)
  
  n_per_page <- 6
  total_pages <- ceiling(length(plot_list) / n_per_page)
  
  for (page in seq_len(total_pages)) {
    from <- (page - 1) * n_per_page + 1
    to <- min(page * n_per_page, length(plot_list))
    do.call(gridExtra::grid.arrange, c(plot_list[from:to], ncol = 2, nrow = 3))
  }
  
  dev.off()
  
  
}


################ expression association #############

expression_boxPlot <- function(all_mutData, all_sampleInfo,
                               path_sample_sheet , path_IDs_maf,
                               path_donorInfo , path_CN, path_FPKM_UQ,
                               candidate, cohort){
  
  donorInfo <- fread(path_donorInfo)
  donorInfo <- donorInfo[!donorInfo$HyperMut_donor,]
  
  complete_count_matrix <- as.data.frame(fread(path_FPKM_UQ))
  complete_cn <- fread(path_CN) # colnames are Tumor sample barcode (the column exists in maf data to match with Donor_ids)
  
  
  cancer_specific_dat <- load_cancer_specific_data(path_donorInfo, cohort, all_mutData, all_sampleInfo)
  cancerSp_mutData = cancer_specific_dat[[1]]
  cancerSp_sampleInfo = cancer_specific_dat[[2]]
  
  cancer_specific_CN_RNA <- load_cohort_specific_CNV_RNA(complete_count_matrix, complete_cn, cancerSp_sampleInfo)
  
  
  cohort_specific_CN <- cancer_specific_CN_RNA$cns_withRNAseqData
  cohort_specific_expression <- cancer_specific_CN_RNA$restricted_count_matrix
  
  candidateDat <- restrict_CN_RNA_forCandidate(cohort_specific_CN, cohort_specific_expression, candidate)
  
  candidate_rna_subset = candidateDat$candidate_rna_subset
  candidate_cn_subset = candidateDat$candidate_cn_subset
  N_with_both_RNA_CN <- candidateDat$N_with_both_RNA_CN
  
  x <- generate_data_glm(candidate, cancerSp_mutData, donorInfo, cohort, candidate_rna_subset, candidate_cn_subset)
  
  final_df <- x$final_df
  
  n_mut = nrow(final_df[which(final_df$MUT == 1),])
  nWT = nrow(final_df[which(final_df$MUT == 0),])
  
  final_df$MUT <- ifelse(final_df$MUT == 1, sprintf('Mutated (N = %s)',n_mut),
                         sprintf('Wild type (N = %s)',nWT))
  
  
  final_df$CNfactor <- ifelse(final_df$CN == 2, 0, ifelse(final_df$CN > 2, 1, -1))
  
  
  
  ggplot(final_df, aes(x = MUT, y = FPKM_UQ)) +
    geom_boxplot(aes(fill = MUT), outlier.shape = NA, fill = NA) +  # Color boxes by MUT and remove outliers
    geom_point(aes(color = CNfactor), position = position_jitter(width = 0.5), alpha = 0.8, size = 2) +  # Add points for FPKM values and color by CN
    scale_y_continuous(trans = 'log10') +
    theme_minimal() +
    labs(y = "FPKM-UQ", x = '') +
    theme(
      axis.text.x = element_text(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      text = element_text(size = 18),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.placement = "outside"
    ) +
    # scale_fill_manual(values = c("Wild type" = "#56B4E9", "Mutated" = "#D55E00")) +  # Customize box colors by group
    scale_color_viridis_c() + scale_color_gradient2(low = "blue", mid = 'skyblue', high = "#e41a1c")  # Color points by CN values with a continuous color scale
  
  
  # save_name = define_element_type(c(candidate))
  # save_name = paste0(save_name[[1]], '_', save_name[[2]])
  # ggsave(paste0("boxPlot", save_name,"_expression.png"),
  #        device = "png", width = 6, height = 5,
  #        bg = 'white',
  #        path = path_save)
}


save_all_expression_boxplots <- function(path_sigHits, path_mutData,
                                         path_sample_sheet, path_IDs_maf,
                                         path_donorInfo, path_CN, path_FPKM_UQ,
                                         output_path, save_name = "expression_boxplots") {
  
  df <- fread(path_sigHits)
  df <- df[which(df$type_of_element != 'enhancers'),]
  df <- data.frame(df[which(!df$cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')),])
  df <- df[!is.na(df$p_value_CN),]
  
  cohorts <- df$cohort
  candidates <- df$PCAWG_ID
  
  all_mutData <- fread(path_mutData)
  all_sampleInfo <- generate_samplesInfo(path_sample_sheet, path_IDs_maf)
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  plot_list <- list()
  
  for (i in 1:length(cohorts)) {
    cohort = cohorts[i]
    candidate = candidates[i]
    
    gene_info = define_element_type2(candidate)
    
    message("Generating plot for ", candidate, " in ", cohort)
    
    p <- tryCatch({
      expression_boxPlot(all_mutData, all_sampleInfo,
                         path_sample_sheet, path_IDs_maf,
                         path_donorInfo, path_CN, path_FPKM_UQ,
                         candidate, cohort)
    }, error = function(e) {
      message("Error in ", candidate, " / ", cohort, ": ", e$message)
      NULL
    })
    
    if (!is.null(p)) {
      title <- paste0( gene_info[[2]][1], " ", gene_info[[1]][1], " in ", cohort)
      p <- p + ggtitle(title)
      plot_list[[length(plot_list) + 1]] <- p
    }
    
  }
  
  # Save plots to PDF in a grid layout (2x3 per page)
  pdf(file = file.path(output_path, paste0(save_name, ".pdf")), width = 12, height = 9)
  n_per_page <- 6
  total_pages <- ceiling(length(plot_list) / n_per_page)
  
  for (page in seq_len(total_pages)) {
    from <- (page - 1) * n_per_page + 1
    to <- min(page * n_per_page, length(plot_list))
    gridExtra::grid.arrange(grobs = plot_list[from:to], ncol = 2, nrow = 3)
  }
  dev.off()
}

save_AUPR_singleCohort_woPCAWGmeyhods <- function(n, PATH_newResults, newRESULT,
                                   element, path_save_plot, based_on, tissue,
                                   
                                   base_method,
                                   path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                                   include_legend = F, mask_threshold = .9
){
  pRes <- create_pRes_MyRes(PATH_newResults, newRES)
  
  preprocessed_pRes <- prepare_pRes_helper(pRes, element, base_method, mask_threshold)
  
  
  dir.create(paste0(path_save_plot),
             showWarnings = F,
             recursive = T)
  
  ann_PCAWG_ID <- fread(path_ann_PCAWG_ID)
  
  
  generate_singlepanel_AUPR_plots(preprocessed_pRes, n, ann_PCAWG_ID, element,
                                  path_save_plot, tissue, based_on, include_legend)
  
}



create_pRes_MyRes <- function(PATH_newResults, newRES){
  
  res <- lapply(PATH_newResults, fread)
  names(res) <- newRES
  class(res) = "pRes"
  
  res <- lapply(res, function(s){
    s <- data.frame(s)
    s <- s[,which(colnames(s) %in% c("PCAWG_ID", "pvals", "p_value"))]
    colnames(s) = c("PCAWG_IDs",  "p_value")
    s
  })
  
  
  res
}


save_heatmap_kTop <- function(path_save_tables_rmHypers, path_saveFig,
                              save_name = 'Fig2ef', k = 5){
  df <- prepare_data_kTop(path_save_tables_rmHypers, k )
  
  # ---- reshape data ----
  mat <- df %>%
    select(method, K_top, elem, Freq) %>%
    pivot_wider(names_from = method, values_from = Freq, values_fill = NA)
  
  row_elem <- mat$elem
  row_k    <- as.factor(mat$K_top)
  
  mat_matrix <- as.matrix(mat %>% select(-K_top, -elem))
  rownames(mat_matrix) <- paste0("K", mat$K_top, "_", seq_len(nrow(mat)))
  
  # ---- reorder columns according to desired order ----
  method_order <- c("iDriver", "DriverPower", "MutSig", "dNdScv", 
                    "ncdDetect", "ActiveDriverWGS", "oncodriveFML (CADD)", "NBR",
                    "oncodriveFML (vest3)", 
                    "LARVA", "ncDriver")
  # keep only columns that exist in mat_matrix
  method_order <- method_order[method_order %in% colnames(mat_matrix)]
  mat_matrix <- mat_matrix[, method_order, drop = FALSE]
  
  # ---- row annotation ----
  
  
  # ---- column annotation ----
  method_factor <- factor(colnames(mat_matrix), levels = colnames(mat_matrix))
  
  
  # ---- color mapping ----
  col_fun <- colorRamp2(
    c(min(mat_matrix, na.rm = TRUE), max(mat_matrix, na.rm = TRUE)),
    c("white", "#c55466")
  )
  
  rownames(mat_matrix) <- sub("_.*", "", rownames(mat_matrix))
  rownames(mat_matrix) <- sub("K", "Top ", rownames(mat_matrix))
  
  # ---- Heatmap ----
  method_colors <- create_method_colours()
  
  ht <- Heatmap(mat_matrix,
                name = "Freq",
                col = col_fun,
                na_col = "grey80",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                # left_annotation = ha_row,
                # top_annotation = ha_col,
                row_split = row_elem,
                column_names_rot = 45,
                column_split = method_factor,
                width = unit(10, "cm"),
                height = unit(10, "cm"),
                column_names_gp = gpar(col = method_colors[colnames(mat_matrix)], fontsize = 13, fontface = 'bold'),
                row_names_gp = gpar(fontsize = 11),
                # heatmap_legend_param = list(title = NULL,
                #                             labels_gp = gpar(fontsize = 10)),
                layer_fun = function(j, i, x, y, w, h, fill) {
                  for (k in seq_along(i)) {
                    val <- mat_matrix[i[k], j[k]]
                    if (!is.na(val)) {
                      grid.text(sprintf("%.0f", val),
                                x[k], y[k],
                                gp = gpar(fontsize = 9, col = "black"))
                    }
                  }
                })
  dir.create(path_saveFig, recursive = T, showWarnings = F)
  # Save as PNG
  png(paste0(path_saveFig, "/", save_name,"_heatmap.png"), width = 2500, height = 3000, res = 500)  # adjust size/resolution
  draw(ht, show_heatmap_legend = FALSE)
  dev.off()
  
}


save_scores_scatterPlots <- function(path_res, path_sigHits, path_save, elemType, save_name = ''){
  
  element_types <- c("gc19_pc.cds", "gc19_pc.3utr", "gc19_pc.5utr", 
                     "gc19_pc.promCore", "enhancers", "gc19_pc.ss") 
  element_colors <- structure(
    # RColorBrewer::brewer.pal(length(element_types), "Set2")[1:length(element_types)],
    c("#A6D854" ,"#FC8D62",  "#8DA0CB", "#E78AC3", "#66C2A5", "#FFD92F", "#E5C494", "#B3B3B3")[1:length(element_types)],
    
    names = element_types
  )
  
  x <- fread(path_res)
  x$averageScores <- x$sum_scores/x$sum_nMut
  sigHits <- fread(path_sigHits)
  sigHits <- sigHits[sigHits$cohort == 'Pancan-no-skin-melanoma-lymph',]
  x$sig <- ifelse(x$PCAWG_ID %in% sigHits$PCAWG_ID, T, F)
  
  # x$genSymbol = unlist(lapply(x$PCAWG_ID, function(s){
  #   strsplit(s, '::')[[1]][3]
  # }
  # ))
  
  elemTypes <- sigHits[,c('PCAWG_ID', 'type_of_element',  "geneSymbol" )]
  elemTypes <- elemTypes[!duplicated(elemTypes),]
  x <- left_join(x, elemTypes)
  
  if (elemType == 'CDS') {
    x <- x[grepl('gc19_pc.cds', x$PCAWG_ID),]
  } else if (elemType == 'NC') {
    x <- x[!grepl('gc19_pc.cds', x$PCAWG_ID),]
  }
  
  
  p <- ggplot(x, aes(x = sum_nMut, y = averageScores)) +
    # non-sig points
    geom_point(data = subset(x, sig == FALSE), alpha = 0.6, color = 'grey') +
    
    # geom_smooth(method = "lm", se = FALSE, color = "#F8D0D0", linetype = "dashed")  +
    scale_x_log10(labels = comma) +
    # scale_y_log10(labels = comma) +
    labs(
      # title = "Scatter Plot of Observed Stats vs Mutations (log scale)",
      x = "Total Number of Mutations",
      y = "Average Score"
    ) +
    theme_minimal()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # legend.title = element_blank(),
          text = element_text(size = 18),
          
          axis.line = element_line(colour = "black"))
  if (elemType == 'CDS') {
    p+
      # sig points  gc19_pc.cds
      geom_point(
        data = subset(x, sig == TRUE & grepl("gc19_pc.cds", PCAWG_ID)),
        color = "#C41E1E", size = 2, alpha = 0.8
      ) +
      geom_text_repel(
        data = subset(x, sig == TRUE & grepl("gc19_pc.cds", PCAWG_ID)),
        aes(label = geneSymbol),
        size = 4,
        max.overlaps = 10,
        color = "#8B0000", 
        fontface = "bold.italic"   # italicize text
      )
  } else if (elemType == 'NC') {
    p +
      # sig non-cds points (use element_colors)
      geom_point(
        data = subset(x, sig == TRUE & (!grepl("gc19_pc.cds", PCAWG_ID)) ),
        aes(color = type_of_element),
        size = 2, alpha = 0.8,
        show.legend = F
      )+
      geom_text_repel(
        data = subset(x, sig == TRUE ),
        aes(label = geneSymbol,
            color = type_of_element),
        size = 4,
        max.overlaps = 7, show.legend = F,
        fontface = "bold.italic"   # italicize text
      )+
      scale_color_manual(values = element_colors, name = "")
  }
  
  
  
  dir.create(path_save, recursive = T, showWarnings = F)
  ggsave(paste0("scatterPlot_scores_nMuts", elemType, save_name,".png"),
         device = "png", width = 6.5, height = 5,#width = 7.5
         bg = 'white',
         path = path_save)
}


nMuts_perCohort_stackedPlot <- function(path_donorInfo, path_save, save_name = '') {
  df <- fread(path_donorInfo) %>% as_tibble()
  df <- df[,c('D_id', 'cohort1', 'freq', 'HyperMut_donor')]  
  colnames(df) <- c('patient_id', 'cohort', 'mutations', 'HyperMut_donor')
  
  ### Prepare for plotting
  df <- df %>%
    group_by(cohort) %>%
    arrange(mutations, .by_group = TRUE) %>%
    mutate(order_in_cohort = row_number()) %>%
    ungroup()
  
  # Count patients per cohort and reorder cohorts
  cohort_sizes <- df %>%
    count(cohort, name = "n_patients") %>%
    arrange(desc(n_patients)) %>%
    filter(n_patients > 17)
  
  df <- df %>% filter(cohort %in% cohort_sizes$cohort)
  df$cohort <- factor(df$cohort, levels = cohort_sizes$cohort)
  
  # Count number of hypermutated donors per cohort
  hyper_counts <- df %>%
    group_by(cohort) %>%
    summarize(hyper_count = sum(HyperMut_donor, na.rm = TRUE)) %>%
    filter(hyper_count > 0)
  
  # Count total patients per cohort
  patient_counts <- df %>%
    group_by(cohort) %>%
    summarize(patient_count = n())
  
  # Merge counts for positioning
  counts <- left_join(patient_counts, hyper_counts, by = "cohort") %>%
    mutate(hyper_count = ifelse(is.na(hyper_count), 0, hyper_count))
  
  breaks_custom <- c(-Inf, 100, 1000, 10000, 20000, 30000, 40000, 
                     50000, 60000, 70000, 80000, 90000, Inf)
  
  labels_custom <- c(
    "<100", "100–1k", "1k–10k", "10k–20k", "20k–30k", 
    "30k–40k", "40k–50k", "50k–60k", "60k–70k", 
    "70k–80k", "80k–90k", ">90k"
  )
  
  # Assign bins with correct ordering
  df$mut_bin <- cut(
    df$mutations,
    breaks = breaks_custom,
    labels = labels_custom,
    include.lowest = TRUE,
    right = FALSE,
    ordered_result = TRUE
  )
  
  # Reverse levels so that >90k (highest bin) comes last (plotted on top)
  df$mut_bin <- factor(df$mut_bin, levels = rev(labels_custom), ordered = TRUE)
  
  # Define colors from light (low) to dark (high)
  colors_custom <- rev(
    
    colors <- c(
      
      '#FFFFFF' , 
      '#F2F2F2'  ,
      '#E0E0E0'  ,
      '#CCCCCC'  ,
      '#B3B3B3'  ,
      "#F8D0D0", # very light red
      "#F4A6A6", # light red
      "#EF7B7B", # soft red
      "#E65454", # medium red
      "#D93030", # stronger red
      "#C41E1E", # vivid red
      
      "#8B0000"
      
    )
    
  #   c(
  #   '#fbe4d1', '#f7c9aa', '#f6a47c', '#f47d57',
  #   '#ed503e', '#eb3d43', '#d92847', '#b91657',
  #   '#921c5b', '#691f55', '#451c47', '#221331'
  # )
    )
  
  # Bar plot with colored fill and text labels
  p <- ggplot(df, aes(x = cohort, y = 1, fill = mut_bin)) +
    geom_bar(stat = "identity", position = "stack") +
    # geom_text(
    #   data = counts %>% filter(hyper_count > 0),
    #   aes(
    #     x = cohort,
    #     y = patient_count + 6,
    #     label = hyper_count
    #   ),
    #   inherit.aes = FALSE,
    #   size = 4.3
    # ) +
    scale_fill_manual(
      values = colors_custom, name = "#Mutations", drop = FALSE
    ) +
    labs(
      x = "",
      y = "Number of Patients"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(colour = "black")
    )
  
  dir.create(path_save, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("stacked_nMutsCohorts", save_name,".png"),
         plot = p,
         device = "png", width = 14, height = 5,
         bg = 'white',
         path = path_save)
}


plot_sumScores_stacked <- function(path_perDonorOnsStats, path_testY, path_save, save_name = '') {
  # sigHits <- fread('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv')
  # sigHits <- sigHits[sigHits$cohort == 'Pancan-no-skin-melanoma-lymph',]
  # sigHits <- sigHits[sigHits$type_of_element == 'gc19_pc.cds',]
  
  
  df <- fread(path_perDonorOnsStats)
  y <- fread(path_testY)
  # y <- y[y$binID %in% sigHits$PCAWG_ID,]
  y <- y[!grepl('lncrna', y$binID),]
  y <- y[!grepl('gc19_pc.ss', y$binID),]
  y <- y[order(y$nSample, decreasing = T),]
  y <- y[1:65,]
  y$geneSymbol <- extract_HugoSymbol(y$binID)
  head(y)
  df <- df[df$PCAWG_ID %in% y$binID,]
  df$geneSymbol <- extract_HugoSymbol(df$PCAWG_ID)
  df <- data.frame(df)
  columns <- c("donor_id", "PCAWG_ID", "sumScores", "nMut", "obs_stat", "type_of_element", "geneSymbol")
  df <- df[,which(colnames(df) %in% columns)]
  
  df <- df %>%
    group_by(geneSymbol) %>%
    arrange(sumScores, .by_group = TRUE) %>%
    mutate(order_in_cohort = row_number()) %>%
    ungroup()
  
  df$geneSymbol <- factor(df$geneSymbol, levels = y$geneSymbol)
  
  # Basic stats
  min_val <- base::min(df$sumScores, na.rm = TRUE)
  q1_val  <- stats::quantile(df$sumScores, 0.25, na.rm = TRUE)
  med_val <- stats::median(df$sumScores, na.rm = TRUE)
  q3_val  <- stats::quantile(df$sumScores, 0.75, na.rm = TRUE)
  max_val <- base::max(df$sumScores, na.rm = TRUE)
  
  # Define custom breaks (in increasing order)
  breaks_custom <- sort(c(
    min_val, -2, -1, -.5, 0, q1_val, .5, med_val, 1.5,
    2, 2.5, 3, q3_val, 6, 15, 20, max_val
  ))
  
  # Labels for bins
  labels_custom <- paste0(
    head(round(breaks_custom, 2), -1), "–", 
    tail(round(breaks_custom, 2), -1)
  )
  # Sort breaks increasing 
  breaks_custom_sorted <- sort(breaks_custom) 
  # Create labels in increasing order 
  labels_custom_sorted <- c( paste0("≤ ", round(breaks_custom_sorted[2], 2)),
                             paste0(round(breaks_custom_sorted[2], 2), "–", round(breaks_custom_sorted[3], 2)), 
                             paste0(round(breaks_custom_sorted[3], 2), "–", round(breaks_custom_sorted[4], 2)), 
                             paste0(round(breaks_custom_sorted[4], 2), "–", round(breaks_custom_sorted[5], 2)), 
                             paste0(round(breaks_custom_sorted[5], 2), "–", round(breaks_custom_sorted[6], 2)), 
                             paste0(round(breaks_custom_sorted[6], 2), "–", round(breaks_custom_sorted[7], 2)), 
                             paste0(round(breaks_custom_sorted[7], 2), "–", round(breaks_custom_sorted[8], 2)), 
                             paste0(round(breaks_custom_sorted[8], 2), "–", round(breaks_custom_sorted[9], 2)), 
                             paste0(round(breaks_custom_sorted[9], 2), "–", round(breaks_custom_sorted[10], 2)), 
                             paste0(round(breaks_custom_sorted[10], 2), "–", round(breaks_custom_sorted[11], 2)), 
                             paste0(round(breaks_custom_sorted[11], 2), "–", round(breaks_custom_sorted[12], 2)), 
                             paste0(round(breaks_custom_sorted[12], 2), "–", round(breaks_custom_sorted[13], 2)), 
                             paste0(round(breaks_custom_sorted[13], 2), "–", round(breaks_custom_sorted[14], 2)), 
                             paste0(round(breaks_custom_sorted[14], 2), "–", round(breaks_custom_sorted[15], 2)), 
                             paste0(round(breaks_custom_sorted[15], 2), "–", round(breaks_custom_sorted[16], 2)), 
                             paste0(round(breaks_custom_sorted[16], 2), "–", round(breaks_custom_sorted[17], 2)) )
  # Assign bins as ordered factor (low → high)
  df$score_bin <- cut(
    df$sumScores,
    breaks = breaks_custom,
    labels = labels_custom,
    include.lowest = TRUE,
    right = TRUE,
    ordered_result = TRUE
  )
  # Make sure levels are from lowest to highest
  df$score_bin <- factor(df$score_bin, 
                         levels = labels_custom_sorted, 
                         ordered = TRUE)
  
  # Reverse the levels so lowest bins are drawn first (at bottom of stack)
  df$score_bin <- forcats::fct_rev(df$score_bin)
  
  colors_custom <- rev(c(
    "#08306B", "#2171B5", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7",
    "#F2F2F2", "#FFFFFF", "#FFF0F0", "#FFE0E0", "#FBDADA", "#F8D0D0",
    "#F4A6A6", "#EF7B7B", "#D93030", "#C41E1E", "#8B0000"
  ))#[1:length(levels(df$score_bin))]
  
  # Define element colors
  element_types <- c("gc19_pc.cds", "gc19_pc.3utr", "gc19_pc.5utr", 
                     "gc19_pc.promCore", "enhancers", "gc19_pc.ss") 
  
  element_colors <- c("#A6D854" ,"#FC8D62",  "#8DA0CB", 
                      "#E78AC3", "#66C2A5", "#FFD92F")
  names(element_colors) <- element_types
  
  # Mapping: geneSymbol → color
  gene_colors <- df %>% 
    select(geneSymbol, type_of_element) %>% 
    distinct() %>% 
    mutate(color = element_colors[type_of_element])
  
  # Ensure geneSymbol order matches `y`
  df$geneSymbol <- factor(df$geneSymbol, levels = y$geneSymbol)
  
  # Create stacked plot (bars stacked by increasing score_bin)
  p <- ggplot(df, aes(x = geneSymbol, y = 1, fill = score_bin)) +
    geom_bar(stat = "identity", position = "stack", show.legend = F) +
    scale_fill_manual(values = colors_custom, drop = FALSE) +
    labs(x = "", y = "# Mutated Samples", fill = "sumScores") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_markdown(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(labels = function(x) {
      # Color each gene symbol
      sapply(x, function(g) {
        col <- gene_colors$color[gene_colors$geneSymbol == g]
        paste0("<span style='color:", col, "'>", g, "</span>")
      })
    })
  
  # Save plot
  dir.create(path_save, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = paste0("stacked_sumScores", save_name,".png"),
    plot = p, device = "png",width = 15, height = 6, # width = 19, height = 6, 
    bg = 'white', path = path_save
  )
}


plot_sumScores_stacked_v2 <- function(path_perDonorOnsStats, path_testY, path_save, save_name = '') {
  df <- fread(path_perDonorOnsStats)
  y <- fread(path_testY)
  y <- y[!grepl('lncrna', y$binID),]
  y <- y[!grepl('gc19_pc.ss', y$binID),]
  y <- y[order(y$nSample, decreasing = TRUE),]
  y <- y[1:65,]
  y$geneSymbol <- extract_HugoSymbol(y$binID)
  
  df <- df[df$PCAWG_ID %in% y$binID,]
  df$geneSymbol <- extract_HugoSymbol(df$PCAWG_ID)
  df <- data.frame(df)
  columns <- c("donor_id", "PCAWG_ID", "sumScores", "nMut", "obs_stat", "type_of_element", "geneSymbol")
  df <- df[, which(colnames(df) %in% columns)]
  
  # order within each geneSymbol so smaller sumScores at bottom
   
  df <- df %>%
    group_by(geneSymbol) %>%
    arrange(desc(sumScores), .by_group = TRUE) %>%
    mutate(order_in_gene = row_number()) %>%
    ungroup()
  
  # keep gene ordering
  df$geneSymbol <- factor(df$geneSymbol, levels = y$geneSymbol)
  
  # stats for midpoint
  med_val <- stats::median(df$sumScores, na.rm = TRUE)
  
  # element colors
  element_types <- c("gc19_pc.cds", "gc19_pc.3utr", "gc19_pc.5utr", 
                     "gc19_pc.promCore", "enhancers", "gc19_pc.ss") 
  element_colors <- c("#A6D854" ,"#FC8D62",  "#8DA0CB", 
                      "#E78AC3", "#66C2A5", "#FFD92F")
  names(element_colors) <- element_types
  
  # map each geneSymbol → color
  gene_colors <- df %>%
    select(geneSymbol, type_of_element) %>%
    distinct() %>%
    mutate(color = element_colors[type_of_element])
  
  # continuous stacked plot (stacking respects order_in_gene)
  p <- ggplot(df, aes(x = geneSymbol, y = 1, fill = sumScores, group = order_in_gene)) +
    geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
    scale_fill_gradient2(low = "#2171B5", mid = "white", high = "#C41E1E", #low = "#08306B", mid = "white", high = "#8B0000",
                         med_val,
                         name = "sumScores") +
    labs(x = "", y = "# Mutated Samples") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_markdown(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(labels = function(x) {
      # color each gene symbol
      sapply(x, function(g) {
        col <- gene_colors$color[gene_colors$geneSymbol == g]
        paste0("<span style='color:", col, "'>", g, "</span>")
      })
    })
  
  # save
  dir.create(path_save, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = paste0("stacked_sumScoresV2", save_name, ".png"),
    plot = p, device = "png", width = 15, height = 6, bg = 'white', path = path_save
  )
}

save_nFPs_heatmap <- function(path_all_methods_nFPs, path_save, elem, save_name = ''){
  # read input
  # df <- fread('../extdata/output_release2.0/simulated/all_methods_nFPs/nFPs_vsPCAWG_dpIdrBase.tsv')
  df <- fread(path_all_methods_nFPs)
  df$cohort <- ifelse(df$cohort == 'Pan_Cancer', 'Pancancer', df$cohort)
  df$cohort <- ifelse(df$cohort == 'Pancan-no-skin-melanoma-lymph',
                      'Pancancer*', df$cohort)
  # filter out excluded methods and cohorts
  # df <- df[method != "oncodriveFML_cadd" &
  #            !cohort %in% c( "Pancan-no-skin-melanoma-lymph")]
  df$method <- ifelse(df$method == "oncodriveFML_cadd", "oncodriveFML (CADD)", df$method)
  df$method <- ifelse(df$method == "oncodriveFML_vest3", "oncodriveFML (vest3)", df$method)
  df <- df[method != "compositeDriver"]
  if (elem != '') {
    df <- df[df$element == elem,]
  }
  
  # make sure all combinations exist (so missing ones can be gray)
  all_methods <- unique(df$method)
  all_cohorts <- unique(df$cohort)
  df_full <- CJ(method = all_methods, cohort = all_cohorts)
  df_full <- merge(df_full, df, by = c("method", "cohort"), all.x = TRUE)
  
  method_colors <- create_method_colours()
  
  # wrap y-axis labels in <span style='color:...'>...</span>
  df_full$method <- factor(df_full$method, levels = names(method_colors))
  
  y_labels_colored <- setNames(
    paste0("<span style='color:", method_colors, "'>", names(method_colors), "</span>"),
    names(method_colors)
  )
  
  p <- ggplot(df_full, aes(x = cohort, y = method, fill = nFPs)) +
    geom_tile(color = "black") +
    geom_text(
      aes(label = ifelse(!is.na(nFPs) & nFPs > 0, nFPs, "")),
      size = 4, color = "black"
    ) +
    scale_fill_gradient(
      low = "white", high = "#d7301f", na.value = "grey80",
      name = "nFPs"
    ) +
    scale_y_discrete(labels = y_labels_colored, position = "left") +  # apply colored labels
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
      axis.text.y = element_markdown(size = 18, face = 'bold'),   # <-- enable markdown coloring
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.position = "none", 
      legend.text = element_text(size = 14),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  # save
  dir.create(path_save, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = paste0("nFPs_heatmap_", elem, "_", save_name, ".png"),
    plot = p, device = "png", width = 14, height = 8, bg = 'white', path = path_save
  )
}

########################### ablation plots ########################### 
prepare_data_perElemType1 <- function(path_tables, elem_type){
  files <- list.files(path_tables, full.names = T)
  
  cds_files <- files[grepl('_CDS', files)]
  NC_files <- files[grepl('_NC', files)]
  
  
  if (elem_type == 'CDS') {
    cancer_types <- cds_files[
      !grepl("Pan_Cancer|Pancan-no-skin-melanoma-lymph", cds_files)
      ]
    
    PanCancer <- cds_files[
      grepl("Pancan-no-skin-melanoma-lymph", cds_files)
      ]
  } else if (elem_type == 'NC'){
    cancer_types <- NC_files[
      !grepl("Pan_Cancer|Pancan-no-skin-melanoma-lymph", NC_files)
      ]
    
    PanCancer <-  NC_files[
      grepl("Pancan-no-skin-melanoma-lymph", NC_files)
      ]
  }
  
  
  PanCancer <- fread(PanCancer)
  
  x <- lapply(cancer_types, fread)
  df <- do.call(rbind, x)
  
  setDT(df)
  
  df_sums <- df[, .(
    nTPs = sum(nTPs, na.rm = TRUE),
    nHits = sum(nHits, na.rm = TRUE)
  ), by = method]
  
  list('Pancancer' = PanCancer, 'Cancer-specific' = df_sums)
  
}

prepare_data_perElemType2 <- function(df_sums, elem_type, cancer_category){
  df_sums <- df_sums[df_sums$method %in% c("iDriver", "Population-level", "Only recurrence", "Only FI"),]
  df_sums$notReported <- df_sums$nHits - df_sums$nTPs
  df_elem <- df_sums[,c("method", "nTPs", "notReported")]
  df_elem$elem_type <- elem_type
  df_elem$cancer_category <- cancer_category
  df_elem
  
}

save_ablation_stackedPlot <- function(path_tables, path_save, save_name = ''){
  
  df <- c()
  for (elem_type in c('CDS', 'NC')) {
    cancerList_df <- prepare_data_perElemType1(path_tables, elem_type)
    df_elem <- rbind(prepare_data_perElemType2(cancerList_df[['Cancer-specific']], elem_type, cancer_category = 'Cancer-specific'),
                     prepare_data_perElemType2(cancerList_df[['Pancancer']], elem_type, cancer_category = 'Pancancer'))
    df <- rbind(df, df_elem)
  }
  
  df$cancer_category <- ifelse(df$cancer_category == 'Pancancer', 'Pancancer*', df$cancer_category)
  
  setDT(df)
  
  # reshape into long format
  df_long <- melt(
    df,
    id.vars = c("method", "elem_type", "cancer_category"),
    measure.vars = c("nTPs", "notReported"),
    variable.name = "type",
    value.name = "count"
  )
  
  # ensure factor ordering
  df_long[, method := factor(method, levels = c("iDriver", "Population-level", "Only recurrence", "Only FI"))]
  df_long[, type := factor(type, levels = c("notReported", "nTPs"))]  # nTPs first → bottom
  
  
  
  p <- ggplot(df_long[df_long$elem_type == "CDS" & df_long$count > 0, ],
         aes(x = method, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = count),
              position = position_stack(vjust = 0.5),
              color = "white", size = 4) +
    facet_grid(~cancer_category ,
               scales = "free_x", space = "free_x") +
    labs(x = "Method", y = "Count", fill = "") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      panel.spacing = unit(1, "lines")
    ) +
    scale_fill_manual(
      values = c("nTPs" ="#d7301f", #,  "#c55466"
                 "notReported" = "grey"),
      breaks = c("nTPs", "notReported")
    )
  dir.create(path_save, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = paste0("Ablation_stackedBar", "_", save_name, ".png"),
    plot = p, device = "png", width = 6, height = 4, bg = 'white', path = path_save
  )
}


save_upsetAblations <- function(path_files, path_save, save_name = '', Just_cds = TRUE) {
  dfs <- lapply(path_files, function(s) {
    df <- data.table::fread(s)
    df <- df[!(df$cohort %in% c('Pan_Cancer', 'Lymph-BNHL', 'Pancan-no-skin-melanoma-lymph')), ]
    if (Just_cds) {
      df <- df[grepl('gc19_pc.cds', df$PCAWG_ID), ]
    }
    df
  })
  
  # --- Combine into gene_list ---
  gene_list <- list(
    iDriver = dfs[[1]]$geneSymbol,
    `Population-level` = dfs[[2]]$geneSymbol,
    `Only-recurrence` = dfs[[3]]$geneSymbol
  )
  
  method_colors <- create_method_colours()
  
  # Open device
  png(file.path(path_save, paste0("UpsetPlot", save_name, ".png")),
      width = ifelse(Just_cds, 2000, 2600), height = 1770, res = 300)
  
  # Force plotting inside device
  print(
    upset(fromList(gene_list), 
          order.by = "freq", 
          mainbar.y.label = "", 
          sets.x.label = "",
          sets.bar.color = method_colors[names(gene_list)],
          text.scale = c(3, 3, 3, 3),
          queries = list(
            list(query = intersects,
                 params = list("iDriver"),
                 color = "darkred",
                 active = TRUE)))
  )
  
  dev.off()
}

get_gene_color <- function(gene, df) {
  gene_rows <- df[df$geneSymbol == gene, ]
  in_CGC_new <- any(gene_rows$in_CGC_new)
  in_onco <- any(gene_rows$in_oncoKB)
  in_pcawg <- any(gene_rows$in_pcawg)
  
  if (!in_CGC_new && !in_onco && !in_pcawg) {
    "red"
  } else if (in_CGC_new) {
    "#238b45"
  } else if (in_onco) {
    "#08306b"
  } else {
    "black"
  }
}


save_oNlyHyper_plot <- function(path_res, path_ann, sig_threshold,
                                justTop100, path_save, save_name = ''){
  
  dir.create(path_save, recursive = T, showWarnings = F)
  
  hits <- get_sigHits_cohort(path_res, sig_threshold)
  
  if (nrow(hits) < 10) {
    print(paste0('Significant elements are:' , nrow(hits), ' We are using top100 pvals'))
    hits <- fread(path_res)
    hits <- hits[order(hits, decreasing = F),]
  }
  # hits<- fread('../extdata/output_release2.0/observed/onlyHyperMuts/betaRhoPancan_ElemTypeCatLengthQuartile/observed_OnlyHyperMuts.tsv')
  # 
  
  ann <- fread(path_ann)
  hits<- left_join(hits, ann, by = c('PCAWG_ID' = 'PCAWG_IDs'))
  hits <- data.frame( hits[order(hits$pvals, decreasing = F),])
  if(justTop100) {
    hits <- hits[1:100,]
  }
  newHits <- hits[(!(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB | hits$in_CGC_literature)),]
  print(newHits$PCAWG_ID)
  
  get_df_elemTypes <- function(hits, elemType){
    if (elemType == 'Coding') {
      hits <- hits[grepl('gc19_pc.cds', hits$PCAWG_ID),]
    } else if (elemType == 'Non-coding'){
      hits <- hits[!grepl('gc19_pc.cds', hits$PCAWG_ID),]
    }
    newHits <- hits[(!(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB | hits$in_CGC_literature)),]
    
    df <- data.frame('CGC' = sum(hits$in_CGC_new), 'OncoKB' = sum(hits$in_oncoKB),
                     'PCAWG' = sum(hits$in_pcawg), 'Not reported' = nrow(newHits))
    
    
    
    
    # Reshape to long format
    df_long <- df %>%
      pivot_longer(cols = everything(), names_to = "category", values_to = "count")
    df_long$category <- ifelse(df_long$category == 'Not.reported', 'Not reported', df_long$category )
    # Assign custom colors based on category
    df_long <- df_long %>%
      mutate(color = case_when(
        category == "Not reported" ~ "red",
        category == "CGC" ~ "#238b45",
        category == "OncoKB" ~ "#08306b",
        category == "PCAWG" ~ "black"
      ))
    
    df_long <- df_long %>%
      mutate(category = factor(category, levels = category[order(-count)]))
    df_long$type <- elemType
    df_long
  }
  
  # Get data
  df_coding <- get_df_elemTypes(hits, 'Coding')
 
  df_noncoding <- get_df_elemTypes(hits, 'Non-coding')
  
  # Function to lighten colors
  lighten <- function(col, factor = 0.5) {
    rgb_val <- col2rgb(col)
    rgb(t(rgb_val + (255 - rgb_val) * factor), maxColorValue = 255)
  }
  
  # Lighten colors for non-coding
  df_noncoding <- df_noncoding %>%
    mutate(color = lighten(color, 0.6))
  
  # Combine
  df_all <- bind_rows(df_coding, df_noncoding)
  
  # Compute totals for ordering
  order_df <- df_all %>%
    group_by(category) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total))
  
  df_all$category <- factor(df_all$category, levels = order_df$category)
  df_all$type <- factor(df_all$type, levels = c("Non-coding", "Coding"))
  
  # Plot
  ggplot(df_all, aes(x = category, y = count, fill = interaction(type, category))) +
    geom_bar(stat = "identity", show.legend = F) +
    scale_fill_manual(
      values = setNames(df_all$color, interaction(df_all$type, df_all$category)),
      labels = paste(df_all$type, df_all$category)
    ) +
    geom_text(aes(label = count),
              position = position_stack(vjust = 0.5),
              color = "white", size = 4) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = "Count", title = "Coding vs Non-coding elements") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey90"),
          panel.spacing = unit(1, "lines")
    )
  
}

