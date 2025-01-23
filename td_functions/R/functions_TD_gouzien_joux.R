###From sequence to frequencies
extract_genes_based_on_functional_annotation <- function(pattern_annotaion, annotaion_data_frame,percentage=1,exclude=FALSE){
  #This function filter the annottaion data frame to keep only genes that have the `pattern_annotation` character in their protein annotation + only genes that are translated into proteins (remove genes coding for ribosomal RNA for example) 
  
  output_df=annotaion_data_frame %>% dplyr::filter(grepl(pattern_annotaion, product)) %>% dplyr::filter(!is.na(protein_id))
  if (exclude){
    output_df=annotaion_data_frame %>% dplyr::filter(!grepl(pattern_annotaion, product)) %>% dplyr::filter(!is.na(protein_id))
    
  }
  output_df=dplyr::slice_sample(output_df,prop=percentage,replace = F)
 
  return(output_df)
}

extract_sequence_from_annotation <- function(row,genome_sequence_file) {
  #This function extract the corresponding sequence from the .fna file given the informations in the .gff file (gene start position, end position, strand location of the gene). If needed (gene on the "-" strand), the sequence returned is the reverse complement of the .fna file. The function is coded to be used with an "apply"-like syntax.
  start = as.integer(row["start"])
  end = as.integer(row["end"])
  strand=row["strand"]
  sequence = subseq(genome_sequence_file, start, end)
  if (strand == "-") {
    sequence = Biostrings::reverseComplement(sequence)
  }
  return(as.character(sequence))
}


add_seq_to_annotation_df <- function(annotation_df,genome_sequence_file){
  #This function apply the previous one to the annotaion data frame
  surexpressed_genes = apply(annotation_df, 1, extract_sequence_from_annotation, genome_sequence_file)
  annotation_df$sequence=surexpressed_genes
  return(annotation_df)
}


slice_gene_from_start_to_stop <- function(row,verbose){
  
  #This function find the first occurrence of the start codon in a sequence, then slice the sequence by 3 starting from the codon start. Then, it inspects the codons list and stop whenever a stop codon is encountered. 
  
  gene_seq=row['sequence']
  gene_name=as.character(row['product'])
  
  
  first_occurence_start_codon=regexpr('ATG', gene_seq)
  start_position <- ifelse(first_occurence_start_codon != -1, first_occurence_start_codon, NA)
  if (is.na(start_position)){
    if (verbose){
    print(sprintf('No start codon found for the "%s" gene',gene_name))
      }
    return(NA)
  }
  
  codons=unlist(sapply(seq(start_position,nchar(gene_seq),3), function(i) substr(gene_seq, i , i+2))) #slicing by 3 starting from start_position
  
  stop_codon_list=c('TAA','TAG','TGA')
  sliced_gene=c()
  for (codon in codons){
    if (codon %in% stop_codon_list){
      sliced_gene=c(sliced_gene,codon)
      return(sliced_gene)
    }
    sliced_gene=c(sliced_gene,codon)
  }
  if (length(tail(sliced_gene))!= 3){
    if (verbose){
    
    print(sprintf('No stop codons found for the "%s " gene',gene_name))
    }
    return(sliced_gene[-length(sliced_gene)]) # if no codon stop encountered, the last codon is removed when its length is != from 3
  }
  return(sliced_gene)
}


add_sliced_codons_to_annotaion_df <- function(anotation_df_with_sequences,verbose=F){
  trimmed_sequences=apply(anotation_df_with_sequences,1,slice_gene_from_start_to_stop,verbose)
  anotation_df_with_sequences$codons=trimmed_sequences
  return(anotation_df_with_sequences)
}



attribute_codons_to_amino_acid <- function(df_with_column_with_codons) {
  # This function gather all the codons related to a given amino acid into a dicitonnary like structure. It is like this : dictionnary[[amino_acid]]=list__of_correponding_amino_acids
  
  all_codons <- unlist(df_with_column_with_codons$codons) 
  valid_codons <- all_codons[all_codons %in% names(Biostrings::GENETIC_CODE)]
  invalid_codons <- all_codons[!all_codons %in% names(Biostrings::GENETIC_CODE)]
  if (length(invalid_codons) > 0) {
    message(sprintf("The following codons are invalid: %s", paste(invalid_codons, collapse = ", ")))
  }
  amino_acid_mapping <- Biostrings::GENETIC_CODE[valid_codons] # create a named character vector where names=codons values=aminoacids
  
  result <- lapply(unique(amino_acid_mapping), function(aa) { # transfer the previous character vector to a dict like structure
    valid_codons[amino_acid_mapping == aa]
  })
  names(result) <- unique(amino_acid_mapping)  # assign names to the dict like structure
  return(result)
}



count_codons_per_amino_acids <- function(dict_codons_to_aa,aminoacid_to_codons){
  # This function counts the number of occurences of codons per amino acid, 
  output_dict <- list()
  
  for (amino_acid in ls(dict_codons_to_aa)) {
    list_codons <- dict_codons_to_aa[[amino_acid]]
    counts <- table(list_codons)
    codon_count_hash <- list()
    
    for (verif_codon in aminoacid_to_codons[[amino_acid]]) {
      if (verif_codon %in% names(counts)) {
        codon_count_hash[[verif_codon]] <- counts[verif_codon]
      } else {
        codon_count_hash[[verif_codon]] <- 0
      }
    }
    
    output_dict[[amino_acid]] <- codon_count_hash
  }  
  return(output_dict)
}


convert_data_for_plotting <- function(per_amino_acid_total_counts){
  #This function takes the output of the count_codons_per_amino_acids function and convert it to a suitable df that will draw the frequency plot. The frequencies are computed, as long as the theoretical uniform probability associated with this frequency and the type of synonymous family the amino_acid belongs to. The df is then ordered and some formats converted to draw a convenient plot afterward. 
  per_amino_acid_freq_df <- data.frame(amino_acid=factor(),codons=factor(),frequency=numeric(),total_per_aa=numeric(),uniform_prob=numeric(),synonymous_fam=numeric())
  
  
  for (amino_acid in names(per_amino_acid_total_counts)){
    codons <- names(per_amino_acid_total_counts[[amino_acid]])
    counts <- unlist(per_amino_acid_total_counts[[amino_acid]])
    total <- sum(counts)
    uniform=1/length(counts)
    for(codon in codons){
      
      count=per_amino_acid_total_counts[[amino_acid]][[codon]]
      freq= count/total
      
      per_amino_acid_freq_df <- rbind(per_amino_acid_freq_df, data.frame(amino_acid=amino_acid,codons=codon, frequency=freq, total_per_aa=total,uniform_prob=uniform,synonymous_fam=length(counts)))
      #per_amino_acid_freq_df <- rbind(per_amino_acid_freq_df, data.frame(amino_acid="",codons="", frequency=NA, uniform_prob=NA))
      
    }
    
  }
  # Re-order +format conversion 
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% 
    arrange(synonymous_fam, amino_acid, codons) %>%
    mutate(
      codons = factor(codons, levels = unique(codons)),
      amino_acid = factor(amino_acid, levels = unique(amino_acid))
    ) 
  return(per_amino_acid_freq_df)
}


# drawing_theoretical_probability_assistant <- function(per_amino_acid_freq_df){
# line_data <- per_amino_acid_freq_df %>% 
#   group_by(synonymous_fam) %>% 
#   summarize(
#     x_start = as.numeric(codons[1]),            # First codon position
#     x_end = as.numeric(codons[n()]),            # Last codon position
#     y = unique(uniform_prob)                    # Uniform probability
#   )
# line_data <- line_data %>% mutate(x_start=x_start-0.5,x_end=x_end+0.5)
# return(line_data)
# }



draw_frequency_plot <- function(per_amino_acid_freq_df,title,radial=FALSE,rainbow=TRUE){
  
  #This function draws the barplot of the frequencies, ordered by ascending synonymous family (AAs with 1 codon in the family, then AAs with 2 codons, then AAs with 3 codons, etc..), groups the codons per amino acids, and display a polar are classical view of the plot. You can also modify the rainbow parameter to display better colors differences between amino acids codons. 
  
  g=ggplot(data=per_amino_acid_freq_df, aes( x = codons, y = frequency, fill = amino_acid,label=codons)) +
    geom_col( position = position_dodge2(preserve = "total")) +
    labs(title = title,x = "codons",y = "frequency",fill = "Amino Acid") +
    geom_text(angle = 90, vjust = 0.5, hjust=-0.1,size=2) + 
    theme_minimal() + 
    theme(axis.text.x = element_blank())+ 
    geom_step(aes( x=as.numeric(codons)-0.5,y=uniform_prob,group=1,color="Expected frequence \n(uniform probability per \namino_acids)"))+
    scale_color_manual(name="",values=c("Expected frequence \n(uniform probability per \namino_acids)"="black"))
  #geom_segment(data=line_data, aes(x=x_start,, xend = x_end,y=y,yend = y), linewidth =0.5,linetype="solid",inherit.aes = FALSE)
  if (!rainbow){
    #We change the order of the colors in the ggplot2 palete to avoid a rainbow look
    n <- length(levels(per_amino_acid_freq_df$amino_acid)) # number of colors
    set.seed(29)
    cols <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1)(n)[order(sample(1:n, n))]
    g=g+scale_fill_manual(values = cols)
  }
  if (radial){
    g=g+coord_radial(rotate.angle = TRUE)
  }
  return (g)
}

from_annotation_to_plot <- function(annotations_df,genome_sequence,pattern_annotation,aminoacid_to_codons,to_exclude=FALSE,percent=1){
  #This function launch all the previous written functions
  ribosomal_genes=extract_genes_based_on_functional_annotation(pattern_annotation,annotations_df,exclude=to_exclude,percentage = percent)
  ribosomal_genes=add_seq_to_annotation_df(ribosomal_genes,genome_sequence)
  ribosomal_genes_scliced=add_sliced_codons_to_annotaion_df(ribosomal_genes,verbose=F)
  repartition_codons=attribute_codons_to_amino_acid(ribosomal_genes_scliced)
  per_amino_acid_total=count_codons_per_amino_acids(repartition_codons,aminoacid_to_codons)
  per_amino_acid_freq_df=convert_data_for_plotting(per_amino_acid_total)
  return(per_amino_acid_freq_df)
}
#   plot=draw_frequency_plot(per_amino_acid_freq_df,radial = radial,rainbow = rainbow,title=title)
#   return(plot)
#   
# }


###CUB computation and plot


compute_effective_number_of_codon_usage <- function(per_amino_acid_freq_df){
  #This functions take the df with the attributed  codon count to each amino acid (per_amino_acid_total) and compute the effective number of codon based on the homozygosity defined in Wright 1990. 
  
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% filter(!amino_acid=="*") #remove stop codons
  per_amino_acid_freq_df =per_amino_acid_freq_df %>% mutate(frequency_2 = frequency^2) #compute the squared frequencies
  
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% group_by(amino_acid) %>%
    mutate(homozygosity=(total_per_aa*sum(frequency_2)-1)/(total_per_aa-1)) #compute homozygosity
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% group_by(synonymous_fam) %>%     mutate(effective_number_of_codon_per_aa=1/(mean(homozygosity))) #compute effective number of codon
  return(per_amino_acid_freq_df)
}

compute_cub <- function(per_amino_acid_freq_df){
  #This function take a df with all the columns needed to compute the cub as defined by wright 1990 and return a unique cub value per df
  size_syn_fam_df=per_amino_acid_freq_df %>% group_by(synonymous_fam) %>% distinct(amino_acid)%>% mutate(syn_fam_size=n()) 
  per_amino_acid_freq_df=merge(size_syn_fam_df,per_amino_acid_freq_df)
  to_compute_cub_df=per_amino_acid_freq_df %>% distinct(synonymous_fam,syn_fam_size,effective_number_of_codon_per_aa) 
  cub= to_compute_cub_df %>% mutate(cub=sum(syn_fam_size/(1/effective_number_of_codon_per_aa))) %>%distinct(cub) %>% pull()
  return(cub)
}

plot_effective_number_of_codon_usage <- function(per_amino_acid_freq_df,title){
  #This function plot the effective number of codon per amino acid
  per_amino_acid_freq_df_unqiue <- per_amino_acid_freq_df %>% distinct(amino_acid,effective_number_of_codon_per_aa)
  per_amino_acid_freq_df_unqiue$amino_acid <- factor(per_amino_acid_freq_df_unqiue$amino_acid, 
                                                     levels = unique(per_amino_acid_freq_df_unqiue$amino_acid))
  
  
  ggplot(data = per_amino_acid_freq_df_unqiue, aes(x=amino_acid,y=effective_number_of_codon_per_aa,fill = amino_acid))+geom_col(show.legend = F)+
    geom_step(aes(x = as.numeric(amino_acid)-0.5, y=synonymous_fam, group=1, color="Expected value if \n uniform distribution of the codons"))+
    scale_color_manual(name="",values=c("Expected value if \n uniform distribution of the codons"="black"))+
    labs(title = title,x = "Amino acids",y = "Effective number of codons")
}



###Corrected CUB by genome compo

get_nucleotide_freq <- function(full_genome){
  #This functions computes A,C,G,T frequencies from the full genome.  
  
  nucleotids=Biostrings::alphabetFrequency(full_genome,baseOnly=T)
  total_counts <- sum(nucleotids)
  freq=nucleotids/total_counts
  A=freq[,"A"]
  C=freq[,"C"]
  G=freq[,"G"]
  T=freq[,'T']
  
  return(c(A,C,G,T))
}


compute_Ka_constant <- function(nucleotides_frequencies,amino_acid_one_letter_code){
  # This function compute the trinucleotide frequency Ka constant as described in November 2002, "Accounting for Background Nucleotide Composition When Measuring Codon Usage Bias"
  Ka=0  
  for (codon in aminoacid_to_codons[[amino_acid_one_letter_code]]){
    nucleotid= unlist(strsplit(codon,''))
    freq_product= prod(nucleotides_frequencies[nucleotid])
    
    Ka=Ka+freq_product
  }
  return(Ka)
}

compute_Ka_apply <- function(row,nucleotides_frequencies){
  amino_acid_one_letter_code=row["amino_acid"]
  Ka=compute_Ka_constant(nucleotides_frequencies,amino_acid_one_letter_code)
  return(Ka)
}

compute_ei <- function(codon,nucleotides_frequencies){
  amino_acid_one_letter_code=Biostrings::GENETIC_CODE[codon]
  Ka=compute_Ka_constant(nucleotides_frequencies,amino_acid_one_letter_code)
  nucleotid= unlist(strsplit(codon,''))
  ei=prod(nucleotides_frequencies[nucleotid])/Ka
  return(ei)
}

compute_ei_apply <- function(row,nucleotides_frequencies){
  codon=row["codons"]
  ei=compute_ei(codon,nucleotides_frequencies )
  return(ei)
}


# compute_chi2 <- function(row,nucleotides_frequencies){
#   codon=row["codons"]
#   na=row["total_per_aa"]
#   ei=compute_ei(codon,nucleotides_frequencies)
#   pi=as.numeric(row['frequency'])
#   departure=na*((ei-pi)*(ei-pi)/ei)
#   if(pi==1){
#     return(1)
#   }
#   return(departure)
# }
# 
# compute_chi2_theoretical <- function(row,nucleotides_frequencies){
#   codon=row["codons"]
#   ei=compute_ei(codon,nucleotides_frequencies)
#   pi=as.numeric(row['uniform_prob'])
#   departure=(ei-pi)*(ei-pi)/ei
#   return(departure)
# }

compute_corrected_effective_number_of_codon_usage <- function(per_amino_acid_freq_df,nucleotide_freq){
  #This function take a frequency df and return the df with the computed correcte by genome composition effective number of codon usage
  
  nucleotide_freq=get_nucleotide_freq(genome_sequence)
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% filter(!amino_acid=="*") #remove stop codons
  per_amino_acid_freq_df$ei=apply(per_amino_acid_freq_df,1,compute_ei_apply,nucleotide_freq)
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% group_by(amino_acid) %>% mutate(chi2=total_per_aa*(sum(((frequency-ei)^2)/ei))) 
  
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% group_by(amino_acid) %>% mutate(F_prim=(chi2+total_per_aa-synonymous_fam)/(synonymous_fam*(total_per_aa-1)))
  per_amino_acid_freq_df <- per_amino_acid_freq_df %>% group_by(synonymous_fam) %>%  mutate(effective_number_of_codon_per_aa=1/(mean(F_prim)))
  return(per_amino_acid_freq_df)
}
