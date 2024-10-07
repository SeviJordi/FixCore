#!/usr/bin/env Rscript

## Load libraries
library("tidyr")
library("dplyr")
library("ape")
library("phylobase")
library("stringr")
library("phangorn")
library("logger")
log_threshold(INFO)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

## Start the dataframe with the isolates to remove
removedf= data.frame(isolate=character(0))


# The algorithm only works if there are more than 3 variants in the VCF
nlines <- readLines(snakemake@input[["vcf"]])[!startsWith(readLines(snakemake@input[["vcf"]]), "#")] %>% length()
if (nlines >= 3){

    ## Get the length of the alignment, which is in the VCF
    msalen = grep("length", readLines(snakemake@input[["vcf"]]), value = T) %>% 
        str_remove(., ".*=") %>% str_remove(., ">") %>% as.integer()
    
    ## Read the VCF file
    vcf <- read.csv(
        snakemake@input[["vcf"]],
        sep = "\t", 
        skip = 3,
        colClasses = c("REF"="character", "ALT"="character") # to avoid reading F and T as booleans
        ) %>%    
        select(
            -c(X.CHROM,QUAL,ID,FILTER,INFO,FORMAT) # Remove unnecessary columns
            ) %>%
        # Change the numbers by the amino acid that each sample has
        mutate(
            across(
                .cols = c(ALT),
                ~str_remove_all( ., ",")
                )
            ) %>%
        mutate(
            across(
                .cols = -c(POS, REF, ALT),
                ~str_sub(ALT, .,.) 
                )
            ) %>% 
        mutate(
            across(
                .cols = -c(POS, REF, ALT),
                ~str_replace_all( ., "^$", REF )
                )
            ) %>%
        select(-c(REF)) %>% # change the reference created by snp-sites for the emboss
        rename(REF = "CONSENSUS")
    
    ## Transform the VCF from matrix to molten table
    vcf_list <- vcf %>%
        select(-ALT) %>%
        pivot_longer(
            cols = -c( POS, REF),
            names_to = "isolate",
            values_to = "ALT"
            ) 
    
    ## First remove all samples that have 15% of gaps
    tempo <- vcf_list %>%  
        # Remove the positions that are common with the REF
        filter(!(REF == "X" & ALT == "X")) %>%
        filter(ALT == "X") %>% 
        group_by(isolate) %>% 
        count() %>%
        filter(n > (0.10 * msalen)) %>% 
        select(isolate)
    
    removedf = rbind(removedf, tempo)
    log_info("ISolates with more than 15% o gaps removed")

    ## Check if there are any samples with more than 3 consecutive variant amino acids
    testnvcf <- vcf_list %>% 
        filter(REF != ALT) %>%
        group_by(isolate) %>%
        mutate(
            GROUP = paste(isolate,cumsum(c(1L, diff(POS) > 1)), sep="-")
            ) %>%
        group_by(isolate, GROUP) %>%
        filter(
            any(REF != "X") & any(ALT != "X") # Remove all groups that are gaps
            ) %>% 
        count() %>%
        filter(n >= 3 ) %>%
        group_by(isolate) %>%
        distinct(isolate)

    # If there are more than 3 variant sites in a row in the VCF, the tree is evaluated
    if (nrow(testnvcf) > 0){
        log_info("There are samples with more than 3 consecutive variant amino acids")

        # Read the fasta
        fasta <- read.FASTA(snakemake@input[["aligned"]])

        # Make a NJ tree and remove the isolates with more than 15% of gaps
        tree <- dist.dna(fasta, pairwise.deletion = T) %>%
            replace(is.na(.), 0) %>%
            replace(is.infinite(.), 1) %>%
            nj() %>%
            drop.tip(pull(removedf, isolate))
        
        L <- list(seq(0, tree$Nnode)) # Create node names

        tr <- makeNodeLabel(tree, "n", nodeList = L)  # Give names to nodes to convert to phylo4
        tr$edge.length[tr$edge.length < 0] <- 0 # replace negative distances by 0
  
        # Convert tree to phylo4 S4 class
        tr4df <- as(tr, "phylo4") %>%
            as("data.frame") %>% # Extract dataframe
            filter(!is.na(edge.length)) # Remove root


        # Standard deviation for terminal branches 
        sdmax_tip <- tr4df %>%  # extract 3 SD from terminal branches
            filter(node.type == "tip") %>% 
            summarise(
                sd = mean(edge.length)+(3*sd(edge.length))
            ) %>% 
            pull(sd)

        log_info(sprintf("Standard deviation for terminal branches: %f", sdmax_tip))

        # Standard deviation for internal branches
        sdmax_int <- tr4df %>% # extract 3 SD from internal branches
            filter(node.type == "internal") %>% 
            summarise(
                sd = mean(edge.length)+(3*sd(edge.length))
                ) %>%
            pull(sd)

        log_info(sprintf("Standard deviation for internal branches: %f", sdmax_int))

        # Search for outliers
        outliers_tips <- tr4df  %>% # select teminal outliers --> those nodes and tips above the 3SD 
            filter((node.type == "tip" & edge.length > sdmax_tip)) %>%
            mutate(Filter = "Test", Outlier = node)
  
        outliers_int <- tr4df  %>% # select internal outliers --> those nodes and tips above the 3SD 
            filter((node.type == "internal" & edge.length > sdmax_int))  
  

        # Evaluation of internal nodes
        nodes <- outliers_int %>% pull(node) # extract all nodes of oulier_int
  
        nodestocheck = data.frame(
            label=character(0),
            node=integer(0),
            ancestor=integer(0), 
            edge.length=numeric(0),
            node.type=character(0),
            Outlier=integer(0),
            Filter=character(0)
        ) 

        # We want to asses if some nodes are inside other nodes 
        # to remove common sites of VCF evaluation
        for (i in nodes){ # for each node
            ntips <- phangorn::Descendants(tr, i, "tips") # we obtain the nodes descendants tips
    
            temp <- (tr4df %>% filter(node %in% ntips[[1]]) %>% # filter the tree df to keep only descendants
                mutate(Outlier = !!i))  # add column with node name
    
            nodestocheck = rbind(nodestocheck, temp)
            }


        # check if any descendant tip of the node has to be tested
        nodestocheck_def <- nodestocheck %>%
            mutate(
            # Evaluate if there is a tip of node within the node that needs to be tested 
            Filter = case_when(
                node.type == "internal" & edge.length > sdmax_int  ~ "Test",
                node.type == "tip" & edge.length > sdmax_tip       ~ "Test",
                TRUE                                               ~ "Pass"
                )
            ) %>%
        full_join(., outliers_tips)  %>% # add outlier tips
        group_by(label) %>% 
        slice(which.max(Outlier)) %>% # remove duplicate tips if in two nodes
        group_by(Outlier) %>% # count tips in nodes
        add_count() 

        if (nrow(nodestocheck_def) != 0){ # If there are any samples to test
            
            # evaluate the consecutive amino acids of each tip of the nodes    
            outliersdef <- nodestocheck_def %>%
            select(Outlier) %>%
            distinct() %>%
            pull() 

            for (outli in outliersdef) { # for each node/tip

                teting_nodedf <- nodestocheck_def %>%
                filter( outli == Outlier)
            
                df5 <- vcf_list %>% 
                    # Merge the VCF with the nodes to test
                    left_join(., teting_nodedf, by=c("isolate"="label")) %>%
                    filter(REF != ALT) %>% # optimizacion
                    # Remove common sites
                    group_by(Outlier) %>%
                    distinct(POS, REF, ALT, .keep_all = T) %>%
                    # Keep only the sites to test
                    filter(Filter == "Test") %>% 
                    group_by(isolate) %>%
                    mutate(
                        # diff(POS) returns the difference between POS X and X-1 (because the default lag is 1, it can be changed)
                        # diff(POS) > 1 returns TRUE/FALSE if between one position and another there is more than 1
                        # the 1L of c(1L, diff(POS) > 1), converts the T/F to binary
                        # and cumsum sums cumulatively and as it is binary consecutive groups remain.
                        GROUP = paste(isolate,cumsum(c(1L, diff(POS) > 1)), sep="-")
                        ) %>%
                    group_by(isolate, GROUP) %>%
                    filter(any(REF != "X") & any(ALT != "X") ) %>% # Remove all groups that are gaps
                    filter(n() >= 3 ) %>%
                    tally(ALT != "X") %>%
                    filter(n >= 3 ) %>%
                    group_by(isolate) %>%
                    distinct(isolate) %>% 
                    filter(isolate != "")
            
                removedf <- rbind(removedf, df5)
            }
        }
    }




}

# Write the isolates to remove
write.table(removedf, snakemake@output[["removed"]],  row.names = F,  col.names=F, quote = F)
log_info(sprintf("Number of isolates to remove: %s", length(removedf$isolate)))

# Subset the alignment
aln <- read.FASTA(snakemake@input[["aligned"]])
aln <- aln[!(names(aln) %in% removedf$isolate)]
write.FASTA(aln, snakemake@output[["filtered"]])

# End of the script

