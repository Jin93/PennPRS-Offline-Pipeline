match <- Filter(length, lapply(1:length(LD), function(i){
  snp.overlap <- intersect(gwas$SNP,rs[[i]])
  if (length(snp.overlap)==0) {
    return()
  } else {
    if (length(snp.overlap)==1){
      ld_blk <- as.matrix(1) # as.matrix(LD[[i]][match(snp.overlap,rs[[i]]),match(snp.overlap,rs[[i]])])
    } else{
      ld_blk <- LD[[i]][match(snp.overlap,rs[[i]]),match(snp.overlap,rs[[i]])]
    }
    rs_blk <- rs[[i]][match(snp.overlap,rs[[i]])]
    gwas <- gwas[match(snp.overlap,gwas$SNP),]
    ld_blk_size <- length(snp.overlap)
    return(list(ld_blk=ld_blk, rs_blk=rs_blk, gwas=gwas, ld_blk_size=ld_blk_size))
  }
}))



pumas_tmp <- PUMAS.II.FUN.I(FunII.GWAS=matched_data$gwas_matched, 
                            FunII.Nt=floor(partitions[2]*min(matched_data$gwas_matched$N)),
                            FunII.LD=matched_data$new_ld_blocks, 
                            FunII.LD.block=matched_data$ld_block_size, 
                            trait_name=trait_name, k=k, chr=chr)

FunII.GWAS=matched_data$gwas_matched
FunII.Nt=floor(partitions[2]*min(matched_data$gwas_matched$N))
FunII.LD=matched_data$new_ld_blocks
FunII.LD.block=matched_data$ld_block_size 

