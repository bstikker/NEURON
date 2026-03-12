# peroperative sequencing functions

library(QDNAseq)
library(DNAcopy)
#required files
#rds file containing a colorlist 
repo_root = normalizePath(file.path(getwd()), mustWork = TRUE)
rscriptslocation = file.path(repo_root, "scripts/")
color_translation = readRDS(paste0(rscriptslocation, "color_translation.rds"))
refgenome = file.path(repo_root, "reference/genome/chm13v2.0.fa.gz")
cnvplot_reffile = file.path(repo_root, "reference/bins_export.bed")
sturgeon_modkit_shellscript = paste0(rscriptslocation, "sturgeon_guppy_modkit_SLURM.sh")
modkitexec = "modkit"

#will be generated
cnvplot_tmpfile = file.path(repo_root, "bins_sample_counts_tmp.bed")
#DNAcopy bins
bins=readRDS(file.path(repo_root, "reference/bins/chm13_v2.0_bins.rds"))
relevant_genes_cnvplot = file.path(repo_root, "scripts/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed")

samtools = "modkit"
#functionchecks, also contain hardcoded paths, can be removed altogether. 
functioncheck = function(){
  if(exists("refgenome")==F){
    print("refgenome not selected, using chm13v2.0")
    refgenome = "/home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa"}
  mega_command = paste0("/home/sturgeon/miniconda3/envs/megalodon/bin/megalodon ~/nanocns/data/functioncheck/ --outputs mods basecalls mappings --mappings-format bam ",
                        "--reference ",refgenome, " --write-mods-text --mod-motif m CG 0 --processes 10 ",
                        "--guppy-server-path /home/sturgeon/nanocns/software/guppy_v5/bin/guppy_basecall_server ",
                        "--guppy-params \"-d /home/sturgeon/nanocns/software/rerio/basecall_models/ --num_callers 2 --ipc_threads 3\"",
                        " --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --devices cuda:0 --overwrite --output-directory ~/nanocns/data/functioncheck/meg_out/ ",
                        "--suppress-progress-bars --suppress-queues-status")
  
  system("rm -r ~/nanocns/data/functioncheck/sturgeonfiles/")
  system("rm -r ~/nanocns/data/functioncheck/meg_out/")
  system("rm ~/nanocns/data/functioncheck/barcoding.tsv")
  
  system(mega_command)
  
  system(paste0("/home/sturgeon/miniconda3/bin/qcat --tsv -k RBK004 -f ",
                "~/nanocns/data/functioncheck/meg_out/basecalls.fastq > ~/nanocns/data/functioncheck//barcoding.tsv"))
  
  system("mkdir ~/nanocns/data/functioncheck/sturgeonfiles")
  system("mv ~/nanocns/data/functioncheck/meg_out/per_read_modified_base_calls.txt ~/nanocns/data/functioncheck/sturgeonfiles/per_read_modified_base_calls.txt")
  system(paste0("~/Carlo/run_sturgeon_tobed_v2.sh ~/nanocns/data/functioncheck/sturgeonfiles/ ~/nanocns/data/functioncheck/sturgeonfiles/"))
  #system2('open', args = "~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.pdf", wait = FALSE)
  res = read.table("~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.csv", sep=",", header=T)
  print(res[1:2])
  system("rm -r ~/nanocns/data/functioncheck/sturgeonfiles/")
  system("rm -r ~/nanocns/data/functioncheck/meg_out/")
  system("rm ~/nanocns/data/functioncheck/barcoding.tsv")
  print("check complete, no errors")
  
}

fix_name=function(x){
  return(gsub(gsub(x = x, pattern = "\\.\\.\\.", replacement = " - "), pattern = "\\.", replacement = " "))}

add_and_plot=function(merged_data, result_file, outfile="x"){
  
  if(exists("color_translation")==F){color_translation=readRDS(file.path(repo_root, "scripts/color_translation.rds"))}
  reslist_2 = t(data.frame(read.table(result_file, sep=",")))
  reslist_2=reslist_2[2:nrow(reslist_2),]
  colnames(reslist_2)=c("class", "score")
  tm = Sys.time()
  tm= paste(format(as.POSIXct(tm), format = "%H:%M"))
  reslist_2=rbind(reslist_2,data.frame(class="TIME", score=tm))
  
  if(nrow(merged_data)==0){
    return(reslist_2)
  }else{
    
    mgd = merge(merged_data, reslist_2, by="class") 
    colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
    
    plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
         xaxt="n") 
    for(i in 1:(nrow(mgd)-1)){
      clr = unlist(unname(color_translation[mgd[i,"class"]]))
      lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
      if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
    }
    maxtics= ifelse(ncol(mgd)>5,ncol(mgd), ncol(mgd)+5)
    axis(side = 1, at = 0:maxtics, labels = 1:(maxtics+1))
    axis(side = 1, at = 0:(ncol(mgd)-2), labels = mgd[mgd$class=="TIME",2:ncol(mgd)], line=1, tick=F, cex.axis=0.7)
    abline(h=0.95, col="red")
    abline(h=0.80, col="orange")
    return(mgd)
  }
}

confidence_over_time_plot <- function(merged_data, color_translation = NULL) {
  # Ensure color_translation is loaded if not provided
  if (is.null(color_translation)) {
    color_translation <- readRDS(file.path(repo_root, "scripts/color_translation.rds"))
  }
  
  mgd <- merged_data
  colnames(mgd)[ncol(mgd)] <- paste0("iteration_", ncol(mgd) - 1)
  
  # Create the base plot
  plot(
    '', xlim = c(0, ifelse(ncol(mgd) > 5, ncol(mgd) + 5, 10)), ylim = c(0, 1),
    main = "Confidence Over Time", xlab = "Iteration", ylab = "Confidence",
    xaxt = "n"
  )
  
  # Add lines and labels for each class
  for (i in 1:(nrow(mgd) - 1)) {
    clr <- unlist(unname(color_translation[mgd[i, "class"]]))
    lines(x = 0:(ncol(mgd) - 2), y = as.numeric(mgd[i, 2:ncol(mgd)]), col = clr)
    if (as.numeric(mgd[i, ncol(mgd)]) > 0.5) {
      text(x = ncol(mgd), y = as.numeric(mgd[i, ncol(mgd)]) - 0.02, labels = mgd[i, "class"])
    }
  }
  
  # Add axes and confidence thresholds
  maxtics <- ifelse(ncol(mgd) > 5, ncol(mgd), ncol(mgd) + 5)
  axis(side = 1, at = 0:maxtics, labels = 1:(maxtics + 1))
  axis(
    side = 1, at = 0:(ncol(mgd) - 2), 
    labels = mgd[mgd$class == "TIME", 2:ncol(mgd)], 
    line = 1, tick = FALSE, cex.axis = 0.7
  )
  abline(h = 0.95, col = "red")
  abline(h = 0.80, col = "orange")
  #return(mgd)
}


#still contains hardcoded paths
wrappert_guppy_R10_guppy6.5 = function(main_folder,fast5, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  guppy_command = paste0( guppy_location, "-r -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                          "-c ~/nanocns/software/guppy_v6.5.7/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_mk1c.cfg --device auto --bam_out ",
                          "--disable_qscore_filtering --min_score_barcode_front 6 ",
                          "--align_ref ", refgenome, " --barcode_kits SQK-RBK114-24 ",
                          "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  
  system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))
  
  #include unclassified reads yes or no?
  if(include_unclassified==T){
    system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_it",iteration,"_unclassified.bam"))}
  
  system(paste(sturgeon_shellscript,main_folder))
  
  
}

wrappert_guppy_R10_guppy6.5_V2 = function(main_folder,fast5_cpy, fast5_folder, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  for(fl in fast5_cpy){
    fullname = paste0(fast5_folder, fl)
    system(paste0("cp ", fullname," ", out_folder))
  }
  
  guppy_command = paste0( guppy_location, " -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                          "-c ~/nanocns/software/guppy_v6.5.7/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_mk1c.cfg --device auto --bam_out ",
                          "--disable_qscore_filtering --min_score_barcode_front 6 ",
                          "--align_ref ", refgenome, " --barcode_kits SQK-RBK114-24 ",
                          "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  
  system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  
  if(length(fast5_cpy)>1){
    mergecom = paste0("samtools merge -@ 10 ", main_folder,"/guppy_output_it",iteration,".bam ",out_folder,"/guppy_out/barcode",barcode,"/*.bam")
    system(mergecom)
  }else{
    system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))}
  
  #include unclassified reads yes or no?
  if(include_unclassified==T){
    if(length(fast5_cpy)>1){
      mergecom_u = paste0("samtools merge -@ 10 ", main_folder,"/guppy_output_it",iteration,"_unclassified.bam ",
                          out_folder,"/guppy_out/unclassified/*.bam")
      system(mergecom_u)
    }else{
      system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_it",iteration,"_unclassified.bam"))}
  }
  
  system(paste(sturgeon_shellscript,main_folder))
  
}

wrappert_guppy_R10_guppy6.5_V2_modkit = function(main_folder,fast5_cpy, fast5_folder, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  for(fl in fast5_cpy){
    fullname = paste0(fast5_folder, fl)
    system(paste0("cp ", fullname," ", out_folder))
  }
  
  guppy_command = paste0( guppy_location, " -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                          "-c ~/nanocns/software/guppy_v6.5.7/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_mk1c.cfg --device auto --bam_out ",
                          "--disable_qscore_filtering --min_score_barcode_front 6 ",
                          "--align_ref ", refgenome, " --barcode_kits SQK-RBK114-24 ",
                          "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  
  system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  
  if(length(fast5_cpy)>1){
    mergecom = paste0("samtools merge -@ 10 ", out_folder,"/merged.bam ",out_folder,"/guppy_out/barcode",barcode,"/*.bam")
    system(mergecom)
  }else{
    system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",out_folder,"/merged.bam"))}
  
  #include unclassified reads yes or no?
  if(include_unclassified==T){
    if(length(fast5_cpy)>1){
      mergecom_u = paste0("samtools merge -@ 10 ", out_folder,"/merged_unclassified.bam ",
                          out_folder,"/guppy_out/unclassified/*.bam")
      system(mergecom_u)
    }else{
      system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",out_folder,"/merged_unclassified.bam"))}
    
    system(paste0("samtools merge -@ 10 ",out_folder, "/merged_both.bam ", out_folder,"/*.bam" ))
    system(paste0("rm ", out_folder,"/merged_unclassified.bam"))
    system(paste0("rm ", out_folder,"/merged.bam"))
    system(paste0("mv ", out_folder,"/merged_both.bam ",out_folder,"/merged.bam"))
  }
  
  system(paste0("samtools sort -@ 10 ",out_folder,"/merged.bam -o ",out_folder,"/merged_srt.bam"))
  system(paste0("rm ", out_folder,"/merged.bam"))
  modkit_extract = (paste0(modkitexec, " extract -t 10 ",out_folder,"/merged_srt.bam ", main_folder, "/iteration_", iteration, "_modkit.txt"))
  system(modkit_extract)
  
  system(paste(sturgeon_modkit_shellscript,main_folder))
  
}

plot_cnv_from_bam_DNAcopy = function(bam, makeplot=T, lines_only=F, binsize=1e6){
  readCounts <- binReadCounts(bins, bamfiles=bam, isNotPassingQualityControls=NA, minMapq=2, isDuplicate=NA, isSecondaryAlignment=NA)
  exportBins(readCounts, file=cnvplot_tmpfile, format="bed", filter=F, logTransform=F)
  reference=read.table(cnvplot_reffile, skip = 1)
  colnames(reference) = c("chrom", "start", "end", "name", "coverage", "orientation")
  sample = read.table(cnvplot_tmpfile, skip = 1)
  system(paste0("rm ",cnvplot_tmpfile))
  colnames(sample) = c("chrom", "start", "end", "name", "coverage", "orientation")
  totreads_title = sum(sample$coverage)
  #print(totreads_title)
  
  sample$coverage = as.numeric(sample$coverage)+0.001
  totreads = sum(sample$coverage)
  sample$refcov = reference$coverage
  sample = sample[sample$chrom!="Y"&sample$chrom!="X",]
  sample$reltoref = sample$coverage/sample$refcov
  ##sample=na.omit(sample)
  blacklistbins = sample[sample$refcov<500,]
  blacklistbins = row.names(blacklistbins)
  sample = sample[!row.names(sample)%in%blacklistbins,]
  #sample=sample[sample$coverage>0,]
  mn = mean(sample$reltoref)
  totreads=sum(sample$coverage)
  sample$logtr = log(sample$reltoref/mn, 2)
  sample$logtr=ifelse(sample$logtr>3,3, ifelse(sample$logtr<(-3),-3, sample$logtr))
  sample$chromname=ifelse(nchar(sample$chrom)==1, paste0("chr0", sample$chrom), paste0("chr", sample$chrom))
  bam_cna = CNA(sample$logtr, chrom=sample$chromname, maploc=sample$start+binsize, data.type="logratio", sampleid = "cnvplot")
  smoothed.CNA.object <- smooth.CNA(bam_cna)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  #plot(segment.smoothed.CNA.object, plot.type="w")
  pointcolorframe = data.frame(pointx=1:length(bam_cna$cnvplot), pointy=bam_cna$cnvplot, pointcolor="black")
  for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]>0.5){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="green"
    }
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]<(-0.5)){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="blue"
    }    }
  if(makeplot){
    plot("", xlim=c(0,length(bam_cna$cnvplot)), ylim=c(-3,3), main=paste0("CNVs nreads=",totreads_title), xaxt="n", xlab="genomic pos",
         ylab="log ratio")
    if(lines_only==F){
      points(x=pointcolorframe$pointx, y=pointcolorframe$pointy, pch=16, cex=0.5, col=pointcolorframe$pointcolor)
    }
    for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
      segments(x0=  segment.smoothed.CNA.object$segRows[x,1], x1=segment.smoothed.CNA.object$segRows[x,2],
               y0=segment.smoothed.CNA.object$output[x,"seg.mean"], y1=segment.smoothed.CNA.object$output[x,"seg.mean"],
               col="red", lwd=2)
    }
    fr = data.frame(chrom=smoothed.CNA.object$chrom, loc=smoothed.CNA.object$maploc, val=smoothed.CNA.object$cnvplot)
    nr = 1
    for(i in 1:(nrow(fr)-1)){
      if(fr$chrom[i]!=fr$chrom[i+1]){
        abline(v=as.numeric(row.names(fr)[i])+0.5, col="grey")
        chrname= gsub(pattern = "chr", replacement = "", x = fr[i,"chrom"])
        ypos = ifelse(nr %% 2 == 1, 2, 1.8)
        nr=nr+1
        text(x = i, adj=c(1,0), labels = chrname, y = ypos)
      }
    }
    if(lines_only==F){
      relevant_genes = read.table(relevant_genes_cnvplot, header=T)
      fr$start=fr$loc-binsize
      for(i in 1:nrow(relevant_genes)){
        genestart = relevant_genes[i,"start"]
        genechrom = relevant_genes[i,"chrom"]
        genechrom = strsplit(genechrom, split="hr")[[1]][2]
        genechrom=ifelse(nchar(genechrom)==1, paste0("chr0", genechrom), paste0("chr", genechrom))
        binline = fr[fr$chrom==genechrom&fr$start<genestart&(fr$start+binsize*2)>genestart,]
        ycoord = binline$val
        text(x=as.numeric(row.names(binline)), y=ycoord-0.1+relevant_genes[i,"yoffset"], cex=0.5, srt=-90,adj=c(0.5,1), labels = relevant_genes[i, "name"], pos=1)
        points(x=as.numeric(row.names(binline)), y=ycoord, col="red")
      }}
  }
  output = list(pointdata=pointcolorframe, segdata=segment.smoothed.CNA.object)
  return(output)
}

plot_cnv_from_live_dir_DNAcopy=function(directory){
 #directory = "~/Documents/Experiments/BS-001/test"
  
  if(dir.exists(paste0(directory,"/merged_bams/"))==F){system(paste0("mkdir ", directory,"/merged_bams/"))}else{
    system(paste0("rm ", directory,"/merged_bams/merged_bam.bam"))
  }
  system(paste0(samtools, " merge -@ 10 -O BAM -o ",directory,"/merged_bams/merged_bam.bam ", directory,"/*.bam" ))
  out=plot_cnv_from_bam_DNAcopy(paste0(directory,"/merged_bams/merged_bam.bam"))
}

