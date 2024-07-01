###### Usage of the GENOVA package for 3D genome analysis as carried on in Rabuffo et al.


library(GENOVA)
library(scales)
library(viridis)
library(ggplot2)

############################################
######         load matrices           #####
############################################

centromeres <- read.delim("/path/to/centromere_file.bed", header = FALSE)
microc_contact_file_2kb <- load_contacts(signal_path = '/path/to/contact_microc_file.mcool',
                                     sample_name = "Micro-C",
                                     balancing = T,
                                     resolution = 2e3,
                                     centromeres = centromeres,
                                     colour = "#528e9e")

hic_contact_file_2kb <- load_contacts(signal_path = '/path/to/contact_hic_file.mcool',
                                     sample_name = "Hi-C",
                                     balancing = T,
                                     resolution = 2e3,
                                     centromeres = centromeres,
                                     colour = "#242466")


############################################
###### APA on loops called by Mustache #####
############################################
loop_file <- read.csv("/path/to/loop_file.tsv", header = T, sep="")

APA_microc_hic  <- APA(list("Micro-C" = microc_contact_file_2kb,
                        'Hi-C' = hic_contact_file_2kb),
                        bedpe = loop_file,
                        dist_thres = c(200e3, Inf),
                        size_bin = 41)

visualise(APA_microc_hic, 
    title = "Micro-C vs Hi-C loops",
    colour_lim = c(0, 40),
    colour_lim_contrast = c(-5, 5),
    metric = "diff",
    contrast = 1,
    raw = TRUE) + scale_fill_continuous(high = "#242466", low = "#ECECFF", limits=c(0, 30))


############################################
####### CSCAN on annotations by type #######
############################################

###### load annotations

sTTSs <- read.csv("/path/to/annotation/file_sTTSs.bed", header = F, sep ="", row.names=NULL)
cTTSs <- read.csv("/path/to/annotation/file_cTTSs.bed", header = F, sep ="", row.names=NULL)
dTSSs <- read.csv("/path/to/annotation/file_dTSSs.bed", header = F, sep ="", row.names=NULL)
sTSSs <- read.csv("/path/to/annotation/file_sTSSs.bed", header = F, sep ="", row.names=NULL)

sTTS <- sTTSs[,c(1,2,3,1,2,3)]
colnames(sTTS) <- c("BIN1_CHR","BIN1_START","BIN1_END",'BIN2_CHR',"BIN2_START","BIN2_END")
cTTS <- cTTSs[,c(1,2,3,1,2,3)]
colnames(cTTS) <- c("BIN1_CHR","BIN1_START","BIN1_END",'BIN2_CHR',"BIN2_START","BIN2_END")
dTSS <- dTSSs[,c(1,2,3,1,2,3)]
colnames(dTSS) <- c("BIN1_CHR","BIN1_START","BIN1_END",'BIN2_CHR',"BIN2_START","BIN2_END")
sTSS <- sTSSs[,c(1,2,3,1,2,3)]
colnames(sTSS) <- c("BIN1_CHR","BIN1_START","BIN1_END",'BIN2_CHR',"BIN2_START","BIN2_END")

###### Analysis and visualization

CScan_dTSS <- CSCAn(
  list("MicroC" = microc_contact_file_2kb),
  bedlist = list(dTSS, dTSS),
  shift = 15000L,
  dist_thres = c(50e3, Inf)
)
visualise(CScan_dTSS, 
    title = "Micro-C at dTSSs",
    metric = "obsexp",
    colour_lim = c(0, 3))

CScan_sTSS <- CSCAn(
  list("MicroC" = microc_contact_file_2kb),
  bedlist = list(sTSS, sTSS),
  shift = 15000L,
  dist_thres = c(50e3, Inf)
)
visualise(CScan_sTSS, 
    title = "Micro-C at sTSSs",
    metric = "obsexp",
    colour_lim = c(0, 3))

CScan_cTTS <- CSCAn(
  list("MicroC" = microc_contact_file_2kb),
  bedlist = list(cTTS, cTTS),
  shift = 15000L,
  dist_thres = c(50e3, Inf)
)

visualise(CScan_cTTS, 
    title = "Micro-C at cTTSs",
    metric = "obsexp",
    colour_lim = c(0, 3))


CScan_sTTS <- CSCAn(
  list("MicroC" = microc_contact_file_2kb),
  bedlist = list(sTTS, sTTS),
  shift = 15000L,
  dist_thres = c(50e3, Inf)
)

visualise(CScan_sTTS, 
    title = "Micro-C at sTTS",
    metric = "obsexp",
    colour_lim = c(0, 3))




############################################
#######        cis trans ratio       #######
############################################


####### load matrices
human <- load_contacts(signal_path = './GSE163625_hg38_wildtype.hic',
                                     sample_name = "Homo sapiens",
                                     balancing = T,
                                     resolution = 5e4)

wheat <- load_contacts(signal_path = './GSM5182739_161010_Chinese_Spring_v1.0_pseudomolecules.hic',
                                     sample_name = "Triticum aestivum",
                                     balancing = T,
                                     resolution = 5e4)

mushroom <- load_contacts(signal_path = './GSM5182736_Agabi_varbisH97_2_HiC.hic',
                                     sample_name = "Agaricus bisporus",
                                     balancing = T,
                                     resolution = 5e4)

tardigrade <- load_contacts(signal_path = './GSM5182729_nHd_3.1_HiC.hic',
                                     sample_name = "Hypsibius dujardini",
                                     balancing = T,
                                     resolution = 5e4)

yeast <- load_contacts(signal_path = './GSM5182737_sacCer3.hic',
                                     sample_name = "Saccharomyces cerevisiae",
                                     balancing = T,
                                     resolution = 5e4)

microc_contact_file_50kb <- load_contacts(signal_path = '/path/to/contact_microc_file.mcool',
                                     sample_name = "Micro-C",
                                     balancing = T,
                                     resolution = 5e4)

hic_contact_file_50kb <- load_contacts(signal_path = '/path/to/contact_hic_file.mcool',
                                     sample_name = "Hi-C",
                                     balancing = T,
                                     resolution = 5e4)

cisChrom_out <- cis_trans(list(tryps, tryps_micro, mushroom, wheat, yeast, tardigrade, human))
barplot(cisChrom_out$cis, names.arg = cisChrom_out$sample, ylim = c(0, 100), ylab="Percentage of cis contacts in Hi-C data", las=2, 
    col = c("#220030", "#010048", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"))




