#!/usr/bin/Rscript
library(treemapify)
library(ggplot2)
library(argparser)

#usage: Rscript TreeMap.R --input /home/wuzhikun/Project/Vigna/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.txt --pdf /home/wuzhikun/Project/Vigna/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.pdf --width 5 --height 4

#usage: Rscript ~/github/NanoHub/script/SVTypeBox.R --input  Samples_summary_SV.xls --outPrefix temp

arg <- arg_parser('Tree Map for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



Length_treeMap <- function(infile, outFile, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)


    data <- read.table(infile, header = FALSE, sep = "\t")
    colnames(data) <- c("Contig", "Length")
    print(head(data))
    data$Length <- data$Length / 1000000

    ### hist count for SV type and number
    # ggplot(data, aes(area = Length, fill = Length, label = Contig)) + 
    #     geom_treemap() +
    #     geom_treemap_text(colour = "black", place = "centre", size = 10) +
    #     labs(fill='Length (Mb)') 

    ggplot(data, aes(area = Length, fill = Length)) + 
        geom_treemap() +
        labs(fill='Length (Mb)') 

 
    ggsave(outFile, width=width, height=height)


    


}


Length_treeMap(argv$input, argv$pdf, argv$width, argv$height)