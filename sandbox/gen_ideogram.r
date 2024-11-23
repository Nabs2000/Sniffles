# # Load necessary libraries
# if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
#     install.packages("BiocManager")
#     BiocManager::install("VariantAnnotation")
# }
# if (!requireNamespace("karyoploteR", quietly = TRUE)) {
#     BiocManager::install("karyoploteR")
# }
# if (!requireNamespace("BSgenome.Hsapiens.UCSC.hs1", quietly = TRUE)) {
#     BiocManager::install("BSgenome.Hsapiens.UCSC.hs1")
# }

library(VariantAnnotation)
library(karyoploteR)
library(BSgenome.Hsapiens.UCSC.hs1)

# # Load your VCF file
vcf_file <- "C:/dev/Sniffles/ASPC-1_results/output.vcf"
vcf <- readVcf(vcf_file)

# Extract structural variant information
info <- info(vcf)
svtypes <- as.character(info$SVTYPE) # Extract SVTYPE field
chromosomes <- as.character(seqnames(rowRanges(vcf))) # Extract chromosome names
positions <- start(rowRanges(vcf)) # Extract start positions
ends <- positions + abs(info$SVLEN) # Extract end positions
strands <- info$STRAND # Extract strands

# Filter valid SVs
valid_indices <- !is.na(svtypes) & !is.na(positions) & !is.na(ends)
chromosomes <- chromosomes[valid_indices]
positions <- positions[valid_indices]
ends <- ends[valid_indices]
svtypes <- svtypes[valid_indices]
strands <- strands[valid_indices]

# Create a GRanges object for SVs
sv_ranges <- GRanges(
    seqnames = chromosomes,
    ranges = IRanges(start = positions, end = ends),
    svtype = svtypes
)

# Define color coding for each SV type
sv_colors <- ifelse(sv_ranges$svtype == "DEL", "red",
    ifelse(sv_ranges$svtype == "INS", "blue",
        ifelse(sv_ranges$svtype == "DUP", "yellow",
            ifelse(sv_ranges$svtype == "INV", "green", "gray")
        )
    )
) # Default to gray for others

# Sort the GRanges object by chromosome and start position
sv_ranges <- sv_ranges[order(seqnames(sv_ranges), start(sv_ranges))]

# Plot the karyotype with properly ordered SVs
output_file <- "C:/dev/Sniffles/ASPC-1_results/ideogram_with_SVs_corrected.png"
png(output_file, width = 1200, height = 600)

kp <- plotKaryotype(genome = "BSgenome.Hsapiens.UCSC.hs1", chromosomes = "all") # Use hs1 genome
# Add vertical lines at 5000 bp intervals for each chromosome
interval <- 5000000 # 5000 bp intervals
for (chr in unique(as.character(seqnames(sv_ranges)))) {
    chr_length <- seqlengths(BSgenome.Hsapiens.UCSC.hs1)[chr] # Get chromosome length
    print("Chromosome:")
    print(chr)
    print("Length:")
    print(chr_length)
    positions <- seq(0, chr_length, by = interval) # Generate positions at 5000 bp intervals
    kpSegments(kp, chr = chr, x0 = positions, x1 = positions, y0 = 0, y1 = 1, col = "black", lty = 2)
}

kpPlotRegions(kp, data = sv_ranges, col = sv_colors) # Add color-coded regions

# Add a title
kpAddMainTitle(kp, "Structural Variants in ASPC-1 Genome (Corrected)")

# Add a legend
legend_labels <- c("Deletion (DEL)", "Insertion (INS)", "Duplication (DUP)", "Inversion (INV)", "Other")
legend_colors <- c("red", "blue", "green", "purple", "gray")

legend("topright", legend = legend_labels, col = legend_colors, pch = 16, cex = 0.8, title = "SV Types")

dev.off()
cat("Karyotype plot with legend saved to:", output_file, "\n")
