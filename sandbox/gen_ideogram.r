# Install karyoploteR
if (!requireNamespace("karyoploteR", quietly = TRUE)) {
    BiocManager::install("karyoploteR")
}

# Install VariantAnnotation for parsing VCF
if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    BiocManager::install("VariantAnnotation")
}

# Load libraries
library(BiocManager)
# remove.packages("BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0")
install("BSgenome.Hsapiens.UCSC.hs1")
library(karyoploteR)
library(VariantAnnotation)

# Load your VCF file
vcf_file <- "C:/dev/Sniffles/ASPC-1_results/output.vcf"
vcf <- readVcf(vcf_file)

# Extract the necessary fields
chrom <- as.character(seqnames(rowRanges(vcf))) # Chromosomes
# Remove "chr" prefix and capitalize "x" or "y"
# chrom <- gsub("^chr", "", chrom, ignore.case = TRUE)
# chrom <- toupper(chrom) # Ensure "x" becomes "X", etc.

start <- start(rowRanges(vcf)) # Start positions
info <- info(vcf) # INFO fields

# Extract SVTYPE and SVLEN
svtype <- info$SVTYPE
svlen <- info$SVLEN

# Create a data frame for the variants
variants <- data.frame(
    chrom = chrom,
    start = start,
    svtype = svtype,
    svlen = svlen,
    stringsAsFactors = FALSE
)

# Save the plot to a PNG file
png("C:/dev/Sniffles/ASPC-1_results/ideogram.png", width = 1200, height = 600)

# Initialize the karyotype plot
kp <- plotKaryotype(genome = "hs1")

# Add structural variants
colors <- c(DEL = "red", INS = "blue", DUP = "orange", INV = "green")
for (type in unique(variants$svtype)) {
    subset_variants <- variants[variants$svtype == type, ]

    # Calculate end positions for each variant
    end_positions <- subset_variants$start + subset_variants$svlen

    # Use `kpSegments()` to overlay lines directly on chromosomes
    kpSegments(
        kp,
        chr = subset_variants$chrom,
        x0 = subset_variants$start, # Start position
        x1 = end_positions, # End position
        y0 = 0, # Align with chromosome
        y1 = 0, # Align with chromosome
        col = colors[type], # Color based on SVTYPE
        lwd = 2 # Line width
    )
}



# Add legend, chromosome names, and title
legend("topright", legend = names(colors), col = colors, pch = 16, cex = 0.8, title = "SV Types")
kpAddChromosomeNames(kp)
kpAddMainTitle(kp, "My Karyotype Plot")

# Close the device
dev.off()

# Retrieve genome information
# genome_data <- getGenomeAndMask(genome = "hs1")

# # Extract chromosome names and lengths
# chromosome_lengths <- as.data.frame(seqinfo(genome_data$genome))

# # Print the lengths of the chromosomes
# print(chromosome_lengths)
