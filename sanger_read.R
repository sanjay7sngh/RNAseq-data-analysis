## analysis of sanger sequence data
## for only one end sequence
## sangerRead https://sangeranalyser.readthedocs.io/en/latest/content/

## Preparing SangerRead AB1 input
inputFilesPath <- system.file("sequencingdata/")
A_chloroticaFFN <- file.path(inputFilesPath,
                             "Allolobophora_chlorotica",
                             "ACHLO",
                             "")

## Creating SangerRead instance from AB1
# using `constructor` function to create SangerRead instance
sangerReadF <- SangerRead(readFeature           = "Forward Read",
                          readFileName          = A_chloroticaFFN)

# using `new` method to create SangerRead instance
sangerReadF <- new("SangerRead",
                   readFeature           = "Forward Read",
                   readFileName          = A_chloroticaFFN)

# Visualizing SangerRead trimmed read
qualityBasePlot(sangerReadF)

# Writing SangerRead FASTA file (AB1)
writeFasta(sangerReadF)

#  Generating SangerRead report (AB1)
generateReport(sangerReadF)