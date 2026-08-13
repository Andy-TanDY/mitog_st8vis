Mitog Str8vis: A Linearised Mitogenome Visualizer

A lightweight, Base R script for generating schematic plots of linear mitogenomes. 
This tool takes annotated mitogenome data (such as protein-coding genes, rRNAs, tRNAs, Transposable Elements (TEs), and Open Reading Frames (ORFs)) and visualizes them as track-based, directional polygons.  

✨ Features
1. Base R Plotting: No heavy ggplot2 dependencies required for the core rendering
2. Directional Polygons: Automatically detects strand direction (forward/reverse) based on Start and Stop coordinates.  
3. Multi-Track Support: Plots elements on specific Y-axis tracks based on your input data. 
4. Highly Customizable: Easily modify color palettes (warna vectors), axis intervals, and species labels directly in the script.  

📦 Prerequisites
Before running the script, ensure you have R installed along with the following libraries:  
Rinstall.packages(c("MetBrewer", "scales", "extrafont"))
Note: The script uses the "Seaford" and "Gill Sans MT" Windows fonts via the extrafont package. You may need to import your system fonts first using font_import().  Feel free to change to your own preferred font. 

📂 Input Data Format
The script reads .csv files containing mitogenome annotations. (see included example -- Ffra_schematic.csv)
Your input CSVs must contain at least the following columns based on the script's logic:  
1. Type: The classification of the element (e.g., tRNA, rRNA, ORF) used for color-coding and the legend.  
2. Start: The starting base-pair position.  
3. Stop: The ending base-pair position.  
4. Y1: The Y-axis track number where the element should be plotted.  

🚀 Usage 
1. Clone the repository and open the R script.
2. Set your working directory: Update the setwd() path at the top of the script to point to your local project folder.
3. Update file paths: The script currently processes specific species sequentially (e.g., Fragum fragum, Tridacna gigas, Cerastoderma edule). 
4. Replace the read.csv() paths with the paths to your own annotation files.  
5. Adjust custom parameters:
    warna: Change the hex color codes to match the number of unique Type elements in your dataset.  
    text(...): Update the species name and plot title for each section.  
    frame_interval: Adjust the X-axis tick intervals (default is 2000).  
6. Run the section of scripts for plotting : The output will generate directly in your R graphical device.  



##########################################

Included here is also concatenated supermatrix of 15 mitochondrial genes for 12 cardiid species, partitioned as follows:

#nexus
begin sets;
    charset ATP6 = 1-1217;
    charset ATP8 = 1218-1331;
    charset COX1 = 1332-2890;
    charset COX2 = 2891-3579;
    charset COX3 = 3580-4410;
    charset NAD1 = 4411-5320;
    charset NAD2 = 5321-6291;
    charset NAD3 = 6292-6648;
    charset NAD4 = 6649-8007;
    charset NAD4L = 8008-8283;
    charset NAD5 = 8284-9965;
    charset NAD6 = 9966-10447;
    charset rrnL = 10448-11515;
    charset rrnS = 11516-12185;
end;

supermatrix was used to reconstruct IQtree and RAxML phylogenies. 
