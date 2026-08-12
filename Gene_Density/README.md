# Gene density distribution 

This script can be used for visualizing gene-density distributions across multiple genomes using a ridgeline (joyplot) representation. The script reads genome-wide **gene-density** values from a directory containing one tab-delimited **text file** per species, estimates the density distribution using a **Gaussian kernel density estimator (KDE)**, orders species by median gene density.

## Requirements
Python 3.8 or newer is recommended.

Required Python packages;
`numpy` `pandas` `matplotlib` `scipy`

To install;
```
pip install numpy pandas matplotlib scipy
```

- An input directory named `gene_density` where there are individual species **gene density** file. Individual **gene density** can be extracted using TBtools as `.txt` file. **The directory should be located in the same working directory from where the script is run.**

Example:[`gene_density`](./gene_density)
```
gene_density/
    species1.txt
    species2.txt
    species3.txt
```


## Usage
```
python3 0.plot-gene-density.py
```

The script automatically searches for `./gene_density/*.txt`

**No input filenames needed to specify.**

# Input files

All **gene density** files that are in `./gene_density` directory should be in `.txt` format.
Each file must contain four columns;
```
scaffold  - Scaffold/chromosome
start     - Window start coordinate
end       - Window end coordinate
density   - Gene density
```

  Example: [`Avul_GeneDensity.txt`](./gene_density/Avul_GeneDensity.txt)
```
  scaffold      start      end      density
  CM147319.1    200000    300000      1.0
  CM147319.1    300000    400000      1.0
  CM147319.1    400000    500000      2.0
  CM147319.1    500000    600000      1.0
  CM147319.1    1100000   1200000     1.0
```

The species name will be taken from the prefix of `.txt` **gene density** file.

## Data Processing

For every input file, the script reads the four columns as;
`scaffold start end density`

The density column is converted to numeric values and invalid or missing density values are removed. Only the files containing more than one valid density observation are included in the final analysis.

### Species Ordering

All the species are ordered according to their median gene density. The script uses: `np.median(species_data[x])` and sorts species in descending order. Thus, species with the highest median **gene density** appears at the top of the plot.


### Kernel Density Estimation

The distribution for each species is estimated using a `Gaussian kernel density estimator`: `gaussian_kde(vals)`. The KDE converts the individual **gene-density** observations into a smooth probability-density distribution.

This provides a visual representation of where **gene-density** values are concentrated across the genome.

### Ridgeline Plot

Each species is represented by one ridge. The ridges are vertically stacked and normalized to a common visual height.
  The x-axis represents: `Gene density (genes per 100 kb)`
  The y-axis represents the `species`

### Median Density

A **vertical line** is drawn for each species representing its median **gene density**. The median is calculated using: `np.median(vals)`
This provides an easy visual reference for comparing the central tendency of gene density among genomes.


# Output

The script generates two output files.
```
Ridgeline plot generated for 5 species.

gene_density_ridgeline.png
gene_density_ridgeline.pdf
```
Script generates high quality, 600 dpi PNG and vector PDF file.

### Gene density units
[!NOTE]
- The **x-axis** is labelled: `Gene density (genes per 100 kb)`
- Therefore, the input density values should be calculated in units of genes per 100 kb.
