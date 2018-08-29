# Time-Curve-Projection
The script was developed for the analyis of gene expression data with two types of transcriptome measurements: (1) A time series of control cells undergoing the process of reprogramming from fibroblasts to induced pluripotent cells; (2) Knockdown of the expression of a number of genes using siRNA measured at one time point. The hypothesis is that the knockdowns cause changes in gene expression, some of which cause the reprogramming to be delayed or accelerated, in addition to changes in gene expression that are unrelated to the reprogramming process. In other words, the question is: If cells are perturbed, what time point of control cells are they most closely related to? To answer this question, the transcriptome data was first subjected to principle component analysis (PCA), such that each data point is represented as a dot in a two-dimensional plane, with the x and y-axis corresponding to principal components 1 and 2 (PC1 and PC2). The script is then used to compare the perturbations to the time series of control cells.
## Limitations
The method only works if the changes associated with the time series are well represented in one of the principal components. It is most intuitive if time is represented on the x-axis. If necessary the x and y-axis can be swapped. For a propper fitting of a polynomial function to the time series data, it is also necessary that there is only one y-axis position for each x of the time series in the PCA space. In some cases it can be useful to rotate all data points (time series data and perturbation data) to obtain a better fit with the polynomial function. The degree of rotation should be such that x-axis corresponds maximally to the progression of time. If the progression of time in PCA space is such that for some positions on x, there are multiple solutions of y, it is not possible to obtain a good fit of the time line data with a polynomial function and the script cannot be used.

## What the script does
The script fits a polynomial function to the data points of the time series of control cells (-t filename with time series). It then projects the data contained in a second file, the experimental perturbations (-d file name) on the polynomial function using the shortest distance (perpendicular to the polynomial function). The output text file contains both the original coordinates (x,y) of each data point and the coordinates of the data point projected onto the polynomial function (x on curve, y on curve). The projection coordinates can be used as an estimation of progress of the perturbed cells relative to the control cells. The distance between the original data point and its projection on the time line can be interpretated as a proxy for gene expression differences that are unrelated to the process of reprogramming (not included in the output text file, but simply calculated using Pythagoras' formula). The formula of the polynomial can be found in the header of this output file. A second output file (PDF, optional) plots the original data, the fitted polynomial curve and the projection of the data points on the curve.

## Usage (projection_on_time_curve.py -h)

projection_on_time_curve.py, Arguments:
 * -h
 * -t   \<input filename time series\>
 * -d   \<input filename gene knockdown data\>
 * -o   \<filename output txt file\>
 * -f   \<fit polynomial order, integer 1..6\>
 * [-p] \<filename output pdf plot\>

General format input files:
  - Tab seperated .txt files
  - Comments and headers should start with a '#'
  
Format time series input file: 
 
  #Day &nbsp;&nbsp;&nbsp;&nbsp; PC1  &nbsp;&nbsp;&nbsp;&nbsp;   PC2

Format gene knockdown data input file:
  
  #Gene &nbsp;&nbsp;&nbsp;&nbsp;  PC1  &nbsp;&nbsp;&nbsp;&nbsp;  PC2


## Example

Executing the following code, with the input files that have been provided, will create an output file output.txt, where the data has been projected on a polynomial of order 2:  

```
python projection_on_time_curve.py -t Input-files\Time_rotated.PCA_2-7.txt -d Input-files\siRNA.Time_rotated.PCA.txt -o output.txt -f 2 -p output.pdf
```
The output file that this would produce can be found in the Output folder.

## Requirements
The script has only been tested on Python 3.6 and it uses the following libraries: Sympy, Numpy, Scipy and Matplotlib.
