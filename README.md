# BUD-ELM-Image-Analysis
Code files relating to BUD-ELM Flask Image Analysis

Each File is used for a separate step of image processing. 

Split_Channels.m is used to isolate the blue channel from the original RGB Image
Material_Image_Analysis.zip is the ilastik image analysis program used to segment the blue channel images into binary masks
BUD_ELM_Processing.m is used to analyze binary mask files for flat surface area, and then map those flat surface areas to a particular Modified Volumetric Power value
Flask Diam Conversions mm_px.csv is used by Final_Mask_Processing to convert pixels to mm for each image file

Blue Channel images and Final Segmented Masks have been provided within this repository. Raw images are available upon request. 
