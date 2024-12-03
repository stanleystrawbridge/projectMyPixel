# projectMyPixel

projectMyPixel is an image processing tool for simple dimensional reduction by the projection of pixels from fluorescent images onto lower dimensional manifolds.


## 2024_chinnaiya_et_al

main.m takes in segmentation and multichannel images of developing chick embryos and outputs anterior-posterior (AP) projects and raw_percentile_traces.

Analysi

### Directory Structure

- **data/sum_stacks**: INPUT sum stacks of raw confocal fluorescence images.
- **data/segmentation**: RGB images generated from images in 'data/sum_stacks' overlaid with roi and ap labels from manual segmentation and annotation.
- **data/segmentation.mat**: INPUT object housing maunual segmentation data.
- **data/projections**: OUTUPT First analysis resulting from dimensional reducation by pixel projection.
- **data/raw_percentile_traces**: OUTUPT Takes in projections from 'data/projections' and calulates the median AP signal trace.
- **src**: Houses classes projectMyPixel.m and percentileTraces.m
-
### Usage

1. **Processing Images**:
   - Process images into RGB, label using the Computer Vision toolbox.
   - Place the processed images in the **max** or **sum** folder, respectively.
   - In the **labels** directory, put the corresponding labeled images.
   - Change the method in **projectMyPixel.m** between either the maximum or sum projection method.

2. **Running the Main Script**:
   - Open **main_ap_patterning_ex_vivo.m**.
   - Provide the following user inputs:
     - **Image Type**: Specify the type of image data (e.g., AP sagittal images).
     - **Projection Method**: Choose the projection method (default is linear).
     - **Data Path**: Set the path where the data is located.
   - Run the script, which consists of three blocks of code:
     1. **Setting up Directories**: Creates necessary directories for storing processed and labeled data.
     2. **Processing Images**: Steps through the images and processes them using the specified projection method.
     3. **Percentile Traces**: Generates percentile traces for further analysis.

### Required Scripts

- **main_ap_patterning_ex_vivo.m**: Main script for processing AP patterning.
- **projectMyPixel.m**: Function for processing and generating processed images.
- **percentileTraces.m**: Function for generating percentile traces.

