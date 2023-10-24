# ProjectMyPixel

ProjectMyPixel is an image processing project designed to analyze patterning in fluorescent images of chick embryos. This project includes scripts to process data and generate projection plots of images.

## Directory Structure

- **test_data/ap_sagittal**: Raw AP sagittal images are located in this directory.
- **ap_sagittal/processed_data**: Processed data files are stored here.
  - **processed_data/images**: Processed images, either max or sum projections, are saved in this folder.
  - **processed_data/labels**: Labeled images corresponding to the processed images are stored here.

## Usage

1. **Processing Images**:
   - Process images are RGB, label using the Computer Vision toolbox.
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

## Required Scripts

- **main_ap_patterning_ex_vivo.m**: Main script for processing AP patterning in biological images.
- **projectMyPixel.m**: Function for processing and generating processed images.
- **percentileTraces.m**: Function for generating percentile traces.


