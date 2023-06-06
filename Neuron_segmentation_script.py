"""
Created on Tue Jun  6 10:07:55 2023

@author: Mohammad pourhosseinzadeh
Huising Lab
"""
#-----------------------------------------------------------------------------#
""" The complete segmentation script below was used to analyze Sst-crextdTomato 
brains by copying and pasting it into the Napari console. The script below was 
used to segment cells within brain images taken with a 4x objective"""
#-----------------------------------------------------------------------------#
import numpy as np
import cv2
import scipy.ndimage.filters as filters
import scipy.ndimage as ndimage
import pandas as pd

image=viewer.layers['image_name'].data

px=11
a1=cv2.GaussianBlur(image,(px,px),0)
a2=cv2.GaussianBlur(image,(px*2+1,px*2+1),0)
a3=a1.astype(int)-a2.astype(int)
aa=(a3>0)*a3
thresh,binary=cv2.threshold(aa.astype('uint16'),0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)

neighborhood_size=5

data_max=filters.maximum_filter(a3,neighborhood_size)
maxima=(a3==data_max)
data_min=filters.minimum_filter(a3,neighborhood_size)
diff=((data_max-data_min)>5) # For real images this number needs to be a small positive number like 5
maxima[diff==0]=0
maxima=maxima*((a3>thresh)*1)

labeled, num_objects = ndimage.label(maxima)
slices = ndimage.find_objects(labeled)
x, y = [], []
for dy,dx in slices:
    x_center = (dx.start + dx.stop - 1)/2
    x.append(x_center)
    y_center = (dy.start + dy.stop - 1)/2    
    y.append(y_center)

m=(a3>thresh)*a3 #for real images, this parameter needs to be a small positive number like 10
mask=np.zeros([len(m[0]),len(m[1])])
x=np.array(x).astype(int)
y=np.array(y).astype(int)

for row in range(0,m.shape[0]):
  for column in range(0,m.shape[1]):
    if m[row,column]>0:
      distance=(((row-y.astype(int))**2)+((column-x.astype(int))**2))**(1/2)
      mask[row,column]=pd.Series(distance).idxmin()+1

mask=mask.astype(int)
viewer.add_labels(mask)

#-----------------------------------------------------------------------------#
"""manually erase areas on edge of brain that are artifacts of the bright 
flourescence near the edge of the brain"""
#-----------------------------------------------------------------------------#

mask_new=viewer.layers['mask'].data
sst_cell_count=len(np.unique(mask_new))

np.savetxt('3_2_23_1205_Dapi_tdTomato_brain_4x_large_image_slide1_brain4_001_post_process_segment.csv',mask_new,delimiter=',') # The mask containing ROIs of all Sst+ cells in the brain is saved as a csv file 

#-----------------------------------------------------------------------------#
"""To load cell mask from csv file copy and past the script below in the 
Napari console"""
#-----------------------------------------------------------------------------#
import pandas as pd
import numpy as np

mask=pd.read_csv('name of file.csv')
mask_array=np.array(mask)
mask_array=mask_array.astype(int)
viewer.add_labels(mask_array)

#-----------------------------------------------------------------------------#
"""To quantify the number of Sst+ neurons in specific brain regions use the 
script below"""
#-----------------------------------------------------------------------------#

# First crop the region of brain you are interested in
# This function can be found in Napari under tools> utilities> crop
# Name the cropped brain regions 'brain region_mask'
cortex=viewer.layers['cortex_mask'].data
hypothalamus=viewer.layers['hypothalamus_mask'].data

cortex_neurons=len(np.unique(cortex))
hypothalamus_neurons=len(np.unique(hypothalamus))
