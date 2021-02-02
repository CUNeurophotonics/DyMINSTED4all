# DyminSTED4all
Open source instructions and code for making your STED microscope a Dymin STED

## Description

Using inexpensive parts and labview code, we implemented Dymin STED Microscopy as detailed from in the paper: [Adaptive Illumination STED Nanoscopy](https://www.pnas.org/content/114/37/9797).  The following files include our code, parts list and instructions to set up.

## Parts list
National Instruments myRIO-1900

## Experimental diagram

## Dymin_FPGA_Main.vi 
The file includes the labview code for the FPGA.

Requirements:  
Labview 2016 (32 bit)  
Labview Real-Time Module  
Labview FPGA Module  
Compile Cloud access (Sign up at https://www.ni.com/en-us/support/documentation/supplemental/14/compile-faster-with-the-labview-fpga-compile-cloud-service.html)  
OR Xilinx Compilation Tools to compile the code on your FPGA.

Run this file before you run the Dymin_RT_Main.vi file to generate the bitfile for your FPGA.  Warning, every time you compile the code after changes it will take several minutes to compile.

## Dymin_RT_Main.vi
Check the connection to Dymin_FPGA_Main.vi (look for bookmark #Dymin_FPGA_Main) and then use this file to run the code after you have compiled the FPGA code.

