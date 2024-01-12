# DyMINSTED4all
Open source instructions and code for making your STED microscope a Dymin STED

## Description
Using inexpensive parts and labview code, we implemented DyMIN STED Microscopy as detailed from in the paper: [Adaptive Illumination STED Nanoscopy](https://www.pnas.org/content/114/37/9797).  The following files include our code, parts list and instructions to set up.

## Parts list
* National Instruments myRIO-1900
* STED microscope: must use an image acquisition program that can provide a pixel clock to send to the DyMIN FPGA circuit, must also use an APD (up to two APDs can be included) with digital pulses for fluorescence detection
* AOM and RF driver for AOM, AOM must be aligned to control the power to the STED microscope using an analog input to the RF driver.
* National Instruments Inc, myRIO-1900 and power supply
* Digilent Inc, Wire Wrap or Protoboard Expansion for NI myRio
* Analog Devices Inc, ADV3221 soldered to an Evaluation Board (50 Ohm side)
* BNC and SMA cables and connectors
* Metal enclosure for circuit (Protocase, Inc.)
* RF function generators for testing
* Oscilloscope for testing
* Soldering supplies, wires, etc.

## DyMIN_FPGA_Main.vi
The file includes the labview code for the FPGA.

## Computer Requirements:  
Labview 2016 (32 bit)  
Labview Real-Time Module  
Labview FPGA Module  
Compile Cloud access (Sign up at https://www.ni.com/en-us/support/documentation/supplemental/14/compile-faster-with-the-labview-fpga-compile-cloud-service.html)  
OR Xilinx Compilation Tools to compile the code on your FPGA.  

Run this file before you run the DyMIN_RT_Main.vi file to generate the bitfile for your FPGA.  Warning, every time you compile the code after changes it will take several minutes to compile.  
## Experimental diagram
See DyMINSTED4ALL Schematic file in Documents folder for experimental schematic.

## AOM RF Driver Analog Voltage Control and Eval Board for Multiplexer
We implemented DyMIN with an AOM and RF driver purchased off of Ebay that is no longer commercially available.  The electrical properties of the analog voltage control of the RF driver (ie. the input impedance and voltage range required to obtain STED powers you want at the sample) are the relevant parameters that will determine how your AOM and RF driver will work with the multiplexer circuit.  With a 50 Ohm input coupling impedance and 0-1V analog modulation values for our RF driver, we were able to get the right range of currents out of the circuit by modifying a few of the resisters on the evaluation board.  Resistors R14 and R21 were replaced with 2k Ohm resistors and resistors R16 and R22 were removed to isolate the 75 Ohm side of the board. See the file in Documents called "Anal. Dev. ADV3221-EVALZ Schematic.pdf" (downloaded from Analog Devices, Inc.) for the evaluation board schematic.  If your RF driver has different analog control properties, you will need to modify these resistors so that the multiplexer circuit can provide the current at the voltage ranges required for your RF driver.

# ADV3221 Multiplexer Information and Truth Table
The multiplexer allows the FPGA's fast DIO lines to be used to quickly select between various Analog output values to be passed on based on a truth table provided by the manufacturer.

# Multiplexer Inputs and Outputs
Digital inputs that control the output:  
A0  
A1  
CK1  
CK2  
CS  

Analog inputs (+/- 3V) that are passed based on the truth table:
IN0  
IN1  
IN2  
IN3  

Analog output (one of IN0,IN1,IN2,IN3 is passed to the OUT based on the truth table):  
OUT=RF Driver analog voltage control input signal  

# Multiplexer Truth Table  
Multiplexer digital input names and truth table:  
0=FALSE, 1=TRUE  

CS A1 A0 CK1 CK2 -> Multiplexer IN# passed  
---------------------------
00000 -> IN0  
00100 -> IN1  
01000 -> IN2  
01100 -> IN3  
1XX00 ->High-Z X means anything  


#  myRio-1900 FPGA Information
Install the myRIO-1900 and the required Labview components on your computer by following the manufacturer's instructions and use the testing features (click on "Configure NI myRIO") when you connect your myRIO-1900 to the computer to test that it is working properly.

Digital inputs from STED microscope to myRIO:  
C/DIO0=PIX CLK  
C/DIO1=APD1 IN  
C/DIO2=APD2 IN  
(If you only have one APD just use APD1 and ignore APD2.  If you have both, the program adds them up and keeps performs DyMIN on the sum of the two APD signals.)

Digital outputs:  
C/DIO5=APD1 OUT  
C/DIO6=APD2 OUT  
C/DIO7=10 MHz OUT (used to sync our STED CLK with the single cycle timed loops (SCTL) of the FPGA)  

Digital outputs to multiplexer from myRio:  
A/DIO7=A0  
A/DIO8=A1  
A/DIO9=CK1=CK2  
A/DIO10=CS  

Analog outputs to multiplexer from myRIO:  
A/AO0=IN0  
C/AO0=IN1  
C/AO1=IN2  
A/AO1=IN3  

##  Instructions for Circuit Building and Testing:
First align your STED laser through your AOM and make sure you get good diffraction into the first order.  Then couple the 1st order diffraction beam into the fiber that couples your STED laser to your STED microscope.

Attach the ADV3221 to the Evaluation board.  Once you have attached the ADV3221 to the board, then wire up the digital and analog connections from the FPGA to the ADV3221 (see above).  Connect the output to an oscilloscope, terminating it with your RF driver's Analog voltage control input coupling impedance so you can check the voltage range output by the multiplexer is what you need for your application by monitoring the output on the scope as you change the input analog voltages and selected output with a testing program.

Compile MultiplerTesting.vi by running the code, which will take several minutes. You will most likely need to modify the resistors on the board to ensure the OUT gives a reasonable output for your RF Driver when the digital inputs you want are selected. You can use this program to test that the multiplexer is working properly.  The program will send the correct digital signals to the multiplexer to select the various output values, and you can also modify the analog output values.  See the comments and instructions in the vi file for further information.

You will likely need to modify some of the resistors on the board for it to control your RF driver correctly and at this point you need to work on that until you get it working properly.  For our application, resistors R14 and R21 were replaced with 2k Ohm resistors and resistors R16 and R22 were removed to isolate the 75 Ohm side of the board

Once your FPGA and multiplexer circuit can provide the RF Driver with the correct signals to switch between the STED powers that you want to select from with your AOM, you will probably want to use the multiplexer testing code to characterize the speed of the transition. First you will  want to characterize the electronic delay between the FPGA and the Multiplexer output.

Then, to characterize the time to change the STED laser power, detect the STED laser after the AOM with a fast photodiode.  Use the Multiplexer testing code and use a trigger signal from the FPGA (in the MultiplexerTesting.vi code DIO's C/DIO2-5 are set to send trigger signals out for this purpose) as the trigger signal to the oscilloscope to measure the time of the transition from digital output to optical power change.  Hopefully it will be <= 1 us.  This time will limit the speed of your DyMIN implementation, so minimizing it is helpful, all else being equal.

Now that you know the multiplexer circuit is working to change your STED power quickly you can move on to implementing DyMIN.  Purchase BNC connectors, SMA and BNC cables, elbows, etc so you can put the FPGA and Eval board inside of a metal enclosure.  We used a Protocase 2U Extruded Aluminum Type C Template with custom holes for the various BNC inputs and outputs. Mount connectors on the box and install your FPGA and Multiplexer circuit in the box.

# DyMIN_FPGA_Main.vi
This file includes the Labview code for the FPGA.

Run this file before you run the DyMIN_RT_Main.vi file to generate the bitfile for your FPGA.  Warning, every time you compile FPGA code after changes it will take several minutes to compile.


# DyMIN_RT_Main.vi
Check the connection to DyMIN_FPGA_Main.vi (look for bookmark #DyMIN_FPGA_Main in bookmark manager) and the myRIO are correct and then run this file to run the FPGA code from the computer.  Every time you change the underlying DyMIN_FPGA_Main.vi file you will need to recompile it by running it directly from the FPGA vi, then you can run it from this RT (RT=Real-Time) file again.

Once your code is compiled and ready to test, the first step of testing is to test it with known inputs, using function generators to provide an APD signal and pixel clock that are electrically equivalent to what you will be using with your STED microscope.  Send the OUT from the multiplexer, with appropriate output coupling resistor, to your oscilloscope and send the APD1 and/or APD2 output to the other channel of your oscilloscope.  As you modify the analog voltages, Threshold values and times for each DyMIN step you should verify that the outputs act as expected.

The next step is to test the counting of your APD signal by your STED image acquisition software.  Connect the pixel clock from your DAQ to the PIX CLK input on your circuit box, connect the APD1 OUT to your usual DAQ input counting digital line, but keep the function generator fake APD signal(s) as inputs to APD1 IN on the circuit.  Now you can make sure that your counts are being counted accurately by the software.

Now you are ready to run the DyMIN on some samples.  After you usual alignment routine, you can use the DyMIN method by just selecting DyMIN ON? to be true and choosing appropriate values for the voltages, time and thresholds for each DyMIN step.

Please contact Emily.Gibson@cuanschutz.edu or stephanie.a.pierce@protonmail.com for questions.
