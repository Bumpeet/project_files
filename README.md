# project_files
This are my MTP files

Step for running the ISOM model
1) Edit the parameters of the PFR model ISOM.mo in the OpenModelica editor.
2) parameters like, length of the reactor, molar flow rate, temperature, pressure, mole fraction, diameter of the rector.
3) Run the simulation till the length of the reactor, i.e., set the time as reactor length.

Steps for tuning
1) Run the external.mos in the command prompt as> omc External.mos
2) These should produce all the object files, XML files and a CSV file.
3) open the original_data.csv present in the same folder for updating the end point mole fraction
4) then open the R file
5) make the changes you want and run the file.
6) for the results check ISOM.PFR_res.csv file
