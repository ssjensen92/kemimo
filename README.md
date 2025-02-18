# kemimo
Three-phase gas-grain astrochemical model. The code is provided "as-is" and not well documented at the moment. 
The code consists of a python wrapper around a fortran core.
Feel free to contact me if you encounter any bugs or have questions.

The code is provided under the GPL v3 license. I ask that you refer to Jensen et al. (2021) [https://ui.adsabs.harvard.edu/link_gateway/2021A&A...649A..66J/doi:10.1051/0004-6361/202040196] if the code of parts hereof are used in published works.

# Chemical network
The chemical network is provided in KIDA format. All files are located in the "data_deuterated_total" directory.
The current gas-phase network provided here is an extension of the network by Majumdar et al.(2017) [https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.4470M/abstract]

# Running the code
To run the "standard" steady-state model, type "python main.py standard" in the kemimo directory. 
The "main.py" script will read the model called "standard" located in the "models" directory. 
In the model directory: 
 - "pathway.py" reads the data directory and holds information on whether to run the model or just pre-process the network. Pre-processing includes reading the chemical network and writing the fortran files for compilation.
 - "main.f90" is the model script. In this file, one must specify the initial abundances and the physical conditions or evolution. A basic template is provided. 
  

# Disclaimer
"kemimo" is provided "as it is", without any warranty. The Authors assume no liability for any damages of any kind (direct or indirect damages, contractual or non-contractual damages, pecuniary or non-pecuniary damages), directly or indirectly derived or arising from the correct or incorrect usage of kemimo, in any possible environment, or arising from the impossibility to use, fully or partially, the software, or any bug or malfunction. Such exclusion of liability expressly includes any damages including the loss of data of any kind (including personal data)
