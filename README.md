# KemiMo
Three-phase gas-grain astrochemical model. The code is provided "as-is" and not well documented at the moment. 
The code consists of a python wrapper around a fortran core.
Feel free to contact me if you encounter any bugs or have questions.

The code is provided under the GPL v3 license. I ask that you refer to Jensen et al. (2021) [https://ui.adsabs.harvard.edu/link_gateway/2021A&A...649A..66J/doi:10.1051/0004-6361/202040196] if the code of parts hereof are used in published works.

# Chemical network
The chemical network is provided in KIDA format. All files are located in the "data_deuterated_total" directory.
The current gas-phase network provided here is an extension of the network by Majumdar et al.(2017) [https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.4470M/abstract]

# Running the code
To run the "standard" steady-state model, type "python main.py standard" in the kemimo directory.
The alternative non-spin KIDA 2024 model can be selected with "python main.py kida2024".
The "main.py" script copies the selected model from the "models" directory into the repository root as a runtime workspace. These copied files and the generated Fortran build products are not tracked.
In the model directory: 
 - "config.nml" contains the chemical-network and preprocessing settings.
 - "pathway.py" controls whether to run the model or only pre-process the network. Pre-processing includes reading the chemical network and writing the Fortran files for compilation.
 - "main.f90" is the model script. In this file, one must specify the initial abundances and the physical conditions or evolution. A basic template is provided. 

Boolean settings in `config.nml` use `1` for true and `0` for false. `h2_spin` controls whether its
ortho/para chemistry is included. Reaction-diffusion competition is enabled
by default and can be disabled with:

```
reaction_diffusion_competition = 0
```

Three-phase models use the thin-ice approximation by default. It relaxes the
mantle tolerances and skips mantle transfer terms only while the mantle is
dynamically insignificant. The normal three-phase ODE remains available with:

```
thin_ice_approximation = 0
```

# Development checks
Run the Python tests with:

```
make test
```

Preprocess and compile the default three-phase spin model without running its
trajectory with:

```
make smoke
```

The model and configuration can be changed on the command line, for example:

```
make smoke MODEL=kida2024 NLAYERS=2 H2_SPIN=0
```

`LAYER_THICKNESS` and `THIN_ICE_APPROXIMATION` can also be overridden for a
compile check:

```
make smoke NLAYERS=2 LAYER_THICKNESS=4 THIN_ICE_APPROXIMATION=0
```

GitHub Actions runs the Python tests and compile checks for the two-phase spin,
three-phase spin, and three-phase no-spin configurations.

# Disclaimer
"kemimo" is provided "as it is", without any warranty. The Authors assume no liability for any damages of any kind (direct or indirect damages, contractual or non-contractual damages, pecuniary or non-pecuniary damages), directly or indirectly derived or arising from the correct or incorrect usage of kemimo, in any possible environment, or arising from the impossibility to use, fully or partially, the software, or any bug or malfunction. Such exclusion of liability expressly includes any damages including the loss of data of any kind (including personal data)
