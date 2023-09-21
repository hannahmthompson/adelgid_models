# Data Driven Models of Hemlock Woolly Adelgid Impacts and Biological Control

Two models of hemlock woolly adeglid, _Adelges tsugae_ Annand (Hemiptera: Adelgidae), an invasive insect pest of hemlock trees, _Tsuga canadensis_ (L.) Carrière, in the eastern US:

* Hemlock-Adelgid model: _T. canadensis_ health and _A. tsugae_ density

* Adelgid-Predator model: _A. tsugae_, and biological control beetle predators _Laricobius nigrinus_ Fender (Coleoptera: Derodontidae), and _Sasajiscymnus tsugae_ (Sasaji and McClure) (Coleoptera: Coccinellidae) densities

The models and data here are presented in a doctoral dissertation from the University of Tennessee Knoxville, and two papers resulting from it.

* Thompson, Hannah M., "Data Driven Models of Hemlock Woolly Adelgid Impacts and Biological Control. " PhD diss., University of Tennessee, 2021.
[https://trace.tennessee.edu/utk_graddiss/6736](https://trace.tennessee.edu/utk_graddiss/6736)

* Thompson, H. M., McAvoy, T. J., Salom, S. M., Wiggins, G. J., Lenhart, S. 2023. Interaction of hemlock woolly adelgid populations with host tree health drives population oscillations. Ecosphere, 14(8).  [https://doi.org/10.1002/ecs2.4637](https://doi.org/10.1002/ecs2.4637)

* Thompson, H. M., Wiggins, G. , Lenhart, S. 2022. Modeling interactions of hemlock woolly adelgid and two biological control predatory beetle species in the Great Smoky Mountains National Park. Biological Control, 173. [https://doi.org/10.1016/j.biocontrol.2022.104986](https://doi.org/10.1016/j.biocontrol.2022.104986)


## Authors

Hannah M. Thompson [0000-0001-8316-3232](https://orcid.org/0000-0001-8316-3232)
University of Tennessee – Knoxville, Department of Mathematics

Thomas J. McAvoy [0000-0003-1157-9399](https://orcid.org/0000-0003-1157-9399)
Virginia Polytechnic Institute and State University, Department of Entomology

Scott M. Salom [0000-0002-8950-3442](https://orcid.org/0000-0002-8950-3442)
Virginia Polytechnic Institute and State University, Department of Entomology

Gregory J. Wiggins [0000-0002-0876-3366](https://orcid.org/0000-0002-0876-3366)
North Carolina Department of Agriculture and Consumer Services, Plant Industry Division

Suzanne Lenhart [0000-0002-6898-5796](https://orcid.org/0000-0002-6898-5796)
University of Tennessee – Knoxville, Department of Mathematics

## Description

_Adelges tsugae_ have spread into more than half of the native range of _T. canadensis_ in the eastern US since their introduction prior to 1951. An _A. tsugae_ infestation can result in _T. canadensis_ death within years, and this damage and mortality has caused significant changes in hemlock forests. We developed two models composed of systems of ordinary differential equations with time dependent parameters to represent seasonality.

The first model (Hemlock-Adelgid model) represents the coupled cycles in _T. canadensis_ health and _A. tsugae_ density.  We use field data from Virginia to develop the model and to perform parameter estimation. The mechanisms we represent in the model, including an _A. tsugae_ density dependent _T. canadensis_ growth rate, a _T. canadensis_ health dependent _A. tsugae_ mortality rate and a density dependent _A. tsugae_ mortality rate, produce the cycles in _T. canadensis_ health and _A. tsugae_ density commonly seen with the _A. tsugae_ system in the eastern United States. We test sets of initial conditions to determine scenarios that will likely lead to _T. canadensis_ mortality and explore the long-term dynamics of the system.

The second model (Adelgid-Predator model) represents the population dynamics of _A. tsugae_ and two introduced biological control organisms, _L. nigrinus_ and _S. tsugae_, released to attempt to control _A. tsugae_ populations and allow _T. canadensis_ survival. We use field data from the Great Smoky Mountains National Park to develop the model and to perform parameter estimation. The model is stage structured, with two classes per species. To represent the seasonality of the system, in addition to time dependent parameters, the model is composed of multiple systems of ordinary differential equations that control the system at different times of the year. We explore the dynamics of the model by running simulations with species absent from the system, and by varying the fixed parameter representing _T. canadensis_ health.

The data used to develop and perform parameter estimation on the Hemlock-Adelgid Model are from a study led by McAvoy, on the impact of chemical control on _A. tsugae_ and the health of _T. canadensis_, which recorded _T. canadensis_ canopy characteristics and _A. tsugae_ densities at a field site in Mountain Lake, Virginia from 2001 to 2015. A preliminary report of the study was presented in 2005.[^1] `ha_data.xlsx` contains the data from that study that were used in this modeling work.

The data used for development and parameter estimation for the Adelgid-Predator model are from a study led by Wiggins and were collected near Elkmont Campground in the Great Smoky Mountains National Park from 2010-2013. The number of larval and adult _L. nigrinus_ and _S. tsugae_ collected through beat sheet sampling were recorded. Part of the study assessed the use of emergence traps, which collect adult _L. nigrinus_ as they emerge from the soil where pupation and aestivation takes place.[^2] `als_data.xlsx` contains the data used in this modeling work.

## Repository Organization 

* `./adelgid_predator/` is the directory containing the Adelgid-Predator model files
    * `als_data.xlsx` contains data used in the Adelgid-Predator model
    * `als_data_plot.m` plots a time series of the adelgid-predator data 
    * `als_data_plot_seasonality.m` plots the predator data over a calendar year to show seasonality 
    * `als_graphical_abstract.pptx` is the graphical abstract for the Adelgid-Predator model
    * `als_model_schematic.pptx` is a schematic diagram of the biological mechanisms represented in the Adelgid-Predator model
    * `als_paraest.m` is the parameter estimation of the Adelgid-Predator model
    * `als_plot.m` plots the data and model simulation results using the optimal parameter values found through the parameter estimation for the Adelgid-Predator model
    * `als_scenarios_hemlock.m` plots the model simulation results using the optimal parameter values found through the parameter estimation, but varying the hemlock carrying capacity parameter value for the Adelgid-Predator model
    * `als_scenarios_predator.m` plots the model simulation results using the optimal parameter values found through the parameter estimation, but varying the the presence or absence of the predators through changing initial conditions for the Adelgid-Predator model
    * `als_seasonality_figure.pptx` is a figure representing the seasonality of the Adelgid-Predator model

* `./hemlock_adelgid/` is the directory containing the Hemlock-Adelgid model files
    * `h_likely_mortality_value.m` uses the hemlock time series data to find a proportion tips alive value where tree mortality is likely 
    * `ha_data.xlsx` contains the data used in the Hemlock-Adelgid model
    * For Group _i=1,2_
        * `ha_group_i_longterm.m` plots the model results over a longer time using the optimal parameter values found through the parameter estimation for Group _i_ for the Hemlock-Adelgid model
        * `ha_group_i_paraest.m` is the parameter estimation for Group _i_ for the Hemlock-Adelgid model
        * `ha_group_i_plot.m` plots the data and model results using the optimal parameter values found through parameter estimation for Group _i_ for the Hemlock-Adelgid model
        * `ha_group_i_scenarios.m` determines if/when model simulation results reach the proportion tips alive value likely to result in hemlock mortality for varying initial conditions for Group _i_ for the Hemlock-Adelgid model
        * `ha_group_i_scenarios_plot.m` plots model simulation results for various initial conditions for Group _i_ for the Hemlock-Adelgid model
    * `ha_scenarios_heatmap.R` plots the results of the varied initial condition scenarios as a heatmap for both groups for the Hemlock-Adelgid model
* `.gitignore` lists the files types for git to not track
* `LICENSE_Apache_2.0.md` is the License for the code in this project
* `LICENSE_Attribution_4.0_International.md` is the license for the data used in this project
* `README.html` is the html version of this document
* `README.md` is this document 
    
## Parameter Estimation

The files `als_paraest.m`, `ha_group_1_paraest.m` and `ha_group_2_paraest.m` perform parameter estimation on the Adelgid-Predator model (`als_paraest.m`) and the Hemlock-Adelgid model (`ha_group_1_paraest.m`, `ha_group_2_paraest.m`) using fmincon and MultiStart in MATLAB. Each file uses a similar structure. 

The models are implemented as several systems of ordinary differential equations where each system applies during some part of the year. Each of these systems has constant parameter values. Initial conditions for each system are taken from state values at the final time of the previous system (or from some other previous state value, as appropriate to represent the biology). Each system is numerically solved with the numerical ordinary differential equations solver ode45.

We construct a constrained optimization problem where we minimize an objective function over the parameter values and initial conditions. We set upper and lower bounds for each of the parameters estimated. The objective function value is the sum of the relative error for each class where we have data: <img src="https://latex.codecogs.com/svg.image?\sum\sqrt\frac{{\sum_{i=1}^{n}(N_{\mathrm{data}}(t_{i})-N_{\mathrm{model}}(t_i))^2}}{{\sum_{i=1}^{n}(N_{\mathrm{data}}(t_i))^2}}" title="\sum\sqrt\frac{{\sum_{i=1}^{n}(N_{\mathrm{data}}(t_{i})-N_{\mathrm{model}}(t_i))^2}}{{\sum_{i=1}^{n}(N_{\mathrm{data}}(t_i))^2}}" /> where <img src="https://latex.codecogs.com/svg.image?t_i" title="t_i" /> is the time of stage _N_ data _i_, and the outside sum is over the different stages where we have data. We use MultiStart to generate starting points and fmincon attempts to find a minimizer from each of the starting points.

There are four main sections of the code:

1. Optimization problem
    * We estimate the values of parameters and initial conditions by solving the optimization problem of minimizing the objective function value. To do this we:
        * save the upper and lower bounds for each parameter we are estimating as a vector,
        * create the optimization problem, specifying:
            * fmincon as the solver,
            * the function we create in section 3 that returns the objective function value,
            * the upper and lower bounds for the parameters,
    * create MultiStart object to run the fmincon solver repeatedly to find multiple local minima of the objective function value,
    * specify the number of start points,
    * run the problem and save the results.
    
2. Model solution with optimal parameter values found
    * Using the results found and saved above, we then solve the model using the best parameter values found, and plot the results.

3. Objective function value function
    * We solve the model and calculate and return the objective function value for a given set of parameter values _z_. The values of _z_ are determined by MultiStart and fmincon.

4. Model
    * We use a function for each system of ordinary differential equations in our model and call these functions in sections 2 and 3 when we are solving our systems.

We can test the code by using test data. We simulate the test data by choosing values for all parameters and solving the model. The test data is the state function values at some time points. We use the same time points where we have real data. Using this test data, run the parameter estimation. The optimal parameter values found should agree with the parameter values used to simulate the test data.

## License

The data used in this project (files `als_data.xlsx` and `ha_data.xlsx`) are licensed under the [Creative Commons Attribution 4.0.](https://creativecommons.org/licenses/by/4.0/) The remaining content of this project is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

[^1]: T. McAvoy, W. T. Mays, S. M. Salom, and L. T. Kok. Impact of imidacloprid on hemlock woolly adelgid _(Adelges tsugae)_ and water quality at Mt. Lake, Virginia. In B. Onken and R. Reardon, editors, _Third Symposium on Hemlock Woolly Adelgid in the Eastern United States_, pages 324–334. Forest Health Technology Enterprise Team, 2005

[^2]: G. J. Wiggins, J. F. Grant, J. R. Rhea, A. E. Mayfield III, A. Hakeem, P. L. Lambdin, and A. B. L. Galloway. Emergence, seasonality, and hybridization of _Laricobius nigrinus_ (Coleoptera: Derodontidae), an introduced predator of hemlock woolly adelgid (Hemiptera: Adelgidae), in the Tennessee Appalachians. _Environmental Entomology_, 45(6):1371–1378, 2016. [https://doi.org/10.1093/ee/nvw128](https://doi.org/10.1093/ee/nvw128)

