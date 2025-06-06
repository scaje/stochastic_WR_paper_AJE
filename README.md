# The Wells-Riley Model Revisited: Randomness, Heterogeneity and Transient Behaviours - AJE
This repository contains the supporting code and data for the following manuscript:
> Edwards, A. J., King, M.-F., Noakes, C. J., Peckham, D., & López-García, M. (2024). The Wells–Riley model revisited: Randomness, heterogeneity, and transient behaviours. Risk Analysis, 44, 2125–2147. https://doi.org/10.1111/risa.14295

# Software
This code is written using Python in Spyder 4.1.4.

# Description
1. To reproduce the results in Case Study 1, use the script 'HCW_visit_AJE.py'. This script relies on the data that can be found in the excel file 'time_in_hospital_room.csv'. The user needs to edit the filepath on line 21 to be the working file path that corresponds to this excel sheet.

2. The file entitled 'time_in_hospital_room.csv' contains the data which corresponds to that used in a previous study (https://doi.org/10.1016/j.jhin.2020.11.016). In this work, only the actual care times are used, not the mock care. The file path to this data should be included in the script 'HCW_visit_AJE.py' on line 21. 

3. To reproduce the results in Case Study 2, use the script entitled 'meal_times_AJE.py'. This script contains data from a previous study (Table 1 in https://doi.org/10.1016/S0195-6663(03)00109-0).
