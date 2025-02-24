# Temperature-estimation-battery-pack
#Steps to change number of sensors to use 
    1-> Change matrix of facobian_measurement_model function
    2 -> Inside data_initialization. Variable R, determine how many RStates wit how many sensors you gonna use
    3 -> Decide in which batteries you want sensors "x_k_1,P_k_1, Breakflag = ekf.measurement_update(States_True[["Ts1", "Ts2", "Ts5", "Ts7"]]")