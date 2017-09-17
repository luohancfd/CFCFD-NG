Example to test the implementation of a slow opening valve in L1d3. The valve 
was implemented to allow modelling of fast acting valves, such as the one used 
in the HDT in Oxford or other facilities.
The valve is modelled in three stages.
(1) before open_time: The valve is modelled as closed.
(2) between open_time and open_time + open_period: Valve opens
(3) after open_time + open_period: Valve is fully open. 

A detailed description of the Valve implementation is available in the Thesis by 
Aaliya Doolabh, A Concept Study for a Piston Driven sCO2 Turbine Test Facility,
Bachelor of Engineering Thesis, School of Mechanical and Mining Engineering, The
University of Queensland. 

During stage (2), the valve is modelled as a contraction, which extends a number 
of points upstream and downstream of the valve, to allow the valve body to be 
modelled with some axial length. The minimum diameter exists at the centre of
the valve. During the opening process the valve area is modelled as:
    A_valve = A_max * (dT / open_period)^2
where dT is the time that has elapsed since the valve has started opening.  


    +----+--               --+----+
            --+--     --+--
                 --+--

               throat = A_valve

                 --+--
            --+--     --+--
    +----+--               --+----+
   -3   -2   -1    0    1    2    3 
         <------------------->
       points that vary diameter
      as part of valve simualtion


A valve is added to a simulation in two steps:
*/ Add code for valve break point
add_valve_point(x_loc, d_mx, n_points)
    x_loc - float - x-coordinate of the valve point
    d_max - float - Maximum diameter of valve when fully open
    n_points - int - Number of points to vary diameter either side of the valve 
                    central node location (e.g. to model valve body)

*/ Define the valve
Valve(x0=5.712, open_time = 0.001, open_period = 0.0, is_open = 0, label="Valve")
    x0 - Float - Define the valve location as an x-position in the tube;
    is_open - Int - (default = 0) - Flag to indicate whether the valve is open (0), 
                                closed(1) or in the process of opening (2);
    open_period - Float - (default = 0.0) - Defines the time it take for the 
                                        valve to open
    open_time - Float - (default = 0.0) - Time after which the valve starts to open
    poly - Vector - (default = [-2.15909091e-24, 1.25227273e-16, 
                    -2.80855519e-9, 3.63299825e-2])
              - Used for automated calculation of open_time. Only used when 
                open_time = 0.0 or not specified. The vector descrbe the
                coefficients of the 3rd order polynomial to define open time as 
                a function of pressure immediately upstream of the valve.
    label - String - A label that will appear in the parameter file


The current case considers a valve separating two slugs of gas. 



