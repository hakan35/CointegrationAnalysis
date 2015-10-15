This is an initial and preliminary implementation of the Johansen et al. (200) methodology based on the functions provided by B. Pfaff in the R package "urca". 
Additionally you can't rely on the available critical values, but rather need to use the methods published by Giles, Godwin (2012). The R code for these calculations can be found on [1].
I will try to make it more user-friendly and also will include some code to set up all the matrices needed to implement the methodology correctly.
Comments, bug reports and suggestions are highly welcome!

- Johannes

[1] http://web.uvic.ca/~dgiles/downloads/johansen/index.html



References:

Giles, D. E. and R. T. Godwin (2012), "Testing for Multivariate Cointegration in the Presence of Structural Breaks: p-Values and Critical Values", Applied Economics Letters, 19, 1561-1565

Johansen, S., R. Mosconi and B. Nielsen (2000), "Cointegration analysis in the presence of structural breaks in the deterministic trend", Econometrics Journal, 3, 216-249. 

 Joyeux, R. (2007), "How to deal with structural breaks in practical cointegration analysis", in B. Bhaskara Rao (ed.), Cointegration for the Applied Economists, (2nd. ed.), Palgrave Macmillan, New York, 195-221.

