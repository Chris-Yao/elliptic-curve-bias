# elliptic-curve-bias
Code for the paper "Applications of Moments of Dirichlet Coefficients in Elliptic Curve Families" (arXiv:2311.17215) by Zoë Batterman, Aditya Jambhale, Steven J. Miller, Akash L. Narayanan, Kishan Sharma, Andrew Yang, and Chris Yao started during the Probability and Number Theory group of the SMALL 2023 REU run by Steven J. Miller. Computes moments of elliptic curve families, simulates random moment data according to the Sato-Tate distribution and plots the data in various forms.

Authors: Akash Narayanan, Aditya Jambhale, and Chris Yao.

To contact the authors, please email chris.yao@berkeley.edu.

## Files

CalculateMoments.cpp - Calculates the moments of elliptic curve families. Outputs data in the form "a(t),b(t),c(t),final.txt" in the "data" folder where the elliptic curve family takes the form y^2 = x^3 + a(t)x^2 + b(t)x + c(t).

CalculateRandomMoments.cpp - Calculates random moments of an elliptic curve family; not used in paper. Outputs data in the form "a(t),b(t),c(t),random.txt" in the "data" folder where the elliptic curve family takes the form y^2 = x^3 + a(t)x^2 + b(t)x + c(t).

sato_tate_simulations.ipynb - Simulates and plots random second moments based on the Sato-Tate distribution.

plotting.ipynb - Plots data.
