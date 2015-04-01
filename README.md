# NE-155

## Homework 4

The files are run with the command "python *filename*", where *filename* is hw4\_6j.py, hw4\_6gs.py, or 
hw4\_6sor.py (it is assumed that the user has Python installed on their system).

## Homework 5

The files are run with the command "python *filename*", where *filename* is hw5\_4j.py, hw5\_4gs.py, 
hw5\_4sor.py, or hw5\_4omega.py.

The output format for these files (except for the calculation of omega\_opt) is:

number of iterations  
solution vector  
error

The error is calculated as the absolute difference between the calculated solution vector and the result 
of using the NumPy solve function.

For the calculation of omega\_opt, the output is a plot of omega vs. number of iterations, called 
hw5\_omega.pdf.

## Homework 6

All of the files pertaining to this assigment have the prefix of "hw6".

The files are run in the same manner as the ones submitted previously, though some of the files request
user input at runtime. This input should be in the format of either integers or floats.

hw6\_2.py creates a plot and saves it as "hw6\_2.pdf". hw6\_3.py prints out the mesh size and the maximum
relative error in the flux solutions; the resultant corresponding plot is a compilation of all 
combinations of mesh sizes and error values and is saved as "hw6\_3.pdf".

The codes for problem 4 all print out the absolute error tolerance used, the chosen mesh size, and the
number of iterations required for convergence. The resultant corresponding plots are compilations of all
combinations of mesh sizes and error values.

hw6\_5.py prints out the dominant eigenvalue (k) as well as the number of iterations required for 
convergence. It saves a plot of the normalized flux solution as "hw6\_5.py".

Problem 6 is solved with built-in NumPy functions as well as a combination of the Power Iteration method
and the SOR method. The code prints out the dominant eigenvalue calculated using the built-in NumPy
functions as well as the norms of the maximum relative and absolute differences between the flux 
solutions. The two normalized flux solutions are plotted and saved as "hw6\_6.pdf".
