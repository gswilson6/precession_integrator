***EARTH'S ROTATION AXIS PRECESSION (intstructions of UNIX terminal)***

PURPOSE: 
This code will calculate the changing position of the 
Earth and Sun in the Earth-Sun orbit system, and how the Earth will
begin to precess over a long period of time. The data gathered by this code can then be run with a plotting language like Mathematica to observe that the Earth orbits the Sun in a nearly circular path over the course of one year.


RUNNING THE CODE:
In order to run this code the user first need to make sure they are
in the correct directory where they stored this program. 

Once the user is in the correct directory they can then proceed to type:

	make integrator

into the command line of the terminal. This will compile the code for the user. Now the user can run the code by entering:

	./integrator 10 100

into the command line. Where the first argument after ./integrator is the number of steps and the second argument is the number of years 


PRINTOUT:
The code will print out to a file called data.txt an array of numbers for each time step in the for loop. The file data.txt will be located in the directory 
where you run the code.
