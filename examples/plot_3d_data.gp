# template for plotting 3D data that has been written by WriteDataGnuplot
# this first plots a solid white sphere, and then the data points.
# the data points may lie on the sphere or not

# USAGE: from the command line, call gnuplot with this file as an additional argument, e.g. type "gnuplot plot_3d_data.gp"

# TIP: if the amount of data is very large and rotation of the plot by mouse is not smooth, rotate the figure while holding the CTRL key pressed.


# set the radius of the solid white sphere. suggestion for data on the sphere: make the radius of the sphere slightly smaller as to avoid overlays in the wrong order
R=0.99
# if the amount of data is very large, you can increase this number. e.g. if plot_only_every=5, then only every fith point will be plotted.
plot_only_every=1
# manually set minimum and maximum range of the values. this range will be resolved by the color space.
set cbrange [-5:5]
# function for scaling the values that determine the color of the plotted points
scale(x)=x

# the name of the file that contains the data that is to be plotted. it must be given as a path name relative to the position of this plot file's position
filename_a='../data/dir_D3_N256_s0_e0.txt'
filename_b='../data/output_demo2.txt'

# no need to adjust the following options
set view equal xyz
set parametric
set isosamples 10,10
set hidden
set urange [-pi/2:pi/2]
set vrange [0:2*pi]
unset key

# the plot command: point locations on the sphere
set term wxt 1
set title "training directions"
splot R*cos(u)*cos(v),R*cos(u)*sin(v),R*sin(u) w l lc rgb "white",\
      filename_a u 1:2:3  pt 5
print "Press ENTER key for next plot."

# pause -1 means that gnuplot will wait for the user to hit the enter key on the console
 pause -1

# the plot command: plot data, in different window
set term wxt 2
R=4.99
set title "interpolation results of linear function"
splot R*cos(u)*cos(v),R*cos(u)*sin(v),R*sin(u) w l lc rgb "white",\
      filename_b u 1:2:3:(scale($4)) every plot_only_every w p ps .2 pt 5 lc palette
print "Press ENTER key for exit."

# pause -1 means that gnuplot will wait for the user to hit the enter key on the console
 pause -1
