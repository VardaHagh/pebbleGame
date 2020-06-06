Instructions:

To compile and use the pebble game on your computer, follow these steps:

1- Open the file "test.f" and go to line 45 where the number of atoms are specified. Enter the number of nodes in your network (for example n = 512). Then save and close.

2- Open terminal and go to the directory where pebble game is located. For example, if it's on your Desktop, you can go to its directory by typing: cd Desktop/pebbleGame

3- Once you're in the pebble game directory, compile the code by running this command in terminal: make

4- You will see a bunch of warnings popping up as the code compiles, but you can ignore them. They cause no harm! :)

5- Once the pebble game is compiled, make sure that your contact list which should be called "points.dat" is inside the directory of pebble game. This file contains information about nodes that are connected to one another. Note that in Fortran, the indexing of arrays and lists starts from 1. So in your contact list, the particle indices should start from 1 (see the example points.dat file I have included in this directory).

6- To run the pebble game on your network and get an output, run this command in terminal:
./pebblegame.e points.dat > pebbleGameOutput.dat
This gets the "points.dat" file, runs the pebble game on it and saves the output in a data file called "pebbleGameOutput.dat". You can save the output in any file that you want.

7- If the pebble game is compiled correctly, you will see that a file named "pebbleGameOutput.dat" is saved in your directory. You can open this file and you should see that it has 4 columns. The first and second columns are the same as your contact list. They show which nodes are connected to each other in a bond. The third column shows the rigid cluster to which the bond belongs. If there is only one rigid cluster spanning the entire network, you will see that the third column shows 1 everywhere. If there are multiple clusters, you will see that there are different numbers in the third column. If you have more than 1 rigid cluster, it means that there should be hinges which are nodes that attach independent rigid clusters. Hinges are the nodes that belong to more than one cluster. The fourth column includes information about whether or not if a bond belongs to an isostatic region or an over-constrained region. If a bond has number 0 in the fourth column, it is isostatic. If it has number 1, it is called "stressed". Regions that are stressed have redundant bonds in them and any bond that is stressed can be a redundant bond, but not all bonds with 1 in their fourth column are redundant (because then in the example network all the bonds would be redundant!).

8- Finally, you can visualize the output of the pebble game by running the python script in the directory titled "visualizeNetwork.py". To run this piece, type the following command in the terminal: python3 visualizeNetwork.py
After running this command, you will see that a picture of the network is saved in the pebbleGame directory. If you don't want the numbers of the nodes to be printed, comment out lines 60 and 61 in the code by adding a # sign (hashtag) in front of the line.