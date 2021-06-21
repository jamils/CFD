using Cxx

cxx""" 
    #include <fstream>
    #include <iostream>
    #include <stdio.h>
    #include "math.h"
    #include <iomanip>
    #include <string>
    #include <sstream>
    #include <algorithm>
    #include <vector>

    using namespace std;
"""
cxx"""
    void inputMesh(double xn[], double yn[], double zn[], double xc[], double yc[], double zc[])
    {
        double readval; // read in value
        vector<double> data; // create an abitrarily long vector
        string line; // string of the line 

        // get filename
        string filename = "Project_Files/Grids/Inlet.33x17.grd";
        const char *FILENAME = filename.c_str();
        // open the instream
        ifstream infile;
        infile.open(FILENAME);

        // check for reading in error
        if (infile.fail()) 
        {
            cerr << "ERROR: Mesh File Cannot Be Opened!\nDid you extract the grids?!?!" << endl;
            exit(1);
        }

        // Read in data from file
        while (infile.good()) 
        {
            while (getline(infile, line)) 
            {
            istringstream streamA(line);
            while (streamA >> readval)
            {
                data.push_back(readval);
            }
            }
        }
        // data[0] is some number thats irrelevant
        int ni = data[1]-1; // 2nd number in the file is the number of cols (xvals i)
        int nj = data[2]-1; // 3rd number in the file is the number of rows (yvals j)
        int nk = data[3]-1; // number of z vals

        // resize output based off of the read input
        xc.resize(nj,ni); // resize (rows, cols) 4 rows by 5 cols : nj=4, ni=5
        yc.resize(nj,ni);
        zc.resize(nj,ni);

        xn.resize(nj+1,ni+1);
        yn.resize(nj+1,ni+1);
        zn.resize(nj+1,ni+1);

        int x_index = 4; // start of x vals
        int y_index = x_index + (nj+1)*(ni+1)*(nk+1);
        int z_index = y_index + (nj+1)*(ni+1)*(nk+1);


        // fill in the nodal values

        for (int row = 0; row < xn.rows(); row++)
        {
            for (int col = 0 ; col < xn.cols(); col++)
            {
            xn(row,col) = data[x_index];
            yn(row,col) = data[y_index];
            zn(row,col) = data[z_index];
            x_index++;
            y_index++;
            z_index++;
            }
        }

        for (int row = 0; row < xc.rows(); row++)
        {
            for (int col = 0 ; col < xc.cols(); col++)
            {
            // take average of all four xvals from corn
            // order: botL + topL + botR + topR 
            xc(row,col) = 0.25*(xn(row,col) + xn(row+1,col) + 
                xn(row,col+1) + xn(row+1,col+1));
            yc(row,col) = 0.25*(yn(row,col) + yn(row+1,col) + 
                yn(row,col+1) + yn(row+1,col+1));
            zc(row,col) = 0.25*(zn(row,col) + zn(row+1,col) + 
                zn(row,col+1) + zn(row+1,col+1));
            }
        }
    }
"""