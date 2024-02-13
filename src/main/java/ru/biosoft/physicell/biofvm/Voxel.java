package ru.biosoft.physicell.biofvm;

/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

/*! \brief Voxels are the basic spatial container for densities, which are networked into meshes. 
 * 
 * Voxels are the basic spatial container for a finite volume method. Voxels are connected to 
 * other voxels into a General_Mesh, here most likely a Cartesian_Mesh. Voxel boundaries are Voxel_Faces. 
 * 
 * A Microenvironment Domain will include a network of Voxels (and Voxel_Faces), 
 * with a vector<double> of densities for each Voxel, along with rate constants, etc. 
 * The Domain may also include a vector<double> of flux coefficients for each Voxel_Face. 
*/
public class Voxel
{
    int meshIndex; /*!< voxel's index in a General_Mesh */
    double volume; /*!< voxel's volume (cubic spatial units) */
    public double[] center; /*!< center of volume */
    boolean isDirichlet;

    private static String tabbing = "\t\t\t\t";
    private static String tabbing2 = "\t\t\t\t\t";

    public Voxel()
    {
        meshIndex = 0;
        volume = 10 * 10 * 10;
        center = new double[3];
        isDirichlet = false;
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<voxel ID=\"" + meshIndex + "\">\n" );
        sb.append( tabbing2 + "<center " + center[0] + "\t" + center[1] + "\t" + center[2] + "/>\n" );
        sb.append( tabbing2 + "<volume>" + volume + "</volume>\n" );
        sb.append( tabbing + "</voxel>" );
        return sb.toString();
    }


    public String stream_output_with_units(String units)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<voxel ID=\"" + meshIndex + "\">\n" );
        sb.append( tabbing2 + "<center " + center + " units=\"" + units + "\" />\n" );
        sb.append( tabbing2 + "<volume units=\"cubic " + units + "\">" + volume + "</volume>\n" );
        sb.append( tabbing + "</voxel>" );
        return sb.toString();
    }
}
