package ru.biosoft.biofvm;

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
public class VoxelFace
{
    int mesh_index;

    double surface_area;
    //        std::vector<double> center; 
    double[] center;
    //        std::vector<double> outward_normal;
    double[] outward_normal;
    //        std::vector<double> inward_normal;
    double[] inward_normal;

    static String tabbing = "\t\t\t\t";
    static String tabbing2 = "\t\t\t\t\t";

    public VoxelFace()
    {
        mesh_index = 0;

        surface_area = 10 * 10;
        center = new double[3];//.assign( 3 , 0.0  ); 
        outward_normal = new double[3];//.assign( 3 , 0.0 ); 
        inward_normal = new double[3];//.assign( 3 , 0.0 ); 
    }

    //    std::ostream& operator<<(std::ostream& os, const Voxel_Face& vf)  
    //    {
    //        static std::string tabbing = "\t\t\t\t"; 
    //        static std::string tabbing2 = "\t\t\t\t\t"; 
    //        os  << tabbing << "<voxel_face ID=\"" << vf.mesh_index << "\">"  << std::endl
    //            << tabbing2 << "<center " << vf.center << " />" << std::endl  
    //            << tabbing2 << "<outward_normal " << vf.outward_normal << " />" << std::endl  
    //            << tabbing2 << "<inward_normal " << vf.inward_normal << " />" << std::endl  
    //            << tabbing2 << "<surface_area>" << vf.surface_area << "</surface_area>" << std::endl  
    //            << tabbing  << "</voxel_face>"; 
    //
    //     return os; 
    //    }
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<voxel_face ID=\"" + mesh_index + "\">\n" );
        sb.append( tabbing2 + "<center " + center + " />\n" );
        sb.append( tabbing2 + "<outward_normal " + outward_normal + " />\n" );
        sb.append( tabbing2 + "<inward_normal " + inward_normal + " />\n" );
        sb.append( tabbing2 + "<surface_area>" + surface_area + "</surface_area>\n" );
        sb.append( tabbing + "</voxel_face>" );
        return sb.toString();
    }

    //
    //    void Voxel_Face::stream_output_with_units( std::ostream& os , std::string units ) const
    //    {
    //        static std::string tabbing = "\t\t\t\t"; 
    //        static std::string tabbing2 = "\t\t\t\t\t"; 
    //        os  << tabbing << "<voxel_face ID=\"" << mesh_index << "\">"  << std::endl
    //            << tabbing2 << "<center units=\"" << units << "\" " << center << " />" << std::endl  
    //            << tabbing2 << "<outward_normal units=\"" << units << "\" " << outward_normal << " />" << std::endl  
    //            << tabbing2 << "<inward_normal units=\"" << units << "\" " << inward_normal << " />" << std::endl  
    //            << tabbing2 << "<surface_area units=\"square " << units << "\">" << surface_area << "</surface_area>" << std::endl  
    //            << tabbing  << "</voxel_face>"; 
    //    }
    public String stream_output_with_units(String units)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<voxel_face ID=\"" + mesh_index + "\">\n" );
        sb.append( tabbing2 + "<center units=\"" + units + "\" " + center + " />\n" );
        sb.append( tabbing2 + "<outward_normal units=\"" + units + "\" " + outward_normal + " />\n" );
        sb.append( tabbing2 + "<inward_normal units=\"" + units + "\" " + inward_normal + " />\n" );
        sb.append( tabbing2 + "<surface_area units=\"square " + units + "\">" + surface_area + "</surface_area>\n" );
        sb.append( tabbing + "</voxel_face>" );
        return sb.toString();
    }
}