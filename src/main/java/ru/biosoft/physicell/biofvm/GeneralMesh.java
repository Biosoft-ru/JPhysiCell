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
public class GeneralMesh
{
    public static final int minXIndex = 0;
    public static final int minYIndex = 1;
    public static final int minZIndex = 2;
    public static final int maxXIndex = 3;
    public static final int maxYIndex = 4;
    public static final int maxZIndex = 5;

    // this stores the indexing of the voxel faces (connect voxel i to voxel j, face stored at k)
    // only for use in a future release
    // std::unordered_map< int,std::unordered_map<int,int> > voxel_face_index_mapping; 

    // [xmin ymin zmin xmax ymax zmax ]
    public double[] boundingBox;

    public Voxel[] voxels;
    VoxelFace[] voxelFaces;

    // each voxel[k] has a list of connected voxels -- helpful for some numerical methods 
    int[][] connectedVoxelIndices;

    boolean cartesianMesh;
    boolean uniformMesh;
    boolean regularMesh;
    boolean useVoxelFaces;

    public String units;

    static String tabbing = "\t";
    static String tabbing2 = "\t\t";
    static String tabbing3 = "\t\t\t";

    public GeneralMesh()
    {
        // x1, x2, y1, y2, z1, z2 
        boundingBox = new double[6];//.assign(6,0.0); 
        boundingBox[minXIndex] = -0.5;
        boundingBox[minYIndex] = -0.5;
        boundingBox[minZIndex] = -0.5;
        boundingBox[maxXIndex] = 0.5;
        boundingBox[maxYIndex] = 0.5;
        boundingBox[maxZIndex] = 0.5;

        voxels = new Voxel[1];//resize( 1 );
        voxelFaces = new VoxelFace[0];////.resize( 0 );

        connectedVoxelIndices = new int[1][0];//resize( 1 );
        //        connected_voxel_indices[0].clear();

        // voxel_face_index_mapping.clear();

        cartesianMesh = false;
        uniformMesh = false;
        regularMesh = false;
        useVoxelFaces = true;
        units = "none";
    }

    public boolean isPositionValid(double x, double y, double z)
    {
        if( x < boundingBox[minXIndex] || x > boundingBox[maxXIndex] )
            return false;
        if( y < boundingBox[minYIndex] || y > boundingBox[maxYIndex] )
            return false;
        if( z < boundingBox[minZIndex] || z > boundingBox[maxZIndex] )
            return false;
        return true;
    }

    public void connect_voxels_faces_only(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!

        // create a new Voxel_Face connecting i to j
        VoxelFace VF1 = new VoxelFace();
        int k = voxelFaces.length;
        VF1.meshIndex = k;
        VF1.surfaceArea = SA;
        VF1.outward_normal = VectorUtil.newDiff( voxels[j].center, voxels[i].center );
        VectorUtil.normalize( VF1.outward_normal );
        VF1.inward_normal = VF1.outward_normal.clone();
        VectorUtil.prod( VF1.inward_normal, -1 );

        // convention: face is oriented from lower index to higher index 
        if( j < i )
        {
            VectorUtil.prod( VF1.outward_normal, -1 );
            VectorUtil.prod( VF1.inward_normal, -1 );
        }

        // add it to the vector of voxel faces 
        VectorUtil.push_back( voxelFaces, VF1 );
    }

    /*! This creates a Voxel_Face from voxels[i] to voxels[j], and another from voxels[j] to 
        voxels[i], both with surface area SA. It also auto-updates connected_voxel_indices[i] 
        and connected_voxel_indices[j]. */
    void connect_voxels(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!
        // create a new Voxel_Face connecting i to j
        VoxelFace VF1 = new VoxelFace();
        int k = voxelFaces.length;
        VF1.meshIndex = k;
        VF1.surfaceArea = SA;
        VF1.outward_normal = VectorUtil.newDiff( voxels[j].center, voxels[i].center );
        VectorUtil.normalize( VF1.outward_normal );
        VF1.inward_normal = VF1.outward_normal.clone();
        VectorUtil.prod( VF1.inward_normal, -1.0 );

        // convention: face is oriented from lower index to higher index 
        if( j < i )
        {
            VectorUtil.prod( VF1.outward_normal, -1.0 );
            VectorUtil.prod( VF1.inward_normal, -1.0 );
        }

        // add it to the vector of voxel faces 
        voxelFaces = VectorUtil.push_back( voxelFaces, VF1 );
        //        voxel_faces.push_back( VF1 ); 
    }

    void connect_voxels_indices_only(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!

        // add j to the list of connected voxels for voxel i 
        connectedVoxelIndices[i] = VectorUtil.push_back( connectedVoxelIndices[i], j );
        connectedVoxelIndices[j] = VectorUtil.push_back( connectedVoxelIndices[j], i );
        //        connected_voxel_indices[i].push_back( j ); 
        //        connected_voxel_indices[j].push_back( i ); 

        return;
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nMesh information: \n" );
        sb.append( "type: general mesh\n" );
        sb.append( "Domain: " );
        sb.append( "[" + boundingBox[0] + "," + boundingBox[3] + "] " + units + " x \n" );
        sb.append( "[" + boundingBox[1] + "," + boundingBox[4] + "] " + units + " x \n" );
        sb.append( "[" + boundingBox[2] + "," + boundingBox[5] + "] " + units + "\n" );
        sb.append( "   voxels: " + voxels.length + "\n" );
        sb.append( "   voxel faces: " + voxelFaces.length + "\n" );
        sb.append( "   volume: " );

        double totalVolume = 0.0;
        for( int i = 0; i < voxels.length; i++ )
        {
            totalVolume += voxels[i].volume;
        }
        sb.append( totalVolume + " cubic " + units + "\n" );

        return sb.toString();
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<mesh type=\"general\" uniform=\"" + uniformMesh + "\" regular=\"" + regularMesh + "\" units=\"" + units
                + "\">\n" );
        sb.append( tabbing2 + "<voxels>\n" );
        for( int i = 0; i < voxels.length; i++ )
        {
            sb.append( voxels[i] + "\n" );
        }
        sb.append( tabbing2 + "</voxels>\n" + tabbing2 + "<voxel_faces>\n" );
        for( int i = 0; i < voxelFaces.length; i++ )
        {
            sb.append( voxelFaces[i] + "\n" );
        }
        sb.append( tabbing2 + "</voxel_faces>\n" );

        sb.append( tabbing2 + "<voxel_connections>\n" );
        for( int i = 0; i < connectedVoxelIndices.length; i++ )
        {
            sb.append( tabbing3 + "<connected_voxel_indices ID=\"" + i + "\">\n" );
            for( int j = 0; j < connectedVoxelIndices[i].length; j++ )
            {
                sb.append( tabbing3 + "\t<index>" + ( connectedVoxelIndices[i] )[j] + "</index>\n" );
            }
            sb.append( tabbing3 + "</connected_voxel_indices>\n" );
        }
        sb.append( tabbing2 + "</voxel_connections>\n" );
        sb.append( tabbing + "</mesh>" );

        // later: output information on connected faces 
        return sb.toString();
    }

}