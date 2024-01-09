package ru.biosoft.physicell.biofvm;

import java.util.Arrays;

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
public class CartesianMesh extends GeneralMesh
{
    public double[] x_coordinates;
    public double[] y_coordinates;
    public double[] z_coordinates;
    public int[][] moore_connected_voxel_indices = new int[0][]; // Keeps the list of voxels in the Moore nighborhood
    public double dx;
    public double dy;
    public double dz;

    double dV;
    double dS;

    double dS_xy;
    double dS_yz;
    double dS_xz;

    static double tolerance = 1e-12;

    public CartesianMesh()
    {
        Cartesian_mesh = true;
        uniform_mesh = true;
        regular_mesh = true;
        use_voxel_faces = false;

        x_coordinates = new double[1];//.assign( 1 , 0.0 ); 
        y_coordinates = new double[1];//.assign( 1 , 0.0 ); 
        z_coordinates = new double[1];//.assign( 1 , 0.0 ); 

        dx = boundingBox[3] - boundingBox[0];
        dy = boundingBox[4] - boundingBox[1];
        dz = boundingBox[5] - boundingBox[2];


        if( Math.abs( dx - dy ) > tolerance || Math.abs( dy - dz ) > tolerance || Math.abs( dx - dz ) > tolerance )
        {
            uniform_mesh = false;
        }

        dV = dx * dy * dz;
        dS = dx * dy;
        dS_xy = dx * dy;
        dS_yz = dy * dz;
        dS_xz = dx * dz;

        voxels = new Voxel[x_coordinates.length * y_coordinates.length * z_coordinates.length];
        for( int i = 0; i < voxels.length; i++ )
        {
            voxels[i] = new Voxel();
            voxels[i].volume = dV;
        }
        voxels[0].center[0] = x_coordinates[0];
        voxels[0].center[1] = y_coordinates[0];
        voxels[0].center[2] = z_coordinates[0];
    }

    public CartesianMesh(int xnodes, int ynodes, int znodes)
    {
        x_coordinates = new double[xnodes];
        y_coordinates = new double[ynodes];
        z_coordinates = new double[znodes];

        dx = 1;
        dy = 1;
        dz = 1;

        dV = dx * dx * dz;
        dS = dx * dy;

        dS_xy = dS;
        dS_yz = dS;
        dS_xz = dS;

        uniform_mesh = true;
        regular_mesh = true;
        use_voxel_faces = false;

        for( int i = 0; i < x_coordinates.length; i++ )
        {
            x_coordinates[i] = i * dx;
        }
        for( int i = 0; i < y_coordinates.length; i++ )
        {
            y_coordinates[i] = i * dy;
        }
        for( int i = 0; i < z_coordinates.length; i++ )
        {
            z_coordinates[i] = i * dz;
        }

        boundingBox[0] = x_coordinates[0] - dx / 2.0;
        boundingBox[3] = x_coordinates[x_coordinates.length - 1] + dx / 2.0;
        boundingBox[1] = y_coordinates[0] - dy / 2.0;
        boundingBox[4] = y_coordinates[y_coordinates.length - 1] + dy / 2.0;
        boundingBox[2] = z_coordinates[0] - dz / 2.0;
        boundingBox[5] = z_coordinates[z_coordinates.length - 1] + dz / 2.0;

        units = "none";

        voxels = new Voxel[x_coordinates.length * y_coordinates.length * z_coordinates.length];
        for( int i = 0; i < voxels.length; i++ )
        {
            voxels[i] = new Voxel();
            voxels[i].volume = dV;
        }

        // initializing and connecting voxels 

        int n = 0;
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length; i++ )
                {
                    voxels[n].center[0] = x_coordinates[i];
                    voxels[n].center[1] = y_coordinates[j];
                    voxels[n].center[2] = z_coordinates[k];
                    voxels[n].mesh_index = n;
                    voxels[n].volume = dV;
                    n++;
                }
            }
        }

        // make connections 
        connected_voxel_indices = VectorUtil.resize( connected_voxel_indices, voxels.length );
        //                        connected_voxel_indices.resize(voxels.length);

        int i_jump = 1;
        int j_jump = x_coordinates.length;
        int k_jump = x_coordinates.length * y_coordinates.length;

        // x-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length - 1; i++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + i_jump, dS_yz );
                }
            }
        }
        // y-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int j = 0; j < y_coordinates.length - 1; j++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + j_jump, dS_xz );
                }
            }
        }
        // z-aligned connections 
        for( int j = 0; j < y_coordinates.length; j++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int k = 0; k < z_coordinates.length - 1; k++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + k_jump, dS_xy );
                }
            }
        }

        if( use_voxel_faces )
        {
            create_voxel_faces();
        }
    }

    void create_moore_neighborhood()
    {
        //                        moore_connected_voxel_indices.resize( voxels.length );
        moore_connected_voxel_indices = VectorUtil.resize( moore_connected_voxel_indices, voxels.length );

        for( int j = 0; j < y_coordinates.length; j++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int k = 0; k < z_coordinates.length; k++ )
                {
                    int center_inex = voxel_index( i, j, k );
                    for( int ii = -1; ii <= 1; ii++ )
                        for( int jj = -1; jj <= 1; jj++ )
                            for( int kk = -1; kk <= 1; kk++ )
                                if( i + ii >= 0 && i + ii < x_coordinates.length && j + jj >= 0 && j + jj < y_coordinates.length
                                        && k + kk >= 0 && k + kk < z_coordinates.length && ! ( ii == 0 && jj == 0 && kk == 0 ) )
                                {
                                    int neighbor_index = voxel_index( i + ii, j + jj, k + kk );
                                    moore_connected_voxel_indices[center_inex] = VectorUtil
                                            .push_back( moore_connected_voxel_indices[center_inex], neighbor_index );
                                }
                }
            }
        }
    }


    void create_voxel_faces()
    {
        // make connections 
        connected_voxel_indices = VectorUtil.resize( connected_voxel_indices, voxels.length );
        //         connected_voxel_indices.resize( voxels.length );

        int i_jump = 1;
        int j_jump = x_coordinates.length;
        int k_jump = x_coordinates.length * y_coordinates.length;

        // x-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length - 1; i++ )
                {
                    int n = voxel_index( i, j, k );
                    connect_voxels_faces_only( n, n + i_jump, dS_yz );
                }
            }
        }
        // y-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int j = 0; j < y_coordinates.length - 1; j++ )
                {
                    int n = voxel_index( i, j, k );
                    connect_voxels_faces_only( n, n + j_jump, dS_xz );
                }
            }
        }
        // z-aligned connections 
        for( int j = 0; j < y_coordinates.length; j++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int k = 0; k < z_coordinates.length - 1; k++ )
                {
                    int n = voxel_index( i, j, k );
                    connect_voxels_faces_only( n, n + k_jump, dS_xy );
                }
            }
        }
    }
    //        int nearest_voxel_face_index( std::vector<double>& position );  


    int voxel_index(int i, int j, int k)
    {
        return ( k * y_coordinates.length + j ) * x_coordinates.length + i;
    }

    int[] cartesian_indices(int n)
    {
        int[] out = new int[3];
        out[0] = out[1] = out[2] = -1;//TODO: maybe not needed;

        // figure out i; 
        int XY = x_coordinates.length * y_coordinates.length;
        out[2] = (int)Math.floor( n / XY );

        // figure out j; 
        out[1] = (int)Math.floor( ( n - out[2] * XY ) / x_coordinates.length );

        // figure out k; 
        out[0] = n - x_coordinates.length * ( out[1] + y_coordinates.length * out[2] );

        return out;
    }

    static double tol = 1e-16;

    void resize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, int x_nodes, int y_nodes,
            int z_nodes)
    {
        x_coordinates = new double[x_nodes];
        y_coordinates = new double[y_nodes];//.assign( y_nodes , 0.0 ); 
        z_coordinates = new double[z_nodes];//.assign( z_nodes , 0.0 ); 

        dx = ( x_end - x_start ) / ( (double)x_nodes );
        if( x_nodes < 2 )
        {
            dx = 1;
        }
        dy = ( y_end - y_start ) / ( (double)y_nodes );
        if( y_nodes < 2 )
        {
            dy = 1;
        }
        dz = ( z_end - z_start ) / ( (double)z_nodes );
        if( z_nodes < 2 )
        {
            dz = 1;
        }

        uniform_mesh = true;
        regular_mesh = true;

        if( Math.abs( dx - dy ) > tol && x_nodes > 1 && y_nodes > 1 )
        {
            uniform_mesh = false;
        }
        if( Math.abs( dy - dz ) > tol && y_nodes > 1 && z_nodes > 1 )
        {
            uniform_mesh = false;
        }
        if( Math.abs( dx - dz ) > tol && x_nodes > 1 && z_nodes > 1 )
        {
            uniform_mesh = false;
        }

        for( int i = 0; i < x_coordinates.length; i++ )
        {
            x_coordinates[i] = x_start + ( i + 0.5 ) * dx;
        }
        for( int i = 0; i < y_coordinates.length; i++ )
        {
            y_coordinates[i] = y_start + ( i + 0.5 ) * dy;
        }
        for( int i = 0; i < z_coordinates.length; i++ )
        {
            z_coordinates[i] = z_start + ( i + 0.5 ) * dz;
        }

        boundingBox[0] = x_start;
        boundingBox[3] = x_end;
        boundingBox[1] = y_start;
        boundingBox[4] = y_end;
        boundingBox[2] = z_start;
        boundingBox[5] = z_end;

        dV = dx * dy * dz;
        dS = dx * dy;

        dS_xy = dx * dy;
        dS_yz = dy * dz;
        dS_xz = dx * dz;

        voxels = new Voxel[x_coordinates.length * y_coordinates.length * z_coordinates.length];
        for( int i = 0; i < voxels.length; i++ )
        {
            voxels[i] = new Voxel();
            voxels[i].volume = dV;
        }

        int n = 0;
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length; i++ )
                {
                    voxels[n].center[0] = x_coordinates[i];
                    voxels[n].center[1] = y_coordinates[j];
                    voxels[n].center[2] = z_coordinates[k];
                    voxels[n].mesh_index = n;
                    voxels[n].volume = dV;

                    n++;
                }
            }
        }

        // make connections 
        connected_voxel_indices = Arrays.copyOf( connected_voxel_indices, voxels.length );
        //            connected_voxel_indices.resize( voxels.length );
        voxel_faces = new VoxelFace[0];//.clear(); 

        for( int i = 0; i < connected_voxel_indices.length; i++ )
        {
            //                connected_voxel_indices[i].clear();//TODO: check
        }

        int i_jump = 1;
        int j_jump = x_coordinates.length;
        int k_jump = x_coordinates.length * y_coordinates.length;

        // x-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length - 1; i++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + i_jump, dS_yz );
                }
            }
        }
        // y-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int j = 0; j < y_coordinates.length - 1; j++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + j_jump, dS_xz );
                }
            }
        }
        // z-aligned connections 
        for( int j = 0; j < y_coordinates.length; j++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int k = 0; k < z_coordinates.length - 1; k++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + k_jump, dS_xy );
                }
            }
        }

        if( use_voxel_faces )
        {
            create_voxel_faces();
        }

        create_moore_neighborhood();
        return;
    }

    public void resize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double dx_new,
            double dy_new,
            double dz_new)
    {
        dx = dx_new;
        dy = dy_new;
        dz = dz_new;

        double eps = 1e-16;
        int x_nodes = (int)Math.ceil( eps + ( x_end - x_start ) / dx );
        int y_nodes = (int)Math.ceil( eps + ( y_end - y_start ) / dy );
        int z_nodes = (int)Math.ceil( eps + ( z_end - z_start ) / dz );

        x_coordinates = new double[x_nodes];
        y_coordinates = new double[y_nodes];
        z_coordinates = new double[z_nodes];

        uniform_mesh = true;
        regular_mesh = true;
        double tol = 1e-16;
        if( Math.abs( dx - dy ) > tol || Math.abs( dy - dz ) > tol || Math.abs( dx - dz ) > tol )
        {
            uniform_mesh = false;
        }

        for( int i = 0; i < x_coordinates.length; i++ )
        {
            x_coordinates[i] = x_start + ( i + 0.5 ) * dx;
        }
        for( int i = 0; i < y_coordinates.length; i++ )
        {
            y_coordinates[i] = y_start + ( i + 0.5 ) * dy;
        }
        for( int i = 0; i < z_coordinates.length; i++ )
        {
            z_coordinates[i] = z_start + ( i + 0.5 ) * dz;
        }

        boundingBox[0] = x_start;
        boundingBox[3] = x_end;
        boundingBox[1] = y_start;
        boundingBox[4] = y_end;
        boundingBox[2] = z_start;
        boundingBox[5] = z_end;

        dV = dx * dy * dz;
        dS = dx * dy;

        dS_xy = dx * dy;
        dS_yz = dy * dz;
        dS_xz = dx * dz;

        voxels = new Voxel[x_coordinates.length * y_coordinates.length * z_coordinates.length];
        for( int i = 0; i < voxels.length; i++ )
        {
            voxels[i] = new Voxel();
            voxels[i].volume = dV;
        }

        int n = 0;
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length; i++ )
                {
                    voxels[n].center[0] = x_coordinates[i];
                    voxels[n].center[1] = y_coordinates[j];
                    voxels[n].center[2] = z_coordinates[k];
                    voxels[n].mesh_index = n;
                    voxels[n].volume = dV;

                    n++;
                }
            }
        }

        // make connections 

        connected_voxel_indices = VectorUtil.resize( connected_voxel_indices, voxels.length );
        //            connected_voxel_indices.resize( voxels.length );
        voxel_faces = new VoxelFace[0];//.clear(); 

        for( int i = 0; i < connected_voxel_indices.length; i++ )
        {
            //connected_voxel_indices[i].clear();//TODO: check
        }

        int i_jump = 1;
        int j_jump = x_coordinates.length;
        int k_jump = x_coordinates.length * y_coordinates.length;

        // x-aligned connections 
        int count = 0;
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int j = 0; j < y_coordinates.length; j++ )
            {
                for( int i = 0; i < x_coordinates.length - 1; i++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + i_jump, dS_yz );
                    count++;
                }
            }
        }

        // y-aligned connections 
        for( int k = 0; k < z_coordinates.length; k++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int j = 0; j < y_coordinates.length - 1; j++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + j_jump, dS_xz );
                }
            }
        }

        // z-aligned connections 
        for( int j = 0; j < y_coordinates.length; j++ )
        {
            for( int i = 0; i < x_coordinates.length; i++ )
            {
                for( int k = 0; k < z_coordinates.length - 1; k++ )
                {
                    int n1 = voxel_index( i, j, k );
                    connect_voxels_indices_only( n1, n1 + k_jump, dS_xy );
                }
            }
        }

        if( use_voxel_faces )
        {
            create_voxel_faces();
        }

        create_moore_neighborhood();
        return;
    }

    void resize(int x_nodes, int y_nodes, int z_nodes)
    {
        resize( 0 - .5, x_nodes - 1 + .5, 0 - .5, y_nodes - 1 + .5, 0 - .5, z_nodes - 1 + .5, x_nodes, y_nodes, z_nodes );

    }

    void resize_uniform(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double dx_new)
    {
        resize( x_start, x_end, y_start, y_end, z_start, z_end, dx_new, dx_new, dx_new );
    }

    public int nearest_voxel_index(double[] position)
    {
        int i = (int)Math.floor( ( position[0] - boundingBox[0] ) / dx );
        int j = (int)Math.floor( ( position[1] - boundingBox[1] ) / dy );
        int k = (int)Math.floor( ( position[2] - boundingBox[2] ) / dz );
        //  add some bounds checking -- truncate to inside the computational domain   
        if( i >= x_coordinates.length )
        {
            i = x_coordinates.length - 1;
        }
        if( i < 0 )
        {
            i = 0;
        }

        if( j >= y_coordinates.length )
        {
            j = y_coordinates.length - 1;
        }
        if( j < 0 )
        {
            j = 0;
        }

        if( k >= z_coordinates.length )
        {
            k = z_coordinates.length - 1;
        }
        if( k < 0 )
        {
            k = 0;
        }

        return ( k * y_coordinates.length + j ) * x_coordinates.length + i;
    }

    int[] nearest_cartesian_indices(double[] position)
    {
        int[] out = new int[3];

        out[0] = (int)Math.floor( ( position[0] - boundingBox[0] ) / dx );
        out[1] = (int)Math.floor( ( position[1] - boundingBox[1] ) / dy );
        out[2] = (int)Math.floor( ( position[2] - boundingBox[2] ) / dz );

        //  add some bounds checking -- truncate to inside the computational domain  

        if( out[0] >= x_coordinates.length )
        {
            out[0] = x_coordinates.length - 1;
        }
        if( out[0] < 0 )
        {
            out[0] = 0;
        }

        if( out[1] >= y_coordinates.length )
        {
            out[1] = y_coordinates.length - 1;
        }
        if( out[1] < 0 )
        {
            out[1] = 0;
        }

        if( out[2] >= z_coordinates.length )
        {
            out[2] = z_coordinates.length - 1;
        }
        if( out[2] < 0 )
        {
            out[2] = 0;
        }

        return out;
    }

    Voxel nearest_voxel(double[] position)
    {
        return voxels[nearest_voxel_index( position )];
    }

    //        void Cartesian_Mesh::display_information( std::ostream& os )
    //        {
    //            os << std::endl << "Mesh information: " << std::endl;
    //            if( uniform_mesh ) 
    //            { os << "type: uniform Cartesian" << std::endl; }
    //            else
    //            {
    //                if( regular_mesh )
    //                { os << "type: regular Cartesian" << std::endl; }
    //                else
    //                { os << "type: general Cartesian" << std::endl; }
    //            }
    //            os << "Domain: " 
    //            << "[" << bounding_box[0] << "," << bounding_box[3] << "] " <<  units << " x " 
    //            << "[" << bounding_box[1] << "," << bounding_box[4] << "] " <<  units << " x " 
    //            << "[" << bounding_box[2] << "," << bounding_box[5] << "] " <<  units << std::endl
    //            << "   resolution: dx = " << dx << " " << units; 
    //            if( !uniform_mesh )
    //            {
    //                os  << ", dy = " << dy << " " << units 
    //                    << ", dz = " << dz << " " << units ; 
    //            }
    //            os << std::endl 
    //            << "   voxels: " << voxels.size() << std::endl
    //            << "   voxel faces: " << voxel_faces.size() << std::endl
    //            << "   volume: " << ( bounding_box[3]-bounding_box[0] )*( bounding_box[4]-bounding_box[1] )*( bounding_box[5]-bounding_box[2] ) 
    //                << " cubic " << units << std::endl;     
    //
    //            return; 
    //        }
    @Override
    public String display_information()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nMesh information\n" );
        if( uniform_mesh )
        {
            sb.append( "type: uniform Cartesian\n" );
        }
        else
        {
            if( regular_mesh )
            {
                sb.append( "type: regular Cartesian\n" );
            }
            else
            {
                sb.append( "type: general Cartesian\n" );
            }
        }
        sb.append( "Domain: " );
        sb.append( "[" + boundingBox[0] + "," + boundingBox[3] + "] " + units + " x " );
        sb.append( "[" + boundingBox[1] + "," + boundingBox[4] + "] " + units + " x " );
        sb.append( "[" + boundingBox[2] + "," + boundingBox[5] + "] " + units + "\n" );
        sb.append( "   resolution: dx = " + dx + " " + units );
        if( !uniform_mesh )
        {
            sb.append( ", dy = " + dy + " " + units );
            sb.append( ", dz = " + dz + " " + units );
        }
        sb.append( "\n" );
        sb.append( "   voxels: " + voxels.length + "\n" );
        sb.append( "   voxel faces: " + voxel_faces.length + "\n" );
        sb.append( "   volume: "
                + ( boundingBox[3] - boundingBox[0] ) * ( boundingBox[4] - boundingBox[1] ) * ( boundingBox[5] - boundingBox[2] ) );
        sb.append( " cubic " + units + "\n" );

        return sb.toString();
    }


    //
    //        void read_from_matlab(String filename )
    //        {
    //            int size_of_each_datum; 
    //            int number_of_data_entries; 
    //            FILE* fp = read_matlab_header( &size_of_each_datum, &number_of_data_entries,  filename ); 
    //
    //            voxel_faces.resize( 0 ); 
    //            
    //            connected_voxel_indices.resize( 1 ); 
    //            connected_voxel_indices[0].clear(); 
    //                
    //            Cartesian_mesh = true; 
    //            uniform_mesh = false; 
    //            regular_mesh = true; 
    //            use_voxel_faces = false; 
    //
    //            // resize the internal data structure 
    //
    //            voxels.resize( number_of_data_entries );
    //            connected_voxel_indices.resize( voxels.length ); 
    //            
    //            x_coordinates.resize( number_of_data_entries );
    //            y_coordinates.resize( number_of_data_entries );
    //            z_coordinates.resize( number_of_data_entries );
    //            
    //            // read in the data
    //            // assumes each column has: x,y,z, dV
    //            
    //            bounding_box[0] = 9e99; 
    //            bounding_box[1] = 9e99; 
    //            bounding_box[2] = 9e99; 
    //
    //            bounding_box[3] = -9e99; 
    //            bounding_box[4] = -9e99; 
    //            bounding_box[5] = -9e99; 
    //         
    //                size_t result;
    //            for( unsigned int i=0; i < number_of_data_entries ; i++ )
    //            {
    //                result = fread( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
    //                result = fread( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
    //                result = fread( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
    //                result = fread( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
    //                
    //                // estimate the bounding box; 
    //                if( voxels[i].center[0] < bounding_box[0] )
    //                { bounding_box[0] = voxels[i].center[0]; }
    //                if( voxels[i].center[0] > bounding_box[3] )
    //                { bounding_box[3] = voxels[i].center[0]; }
    //
    //                if( voxels[i].center[1] < bounding_box[1] )
    //                { bounding_box[1] = voxels[i].center[1]; }
    //                if( voxels[i].center[1] > bounding_box[4] )
    //                { bounding_box[4] = voxels[i].center[1]; }
    //
    //                if( voxels[i].center[2] < bounding_box[2] )
    //                { bounding_box[2] = voxels[i].center[2]; }
    //                if( voxels[i].center[2] > bounding_box[5] )
    //                { bounding_box[5] = voxels[i].center[2]; }
    //            } 
    //            
    //            // figure out dx, dy, dz 
    //
    //            double xmin = bounding_box[0]; // voxels[0].center[0]; 
    //            double ymin = bounding_box[1]; // voxels[0].center[1]; 
    //            double zmin = bounding_box[2]; // voxels[0].center[2]; 
    //
    //            // int n = voxels.size(); 
    //            double xmax = bounding_box[3]; // voxels[n-1].center[0]; 
    //            double ymax = bounding_box[4]; // voxels[n-1].center[1]; 
    //            double zmax = bounding_box[5]; // voxels[n-1].center[2]; 
    //
    //            // figure out number of x nodes  
    //            int xnodes = 0; 
    //            while( fabs( voxels[xnodes].center[0] - xmax ) > 1e-15 )
    //            { xnodes++; }
    //            xnodes++; 
    //
    //            // figure out number of y nodes 
    //            int ynodes = 0; 
    //
    //            while( fabs( voxels[ynodes*xnodes].center[1] - ymax ) > 1e-15 )
    //            { ynodes += 1; }
    //            ynodes++;
    //
    //            // figure out number of z nodes 
    //
    //            int znodes = 0; 
    //
    //            while( fabs( voxels[ynodes*xnodes*znodes].center[2] - zmax ) > 1e-15 )
    //            { znodes += 1; }
    //            znodes++;
    //
    //            // figure out differentials
    //
    //            dx = ( xmax - xmin ) / ( (double) xnodes - 1.0  ); 
    //            dy = ( ymax - ymin ) / ( (double) ynodes - 1.0  ); 
    //            dz = ( zmax - zmin ) / ( (double) znodes - 1.0  ); 
    //            
    //            dV = dx*dy*dz; 
    //            dS = dx*dy;
    //
    //            dS_xy = dx*dy;
    //            dS_yz = dy*dz; 
    //            dS_xz = dx*dz;
    //
    //            uniform_mesh = true; 
    //            double tol = 1e-16; 
    //            if( fabs( dx - dy ) > tol || fabs( dy - dz ) > tol || fabs( dx - dz ) > tol )
    //            { uniform_mesh = false; }
    //            
    //            // correct the bounding box 
    //            
    //            double half_step = dx * 0.5; 
    //            
    //            bounding_box[0] -= half_step; 
    //            bounding_box[3] += half_step;
    //            
    //            half_step = dy * 0.5; 
    //            bounding_box[1] -= half_step; 
    //            bounding_box[4] += half_step;
    //            
    //            half_step = dz * 0.5; 
    //            bounding_box[2] -= half_step; 
    //            bounding_box[5] += half_step;
    //         
    //            // write out the x,y,z coordinates; 
    //            x_coordinates.resize( xnodes ); 
    //            y_coordinates.resize( ynodes ); 
    //            z_coordinates.resize( znodes ); 
    //
    //            for(  int i=0; i < x_coordinates.length ; i++ )
    //            { x_coordinates[i] = xmin + i*dx ;   }
    //
    //            for(  int i=0; i < y_coordinates.length ; i++ )
    //            { y_coordinates[i] = ymin + i*dy ; }
    //
    //            for(  int i=0; i < z_coordinates.length ; i++ )
    //            { z_coordinates[i] = zmin + i*dz ; }
    //
    //            dV = dx*dy*dz; 
    //            
    //            units = "none";  
    //            
    //            // still need to figure out connected indices later 
    //            std::cout << "Warning: Cartesian_Mesh::read_from_matlab is incomplete. No connection information read." << std::endl; 
    //
    //            fclose( fp) ; 
    //            return; 
    //        }
}
