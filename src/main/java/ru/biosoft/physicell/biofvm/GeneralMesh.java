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
    public static final int mesh_min_x_index = 0;
    public static final int mesh_min_y_index = 1;
    public static final int mesh_min_z_index = 2;
    public static final int mesh_max_x_index = 3;
    public static final int mesh_max_y_index = 4;
    public static final int mesh_max_z_index = 5;

    // this stores the indexing of the voxel faces (connect voxel i to voxel j, face stored at k)
    // only for use in a future release
    // std::unordered_map< int,std::unordered_map<int,int> > voxel_face_index_mapping; 

    // [xmin ymin zmin xmax ymax zmax ]
    //    std::vector<double>bounding_box;
    public double[] boundingBox;

    //    std::vector<Voxel> voxels;
    //    std::vector<Voxel_Face> voxel_faces;
    public Voxel[] voxels;
    VoxelFace[] voxel_faces;

    // each voxel[k] has a list of connected voxels -- helpful for some numerical methods 
    //    std::vector<std::vector<int>>connected_voxel_indices;
    int[][] connected_voxel_indices;

    boolean Cartesian_mesh;
    boolean uniform_mesh;
    boolean regular_mesh;
    boolean use_voxel_faces;

    public String units;

    static String tabbing = "\t";
    static String tabbing2 = "\t\t";
    static String tabbing3 = "\t\t\t";


    public GeneralMesh()
    {
        // x1, x2, y1, y2, z1, z2 
        boundingBox = new double[6];//.assign(6,0.0); 
        boundingBox[mesh_min_x_index] = -0.5;
        boundingBox[mesh_min_y_index] = -0.5;
        boundingBox[mesh_min_z_index] = -0.5;
        boundingBox[mesh_max_x_index] = 0.5;
        boundingBox[mesh_max_y_index] = 0.5;
        boundingBox[mesh_max_z_index] = 0.5;

        voxels = new Voxel[1];//resize( 1 );
        voxel_faces = new VoxelFace[0];////.resize( 0 );

        connected_voxel_indices = new int[1][0];//resize( 1 );
        //        connected_voxel_indices[0].clear();

        // voxel_face_index_mapping.clear();

        Cartesian_mesh = false;
        uniform_mesh = false;
        regular_mesh = false;
        use_voxel_faces = true;
        units = "none";
    }

    public boolean isPositionValid(double x, double y, double z)
    {
        if( x < boundingBox[mesh_min_x_index] || x > boundingBox[mesh_max_x_index] )
            return false;
        if( y < boundingBox[mesh_min_y_index] || y > boundingBox[mesh_max_y_index] )
            return false;
        if( z < boundingBox[mesh_min_z_index] || z > boundingBox[mesh_max_z_index] )
            return false;
        return true;
    }

    /* the following help manage the voxel faces */

    // returns the index of the voxel face connecting from voxels[i] to voxels[j] 
    //    int voxel_face_index(int i, int j)
    //    {
    //
    //    }

    public void connect_voxels_faces_only(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!

        // create a new Voxel_Face connecting i to j

        VoxelFace VF1 = new VoxelFace();
        int k = voxel_faces.length;
        VF1.mesh_index = k;
        VF1.surface_area = SA;
        VF1.outward_normal = VectorUtil.newDiff( voxels[j].center, voxels[i].center );
        VectorUtil.normalize( VF1.outward_normal );
        VF1.inward_normal = VF1.outward_normal.clone();
        //        VF1.inward_normal *= -1.0;
        VectorUtil.prod( VF1.inward_normal, -1 );

        // convention: face is oriented from lower index to higher index 
        if( j < i )
        {
            //            VF1.outward_normal *= -1.0;
            //            VF1.inward_normal *= -1.0;
            VectorUtil.prod( VF1.outward_normal, -1 );
            VectorUtil.prod( VF1.inward_normal, -1 );
        }

        // add it to the vector of voxel faces 
        VectorUtil.push_back( voxel_faces, VF1 );
        //        voxel_faces.push_back( VF1 );
    }

    // returns the Voxel_Face connecting voxels[i] to voxels[j] 
    //    Voxel_Face voxel_face(int i, int j)
    //    {
    //
    //    }
    // returns the normal vector from voxels[i] to voxels[j] 
    //    double[] outward_normal(int i, int j)
    //    {
    //
    //    }

    /*! This creates a Voxel_Face from voxels[i] to voxels[j], and another from voxels[j] to 
        voxels[i], both with surface area SA. It also auto-updates connected_voxel_indices[i] 
        and connected_voxel_indices[j]. */
    void connect_voxels(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!

        // create a new Voxel_Face connecting i to j

        VoxelFace VF1 = new VoxelFace();
        int k = voxel_faces.length;
        VF1.mesh_index = k;
        VF1.surface_area = SA;
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
        VectorUtil.push_back( voxel_faces, VF1 );
        //        voxel_faces.push_back( VF1 ); 
    }

    void connect_voxels_indices_only(int i, int j, double SA)
    {
        // check to see if the voxels are connected -- implement later!

        // add j to the list of connected voxels for voxel i 
        VectorUtil.push_back( connected_voxel_indices[i], j );
        VectorUtil.push_back( connected_voxel_indices[j], i );
        //        connected_voxel_indices[i].push_back( j ); 
        //        connected_voxel_indices[j].push_back( i ); 

        return;
    }

    /*! This removes all connections between voxels[i] and voxels[j], and deletes the associated 
        Voxel_Face(s). */
    void disconnect_voxels(int i, int j)
    {

    }

    void clear_voxel_face_index_mapping()
    {

    }


    //    void display_information( std::ostream& os);

    public void write_to_matlab(String filename)
    {
        //        unsigned int number_of_data_entries = voxels.size();
        //        unsigned int size_of_each_datum = 3 + 1; // x,y,z, volume 
        //
        //        FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "mesh" );  
        //
        //        // storing data as cols 
        //        for( unsigned int i=0; i < number_of_data_entries ; i++ )
        //        {
        //            fwrite( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
        //            fwrite( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
        //            fwrite( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
        //            fwrite( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
        //        }
        //
        //        fclose( fp ); 
    }

    public void read_from_matlab(String filename)
    {
        //        unsigned int size_of_each_datum; 
        //        unsigned int number_of_data_entries; 
        //        FILE* fp = read_matlab_header( &size_of_each_datum, &number_of_data_entries,  filename ); 
        //
        //        voxel_faces.resize( 0 ); 
        //        
        //        connected_voxel_indices.resize( 1 ); 
        //        connected_voxel_indices[0].clear(); 
        //        
        //        Cartesian_mesh = false; 
        //        uniform_mesh = false; 
        //        regular_mesh = false; 
        //        use_voxel_faces = false; 
        //
        //        // resize the internal data structure 
        //
        //        voxels.resize( number_of_data_entries );
        //        connected_voxel_indices.resize( voxels.size() ); 
        //        
        //        // read in the data
        //        // assumes each column has: x,y,z, dV
        //        
        //        bounding_box[0] = 9e99; 
        //        bounding_box[1] = 9e99; 
        //        bounding_box[2] = 9e99; 
        //
        //        bounding_box[3] = -9e99; 
        //        bounding_box[4] = -9e99; 
        //        bounding_box[5] = -9e99; 
        //     
        //            size_t result;
        //        for( unsigned int i=0; i < number_of_data_entries ; i++ )
        //        {
        //            result = fread( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
        //            result = fread( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
        //            result = fread( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
        //            result = fread( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
        //            
        //            // estimate the bounding box; 
        //            if( voxels[i].center[0] < bounding_box[0] )
        //            { bounding_box[0] = voxels[i].center[0]; }
        //            if( voxels[i].center[0] > bounding_box[3] )
        //            { bounding_box[3] = voxels[i].center[0]; }
        //
        //            if( voxels[i].center[1] < bounding_box[1] )
        //            { bounding_box[1] = voxels[i].center[1]; }
        //            if( voxels[i].center[1] > bounding_box[4] )
        //            { bounding_box[4] = voxels[i].center[1]; }
        //
        //            if( voxels[i].center[2] < bounding_box[2] )
        //            { bounding_box[2] = voxels[i].center[2]; }
        //            if( voxels[i].center[2] > bounding_box[5] )
        //            { bounding_box[5] = voxels[i].center[2]; }
        //        } 
        //        
        //        std::cout << "Warning: General_Mesh::read_from_matlab is incomplete. No connection information read." << std::endl; 
        //
        //        fclose( fp) ; 
        //        return; 
    }

    public String display_information()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nMesh information: \n" );
        sb.append( "type: general mesh\n" );
        sb.append( "Domain: " );
        sb.append( "[" + boundingBox[0] + "," + boundingBox[3] + "] " + units + " x \n" );
        sb.append( "[" + boundingBox[1] + "," + boundingBox[4] + "] " + units + " x \n" );
        sb.append( "[" + boundingBox[2] + "," + boundingBox[5] + "] " + units + "\n" );
        sb.append( "   voxels: " + voxels.length + "\n" );
        sb.append( "   voxel faces: " + voxel_faces.length + "\n" );
        sb.append( "   volume: " );

        double total_volume = 0.0;
        for( int i = 0; i < voxels.length; i++ )
        {
            total_volume += voxels[i].volume;
        }
        sb.append( total_volume + " cubic " + units + "\n" );

        return sb.toString();
    }

    //    std::ostream& operator<<(std::ostream& os, const General_Mesh& mesh)  
    //    {
    //        std::boolalpha( os ); 
    //        static std::string tabbing = "\t"; 
    //        static std::string tabbing2 = "\t\t"; 
    //        static std::string tabbing3 = "\t\t\t"; 
    //        os  << tabbing << "<mesh type=\"general\" uniform=\"" << mesh.uniform_mesh << "\" regular=\"" << mesh.regular_mesh << "\" units=\"" << mesh.units << "\">"  << std::endl
    //            << tabbing2 << "<voxels>" << std::endl;
    //        for( unsigned int i=0; i < mesh.voxels.size() ; i++ )
    //        { os << mesh.voxels[i] << std::endl; } 
    //        os  << tabbing2 << "</voxels>" << std::endl 
    //            << tabbing2 << "<voxel_faces>" << std::endl;
    //        for( unsigned int i=0; i < mesh.voxel_faces.size() ; i++ )
    //        { os << mesh.voxel_faces[i] << std::endl; } 
    //        os  << tabbing2 << "</voxel_faces>" << std::endl; 
    //        
    //        os  << tabbing2 << "<voxel_connections>" << std::endl;
    //        for( unsigned int i=0 ; i < mesh.connected_voxel_indices.size() ; i++ )
    //        {
    //            os << tabbing3 << "<connected_voxel_indices ID=\"" << i << "\">" << std::endl; 
    //            for( unsigned int j=0; j < mesh.connected_voxel_indices[i].size() ; j++ )
    //            {
    //                os  << tabbing3 << "\t<index>" << (mesh.connected_voxel_indices[i])[j] << "</index>" << std::endl; 
    //            }
    //            os << tabbing3 << "</connected_voxel_indices>" << std::endl; 
    //        }
    //        
    //        os << tabbing2 << "</voxel_connections>" << std::endl; 
    //        
    //        
    //        os  << tabbing  << "</mesh>"; 
    //            
    //            // later: output information on connected faces 
    //            
    //            // later: 
    //
    //     return os; 
    //    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( tabbing + "<mesh type=\"general\" uniform=\"" + uniform_mesh + "\" regular=\"" + regular_mesh + "\" units=\"" + units
                + "\">\n" );
        sb.append( tabbing2 + "<voxels>\n" );
        for( int i = 0; i < voxels.length; i++ )
        {
            sb.append( voxels[i] + "\n" );
        }
        sb.append( tabbing2 + "</voxels>\n" + tabbing2 + "<voxel_faces>\n" );
        for( int i = 0; i < voxel_faces.length; i++ )
        {
            sb.append( voxel_faces[i] + "\n" );
        }
        sb.append( tabbing2 + "</voxel_faces>\n" );

        sb.append( tabbing2 + "<voxel_connections>\n" );
        for( int i = 0; i < connected_voxel_indices.length; i++ )
        {
            sb.append( tabbing3 + "<connected_voxel_indices ID=\"" + i + "\">\n" );
            for( int j = 0; j < connected_voxel_indices[i].length; j++ )
            {
                sb.append( tabbing3 + "\t<index>" + ( connected_voxel_indices[i] )[j] + "</index>\n" );
            }
            sb.append( tabbing3 + "</connected_voxel_indices>\n" );
        }
        sb.append( tabbing2 + "</voxel_connections>\n" );
        sb.append( tabbing + "</mesh>" );

        // later: output information on connected faces 
        return sb.toString();
    }

}