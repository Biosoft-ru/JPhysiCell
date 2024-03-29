package ru.biosoft.physicell.core;

import java.util.HashMap;
import java.util.Map;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
public class PhysiCellConstants
{

    static double cell_removal_threshold_volume = 20; // 20 cubic microns -- about 1% of typical cell

    static int keep_pushed_out_cells_in_outer_voxel = 1;
    static int solid_boundary = 2;
    static int default_boundary_condition_for_pushed_out_agents = keep_pushed_out_cells_in_outer_voxel;

    public static int deterministic_necrosis = 0;
    static int stochastic_necrosis = 1;

    public static int mesh_min_x_index = 0;
    public static int mesh_min_y_index = 1;
    public static int mesh_min_z_index = 2;
    public static int mesh_max_x_index = 3;
    public static int mesh_max_y_index = 4;
    public static int mesh_max_z_index = 5;
    static int custom_cycle_model = 9999;

    static int mesh_lx_face_index = 0;
    static int mesh_ly_face_index = 1;
    static int mesh_lz_face_index = 2;
    static int mesh_ux_face_index = 3;
    static int mesh_uy_face_index = 4;
    static int mesh_uz_face_index = 5;

    // currently recognized cell cycle models

    public static int advanced_Ki67_cycle_model = 0;
    public static int basic_Ki67_cycle_model = 1;
    public static int flow_cytometry_cycle_model = 2;
    public static int live_apoptotic_cycle_model = 3;
    public static int total_cells_cycle_model = 4;
    public static int live_cells_cycle_model = 5;
    public static int flow_cytometry_separated_cycle_model = 6;
    public static int cycling_quiescent_model = 7;

    // currently recognized death models

    public static int apoptosis_death_model = 100;
    public static int necrosis_death_model = 101;
    public static int autophagy_death_model = 102;



    // currently recognized cell cycle and death phases // cycle phases

    public static int Ki67_positive_premitotic = 0;
    public static int Ki67_positive_postmitotic = 1;
    public static int Ki67_positive = 2;
    public static int Ki67_negative = 3;
    public static int G0G1_phase = 4;
    public static int G0_phase = 5;
    public static int G1_phase = 6;
    public static int G1a_phase = 7;
    public static int G1b_phase = 8;
    public static int G1c_phase = 9;
    public static int S_phase = 10;
    public static int G2M_phase = 11;
    public static int G2_phase = 12;
    public static int M_phase = 13;
    public static int live = 14;


    public static int G1pm_phase = 15;
    public static int G1ps_phase = 16;


    public static int cycling = 17;
    public static int quiescent = 18;


    static int custom_phase = 9999; // death phases

    public static int apoptotic = 100;
    public static int necrotic_swelling = 101;
    public static int necrotic_lysed = 102;
    public static int necrotic = 103;
    public static int debris = 104;

    private static Map<String, Integer> cycle_model_codes = new HashMap<>();
    {
        //        {
            cycle_model_codes.put( "Ki67 (advanced)", advanced_Ki67_cycle_model );
            cycle_model_codes.put( "Ki67 (basic)", basic_Ki67_cycle_model );
            cycle_model_codes.put( "Flow cytometry model (basic)", flow_cytometry_cycle_model );
            // { ,PhysiCell_constants::live_apoptotic_cycle_model}, // not implemented 
            // { ,PhysiCell_constants::total_cells_cycle_model}, // not implemented 
            cycle_model_codes.put( "Live", live_cells_cycle_model );
            cycle_model_codes.put( "Flow cytometry model (separated)", flow_cytometry_separated_cycle_model );
            cycle_model_codes.put( "Cycling-Quiescent model", cycling_quiescent_model );

            // currently recognized death models 
            cycle_model_codes.put( "Apoptosis", apoptosis_death_model );
            cycle_model_codes.put( "Necrosis", necrosis_death_model );
            // { ,PhysiCell_constants::autophagy_death_model}, // not implemented 

            cycle_model_codes.put( "ki67 (advanced)", advanced_Ki67_cycle_model );
            cycle_model_codes.put( "ki67 (basic)", basic_Ki67_cycle_model );
            cycle_model_codes.put( "flow cytometry model (basic)", flow_cytometry_cycle_model );
            cycle_model_codes.put( "live", live_cells_cycle_model );
            cycle_model_codes.put( "flow cytometry model (separated)", flow_cytometry_separated_cycle_model );
            cycle_model_codes.put( "cycling-quiescent model", cycling_quiescent_model );
            cycle_model_codes.put( "apoptosis", apoptosis_death_model );
            cycle_model_codes.put( "necrosis", necrosis_death_model );
            //        }
    }

    public static int find_cycle_model_code(String model_name)
    {
        if( cycle_model_codes.containsKey( model_name ) )
            return cycle_model_codes.get( model_name );
        return -1;
    }
}
